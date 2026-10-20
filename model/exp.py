"""Model-improvement experiments. All predictions are returned at 20 bp so every variant is scored on the same bins.
Options: --res 20|50  --arch dilated|unet  --soft 0|1  --spe --epochs --seed --tracks --spec"""
import sys, os, time, json, argparse, numpy as np, torch, torch.nn as nn, torch.nn.functional as F
sys.path[:0] = ['.', 'boundary_model']
import pfprep as P, seg_features as SF, segmodel as SM
from bm.models import ResBlock

# ---------------- resolution ----------------
def to_res(a, res):
    """(C, n20) array at 20 bp -> (C, n_res) exact average (via 10-bp expansion)."""
    if res == 20: return a
    f = res // 10; a10 = np.repeat(a, 2, axis=-1); n = a10.shape[-1] // f
    return a10[..., :n * f].reshape(*a.shape[:-1], n, f).mean(-1).astype(a.dtype)
def lab_to_res(l, res):
    if res == 20: return l
    f = res // 10; l10 = np.repeat(l, 2, axis=-1); n = l10.shape[-1] // f
    return l10[..., :n * f].reshape(*l.shape[:-1], n, f)[..., f // 2]
def pred_to20(p, res, n20):
    if res == 20: return p
    idx = np.minimum((np.arange(n20) * 20 + 10) // res, p.shape[-1] - 1); return p[..., idx]

# ---------------- stranded RNA evidence ----------------
def stranded(spec=None):
    """spec tokens: 'strand' -> 9 inferred channels per strand (exon, intron, polyA x 3 stages); 'daft' -> 3 DAFT-seq stage channels per strand.
    Returned rows: + block (all + channels) then - block in the same order, so the flip augmentation swaps the two blocks."""
    spec = {'strand'} if spec is None else set(spec)
    blocks = []
    if 'strand' in spec: blocks.append(_stranded_inferred())
    if 'daft' in spec: blocks.append(_stranded_daft())
    ks = set.intersection(*[set(b) for b in blocks]); out = {}
    for k in ks:
        n = min(b[k].shape[1] for b in blocks); h = [b[k].shape[0] // 2 for b in blocks]
        out[k] = np.concatenate([b[k][:hh, :n] for b, hh in zip(blocks, h)] + [b[k][hh:, :n] for b, hh in zip(blocks, h)])
    return out
def _stranded_daft():
    z = np.load("data/daft_stage20.npz"); S3 = P.STAGES_REAL; raw = {}
    for c in P.CHROMS: raw[c] = np.concatenate([z[f"{s}__{c}"][0:1].astype(np.float32) for s in S3] + [z[f"{s}__{c}"][1:2].astype(np.float32) for s in S3])
    tr = np.concatenate([raw[c] for c in P.TRAIN], 1); mu = 0.5 * (tr[:3].mean(1, keepdims=True) + tr[3:].mean(1, keepdims=True)); sd = 0.5 * (tr[:3].std(1, keepdims=True) + tr[3:].std(1, keepdims=True)) + 1e-6
    mu, sd = np.concatenate([mu, mu]), np.concatenate([sd, sd])
    return {f"pool|{c}": ((raw[c] - mu) / sd).astype(np.float16) for c in P.CHROMS}
def _stranded_inferred():
    """exon/intron coverage of motif-stranded spliced reads + poly(A) sites; log1p, z-scored on training chromosomes. -> key -> (6, n): +exon,+intron,+polyA,-exon,-intron,-polyA"""
    out = {}
    for s, S_ in (("ring", "Ring"), ("troph", "Trophozoite"), ("schizont", "Schizont")):
        z = np.load(f"data/stranded_final_{S_}.npz"); raw = {}
        for c in P.CHROMS:
            e, i, a = (np.log1p(z[f"{c}__{t}"].astype(np.float32)) for t in ("exon", "intron", "polyA"))
            raw[c] = np.stack([e[0], i[0], a[0], e[1], i[1], a[1]])
        tr = np.concatenate([raw[c] for c in P.TRAIN], 1); mu = tr.mean(1, keepdims=True); sd = tr.std(1, keepdims=True) + 1e-6
        mu[3:], sd[3:] = mu[:3], sd[:3]  # same scaling for both strands keeps the pair exchangeable under flip
        m2 = 0.5 * (tr[:3].mean(1, keepdims=True) + tr[3:].mean(1, keepdims=True)); s2 = 0.5 * (tr[:3].std(1, keepdims=True) + tr[3:].std(1, keepdims=True)) + 1e-6
        mu = np.concatenate([m2, m2]); sd = np.concatenate([s2, s2])
        for c in P.CHROMS: out[f"{s}|{c}"] = ((raw[c] - mu) / sd).astype(np.float16)
    if P.STAGES == ["pool"]:   # + block = ring+,troph+,schizont+ ; - block in the same order -> flip swaps blocks
        return {f"pool|{c}": np.concatenate([out[f"{s}|{c}"][:3] for s in P.STAGES_REAL] + [out[f"{s}|{c}"][3:] for s in P.STAGES_REAL]) for c in P.CHROMS}
    return out

# ---------------- soft targets ----------------
def soft_targets(L, sigma_bins):
    """L int8 (2, n): 0 IG, 1 TSS, 2 body, 3 TTS, -1 ignore -> float16 (2, 4, n) soft class targets."""
    out = np.zeros((2, 4, L.shape[1]), np.float32); x = np.arange(L.shape[1])
    for s in range(2):
        l = L[s]; direction = 1 if s == 0 else -1
        base = np.where(l == 2, 2, 0)
        g = {}
        for ci in (1, 3):
            m = (l == ci).astype(np.int8); d = np.diff(np.r_[0, m, 0]); st, en = np.where(d == 1)[0], np.where(d == -1)[0]
            cen = ((st + en - 1) / 2).astype(int)
            if len(cen) == 0: g[ci] = np.zeros(len(l)); continue
            j = np.clip(np.searchsorted(cen, x), 1, len(cen) - 1) if len(cen) > 1 else np.zeros(len(l), int)
            dist = np.minimum(np.abs(x - cen[j - 1]), np.abs(x - cen[j])) if len(cen) > 1 else np.abs(x - cen[0])
            near = cen[j - 1] if len(cen) > 1 else np.full(len(l), cen[0])
            near = np.where(np.abs(x - cen[j]) < np.abs(x - near), cen[j], near) if len(cen) > 1 else near
            g[ci] = np.exp(-0.5 * (dist / sigma_bins) ** 2)
            inrun = l == ci; after = (x - near) * direction > 0
            base = np.where(inrun, np.where(after, 2 if ci == 1 else 0, 0 if ci == 1 else 2), base)
        gt, ge = g[1], g[3]; tot = gt + ge; sc = np.where(tot > 1, 1 / np.maximum(tot, 1e-9), 1); gt, ge = gt * sc, ge * sc
        out[s, 1] = gt; out[s, 3] = ge; rest = 1 - gt - ge
        out[s, 0] += np.where(base == 0, rest, 0); out[s, 2] += np.where(base == 2, rest, 0)
    return out.astype(np.float16)

# ---------------- U-Net ----------------
class Down(nn.Module):
    def __init__(s, ci, co): super().__init__(); s.c = nn.Conv1d(ci, co, 3, padding=1); s.r = ResBlock(co, 1)
    def forward(s, h): return s.r(F.gelu(s.c(F.avg_pool1d(h, 2))))
class Up(nn.Module):
    def __init__(s, ci, cs, co): super().__init__(); s.c = nn.Conv1d(ci + cs, co, 3, padding=1); s.r = ResBlock(co, 1)
    def forward(s, h, skip): return s.r(F.gelu(s.c(torch.cat([F.interpolate(h, size=skip.shape[-1], mode='nearest'), skip], 1))))
class UNet(nn.Module):
    def __init__(s, n_in, chs=(64, 80, 96, 128), n_mid=6, p=0.1):
        super().__init__()
        s.stem = nn.Sequential(nn.Conv1d(n_in, chs[0], 11, padding=5), nn.GELU()); s.r0 = nn.Sequential(ResBlock(chs[0], 1, p=p), ResBlock(chs[0], 2, p=p))
        s.down = nn.ModuleList(Down(chs[i], chs[i + 1]) for i in range(len(chs) - 1))
        s.mid = nn.ModuleList(ResBlock(chs[-1], 2 ** i, p=p) for i in range(n_mid))
        s.up = nn.ModuleList(Up(chs[i + 1], chs[i], chs[i]) for i in reversed(range(len(chs) - 1)))
        s.head = nn.Conv1d(chs[0], 8, 1)
    def forward(s, x):
        h = s.r0(s.stem(x)); skips = [h]
        for d in s.down: h = d(h); skips.append(h)
        for m in s.mid: h = m(h)
        for u, sk in zip(s.up, reversed(skips[:-1])): h = u(h, sk)
        return s.head(h)

# ---------------- multi-resolution U-Net (Basenji/Enformer-style) ----------------
class AlibiAttn(nn.Module):
    """Pre-norm transformer block on the coarsest level; ALiBi distance bias (translation-invariant, so sliding windows are consistent)."""
    def __init__(s, d, heads=4, p=0.1):
        super().__init__(); s.h = heads; s.n1 = nn.LayerNorm(d); s.att = nn.MultiheadAttention(d, heads, dropout=p, batch_first=True)
        s.n2 = nn.LayerNorm(d); s.ff = nn.Sequential(nn.Linear(d, 2 * d), nn.GELU(), nn.Dropout(p), nn.Linear(2 * d, d)); s.dr = nn.Dropout(p)
        s.register_buffer('slopes', torch.tensor([2 ** (-(i + 1) * 8 / heads) for i in range(heads)]), persistent=False)
    def forward(s, h):  # h (B, C, L)
        x = h.transpose(1, 2); B, L, _ = x.shape
        pos = torch.arange(L, device=x.device); dist = (pos[None] - pos[:, None]).abs().float()
        bias = (-s.slopes[:, None, None] * dist[None]).repeat(B, 1, 1)
        y = s.n1(x); x = x + s.dr(s.att(y, y, y, attn_mask=bias, need_weights=False)[0]); x = x + s.dr(s.ff(s.n2(x)))
        return x.transpose(1, 2)
class Level(nn.Module):
    def __init__(s, ci, co, f, p):
        super().__init__(); s.f = f; s.c = nn.Conv1d(ci, co, 5, padding=2); s.r = nn.Sequential(ResBlock(co, 1, p=p), ResBlock(co, 2, p=p))
    def forward(s, h): return s.r(F.gelu(s.c(F.avg_pool1d(h, s.f, ceil_mode=True))))
class UpL(nn.Module):
    def __init__(s, ci, cs, co, p):
        super().__init__(); s.c = nn.Conv1d(ci + cs, co, 5, padding=2); s.r = nn.Sequential(ResBlock(co, 1, p=p), ResBlock(co, 2, p=p))
    def forward(s, h, skip):
        return s.r(F.gelu(s.c(torch.cat([F.interpolate(h, size=skip.shape[-1], mode='linear', align_corners=False), skip], 1))))
class UNet2(nn.Module):
    """20 bp in/out. pools e.g. (5, 2): 20 -> 100 -> 200 bp. Dilated tower (+ optional attention) at the coarsest level."""
    def __init__(s, n_in, pools=(5, 2), chs=(64, 96, 128), n_mid=6, n_attn=0, p=0.1):
        super().__init__(); assert len(chs) == len(pools) + 1
        s.stem = nn.Sequential(nn.Conv1d(n_in, chs[0], 11, padding=5), nn.GELU()); s.r0 = nn.Sequential(*[ResBlock(chs[0], d, p=p) for d in (1, 2, 4)])
        s.down = nn.ModuleList(Level(chs[i], chs[i + 1], f, p) for i, f in enumerate(pools))
        s.mid = nn.ModuleList(ResBlock(chs[-1], 2 ** i, p=p) for i in range(n_mid)); s.attn = nn.ModuleList(AlibiAttn(chs[-1], 4, p) for _ in range(n_attn))
        s.up = nn.ModuleList(UpL(chs[i + 1], chs[i], chs[i], p) for i in reversed(range(len(pools))))
        s.head = nn.Sequential(nn.Conv1d(chs[0], chs[0], 1), nn.GELU(), nn.Conv1d(chs[0], 8, 1))
    def forward(s, x):
        h = s.r0(s.stem(x)); sk = [h]
        for d in s.down: h = d(h); sk.append(h)
        for m in s.mid: h = m(h)
        for a in s.attn: h = a(h)
        for u, k in zip(s.up, reversed(sk[:-1])): h = u(h, k)
        return s.head(h)

# ---------------- training ----------------
class Win(torch.utils.data.Dataset):
    def __init__(s, X, Lab, Q, keys, win, step):
        s.X, s.Lab, s.Q, s.win = X, Lab, Q, win; s.idx = [(k, t) for k in keys for t in range(0, X[k].shape[1] - win + 1, step)]
    def __len__(s): return len(s.idx)
    def __getitem__(s, i):
        k, t = s.idx[i]; sl = slice(t, t + s.win); lab = s.Lab[k][:, sl].astype(np.int64)
        q = torch.from_numpy(s.Q[k][:, :, sl].astype(np.float32)) if s.Q is not None else torch.zeros(1)
        return torch.from_numpy(s.X[k][:, sl].astype(np.float32)), torch.from_numpy(lab), torch.from_numpy((lab >= 0).astype(np.float32)), q
NSTR = 0  # number of stranded evidence channels per strand, stored as the last 2*NSTR input rows (+ block then - block)
def flip_all(x, y, m, q, nu, kg, oh):
    x, y, m = SM.flip(x, y, m, nu, kg, oh)
    if NSTR:
        a, b = x[:, -2 * NSTR:-NSTR].clone(), x[:, -NSTR:].clone(); x[:, -2 * NSTR:-NSTR], x[:, -NSTR:] = b, a
    if q.dim() > 2: q = q.flip(-1)[:, [1, 0]]
    return x, y, m, q
def loss_fn(logits, y, m, q, w, soft):
    B, _, L = logits.shape; lg = logits.view(B, 2, 4, L)
    if not soft: return SM.seg_loss(logits, y, m, w)
    lp = F.log_softmax(lg, 2); wq = q * w.view(1, 1, 4, 1)
    l = -(wq * lp).sum(2); return (l * m).sum() / ((wq.sum(2)) * m).sum().clamp_min(1.0)
def train(net, X, Lab, Q, tr, va, nu, kg, oh, epochs, spe, win, crop, lr, batch, seed, patience, log, soft, wd=1e-4, chdrop=0.0, noise=0.0):
    torch.manual_seed(seed); np.random.seed(seed)
    allab = np.concatenate([Lab[k].ravel() for k in tr]); cnt = np.array([(allab == c).sum() for c in range(4)], float)
    w = torch.tensor((cnt.sum() / cnt) ** 0.5, dtype=torch.float32); w = w / w.mean()
    opt = torch.optim.AdamW(net.parameters(), lr=lr, weight_decay=wd)
    steps = epochs * int(np.ceil(spe / batch)); sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=lr, total_steps=steps, pct_start=0.1)
    trd = Win(X, Lab, Q, tr, win, win // 4)
    vad = Win(X, Lab, Q, va, win, win // 2); vsub = torch.utils.data.Subset(vad, np.random.default_rng(1).choice(len(vad), min(800, len(vad)), replace=False))
    vl_ = torch.utils.data.DataLoader(vsub, batch); cs = slice(crop, win - crop); best, state, bad, hist = np.inf, None, 0, []
    for ep in range(epochs):
        net.train(); tl = []
        dl = torch.utils.data.DataLoader(trd, batch, sampler=torch.utils.data.RandomSampler(trd, num_samples=min(spe, len(trd))))
        for x, y, m, q in dl:
            if np.random.rand() < 0.5: x, y, m, q = flip_all(x, y, m, q, nu, kg, oh)
            if chdrop > 0:
                keep = (torch.rand(x.shape[0], x.shape[1], 1) > chdrop).float()
                if NSTR: keep[:, -NSTR:] = keep[:, -2 * NSTR:-NSTR]   # drop +/- evidence pairs together
                x = x * keep
            if noise > 0: x = x + noise * torch.randn_like(x)
            loss = loss_fn(net(x)[..., cs], y[..., cs], m[..., cs], q[..., cs] if soft else q, w, soft)
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 1.0); opt.step(); sched.step(); tl.append(loss.item())
        net.eval()
        with torch.no_grad(): v = float(np.mean([SM.seg_loss(net(x)[..., cs], y[..., cs], m[..., cs], w).item() for x, y, m, q in vl_]))
        hist.append((ep, float(np.mean(tl)), v)); log(f"ep{ep:02d} train {np.mean(tl):.4f} val_hardCE {v:.4f}")
        if v < best - 1e-4: best, bad, state = v, 0, {k: t.detach().clone() for k, t in net.state_dict().items()}
        else:
            bad += 1
            if bad >= patience: break
    net.load_state_dict(state); net.eval(); return net, hist, w

def make_net(a, n_in):
    g = (lambda k: a[k]) if isinstance(a, dict) else (lambda k: getattr(a, k))
    arch = g('arch')
    if arch == 'dilated': return SM.SegNet(n_in, n_layers=8, p=g('dropout'))
    if arch == 'unet': return UNet(n_in, p=g('dropout'))
    return UNet2(n_in, tuple(int(x) for x in g('pools').split(',')), tuple(int(x) for x in g('chs').split(',')), g('n_mid'), g('n_attn'), g('dropout'))

@torch.no_grad()
def predict(net, X, keys, nu, kg, oh, win, crop):
    out = {}
    for k in keys:
        x = X[k]; nb = x.shape[1]; acc = np.zeros((8, nb), np.float32); cnt = np.zeros(nb, np.float32)
        step = win - 2 * crop; starts = list(range(0, max(1, nb - win + 1), step))
        if starts[-1] + win < nb: starts.append(max(0, nb - win))
        for s in starts:
            xb = torch.from_numpy(x[:, s:s + win].astype(np.float32))[None]; p = torch.softmax(net(xb).view(1, 2, 4, -1), 2)
            d = torch.zeros(1, 2, xb.shape[-1], dtype=torch.long); xf, _, _, _ = flip_all(xb, d, d.float(), torch.zeros(1, 1), nu, kg, oh)
            p = 0.5 * (p + torch.softmax(net(xf).view(1, 2, 4, -1), 2).flip(-1)[:, [1, 0]])
            p = p.reshape(8, -1).numpy(); L = p.shape[1]; lo = 0 if s == 0 else crop; hi = L if s + win >= nb else L - crop
            acc[:, s + lo:s + hi] += p[:, lo:hi]; cnt[s + lo:s + hi] += 1
        out[k] = (acc / np.maximum(cnt, 1)).reshape(2, 4, nb).astype(np.float16)
    return out

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    for k, d in dict(name=None, spec='chrom,onehot,glm,rna', tracks='', res=20, arch='dilated', soft=0, sigma_bp=60, spe=900, epochs=20, seed=0,
                     threads=4, patience=100, lr=2e-3, dropout=0.1, wd=1e-4, chdrop=0.0, noise=0.0, pools='5,2', chs='64,96,128', n_mid=6, n_attn=0, crop_bp=3200, test='1,2,3,4,5', val='6,7', out='exp', win_bp=40960).items():
        ap.add_argument('--' + k, type=type(d) if d is not None else str, default=d)
    a = ap.parse_args(); torch.set_num_threads(a.threads); os.makedirs(a.out, exist_ok=True)
    if os.path.exists(f"{a.out}/pred_{a.name}.npz"): print('exists, skipping', a.name); sys.exit(0)
    T0 = time.time(); log = lambda *s: print(f"[{time.time()-T0:6.0f}s]", *s, flush=True)
    TRK = [x for x in a.tracks.split(',') if x] or None
    spec = set(a.spec.split(',')); X20, nu, kg, oh = SF.assemble(spec - {'strand'}, TRK)
    if 'strand' in spec:
        ST = stranded(); NSTR = next(iter(ST.values())).shape[0] // 2; globals()['NSTR'] = NSTR
        for k in X20: n = min(X20[k].shape[1], ST[k].shape[1]); X20[k] = np.concatenate([X20[k][:, :n], ST[k][:, :n]])
    L20 = {f"{s}|{c}": v for s in P.STAGES for c, v in np.load(f"data/{os.environ.get('LABELS', 'L4')}_{s}.npz").items()}
    for k in X20: n = min(X20[k].shape[1], L20[k].shape[1]); X20[k] = X20[k][:, :n]; L20[k] = L20[k][:, :n]
    X = {k: to_res(v, a.res) for k, v in X20.items()}; Lab = {k: lab_to_res(v, a.res) for k, v in L20.items()}
    Q = {k: soft_targets(v, a.sigma_bp / a.res) for k, v in Lab.items()} if a.soft else None
    te = [int(i) for i in a.test.split(',')]; va = [int(i) for i in a.val.split(',')]; num = lambda k: int(k.split('|')[1][6:8])
    Kf = lambda f: [k for k in X if f(num(k))]
    tr_k, va_k, te_k = Kf(lambda n: n not in te and n not in va), Kf(lambda n: n in va), Kf(lambda n: n in te)
    win = a.win_bp // a.res; win -= win % 8; crop = max(8, (a.crop_bp // a.res) // 8 * 8)
    n_in = X[tr_k[0]].shape[0]
    net = make_net(a, n_in)
    log(a.name, 'res', a.res, 'arch', a.arch, 'soft', a.soft, 'in_ch', n_in, 'win', win, 'params', sum(p.numel() for p in net.parameters()))
    net, hist, w = train(net, X, Lab, Q, tr_k, va_k, nu, kg, oh, a.epochs, a.spe, win, crop, a.lr, 24, a.seed, a.patience, log, a.soft, a.wd, a.chdrop, a.noise)
    torch.save(dict(state=net.state_dict(), args=vars(a), n_in=n_in, nu=nu, kg=kg, oh=oh, w=w), f"{a.out}/{a.name}.pt")
    Pd = predict(net, X, va_k + te_k, nu, kg, oh, win, crop)
    Pd = {k: pred_to20(v, a.res, L20[k].shape[1]) for k, v in Pd.items()}
    np.savez_compressed(f"{a.out}/pred_{a.name}.npz", **{k.replace('|', '__'): v for k, v in Pd.items()})
    json.dump(hist, open(f"{a.out}/{a.name}_hist.json", "w")); log('done', a.name)
