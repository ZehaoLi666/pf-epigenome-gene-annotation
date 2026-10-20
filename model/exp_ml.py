"""Multi-label boundary model (spec v3): per strand three independent sigmoid channels [TSS, TTS, body]; flank = all 0.
Reuses the inputs, architectures and strand-flip augmentation of exp.py; labels from data/Lml_pool.npz (labels_ml.py).
Run with POOL=1. Weighted BCE: weight = label weight x (sqrt(neg/pos) for positives), per channel; masked bins (weight 0) ignored."""
import sys, os, time, json, argparse, numpy as np, torch, torch.nn as nn, torch.nn.functional as F
sys.path[:0] = ['.', 'boundary_model']
assert os.environ.get('POOL'), 'run with POOL=1'
import pfprep as P, seg_features as SF, segmodel as SM, exp as X_
def load_labels(res):
    z = np.load(os.environ.get('ML_LABELS', 'data/Lml_pool.npz')); Y, W = {}, {}
    for c in P.CHROMS:
        Y[f'pool|{c}'] = X_.lab_to_res(z[f'{c}__Y'], res); W[f'pool|{c}'] = X_.lab_to_res(z[f'{c}__W'].astype(np.float32), res)
    return Y, W
class Win(torch.utils.data.Dataset):
    def __init__(s, X, Y, W, keys, win, step):
        s.X, s.Y, s.W, s.win = X, Y, W, win
        s.idx = [(k, t) for k in keys for t in range(0, X[k].shape[1] - win + 1, step) if W[k][..., t:t + win].sum() > 0]
    def __len__(s): return len(s.idx)
    def __getitem__(s, i):
        k, t = s.idx[i]; sl = slice(t, t + s.win)
        return (torch.from_numpy(s.X[k][:, sl].astype(np.float32)), torch.from_numpy(s.Y[k][..., sl].astype(np.float32)),
                torch.from_numpy(s.W[k][..., sl].astype(np.float32)))
def flip3(x, y, w, nu, kg, oh):
    x, y, w = SM.flip(x, y, w, nu, kg, oh)      # flips last axis, swaps strand axis 1 (works for (B,2,3,L))
    if X_.NSTR:
        n = X_.NSTR; a, b = x[:, -2 * n:-n].clone(), x[:, -n:].clone(); x[:, -2 * n:-n], x[:, -n:] = b, a
    return x, y, w
def logits3(net, x): o = net(x); B, _, L = o.shape; return o.view(B, 2, 4, L)[:, :, :3]
def bce(lg, y, w, pw):
    wt = w * (1 + (pw.view(1, 1, 3, 1) - 1) * y)   # = pw for y=1, 1 for y=0; interpolates for soft targets
    l = F.binary_cross_entropy_with_logits(lg, y, reduction='none'); return (l * wt).sum() / wt.sum().clamp_min(1.0)
def train(net, X, Y, W, tr, va, nu, kg, oh, a, win, crop, log):
    torch.manual_seed(a.seed); np.random.seed(a.seed)
    pos = np.array([sum(float((W[k][:, c] * Y[k][:, c].astype(np.float32)).sum()) for k in tr) for c in range(3)])
    neg = np.array([sum(float((W[k][:, c] * (1 - Y[k][:, c].astype(np.float32))).sum()) for k in tr) for c in range(3)])
    pw = torch.tensor(np.sqrt(neg / np.maximum(pos, 1)), dtype=torch.float32); log('pos_weight', pw.numpy().round(2).tolist())
    opt = torch.optim.AdamW(net.parameters(), lr=a.lr, weight_decay=a.wd)
    steps = a.epochs * int(np.ceil(a.spe / 24)); sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=a.lr, total_steps=steps, pct_start=0.1)
    trd = Win(X, Y, W, tr, win, win // 4); vad = Win(X, Y, W, va, win, win // 2)
    vsub = torch.utils.data.Subset(vad, np.random.default_rng(1).choice(len(vad), min(800, len(vad)), replace=False))
    vl = torch.utils.data.DataLoader(vsub, 24); cs = slice(crop, win - crop); best, state, bad, hist = np.inf, None, 0, []
    for ep in range(a.epochs):
        net.train(); tl = []
        dl = torch.utils.data.DataLoader(trd, 24, sampler=torch.utils.data.RandomSampler(trd, num_samples=min(a.spe, len(trd))))
        for x, y, w in dl:
            if np.random.rand() < 0.5: x, y, w = flip3(x, y, w, nu, kg, oh)
            if a.chdrop > 0: x = x * (torch.rand(x.shape[0], x.shape[1], 1) > a.chdrop).float()
            loss = bce(logits3(net, x)[..., cs], y[..., cs], w[..., cs], pw)
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 1.0); opt.step(); sched.step(); tl.append(loss.item())
        net.eval()
        with torch.no_grad(): v = float(np.mean([bce(logits3(net, x)[..., cs], y[..., cs], w[..., cs], pw).item() for x, y, w in vl]))
        hist.append((ep, float(np.mean(tl)), v)); log(f'ep{ep:02d} train {np.mean(tl):.4f} val {v:.4f}')
        if v < best - 1e-4: best, bad, state = v, 0, {k: t.detach().clone() for k, t in net.state_dict().items()}
        else:
            bad += 1
            if bad >= a.patience: break
    net.load_state_dict(state); net.eval(); return net, hist, pw
@torch.no_grad()
def predict(net, X, keys, nu, kg, oh, win, crop):
    out = {}
    for k in keys:
        x = X[k]; nb = x.shape[1]; acc = np.zeros((6, nb), np.float32); cnt = np.zeros(nb, np.float32)
        step = win - 2 * crop; starts = list(range(0, max(1, nb - win + 1), step))
        if starts[-1] + win < nb: starts.append(max(0, nb - win))
        for s in starts:
            xb = torch.from_numpy(x[:, s:s + win].astype(np.float32))[None]; p = torch.sigmoid(logits3(net, xb))
            d = torch.zeros(1, 2, 3, xb.shape[-1]); xf, _, _ = flip3(xb, d, d, nu, kg, oh)
            p = 0.5 * (p + torch.sigmoid(logits3(net, xf)).flip(-1)[:, [1, 0]])
            p = p.reshape(6, -1).numpy(); L = p.shape[1]; lo = 0 if s == 0 else crop; hi = L if s + win >= nb else L - crop
            acc[:, s + lo:s + hi] += p[:, lo:hi]; cnt[s + lo:s + hi] += 1
        out[k] = (acc / np.maximum(cnt, 1)).reshape(2, 3, nb).astype(np.float16)
    return out
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    for k, d in dict(name=None, spec='chrom,onehot,glm,rna,strand', tracks='', res=50, arch='dilated', spe=900, epochs=20, seed=0, threads=4, patience=100,
                     lr=2e-3, dropout=0.1, wd=1e-4, chdrop=0.0, pools='5,2', chs='64,96,128', n_mid=6, n_attn=0, crop_bp=3200, test='1,2,3,4,5', val='6,7',
                     out='exp_ml', win_bp=40960).items():
        ap.add_argument('--' + k, type=type(d) if d is not None else str, default=d)
    a = ap.parse_args(); torch.set_num_threads(a.threads); os.makedirs(a.out, exist_ok=True)
    if os.path.exists(f'{a.out}/pred_{a.name}.npz'): print('exists, skipping', a.name); sys.exit(0)
    T0 = time.time(); log = lambda *s: print(f'[{time.time()-T0:6.0f}s]', *s, flush=True)
    TRK = [x for x in a.tracks.split(',') if x] or None
    spec = set(a.spec.split(',')); X20, nu, kg, oh = SF.assemble(spec - {'strand', 'daft'}, TRK)
    if spec & {'strand', 'daft'}:
        ST = X_.stranded(spec); X_.NSTR = next(iter(ST.values())).shape[0] // 2
        for k in X20: n = min(X20[k].shape[1], ST[k].shape[1]); X20[k] = np.concatenate([X20[k][:, :n], ST[k][:, :n]])
    Y, W = load_labels(20)
    n20 = {k: min(X20[k].shape[1], Y[k].shape[-1]) for k in X20}
    X = {k: X_.to_res(X20[k][:, :n20[k]], a.res) for k in X20}
    Yr = {k: X_.lab_to_res(Y[k][..., :n20[k]], a.res) for k in X20}; Wr = {k: X_.lab_to_res(W[k][..., :n20[k]], a.res) for k in X20}
    te = [int(i) for i in a.test.split(',')]; va = [int(i) for i in a.val.split(',')]; num = lambda k: int(k.split('|')[1][6:8])
    Kf = lambda f: [k for k in X if f(num(k))]
    tr_k, va_k, te_k = Kf(lambda n: n not in te and n not in va), Kf(lambda n: n in va), Kf(lambda n: n in te)
    win = a.win_bp // a.res; win -= win % 8; crop = max(8, (a.crop_bp // a.res) // 8 * 8); n_in = X[tr_k[0]].shape[0]
    net = X_.make_net(a, n_in); log(a.name, 'res', a.res, 'in_ch', n_in, 'win', win, 'params', sum(p.numel() for p in net.parameters()))
    net, hist, pw = train(net, X, Yr, Wr, tr_k, va_k, nu, kg, oh, a, win, crop, log)
    torch.save(dict(state=net.state_dict(), args=vars(a), n_in=n_in, nu=nu, kg=kg, oh=oh, pw=pw, nstr=X_.NSTR), f'{a.out}/{a.name}.pt')
    Pd = predict(net, X, va_k + te_k, nu, kg, oh, win, crop)
    Pd = {k: X_.pred_to20(v, a.res, n20[k]) for k, v in Pd.items()}
    np.savez_compressed(f'{a.out}/pred_{a.name}.npz', **{k.replace('|', '__'): v for k, v in Pd.items()})
    json.dump(hist, open(f'{a.out}/{a.name}_hist.json', 'w')); log('done', a.name)
