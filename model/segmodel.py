"""4-class stranded gene segmentation (intergenic / TSS / body / TTS on each strand).
Input channels: [unstranded tracks (chromatin, optional mRNA/GRO)] + [stranded DNA: onehot fwd (4) + gLM fwd (k) + gLM rc (k)].
Output: 8 logits = strand(+,-) x class(0 intergenic,1 TSS,2 body,3 TTS)."""
import numpy as np, torch, torch.nn as nn, torch.nn.functional as F
import sys; sys.path.insert(0, 'boundary_model')
from bm.models import ResBlock

class SegNet(nn.Module):
    def __init__(self, n_in, ch=64, n_layers=8, p=0.1):
        super().__init__()
        self.stem = nn.Sequential(nn.Conv1d(n_in, ch, 11, padding=5), nn.GELU())
        self.blocks = nn.ModuleList(ResBlock(ch, 2 ** i, p=p) for i in range(n_layers))
        self.mix = nn.Conv1d(2 * ch, ch, 1)
        self.head = nn.Conv1d(ch, 8, 1)
        self.receptive_field_bins = 11 + sum(2 * 2 ** i for i in range(n_layers))
    def forward(self, x):
        h0 = self.stem(x); h = h0
        for b in self.blocks: h = b(h)
        return self.head(F.gelu(self.mix(torch.cat([h, h0], 1))))  # (B, 8, L)

def flip(x, y, m, n_unstr, k_glm, has_onehot):
    """Reverse window = view from the opposite strand. Unstranded channels just reverse; one-hot is
    reverse-complemented (A<->T, C<->G); gLM fwd/rc blocks swap; labels swap strands."""
    x = x.flip(-1).clone()
    o = n_unstr
    if has_onehot:
        x[:, o:o + 4] = x[:, o:o + 4][:, [3, 2, 1, 0]]; o += 4
    if k_glm:
        a, b = x[:, o:o + k_glm].clone(), x[:, o + k_glm:o + 2 * k_glm].clone()
        x[:, o:o + k_glm], x[:, o + k_glm:o + 2 * k_glm] = b, a
    y = y.flip(-1)[:, [1, 0]]; m = m.flip(-1)[:, [1, 0]]
    return x, y, m

def seg_loss(logits, y, m, w):
    B, _, L = logits.shape
    lg = logits.view(B, 2, 4, L).permute(0, 1, 3, 2).reshape(-1, 4)
    yy = y.reshape(-1).clamp_min(0); mm = m.reshape(-1)
    l = F.cross_entropy(lg, yy, weight=w, reduction='none')
    return (l * mm).sum() / (w[yy] * mm).sum().clamp_min(1.0)

class Win(torch.utils.data.Dataset):
    def __init__(self, X, Lab, keys, win, step):
        self.X, self.Lab, self.win = X, Lab, win
        self.idx = [(k, s) for k in keys for s in range(0, X[k].shape[1] - win + 1, step)]
    def __len__(self): return len(self.idx)
    def __getitem__(self, i):
        k, s = self.idx[i]; sl = slice(s, s + self.win)
        lab = self.Lab[k][:, sl].astype(np.int64)
        return torch.from_numpy(self.X[k][:, sl].astype(np.float32)), torch.from_numpy(lab), torch.from_numpy((lab >= 0).astype(np.float32))

def train(X, Lab, tr, va, n_unstr, k_glm, has_onehot, epochs=20, win=2048, crop=160, spe=900, lr=2e-3, batch=24, seed=0, patience=4, log=print):
    torch.manual_seed(seed); np.random.seed(seed)
    n_in = next(iter(X.values())).shape[0]
    net = SegNet(n_in, n_layers=8)
    allab = np.concatenate([Lab[k].ravel() for k in tr]); cnt = np.array([(allab == c).sum() for c in range(4)], float)
    w = torch.tensor((cnt.sum() / cnt) ** 0.5, dtype=torch.float32); w = w / w.mean()
    opt = torch.optim.AdamW(net.parameters(), lr=lr, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)
    trd = Win(X, Lab, tr, win, win // 4)
    samp = torch.utils.data.RandomSampler(trd, replacement=False, num_samples=min(spe, len(trd)))
    dl = torch.utils.data.DataLoader(trd, batch, sampler=samp)
    vad = Win(X, Lab, va, win, win); vsub = torch.utils.data.Subset(vad, np.random.default_rng(1).choice(len(vad), min(300, len(vad)), replace=False))
    vl_ = torch.utils.data.DataLoader(vsub, batch)
    cs = slice(crop, win - crop); best, state, bad, hist = np.inf, None, 0, []
    for ep in range(epochs):
        net.train(); tl = []
        for x, y, m in dl:
            if np.random.rand() < 0.5: x, y, m = flip(x, y, m, n_unstr, k_glm, has_onehot)
            loss = seg_loss(net(x)[..., cs], y[..., cs], m[..., cs], w)
            opt.zero_grad(); loss.backward(); opt.step(); tl.append(loss.item())
        sched.step(); net.eval()
        with torch.no_grad(): v = float(np.mean([seg_loss(net(x)[..., cs], y[..., cs], m[..., cs], w).item() for x, y, m in vl_]))
        hist.append((ep, float(np.mean(tl)), v)); log(f"ep{ep:02d} train {np.mean(tl):.4f} val {v:.4f}")
        if v < best - 1e-4: best, bad, state = v, 0, {k: t.detach().clone() for k, t in net.state_dict().items()}
        else:
            bad += 1
            if bad >= patience: break
    net.load_state_dict(state); net.eval(); return net, hist, w

@torch.no_grad()
def predict(net, X, keys, n_unstr, k_glm, has_onehot, win=2048, crop=160, flip_avg=True):
    """Returns dict key -> float16 (2, 4, n_bins) softmax probabilities."""
    out = {}
    for k in keys:
        x = X[k]; nb = x.shape[1]; acc = np.zeros((8, nb), np.float32); cnt = np.zeros(nb, np.float32)
        step = win - 2 * crop; starts = list(range(0, max(1, nb - win + 1), step))
        if starts[-1] + win < nb: starts.append(max(0, nb - win))
        for s in starts:
            xb = torch.from_numpy(x[:, s:s + win].astype(np.float32))[None]
            p = torch.softmax(net(xb).view(1, 2, 4, -1), 2)
            if flip_avg:
                dummy = torch.zeros(1, 2, xb.shape[-1], dtype=torch.long)
                xf, _, _ = flip(xb, dummy, dummy.float(), n_unstr, k_glm, has_onehot)
                pf = torch.softmax(net(xf).view(1, 2, 4, -1), 2).flip(-1)[:, [1, 0]]
                p = 0.5 * (p + pf)
            p = p.reshape(8, -1).numpy(); L = p.shape[1]
            lo = 0 if s == 0 else crop; hi = L if s + win >= nb else L - crop
            acc[:, s + lo:s + hi] += p[:, lo:hi]; cnt[s + lo:s + hi] += 1
        out[k] = (acc / np.maximum(cnt, 1)).reshape(2, 4, nb).astype(np.float16)
    return out
