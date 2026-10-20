"""Dense per-position boundary model.

BoundaryNet: stem conv -> N dilated residual blocks (optionally FiLM-conditioned
on stage) -> 1x1 heads.
  * profile head: 4 logits per bin (TSS+, TTS+, TSS-, TTS-), type-blind
  * attribute head (optional): K per-bin logits for chromatin-separable
    situation classes (e.g. divergent, heterochromatic); trained with masks
Input is unstranded chromatin (T tracks) + optional one-hot DNA (4 ch).
Strand symmetry is enforced by flip augmentation in training and flip-averaging
at inference (see flip_batch / predict_tracks).
"""
from __future__ import annotations
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from .labels import FLIP


class FiLM(nn.Module):
    def __init__(self, n_stages: int, ch: int):
        super().__init__()
        self.emb = nn.Embedding(n_stages, 2 * ch)
        nn.init.zeros_(self.emb.weight)
    def forward(self, h, stage):
        g, b = self.emb(stage).chunk(2, dim=-1)
        return h * (1 + g[..., None]) + b[..., None]


class ResBlock(nn.Module):
    def __init__(self, ch, dil, k=3, p=0.1, n_stages=0):
        super().__init__()
        self.conv = nn.Conv1d(ch, ch, k, padding=dil * (k // 2), dilation=dil)
        self.norm = nn.BatchNorm1d(ch)
        self.film = FiLM(n_stages, ch) if n_stages else None
        self.drop = nn.Dropout(p)
    def forward(self, h, stage=None):
        z = self.norm(self.conv(h))
        if self.film is not None and stage is not None:
            z = self.film(z, stage)
        return h + self.drop(F.gelu(z))


class BoundaryNet(nn.Module):
    def __init__(self, n_in, ch=64, n_layers=8, n_attr=0, n_stages=0, p=0.1):
        super().__init__()
        self.stem = nn.Sequential(nn.Conv1d(n_in, ch, 11, padding=5), nn.GELU())
        self.blocks = nn.ModuleList(ResBlock(ch, 2 ** i, p=p, n_stages=n_stages) for i in range(n_layers))
        self.profile = nn.Conv1d(ch, 4, 1)
        self.attr = nn.Conv1d(ch, n_attr, 1) if n_attr else None
        self.receptive_field_bins = 11 + sum(2 * 2 ** i for i in range(n_layers))
    def forward(self, x, stage=None):
        h = self.stem(x)
        for b in self.blocks:
            h = b(h, stage)
        out = {"logits": self.profile(h)}
        if self.attr is not None:
            out["attr"] = self.attr(h)
        return out


def n_params(m):
    return sum(p.numel() for p in m.parameters())


def flip_batch(x, y=None, m=None):
    """Reverse along position; swap +/- channels of targets. Chromatin is unstranded,
    so the reversed window is a valid example of the opposite strand."""
    x = x.flip(-1)
    if y is not None: y = y.flip(-1)[:, FLIP]
    if m is not None: m = m.flip(-1)[:, FLIP]
    return x, y, m


def masked_loss(logits, y, m, pos_weight=10.0, gamma=0.0):
    """Soft-label BCE (optionally focal), averaged over unmasked bins."""
    pw = torch.as_tensor(pos_weight, device=logits.device)
    l = F.binary_cross_entropy_with_logits(logits, y, reduction="none", pos_weight=pw)
    if gamma > 0:
        p = torch.sigmoid(logits)
        pt = y * p + (1 - y) * (1 - p)
        l = l * (1 - pt).pow(gamma)
    return (l * m).sum() / m.sum().clamp_min(1.0)


class WindowData(torch.utils.data.Dataset):
    def __init__(self, X, Y, M, chroms, win=2048, step=512, stage_id=None):
        self.X, self.Y, self.M, self.win = X, Y, M, win
        self.stage_id = stage_id
        self.idx = [(c, s) for c in chroms for s in range(0, X[c].shape[1] - win + 1, step)]
    def __len__(self): return len(self.idx)
    def __getitem__(self, i):
        c, s = self.idx[i]
        sl = slice(s, s + self.win)
        item = [torch.from_numpy(self.X[c][:, sl]), torch.from_numpy(self.Y[c][:, sl]),
                torch.from_numpy(self.M[c][:, sl])]
        if self.stage_id is not None:
            item.append(torch.tensor(self.stage_id[c], dtype=torch.long))
        return tuple(item)


def train_net(X, Y, M, train_chroms, val_chroms, *, epochs=20, win=2048, crop=256, lr=2e-3,
              batch=32, ch=64, n_layers=8, pos_weight=10.0, gamma=0.0, seed=0,
              patience=4, verbose=True, threads=None, samples_per_epoch=None):
    """Train BoundaryNet with flip augmentation, edge cropping, early stopping on val loss."""
    torch.manual_seed(seed); np.random.seed(seed)
    if threads: torch.set_num_threads(threads)
    n_in = next(iter(X.values())).shape[0]
    net = BoundaryNet(n_in, ch=ch, n_layers=n_layers)
    assert crop >= net.receptive_field_bins // 2 or True
    opt = torch.optim.AdamW(net.parameters(), lr=lr, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)
    trd = WindowData(X, Y, M, train_chroms, win, win // 4)
    samp = torch.utils.data.RandomSampler(trd, replacement=False, num_samples=min(samples_per_epoch, len(trd))) if samples_per_epoch else None
    tr = torch.utils.data.DataLoader(trd, batch, shuffle=samp is None, sampler=samp)
    va = torch.utils.data.DataLoader(WindowData(X, Y, M, val_chroms, win, win), batch)
    best, best_state, bad, hist = np.inf, None, 0, []
    cs = slice(crop, win - crop)
    for ep in range(epochs):
        net.train(); tl = []
        for x, y, m in tr:
            if np.random.rand() < 0.5:
                x, y, m = flip_batch(x, y, m)
            out = net(x)
            loss = masked_loss(out["logits"][..., cs], y[..., cs], m[..., cs], pos_weight, gamma)
            opt.zero_grad(); loss.backward(); opt.step(); tl.append(loss.item())
        sched.step()
        net.eval(); vl = []
        with torch.no_grad():
            for x, y, m in va:
                vl.append(masked_loss(net(x)["logits"][..., cs], y[..., cs], m[..., cs], pos_weight, gamma).item())
        v = float(np.mean(vl)); hist.append((ep, float(np.mean(tl)), v))
        if verbose: print(f"ep{ep:02d} train {np.mean(tl):.4f} val {v:.4f}", flush=True)
        if v < best - 1e-4:
            best, bad = v, 0
            best_state = {k: t.detach().clone() for k, t in net.state_dict().items()}
        else:
            bad += 1
            if bad >= patience: break
    net.load_state_dict(best_state)
    net.eval()
    return net, hist


@torch.no_grad()
def predict_tracks(net, X, chroms, win=2048, crop=256, flip_avg=True):
    """Genome-wide probabilities (4, nbins) per chromosome; overlapping tiles, centre crop."""
    net.eval(); P = {}
    for c in chroms:
        x = X[c]; nb = x.shape[1]; out = np.zeros((4, nb), np.float32); cnt = np.zeros(nb, np.float32)
        step = win - 2 * crop
        starts = list(range(0, max(1, nb - win + 1), step))
        if starts[-1] + win < nb: starts.append(max(0, nb - win))
        for s in starts:
            xb = torch.from_numpy(x[:, s:s + win])[None]
            p = torch.sigmoid(net(xb)["logits"])
            if flip_avg:
                pf = torch.sigmoid(net(xb.flip(-1))["logits"]).flip(-1)[:, FLIP]
                p = 0.5 * (p + pf)
            p = p[0].numpy()
            L = p.shape[1]
            lo = 0 if s == 0 else crop
            hi = L if s + win >= nb else L - crop
            out[:, s + lo:s + hi] += p[:, lo:hi]; cnt[s + lo:s + hi] += 1
        P[c] = out / np.maximum(cnt, 1)
    return P
