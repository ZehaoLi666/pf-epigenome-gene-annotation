"""Window-feature baselines comparable to the manuscript's Random Forest.

Features: mean of each track in 20 x 100-bp sub-bins spanning +/-1 kb around a
position, ORIENTED to the strand being scored (reverse the window for '-').
One model per end type (TSS, TTS), applied to both orientations.
"""
from __future__ import annotations
import numpy as np
from .tracks import RES

SUB = 100 // RES          # bins per sub-bin
HALF = 1000 // RES        # half window in bins
OFFS = np.arange(-HALF, HALF, SUB)


def window_features(x: np.ndarray, pos: np.ndarray, orient: int) -> np.ndarray:
    T, nb = x.shape
    cs = np.concatenate([np.zeros((T, 1), x.dtype), np.cumsum(x, 1)], 1)
    s = np.clip(pos[:, None] + OFFS[None, :], 0, nb)
    e = np.clip(s + SUB, 0, nb)
    f = (cs[:, e] - cs[:, s]) / np.maximum(e - s, 1)          # (T, n, 20)
    if orient < 0:
        f = f[:, :, ::-1]
    return f.transpose(1, 0, 2).reshape(len(pos), -1).astype(np.float32)


def sample_training(X, ref, chroms, end_type, rng, neg_ratio=1.0, hard_frac=0.5, tiers=("gold", "silver")):
    Xs, ys = [], []
    r = ref[(ref.end_type == end_type) & ref.tier.isin(tiers) & ref.chrom.isin(chroms)]
    for c, g in r.groupby("chrom"):
        nb = X[c].shape[1]
        bset = (ref[(ref.chrom == c)].pos.values // RES)
        for strand, orient in (("+", 1), ("-", -1)):
            p = (g[g.strand == strand].pos.values // RES).astype(int)
            if len(p) == 0: continue
            n_neg = int(len(p) * neg_ratio)
            n_hard = int(n_neg * hard_frac)
            hard = p[rng.integers(0, len(p), n_hard)] + rng.choice([-1, 1], n_hard) * rng.integers(20, 61, n_hard)
            rand = rng.integers(HALF, nb - HALF, n_neg - n_hard)
            neg = np.concatenate([hard, rand]).clip(HALF, nb - HALF - 1)
            # drop negatives within 150 bp of any boundary
            d = np.abs(neg[:, None] - bset[None, :]).min(1) if len(bset) else np.full(len(neg), 1e9)
            neg = neg[d > 15]
            Xs += [window_features(X[c], p, orient), window_features(X[c], neg, orient)]
            ys += [np.ones(len(p)), np.zeros(len(neg))]
    return np.concatenate(Xs), np.concatenate(ys)


def make_model(kind: str, seed=0):
    if kind == "logistic":
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler
        from sklearn.linear_model import LogisticRegression
        return make_pipeline(StandardScaler(), LogisticRegression(C=0.1, max_iter=2000))
    if kind == "rf":
        from sklearn.ensemble import RandomForestClassifier
        return RandomForestClassifier(n_estimators=500, min_samples_leaf=2, n_jobs=-1, random_state=seed)
    if kind == "lgbm":
        import lightgbm as lgb
        return lgb.LGBMClassifier(n_estimators=400, learning_rate=0.03, num_leaves=31,
                                  subsample=0.8, subsample_freq=1, colsample_bytree=0.8,
                                  random_state=seed, verbose=-1)
    raise ValueError(kind)


def fit_baseline(kind, X, ref, train_chroms, seed=0):
    rng = np.random.default_rng(seed)
    models = {}
    for et in ("TSS", "TTS"):
        F, y = sample_training(X, ref, train_chroms, et, rng)
        models[et] = make_model(kind, seed).fit(F, y)
    return models


def predict_baseline(models, X, chroms, step_bins=5):
    """Score every `step_bins` bins in both orientations -> (4, nbins) tracks."""
    P = {}
    for c in chroms:
        nb = X[c].shape[1]
        pos = np.arange(0, nb, step_bins)
        out = np.zeros((4, nb), np.float32)
        for ch, (et, orient) in enumerate([("TSS", 1), ("TTS", 1), ("TSS", -1), ("TTS", -1)]):
            s = models[et].predict_proba(window_features(X[c], pos, orient))[:, 1]
            out[ch] = np.interp(np.arange(nb), pos, s)
        P[c] = out
    return P
