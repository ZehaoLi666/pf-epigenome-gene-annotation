"""Evaluation contract: every model row reports the same numbers.

  AUPRC (per channel, masked bins, target >= 0.5)
  positional error: median |pred - true| of one-to-one matches within tolerance
  recall@{50,100,200,500} bp, precision@tol
  empirical FDR at the operating threshold from an independent-track-shift null
  optional stratification by any reference column (tier, expr_quantile, situation)
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score
from .labels import CHANNELS, channel_of
from .decode import calls_table
from .tracks import RES


def auprc(P, Y, M, chroms):
    res = {}
    for ch, name in enumerate(CHANNELS):
        p = np.concatenate([P[c][ch][M[c][ch] > 0] for c in chroms])
        y = np.concatenate([Y[c][ch][M[c][ch] > 0] for c in chroms]) >= 0.5
        res[name] = average_precision_score(y, p) if y.any() else np.nan
    return res


def match(pred: np.ndarray, true: np.ndarray, tol: int):
    """Greedy one-to-one nearest matching. Returns (pred_idx, true_idx, signed_err)."""
    if len(pred) == 0 or len(true) == 0: return np.array([], int), np.array([], int), np.array([])
    d = pred[:, None] - true[None, :]
    ii, jj = np.where(np.abs(d) <= tol)
    order = np.argsort(np.abs(d[ii, jj]), kind="stable")
    up, ut, pi, ti = set(), set(), [], []
    for k in order:
        a, b = ii[k], jj[k]
        if a in up or b in ut: continue
        up.add(a); ut.add(b); pi.append(a); ti.append(b)
    pi, ti = np.array(pi, int), np.array(ti, int)
    return pi, ti, (pred[pi] - true[ti]) if len(pi) else np.array([])


def positional_metrics(calls: pd.DataFrame, ref: pd.DataFrame, chroms, tols=(50, 100, 200, 500),
                       strat_col: str | None = None):
    rows = []
    ref = ref[ref.chrom.isin(chroms)]
    calls = calls[calls.chrom.isin(chroms)]
    groups = [("all", ref)] if strat_col is None else list(ref.groupby(strat_col))
    for gname, rg in groups:
        for et in ("TSS", "TTS"):
            n_true = n_pred = 0; errs = []; hits = {t: 0 for t in tols}; tp_pred = 0
            for (c, st), r in rg[rg.end_type == et].groupby(["chrom", "strand"]):
                pr = calls[(calls.chrom == c) & (calls.strand == st) & (calls.end_type == et)].pos.values
                tr = r.pos.values
                n_true += len(tr)
                for t in tols:
                    hits[t] += len(match(pr, tr, t)[0])
                _, _, e = match(pr, tr, max(tols)); errs += list(e)
            if strat_col is None:
                n_pred = int(((calls.end_type == et)).sum())
                tp_pred = hits[max(tols)]
            row = {"stratum": gname, "end_type": et, "n_true": n_true,
                   "median_abs_err_bp": float(np.median(np.abs(errs))) if errs else np.nan}
            row.update({f"recall@{t}": hits[t] / max(n_true, 1) for t in tols})
            if strat_col is None:
                row["n_calls"] = n_pred
                row[f"precision@{max(tols)}"] = tp_pred / max(n_pred, 1)
            rows.append(row)
    return pd.DataFrame(rows)


def shift_null_fdr(predict_fn, X, chroms, thr, n_rep=3, min_shift_bp=100_000, seed=0):
    """Independent circular shift of each track destroys cross-track co-occurrence
    while preserving each track's own autocorrelation. Calls surviving on shifted
    input estimate the false-call rate at `thr`.  FDR ~= null calls / real calls."""
    rng = np.random.default_rng(seed)
    real = calls_table(predict_fn(X, chroms), thr)
    nulls = []
    for _ in range(n_rep):
        Xs = {}
        for c in chroms:
            x = X[c].copy(); nb = x.shape[1]
            lo = min(min_shift_bp // RES, nb // 3)
            for t in range(x.shape[0]):
                x[t] = np.roll(x[t], int(rng.integers(lo, nb - lo)))
            Xs[c] = x
        nulls.append(len(calls_table(predict_fn(Xs, chroms), thr)))
    fdr = float(np.mean(nulls) / max(len(real), 1))
    return {"real_calls": len(real), "null_calls_mean": float(np.mean(nulls)), "fdr_est": min(fdr, 1.0)}
