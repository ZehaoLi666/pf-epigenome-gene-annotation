"""From probability tracks to boundary calls and transcription-unit models.

call_peaks     local maxima above threshold with non-maximum suppression
fit_length_prior   log-length KDE from reference units
decode_units   gene-grammar decoding: choose non-overlapping TSS->TTS units per
               strand maximising sum(logit TSS + logit TTS + log length prior)
               (weighted interval scheduling = exact semi-Markov decoding with
               an intergenic state of zero cost). Several units per intergenic
               interval are allowed, so a lncRNA and a gene can coexist.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from .tracks import RES
from .labels import CHANNELS


def call_peaks(p: np.ndarray, thr=0.5, min_dist_bins=15):
    cand = np.where((p >= thr) & (p >= np.roll(p, 1)) & (p >= np.roll(p, -1)))[0]
    if len(cand) == 0: return np.array([], int), np.array([])
    order = cand[np.argsort(-p[cand])]
    taken = np.zeros(len(p), bool); keep = []
    for i in order:
        if taken[i]: continue
        keep.append(i); taken[max(0, i - min_dist_bins):i + min_dist_bins + 1] = True
    keep = np.sort(np.array(keep))
    return keep, p[keep]


def calls_table(P: dict, thr=0.5, min_dist_bins=15, res=RES) -> pd.DataFrame:
    rows = []
    for c, p in P.items():
        for ch, name in enumerate(CHANNELS):
            idx, sc = call_peaks(p[ch], thr, min_dist_bins)
            et, st = name[:3], name[3]
            rows += [(c, int(i) * res + res // 2, st, et, float(s)) for i, s in zip(idx, sc)]
    return pd.DataFrame(rows, columns=["chrom", "pos", "strand", "end_type", "score"])


def fit_length_prior(lengths_bp):
    from scipy.stats import gaussian_kde
    k = gaussian_kde(np.log10(np.asarray(lengths_bp, float).clip(50)))
    return lambda L: np.log(k(np.log10(np.asarray(L, float).clip(50))) + 1e-12)


def _logit(p): p = np.clip(p, 1e-6, 1 - 1e-6); return np.log(p / (1 - p))


def decode_units(calls: pd.DataFrame, len_logprior, min_len=200, max_len=20000, bias=0.0):
    out = []
    for (c, st), g in calls.groupby(["chrom", "strand"]):
        tss = g[g.end_type == "TSS"][["pos", "score"]].values
        tts = g[g.end_type == "TTS"][["pos", "score"]].values
        if len(tss) == 0 or len(tts) == 0: continue
        segs = []
        tpos = tts[:, 0]
        for a, sa in tss:
            if st == "+": j = np.where((tpos - a >= min_len) & (tpos - a <= max_len))[0]
            else:         j = np.where((a - tpos >= min_len) & (a - tpos <= max_len))[0]
            if len(j) == 0: continue
            L = np.abs(tpos[j] - a)
            s = _logit(sa) + _logit(tts[j, 1]) + len_logprior(L) + bias
            for jj, ss in zip(j, s):
                if ss > 0:
                    lo, hi = (a, tpos[jj]) if st == "+" else (tpos[jj], a)
                    segs.append((lo, hi, ss, a, tpos[jj]))
        if not segs: continue
        segs.sort(key=lambda z: z[1])
        ends = np.array([s[1] for s in segs])
        best = np.zeros(len(segs) + 1); choice = [None] * (len(segs) + 1)
        for i, (lo, hi, sc, a, b) in enumerate(segs, 1):
            k = np.searchsorted(ends, lo, side="left")  # segments ending before lo
            take = best[k] + sc
            if take > best[i - 1]: best[i], choice[i] = take, ("t", k)
            else: best[i], choice[i] = best[i - 1], ("s", i - 1)
        i = len(segs)
        while i > 0:
            kind, k = choice[i]
            if kind == "t":
                lo, hi, sc, a, b = segs[i - 1]
                out.append((c, int(lo), int(hi), st, int(a), int(b), float(sc))); i = k
            else:
                i = k
    return pd.DataFrame(out, columns=["chrom", "start", "end", "strand", "tss", "tts", "score"])
