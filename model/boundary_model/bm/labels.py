"""Stranded soft label tracks from a boundary table.

Output channels (fixed order):  0 TSS+, 1 TTS+, 2 TSS-, 3 TTS-
Soft label = Gaussian bump, sigma from the boundary's read-end IQR (clipped),
max-combined where bumps overlap. Ignore mask = 1 where loss is computed.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from .tracks import RES

CHANNELS = ["TSS+", "TTS+", "TSS-", "TTS-"]
FLIP = [2, 3, 0, 1]  # channel permutation under strand flip


def channel_of(end_type: str, strand: str) -> int:
    return {("TSS", "+"): 0, ("TTS", "+"): 1, ("TSS", "-"): 2, ("TTS", "-"): 3}[(end_type, strand)]


def sigma_bins(iqr_bp, res=RES, lo=20.0, hi=150.0):
    iqr = np.nan_to_num(np.asarray(iqr_bp, float), nan=50.0)
    return np.clip(iqr / 1.349, lo, hi) / res  # IQR -> sd, clip to [20,150] bp


def build_labels(ref: pd.DataFrame, chrom_sizes: dict[str, int], res: int = RES,
                 tiers=("gold", "silver"), exclusion: pd.DataFrame | None = None,
                 ignore_radius_bp: int = 300) -> tuple[dict, dict]:
    """Return (Y, M): Y[c] (4, nbins) soft labels; M[c] (4, nbins) loss mask.

    Masked out (M=0):
      * boundaries of lower tiers not used for training (bronze) +/- ignore_radius
      * exclusion regions (e.g. CAGE clusters / catalogued lncRNA ends not in the
        reference, unfiltered long-read ends) -> treated as unlabeled, not negative
    """
    Y = {c: np.zeros((4, L // res), np.float32) for c, L in chrom_sizes.items()}
    M = {c: np.ones((4, L // res), np.float32) for c, L in chrom_sizes.items()}
    r = ref.copy()
    r["sig"] = sigma_bins(r["iqr_bp"] if "iqr_bp" in r else np.full(len(r), np.nan), res)
    rad = ignore_radius_bp // res
    for row in r.itertuples(index=False):
        c = row.chrom
        if c not in Y: continue
        ch = channel_of(row.end_type, row.strand)
        b = int(row.pos) // res
        nb = Y[c].shape[1]
        if row.tier not in tiers:
            M[c][ch, max(0, b - rad):min(nb, b + rad + 1)] = 0
            continue
        w = int(np.ceil(4 * row.sig))
        lo, hi = max(0, b - w), min(nb, b + w + 1)
        x = np.arange(lo, hi) - b
        Y[c][ch, lo:hi] = np.maximum(Y[c][ch, lo:hi], np.exp(-0.5 * (x / row.sig) ** 2))
    if exclusion is not None:
        for row in exclusion.itertuples(index=False):
            c = row.chrom
            if c not in M: continue
            lo, hi = int(row.start) // res, int(row.end) // res + 1
            chans = range(4) if getattr(row, "strand", ".") not in "+-" else \
                ([0, 1] if row.strand == "+" else [2, 3])
            for ch in chans:
                keep = Y[c][ch, lo:hi] > 0.05  # never mask a labelled positive
                M[c][ch, lo:hi] = np.where(keep, 1.0, 0.0)
    return Y, M
