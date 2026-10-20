"""Chromatin readouts around TSS / ATG / TTS anchors. Anchors come from long-read truth (confirmed genes) or from the chromatin-FREE model
(other genes), so chromatin is never both the predictor and the readout.
Signals: per-stage normalised chromatin (data/X_chrom_<stage>.npz, 25 tracks, 20 bp), z-scored genome-wide per track; oriented to the gene.
Outputs (prefix): _meta.npz (mean profiles), _tables.tsv."""
import sys, json, numpy as np, pandas as pd
RES = 20; S = ['ring', 'troph', 'schizont']
TR = json.load(open('data/meta.json'))['chrom_tracks']
X = {s: {c: v.astype(np.float32) for c, v in np.load(f'data/X_chrom_{s}.npz').items()} for s in S}
for s in S:
    allv = np.concatenate([X[s][c] for c in X[s]], 1); mu, sd = allv.mean(1, keepdims=True), allv.std(1, keepdims=True) + 1e-6
    for c in X[s]: X[s][c] = (X[s][c] - mu) / sd
def prof(s, rows, col, half=75):
    """rows: DataFrame with chrom, strand and anchor column; returns (n_tracks, 2*half+1) mean oriented profile and n."""
    acc = []; 
    for c, st, p in zip(rows.chrom, rows.strand, rows[col]):
        if pd.isna(p) or c not in X[s]: continue
        b = int(p) // RES; a = X[s][c]
        if b - half < 0 or b + half + 1 > a.shape[1]: continue
        v = a[:, b - half:b + half + 1]; acc.append(v if st == '+' else v[:, ::-1])
    return (np.mean(acc, 0) if acc else np.full((len(TR), 2 * half + 1), np.nan)), len(acc)
