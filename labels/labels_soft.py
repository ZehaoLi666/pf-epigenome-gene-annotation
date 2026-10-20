"""Gaussian soft end labels from labels v2 (Lml2_pool.npz): TSS/TTS channels set to exp(-d^2/2 sigma^2) around every positive end
cluster within +-200 bp (weight = end weight), body channel unchanged. Usage: python labels_soft.py <sigma_bp> <out.npz>"""
import sys, numpy as np, pandas as pd
sig = float(sys.argv[1]) / 20; HALF = 10; out_f = sys.argv[2]
import os
Z = np.load(os.environ.get('SOFT_BASE', 'data/Lml2_pool.npz')); CH = sorted({k.split('__')[0] for k in Z.files})
PP = pd.concat([pd.read_csv(f, sep='\t')[['chrom', 'strand', 'pos', 'end_type', 'weight']] for f in os.environ.get('SOFT_POS', 'data/ml_end_positives.tsv,data/ml2_added_end_positives.tsv').split(',')], ignore_index=True)
out = {}
for c in CH:
    Y = Z[f'{c}__Y'].astype(np.float16); W = Z[f'{c}__W'].astype(np.float32); n = Y.shape[-1]
    for ch, et in ((0, 'TSS'), (1, 'TTS')):
        for k, st in ((0, '+'), (1, '-')):
            q = PP[(PP.chrom == c) & (PP.strand == st) & (PP.end_type == et)]; pos_bins = Y[k, ch] >= 1
            soft = np.zeros(n, np.float32); wmax = np.zeros(n, np.float32)
            for b0, w0 in zip(q.pos.values // 20, q.weight.values):
                lo, hi = max(0, b0 - HALF), min(n, b0 + HALF + 1); d = np.arange(lo, hi) - b0
                soft[lo:hi] = np.maximum(soft[lo:hi], np.exp(-d ** 2 / (2 * sig ** 2))); wmax[lo:hi] = np.maximum(wmax[lo:hi], w0)
            region = wmax > 0
            Y[k, ch] = np.where(region, soft, np.where(pos_bins, 0, Y[k, ch])).astype(np.float16)
            W[k, ch] = np.where(region, np.maximum(W[k, ch], wmax), np.where(pos_bins, 0.0, W[k, ch]))
    out[f'{c}__Y'] = Y; out[f'{c}__W'] = W.astype(np.float16)
np.savez_compressed(out_f, **out); print('ok', out_f)
