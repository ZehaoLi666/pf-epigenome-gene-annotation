"""Independent CAGE test of predicted TSSs (CAGE was not used in v3 labels or inputs).
CAGE score at a position = max over stages (ring, trophozoite, schizont) of the same-strand CTSS CPM summed within +-100 bp.
Groups: genes with a v3 long-read full-CDS model (truth available) and remaining genes. Positions compared per gene:
  predicted TSS, annotated v68 TSS, ATG, and 3 random positions drawn uniformly in the same TSS search window (control).
Also: distance from the predicted TSS to the strongest CAGE position in the window, and for truth genes the distance to the long-read TSS."""
from pfconfig import LAB_DIR  # paths: config/paths.yaml
import sys, numpy as np, pandas as pd, pyBigWig
LAB = LAB_DIR; inp, out = sys.argv[1], sys.argv[2]
G = pd.read_csv(inp, sep='\t').merge(pd.read_csv('gene_truth_v3.tsv', sep='\t'), on='gene_id', how='left'); G['lr_group'] = G.lr_group.fillna('no HC model')
BW = {(s, st): pyBigWig.open(f"{LAB}/lncRNA project tracks/cage_tss/CTSS_{s}_{st}_CPM.bw") for s in ('ring', 'trophozoite', 'schizont') for st in ('pos', 'neg')}
def track(c, strand, a, b):
    k = 'pos' if strand == '+' else 'neg'
    return np.max([np.nan_to_num(np.array(BW[(s, k)].values(c, max(0, a), b), float)) for s in ('ring', 'trophozoite', 'schizont')], axis=0)
def score(c, strand, p, w=100):
    if pd.isna(p): return np.nan
    p = int(p); v = track(c, strand, p - w, p + w + 1); return float(v.sum())
rng = np.random.default_rng(0); rows = []; UTR = G.utr5_pred.dropna().values.astype(int)
for r in G.itertuples():
    d = dict(gene_id=r.gene_id)
    for nm in ('pred_TSS', 'annot_tss', 'atg'): d[f'cage_{nm}'] = score(r.chrom, r.strand, getattr(r, nm))
    lo, hi = int(r.tss_win_lo), int(r.tss_win_hi)
    if hi - lo > 50:
        d['cage_random'] = float(np.mean([score(r.chrom, r.strand, x) for x in rng.integers(lo, hi, 3)]))
        d['cage_random_pos1'] = score(r.chrom, r.strand, int(rng.integers(lo, hi)))
        v = track(r.chrom, r.strand, lo - 100, hi + 100); sm = np.convolve(v, np.ones(201), 'valid')       # +-100 bp sums at lo..hi
        d['cage_window_max_pos'] = lo + int(np.argmax(sm)) if sm.max() > 0 else np.nan; d['cage_window_max'] = float(sm.max())
        for nm in ('pred_TSS', 'annot_tss', 'atg', 'lr_TSS'):
            p = getattr(r, nm, np.nan)
            if pd.notna(p) and lo <= int(p) < lo + len(sm) and sm.max() > 0: d[f'pct_{nm}'] = float((sm < sm[int(p) - lo]).mean() + 0.5 * (sm == sm[int(p) - lo]).mean())
        u = UTR[rng.integers(0, len(UTR))]; q = int(r.atg - u) if r.strand == '+' else int(r.atg + u)       # distance-matched control (another gene's predicted 5'UTR)
        if lo <= q < lo + len(sm) and sm.max() > 0: d['pct_dist_matched'] = float((sm < sm[q - lo]).mean() + 0.5 * (sm == sm[q - lo]).mean())
        x = int(rng.integers(0, len(sm))); d['pct_random'] = float((sm < sm[x]).mean() + 0.5 * (sm == sm[x]).mean()) if sm.max() > 0 else np.nan
    rows.append(d)
C = pd.DataFrame(rows); G = G.merge(C, on='gene_id'); G['d_pred_to_cage_peak'] = (G.pred_TSS - G.cage_window_max_pos).abs()
G.to_csv(out, sep='\t', index=False); print(len(G))
