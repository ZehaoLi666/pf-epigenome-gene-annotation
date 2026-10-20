"""Local scoring of cross-validated per-gene calls against independent truths (TSS consensus; Siegel 2014 poly(A)) and long-read ends,
overall and in gene subsets (lowest expression quintile; heterochromatic / clonally variant genes). Paired sign test vs a reference model."""
import numpy as np, pandas as pd
from scipy.stats import binomtest
def score(G, W, groups, ref, tol=200, nb=1000, seed=0):
    rng = np.random.default_rng(seed); X = G.merge(W, on='gene_id').merge(groups, on='gene_id', how='left'); rows = []
    for et, tcol, tl in (('TSS', 'cons_TSS', 'TSS consensus (independent)'), ('TTS', 'siegel_PAS', 'Siegel 2014 poly(A) (independent)'), ('TSS', 'lr_TSS', 'Long-read TSS'), ('TTS', 'lr_TTS', 'Long-read TTS')):
        for gname in ['all'] + [c for c in groups.columns if c != 'gene_id']:
            Q = X[X[tcol].notna() & (True if gname == 'all' else X[gname].fillna(False).astype(bool))]
            R = Q[Q.model == ref].set_index('gene_id')
            for m, q in Q.groupby('model'):
                q = q.set_index('gene_id'); h = ((q[et] - q[tcol]).abs() <= tol).values
                if len(h) == 0: continue
                bs = [h[rng.integers(0, len(h), len(h))].mean() for _ in range(nb)]
                common = q.index.intersection(R.index); a = (q.loc[common, et] - q.loc[common, tcol]).abs() <= tol; b = (R.loc[common, et] - R.loc[common, tcol]).abs() <= tol
                w, l = int((a & ~b).sum()), int((~a & b).sum())
                rows.append(dict(model=m, folds=q.folds.iloc[0], truth=tl, end=et, subset=gname, n=len(h), acc200=h.mean(), lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5),
                                 coverage=q[et].notna().mean(), delta_vs_ref=h.mean() - b.mean() if len(common) == len(h) else np.nan, better=w, worse=l, p_vs_ref=binomtest(w, w + l).pvalue if w + l else 1.0))
    return pd.DataFrame(rows)
