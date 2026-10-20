"""Score ablation per-gene calls (chr 1-5) against independent truths: TSS consensus (>=2 of 3 5'-end datasets) and Siegel 2014 major poly(A) site.
Accuracy@200 with no-call = miss; bootstrap CI over genes; paired McNemar-type sign test vs the full model (A_full)."""
import numpy as np, pandas as pd
from scipy.stats import binomtest
W = pd.read_csv('eval/gene_windows_truth.tsv', sep='\t'); Gc = pd.read_csv('res/abl_gene_calls.tsv', sep='\t').merge(W, on='gene_id')
rng = np.random.default_rng(0); rows = []
for (et, tcol), tl in ((('TSS', 'cons_TSS'), 'TSS consensus (independent)'), (('TSS', 'lr_TSS'), 'Long-read TSS'), (('TTS', 'siegel_PAS'), 'Siegel 2014 poly(A) (independent)'), (('TTS', 'lr_TTS'), 'Long-read TTS')):
    ref = Gc[Gc.model == 'A_full'].set_index('gene_id')
    for m, q in Gc.groupby('model'):
        q = q[q[tcol].notna()].set_index('gene_id'); h = ((q[et] - q[tcol]).abs() <= 200).astype(float)
        bs = [h.values[rng.integers(0, len(h), len(h))].mean() for _ in range(1000)]
        hr = ((ref.loc[q.index, et] - ref.loc[q.index, tcol]).abs() <= 200)
        w, l = int((h.astype(bool) & ~hr).sum()), int((~h.astype(bool) & hr).sum())
        rows.append(dict(model=m, end=et, truth=tl, n=len(h), acc200=h.mean(), lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5), coverage=q[et].notna().mean(),
                         median_err=float((q[et] - q[tcol]).abs().median()), better_than_full=w, worse_than_full=l, p_vs_full=binomtest(w, w + l).pvalue if w + l else 1.0))
pd.DataFrame(rows).to_csv('res/abl_independent_truth.tsv', sep='\t', index=False); print('ok')
