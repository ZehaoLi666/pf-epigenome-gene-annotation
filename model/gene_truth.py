"""Per-gene long-read truth (v3): dominant TSS / TTS of v3-confirmed full-CDS models, FL-weighted (ends merged <=50 bp), and gene group:
'confirmed' (has v3-confirmed full-CDS model), 'truncated-only' (only v3-confirmed truncated models), 'HC model, unconfirmed', 'no long-read model'."""
import numpy as np, pandas as pd
T = pd.read_csv('../v3/models_v3.tsv.gz', sep='\t', low_memory=False); T = T[T.unit_class == 'mRNA']
S = ['ring', 'troph', 'schizont']; T['FL'] = T[[f'FL_{s}' for s in S]].sum(axis=1) + 0.01
T['tss'] = np.where(T.strand == '+', T.start, T.end - 1); T['tes'] = np.where(T.strand == '+', T.end - 1, T.start)
def dom(p, w):
    o = np.argsort(p); p, w = p[o], w[o]; cl = [[p[0], w[0], [p[0]], [w[0]]]]
    for x, y in zip(p[1:], w[1:]):
        if x - cl[-1][2][-1] <= 50: cl[-1][2].append(x); cl[-1][3].append(y); cl[-1][1] += y
        else: cl.append([x, y, [x], [y]])
    b = max(cl, key=lambda c: c[1]); return int(np.average(b[2], weights=b[3])), len(cl)
rows = []
F = T[T.full_cds & T.conf3_any]
for loc, d in F.groupby('locus'):
    t, n5 = dom(d.tss.values, d.FL.values); e, n3 = dom(d.tes.values, d.FL.values)
    rows.append(dict(gene_id=loc, lr_TSS=t, lr_TTS=e, lr_n_TSS_clusters=n5, lr_n_TTS_clusters=n3))
R = pd.DataFrame(rows)
grp = pd.Series('no long-read model', index=pd.Index(T.locus.unique()))
grp[T[~T.conf3_any].locus.unique()] = 'HC model, unconfirmed'
grp[T[(T.trunc5 | T.trunc3) & T.conf3_any].locus.unique()] = 'truncated-only'
grp[F.locus.unique()] = 'confirmed'
R = grp.rename('lr_group').rename_axis('gene_id').reset_index().merge(R, on='gene_id', how='left')
R.to_csv('gene_truth_v3.tsv', sep='\t', index=False); print(R.lr_group.value_counts().to_dict())
