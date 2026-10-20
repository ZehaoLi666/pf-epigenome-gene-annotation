"""Session-derived computation of siegel_pas_tts_benchmark.tsv.

Recovered from the record of the interactive analysis session in which siegel_pas_tts_benchmark.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: gene_boundaries_v3_all_genes.tsv, 12864_2013_7013_MOESM10_ESM.xlsx, benchmark_method_predictions.tsv
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import result  # paths: config/paths.yaml
import numpy as np
import pandas as pd
import openpyxl
from scipy.stats import binomtest

G = pd.read_csv(result('gene_boundaries_v3_all_genes.tsv'), sep='\t')

wb = openpyxl.load_workbook(result('12864_2013_7013_MOESM10_ESM.xlsx'), read_only=True)
ws = wb[wb.sheetnames[0]]
rows = list(ws.iter_rows(values_only=True))
h = [i for i, r in enumerate(rows) if r and r[0] == 'gene ID'][0]
SG = pd.DataFrame([r[:11] for r in rows[h + 1:] if r[0]], columns=['gene_id', 'strand', 'orf3', 'orf5', 'n_clusters', 'd_closest', 'pos_closest', 'd_furthest', 'pos_furthest', 'd_major', 'pos_major'])
for c_ in ['orf3', 'orf5', 'n_clusters', 'pos_closest', 'pos_furthest', 'pos_major', 'd_major']:
    SG[c_] = pd.to_numeric(SG[c_], errors='coerce')

SG['PAS_major'] = SG.pos_major - 1
SG['PAS_closest'] = SG.pos_closest - 1
SG['PAS_furthest'] = SG.pos_furthest - 1

Mx = pd.read_csv(result('benchmark_method_predictions.tsv'), sep='\t')

cols = {'Ours_chromatin__TTS': 'Ours (chromatin)', 'Ours_nochromatin__TTS': 'Ours (no chromatin)', 'PlasmoDB_v71__TTS': 'PlasmoDB v71', 'FixedOffset__TTS': 'Fixed 3′UTR length',
        'StringTie_short__TTS': 'StringTie, short reads', 'StringTie_long__TTS': 'StringTie, long reads', 'LongRead_v3__TTS': 'Long-read TTS (training truth)'}
E = SG[SG.n_clusters > 0][['gene_id', 'n_clusters', 'PAS_major', 'PAS_closest', 'PAS_furthest']].merge(Mx[['gene_id'] + list(cols)], on='gene_id').merge(G[['gene_id', 'strand', 'stop', 'lr_TTS']], on='gene_id')

rng = np.random.default_rng(2)
rows = []

def score(sub, tag):
    for k, lab in cols.items():
        q = sub[sub[k].notna()]
        d = (q[k] - q.PAS_major).abs()
        da = np.minimum.reduce([(q[k] - q[c_]).abs() for c_ in ('PAS_major', 'PAS_closest', 'PAS_furthest')])
        hit = ((sub[k] - sub.PAS_major).abs() <= 200).values.astype(float)
        bs = [hit[rng.integers(0, len(hit), len(hit))].mean() for _ in range(500)]
        rows.append(dict(subset=tag, method=lab, n_genes=len(sub), coverage=q.shape[0] / len(sub), median_err_major=d.median(), **{f'acc@{t}': float(((sub[k] - sub.PAS_major).abs() <= t).mean()) for t in (50, 100, 200, 500)},
                         acc200_lo=np.percentile(bs, 2.5), acc200_hi=np.percentile(bs, 97.5), acc200_any_PAS=float((da <= 200).sum() / len(sub)),
                         signed_median=float(np.median((q[k] - q.PAS_major) * np.where(q.strand == '+', 1, -1)))))

score(E, 'All genes with Siegel PAS')
score(E[E.lr_TTS.isna()], 'No long-read TTS')
score(E[E.n_clusters == 1], 'Single PAS cluster')
R = pd.DataFrame(rows)
R.to_csv('siegel_pas_tts_benchmark.tsv', sep='\t', index=False)
