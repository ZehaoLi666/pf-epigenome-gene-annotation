"""Session-derived computation of fig3F_chromatin_by_expression_longread.tsv.

Recovered from the record of the interactive analysis session in which fig3F_chromatin_by_expression_longread.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: gene_daft_expression_truth.tsv.gz, fig5A_input_ablation_cv.tsv, gene_windows_truth.tsv, cv_gene_calls.tsv.gz, chappell2020_5UTRseq_tagclusters.tsv, gene_boundaries_v3_all_genes.tsv, hc_reference.bed.gz, gene_transcript_end_catalogue.tsv, ablation_cv_input_classes.tsv, benchmark_method_predictions.tsv, gene_vs_readend_v48.tsv, transcript_stage_confirmation.tsv.gz, shaw2021_5CAPture_5CN_clusters.tsv
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import result  # paths: config/paths.yaml
import pandas as pd
import numpy as np
from scipy import stats

GE = pd.read_csv(result('gene_daft_expression_truth.tsv.gz'), sep='\t').set_index('gene_id')
WT = pd.read_csv(result('gene_windows_truth.tsv'), sep='\t').set_index('gene_id')
CG = pd.read_csv(result('cv_gene_calls.tsv.gz'), sep='\t')

CGi = {m: q.set_index('gene_id') for m, q in CG.groupby('model')}

def hits(model, tcol, end):
    t = WT[tcol].dropna()
    c = CGi[model][end].reindex(t.index)
    return ((c - t).abs() <= 200).astype(float)

rng = np.random.default_rng(0)
NB = 2000
QS = ['Q1 (lowest)', 'Q2', 'Q3', 'Q4', 'Q5 (highest)']

rows = []
for tcol, end, tl in [('lr_TSS', 'TSS', 'Long-read TSS'), ('lr_TTS', 'TTS', 'Long-read TTS')]:
    hm = {m: hits(m, tcol, end) for m in ('A_full', 'L_nochrom')}
    q = GE.Qd.reindex(hm['A_full'].index)
    for k in QS:
        r = dict(truth=tl, stratum=k, n=int((q == k).sum()),
                 daft_median=GE.DAFT_max.reindex(hm['A_full'].index)[q == k].median())
        for m, hh in hm.items():
            v = hh[q == k].values
            bs = v[rng.integers(0, len(v), (NB, len(v)))].mean(1)
            r.update({f'acc_{m}': v.mean(), f'lo_{m}': np.percentile(bs, 2.5), f'hi_{m}': np.percentile(bs, 97.5)})
        d = (hm['L_nochrom'] - hm['A_full'])[q == k]
        w, l = int((d > 0).sum()), int((d < 0).sum())
        r.update(delta=d.mean(), better=w, worse=l,
                 p_sign=stats.binomtest(w, w + l).pvalue if w + l else 1.0)
        rows.append(r)

TF3 = pd.DataFrame(rows)
TF3.to_csv('fig3F_chromatin_by_expression_longread.tsv', sep='\t', index=False)
