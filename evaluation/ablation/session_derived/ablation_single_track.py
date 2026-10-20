"""Session-derived computation of ablation_single_track.tsv.

Recovered from the record of the interactive analysis session in which ablation_single_track.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: stage_TPM_9stages.tsv, gene_heterochromatin_groups.tsv, cv_score.py, gene_windows_truth.tsv, cv_gene_calls.tsv.gz, cv_point_metrics.tsv
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import result  # paths: config/paths.yaml
import numpy as np
import pandas as pd
from scipy.stats import binomtest

exec(open(result('cv_score.py')).read())

W = pd.read_csv(result('gene_windows_truth.tsv'), sep='\t')
G = pd.read_csv(result('cv_gene_calls.tsv.gz'), sep='\t')
PT = pd.read_csv(result('cv_point_metrics.tsv'), sep='\t')
tpm = pd.read_csv(result('stage_TPM_9stages.tsv'), sep='\t').set_index('Geneid')
het = pd.read_csv(result('gene_heterochromatin_groups.tsv'), sep='\t')

grp = W[['gene_id']].copy()
e = grp.gene_id.map(tpm[['Ring', 'Trophozoite', 'Schizont']].max(axis=1))
grp['low_expression_Q1'] = e <= e.quantile(0.2)
grp['heterochromatic'] = grp.gene_id.isin(het[het.group.isin(['var/rifin/stevor/Pfmc-2TM', 'other HP1+H3K9me3-high'])].gene_id)

S2 = score(G[G.model.isin(['A_full2']) | G.model.str.startswith('T')], W, grp, 'A_full2')

def bh(p):
    p = np.asarray(p, float); o = np.argsort(p); r = np.empty(len(p)); m = len(p)
    q = p[o] * m / np.arange(1, m + 1); q = np.minimum.accumulate(q[::-1])[::-1]; r[o] = np.minimum(q, 1); return r

tk = S2[S2.model.str.startswith(('TL_', 'TA_')) & (S2.subset == 'all') & S2.truth.str.contains('indep')].copy()
tk['q_BH'] = bh(tk.p_vs_ref)

tr_out = S2[S2.model.str.startswith(('TL_', 'TA_', 'TS_')) | (S2.model == 'A_full2')].merge(tk[['model', 'truth', 'q_BH']], on=['model', 'truth'], how='left')
tr_out.to_csv('ablation_single_track.tsv', sep='\t', index=False)
