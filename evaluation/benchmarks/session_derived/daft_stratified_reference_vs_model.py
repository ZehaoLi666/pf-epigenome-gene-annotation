"""Session-derived computation of daft_stratified_reference_vs_model.tsv.

Recovered from the record of the interactive analysis session in which daft_stratified_reference_vs_model.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: stage_TPM_9stages.tsv, gene_heterochromatin_groups.tsv, gene_boundaries_v3_all_genes.tsv, benchmark_method_predictions.tsv, gene_expression_strata_calls.tsv.gz
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import result  # paths: config/paths.yaml
import os, numpy as np, pandas as pd
from scipy.stats import spearmanr

# Load data
B = pd.read_csv(result('gene_boundaries_v3_all_genes.tsv'), sep='\t')
MX = pd.read_csv(result('benchmark_method_predictions.tsv'), sep='\t')
D = pd.read_csv(result('gene_expression_strata_calls.tsv.gz'), sep='\t')

# Load DAFT-seq gene expression (pre-computed)
E = pd.read_csv('daft_gene_expression_v68.tsv.gz', sep='\t')
H = [c for c in E.columns if c.startswith('h')]
E['DAFT_max'] = E[H].max(axis=1)

# Merge everything
X = D.merge(B[['gene_id', 'final_src', 'final_TSS', 'final_TTS']], on='gene_id').merge(E[['gene_id', 'DAFT_max', 'antisense_max', 'cds_len'] + H], on='gene_id', how='left')

QL = ['Q1 (lowest)', 'Q2', 'Q3', 'Q4', 'Q5 (highest)']
X['Qd'] = pd.qcut(X.DAFT_max.rank(method='first'), 5, labels=QL)
X['final__TSS'] = X.final_TSS
X['final__TTS'] = X.final_TTS

rng = np.random.default_rng(0)
rows = []
GRP = {'Reference genes (long-read)': X[X.final_src == 'long-read'], 'Model-predicted genes': X[X.final_src == 'model'], 'No prediction': X[X.final_src == 'none']}
METH = [('final', 'Final annotation (reference / model)'), ('Ours_chromatin', 'Model (out-of-fold)'), ('PlasmoDB_v71', 'PlasmoDB v71'), ('FixedOffset', 'Fixed UTR length')]
for gname, G_ in GRP.items():
    for et, tcol, tl in (('TSS', 'Consensus_5end__TSS', 'TSS consensus'), ('TTS', 'PAS_major', 'Siegel 2014 poly(A)')):
        for qv in QL + ['all']:
            q = G_[G_[tcol].notna() & ((G_.Qd == qv) if qv != 'all' else True)]
            if len(q) == 0: continue
            for m, ml in METH:
                h = ((q[f'{m}__{et}'] - q[tcol]).abs() <= 200).values; bs = [h[rng.integers(0, len(h), len(h))].mean() for _ in range(1000)]
                rows.append(dict(group=gname, truth=tl, stratum=qv, method=ml, n=len(q), DAFT_median=q.DAFT_max.median(), acc200=h.mean(), lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5),
                                 acc50=((q[f'{m}__{et}'] - q[tcol]).abs() <= 50).mean(), median_err=(q[f'{m}__{et}'] - q[tcol]).abs().median(), has_call=q[f'{m}__{et}'].notna().mean()))
R = pd.DataFrame(rows)
R.to_csv('daft_stratified_reference_vs_model.tsv', sep='\t', index=False)
