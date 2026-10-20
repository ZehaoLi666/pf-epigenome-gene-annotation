"""Session-derived computation of gene_daft_expression_truth.tsv.gz.

Recovered from the record of the interactive analysis session in which gene_daft_expression_truth.tsv.gz was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: stage_TPM_9stages.tsv, gene_heterochromatin_groups.tsv, gene_boundaries_v3_all_genes.tsv, benchmark_method_predictions.tsv, gene_expression_strata_calls.tsv.gz
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import result  # paths: config/paths.yaml
import os, numpy as np, pandas as pd

# Load main data tables
B = pd.read_csv(result('gene_boundaries_v3_all_genes.tsv'), sep='\t')
MX = pd.read_csv(result('benchmark_method_predictions.tsv'), sep='\t')
tpm = pd.read_csv(result('stage_TPM_9stages.tsv'), sep='\t')

# Load pre-computed gene expression strata calls (checkpoint)
D = pd.read_csv(result('gene_expression_strata_calls.tsv.gz'), sep='\t')

# Load DAFT-seq gene expression table
E = pd.read_csv('daft_gene_expression_v68.tsv.gz', sep='\t')
H = [c for c in E.columns if c.startswith('h')]
E['DAFT_max'] = E[H].max(axis=1)

# Merge with boundary table for final_src, final_TSS, final_TTS
X = D.merge(B[['gene_id', 'final_src', 'final_TSS', 'final_TTS']], on='gene_id').merge(
    E[['gene_id', 'DAFT_max', 'antisense_max', 'cds_len'] + H], on='gene_id', how='left')

T = tpm.set_index('Geneid')
D2 = B[['gene_id', 'chrom', 'strand', 'atg', 'stop', 'lr_group', 'lr_TSS', 'lr_TTS']].copy()
D2['IDC_TPM'] = D2.gene_id.map(T[['Ring', 'Trophozoite', 'Schizont']].max(axis=1))

QL = ['Q1 (lowest)', 'Q2', 'Q3', 'Q4', 'Q5 (highest)']
X = X.merge(D2[['gene_id', 'IDC_TPM']], on='gene_id', how='left', suffixes=('', '_d2'))
if 'IDC_TPM_d2' in X.columns:
    X['IDC_TPM'] = X['IDC_TPM'].fillna(X['IDC_TPM_d2'])
    X = X.drop(columns=['IDC_TPM_d2'])

X['Qd'] = pd.qcut(X.DAFT_max.rank(method='first'), 5, labels=QL)
X['final__TSS'] = X.final_TSS
X['final__TTS'] = X.final_TTS

# Merge independent truth columns if not already present
from scipy.stats import spearmanr

TR = pd.read_csv(result('benchmark_method_predictions.tsv'), sep='\t')
if 'Consensus_5end__TSS' not in X.columns:
    # Need to bring in Consensus_5end__TSS from benchmark predictions or truth sets
    if 'Consensus_5end__TSS' in TR.columns:
        X = X.merge(TR[['gene_id', 'Consensus_5end__TSS']], on='gene_id', how='left')

# PAS_major from Siegel data - check if already in D
if 'PAS_major' not in X.columns and 'siegel_PAS' in X.columns:
    X['PAS_major'] = X['siegel_PAS']

# Bring in Ours_chromatin and PlasmoDB_v71 columns from MX if not present
for col in ['Ours_chromatin__TSS', 'Ours_chromatin__TTS', 'PlasmoDB_v71__TSS', 'PlasmoDB_v71__TTS', 'FixedOffset__TSS', 'FixedOffset__TTS']:
    if col not in X.columns and col in MX.columns:
        X = X.merge(MX[['gene_id', col]], on='gene_id', how='left')

X[['gene_id', 'final_src', 'DAFT_max', 'antisense_max'] + H + ['IDC_TPM', 'Qd', 'final_TSS', 'final_TTS', 'Ours_chromatin__TSS', 'Ours_chromatin__TTS', 'Consensus_5end__TSS', 'PAS_major']].to_csv('gene_daft_expression_truth.tsv.gz', sep='\t', index=False)
