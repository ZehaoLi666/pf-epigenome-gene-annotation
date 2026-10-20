"""Session-derived computation of ablation_cv_input_classes.tsv.

Recovered from the record of the interactive analysis session in which ablation_cv_input_classes.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: stage_TPM_9stages.tsv, gene_heterochromatin_groups.tsv, cv_score.py, gene_windows_truth.tsv, cv_gene_calls.tsv.gz
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import result  # paths: config/paths.yaml
import os, numpy as np, pandas as pd
from scipy.stats import binomtest

# Load cv_score.py helper
exec(open(result('cv_score.py')).read())

W = pd.read_csv(result('gene_windows_truth.tsv'), sep='\t')
G = pd.read_csv(result('cv_gene_calls.tsv.gz'), sep='\t')
tpm = pd.read_csv(result('stage_TPM_9stages.tsv'), sep='\t').set_index('Geneid')
het = pd.read_csv(result('gene_heterochromatin_groups.tsv'), sep='\t')

grp = W[['gene_id']].copy()
e = grp.gene_id.map(tpm[['Ring', 'Trophozoite', 'Schizont']].max(axis=1))
grp['low_expression_Q1'] = e <= e.quantile(0.2)
grp['heterochromatic'] = grp.gene_id.isin(het[het.group.isin(['var/rifin/stevor/Pfmc-2TM', 'other HP1+H3K9me3-high'])].gene_id)

P1 = ['A_full', 'L_nochrom', 'L_noonehot', 'L_noglm', 'L_nogro', 'L_nopoly', 'L_nodaft', 'L_nostrand', 'N_noshortread', 'S_chrom', 'S_dna', 'S_rna', 'S_chromdna']

S = score(G[G.folds == 'ABC'], W, grp, 'A_full')

S.to_csv('ablation_cv_input_classes.tsv', sep='\t', index=False)
