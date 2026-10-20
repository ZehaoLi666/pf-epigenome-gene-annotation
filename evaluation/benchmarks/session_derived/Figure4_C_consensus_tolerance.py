"""Session-derived computation of Figure4_C_consensus_tolerance.tsv.

Recovered from the record of the interactive analysis session in which Figure4_C_consensus_tolerance.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: benchmark_method_predictions.tsv, benchmark_truth_sets.tsv, gene_boundaries_v3_all_genes.tsv
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import result  # paths: config/paths.yaml
import pandas as pd
import numpy as np

M = pd.read_csv(result('benchmark_method_predictions.tsv'), sep='\t')
T = pd.read_csv(result('benchmark_truth_sets.tsv'), sep='\t')
G = pd.read_csv(result('gene_boundaries_v3_all_genes.tsv'), sep='\t')[['gene_id', 'chrom', 'strand', 'atg', 'stop', 'lr_group']]

X = M.merge(T, on='gene_id').merge(G, on='gene_id', how='left')

cons = X['Consensus_5end__TSS'].notna()
Q = X[cons].copy()
NC = len(Q)

METH = ['Ours_chromatin', 'FixedOffset', 'PlasmoDB_v71', 'ChromatinPeak', 'ChromHMM', 'StringTie_short', 'StringTie_long']
LAB = {'Ours_chromatin': 'Our model', 'FixedOffset': 'Fixed UTR length', 'PlasmoDB_v71': 'PlasmoDB v71',
       'ChromatinPeak': 'Chromatin peaks', 'ChromHMM': 'ChromHMM promoter', 'StringTie_short': 'StringTie, short reads',
       'StringTie_long': 'StringTie, long reads'}

TOLS = [25, 50, 100, 150, 200, 300, 500]

def hit(df, m, tc, tol=200, end='TSS'):
    return ((df[f'{m}__{end}'] - df[tc]).abs() <= tol).values

TC = pd.DataFrame([dict(method=LAB[m], key=m, tolerance_bp=t, acc=hit(Q, m, 'Consensus_5end__TSS', t).mean()) for m in METH for t in TOLS])

TC.to_csv('Figure4_C_consensus_tolerance.tsv', sep='\t', index=False, float_format='%.4g')
