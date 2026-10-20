"""Session-derived computation of Figure4_D_per_truth_set_acc200.tsv.

Recovered from the record of the interactive analysis session in which Figure4_D_per_truth_set_acc200.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: benchmark_method_predictions.tsv, benchmark_truth_sets.tsv, gene_boundaries_v3_all_genes.tsv, benchmark_acc200_bootstrap.tsv
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import result  # paths: config/paths.yaml
import pandas as pd
import numpy as np

M = pd.read_csv(result('benchmark_method_predictions.tsv'), sep='\t')
T = pd.read_csv(result('benchmark_truth_sets.tsv'), sep='\t')
G = pd.read_csv(result('gene_boundaries_v3_all_genes.tsv'), sep='\t')[['gene_id', 'chrom', 'strand', 'atg', 'stop', 'lr_group']]
BOOT = pd.read_csv(result('benchmark_acc200_bootstrap.tsv'), sep='\t')

X = M.merge(T, on='gene_id').merge(G, on='gene_id', how='left')

METH = ['Ours_chromatin', 'FixedOffset', 'PlasmoDB_v71', 'ChromatinPeak', 'ChromHMM', 'StringTie_short', 'StringTie_long']
LAB = {'Ours_chromatin': 'Our model', 'FixedOffset': 'Fixed UTR length', 'PlasmoDB_v71': 'PlasmoDB v71',
       'ChromatinPeak': 'Chromatin peaks', 'ChromHMM': 'ChromHMM promoter', 'StringTie_short': 'StringTie, short reads',
       'StringTie_long': 'StringTie, long reads'}
TS = ['Adjalley2016_capseq', 'Chappell2020_5UTRseq', 'Shaw2021_5CAPture']

def hit(df, m, tc, tol=200, end='TSS'):
    return ((df[f'{m}__{end}'] - df[tc]).abs() <= tol).values

rowsD = []
for t in TS:
    q = X[X[t + '__TSS'].notna()]
    for m in METH:
        b = BOOT[(BOOT.truth == t) & (BOOT.method == m)].iloc[0]
        a = hit(q, m, t + '__TSS').mean()
        assert abs(a - b.acc200) < 1e-9
        rowsD.append(dict(truth=t, method=LAB[m], key=m, n=len(q), acc200=a, ci_lo=b.lo, ci_hi=b.hi))
TD = pd.DataFrame(rowsD)

TD.to_csv('Figure4_D_per_truth_set_acc200.tsv', sep='\t', index=False, float_format='%.4g')
