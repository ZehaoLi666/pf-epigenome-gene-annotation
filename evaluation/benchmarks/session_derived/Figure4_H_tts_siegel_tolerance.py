"""Session-derived computation of Figure4_H_tts_siegel_tolerance.tsv.

Recovered from the record of the interactive analysis session in which Figure4_H_tts_siegel_tolerance.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: benchmark_method_predictions.tsv, benchmark_truth_sets.tsv, siegel2014_pas_per_gene.tsv
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import result  # paths: config/paths.yaml
import pandas as pd
import numpy as np

M = pd.read_csv(result('benchmark_method_predictions.tsv'), sep='\t')
T = pd.read_csv(result('benchmark_truth_sets.tsv'), sep='\t')
SP = pd.read_csv(result('siegel2014_pas_per_gene.tsv'), sep='\t')

TM = ['Ours_chromatin', 'FixedOffset', 'PlasmoDB_v71', 'StringTie_short', 'StringTie_long']
TLB = {'Ours_chromatin': 'Our model', 'FixedOffset': 'Fixed UTR length', 'PlasmoDB_v71': 'PlasmoDB v71',
       'StringTie_short': 'StringTie, short reads', 'StringTie_long': 'StringTie, long reads'}

TOLS = [25, 50, 100, 150, 200, 300, 500]

def hit(df, m, tc, tol=200, end='TSS'):
    return ((df[f'{m}__{end}'] - df[tc]).abs() <= tol).values

S = SP[SP.n_clusters > 0].merge(M, on='gene_id').merge(T[['gene_id', 'LongRead__TTS']], on='gene_id')

rowsH = []
for m in TM:
    for t in TOLS:
        rowsH.append(dict(method=TLB[m], key=m, tolerance_bp=t, acc=hit(S, m, 'PAS_major', t, 'TTS').mean(), n=len(S)))
TH = pd.DataFrame(rowsH)

TH.to_csv('Figure4_H_tts_siegel_tolerance.tsv', sep='\t', index=False, float_format='%.4g')
