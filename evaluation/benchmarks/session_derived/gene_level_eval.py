"""Session-derived computation of gene_level_eval.tsv.

Recovered from the record of the interactive analysis session in which gene_level_eval.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: gene_boundaries_v3_all_genes.tsv, gene_end_classes.tsv, v71_vs_reference_ends.tsv, gene_transcript_end_catalogue.tsv, gene_vs_readend_v48.tsv, transcript_stage_confirmation.tsv.gz
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import result  # paths: config/paths.yaml
import pandas as pd
import numpy as np

BOUND = result('gene_boundaries_v3_all_genes.tsv')
CGPATH = result('gene_end_classes.tsv')

RB = pd.read_csv(BOUND, sep="\t").set_index('gene_id')
S = pd.read_csv(result('v71_vs_reference_ends.tsv'), sep="\t").set_index('gene_id')
C = RB[RB.lr_group == 'confirmed'].copy()
sg = np.where(C.strand == '+', 1, -1)
C['dTSS'] = (C.model_TSS - C.lr_TSS) * sg
C['dTTS'] = (C.model_TTS - C.lr_TTS) * sg

C = C.join(S[['ref_TSS_dists', 'ref_TTS_dists']])
C['u5_lr'] = (C.atg - C.lr_TSS) * sg
C['u5_m'] = (C.atg - C.model_TSS) * sg
C['u3_lr'] = (C.lr_TTS - C.stop) * sg
C['u3_m'] = (C.model_TTS - C.stop) * sg


def parse(x):
    return [float(v) for v in str(x).split(';')] if isinstance(x, str) and x else []


def cat(row, end):
    m = row['u5_m'] if end == 'TSS' else row['u3_m']
    lr = row['u5_lr'] if end == 'TSS' else row['u3_lr']
    refs = parse(row['ref_TSS_dists'] if end == 'TSS' else row['ref_TTS_dists'])
    if pd.isna(m):
        return 'no call'
    if abs(m - lr) <= 200:
        return 'main end'
    if any(abs(m - r) <= 200 for r in refs):
        return 'alternative end'
    return 'longer' if m > lr else 'shorter'


for end in ['TSS', 'TTS']:
    C['cat_' + end] = C.apply(lambda r: cat(r, end), axis=1)

P = pd.read_csv("cov/gene_longread_coverage.tsv", sep="\t", index_col=0)
C['reads'] = P.ALL_reads.reindex(C.index)
C['expr_q'] = pd.qcut(C.reads, 4, labels=['Q1 (low)', 'Q2', 'Q3', 'Q4 (high)'])
C['n_tss'] = C.ref_TSS_dists.apply(lambda x: len(parse(x)))
C['ntss_bin'] = pd.cut(C.n_tss, [0, 1, 2, 3, 100], labels=['1', '2', '3', '≥4'])
C['u5_bin'] = pd.cut(C.u5_lr, [-1e9, 200, 500, 1000, 1e9], labels=['<200', '200–500', '500–1,000', '>1,000'])
C['u3_bin'] = pd.cut(C.u3_lr, [-1e9, 200, 500, 1000, 1e9], labels=['<200', '200–500', '500–1,000', '>1,000'])

import os
os.makedirs("fig3", exist_ok=True)
C.to_csv("fig3/gene_level_eval.tsv", sep="\t")
