"""Session-derived computation of funnel.tsv.

Recovered from the record of the interactive analysis session in which funnel.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: gene_boundaries_v3_all_genes.tsv, gene_transcript_end_catalogue.tsv, gene_vs_readend_v48.tsv
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import result  # paths: config/paths.yaml
import pandas as pd
import numpy as np
import sys
sys.path.insert(0, "fig1")
import tracks

# Load transcript catalogue
tsc = pd.read_csv(result('gene_transcript_end_catalogue.tsv'), sep="\t", low_memory=False)

# Load GFF and get CDS coordinates
g = tracks.load_gff()
cds = g[g.type == 'CDS'].copy()
cds['gene2'] = cds.attr.str.extract(r'Parent=([^;.]+)')[0]
env = cds.groupby('gene2').agg(cds_lo=('start', 'min'), cds_hi=('end', 'max'))
env['cds_lo'] -= 1

T = tsc.join(env, on='locus')
mr = T.unit_class == 'mRNA'
T['full_cds'] = mr & T.cds_lo.notna() & (T.start <= T.cds_lo) & (T.end >= T.cds_hi)

S3 = ['ring', 'troph', 'schizont']
det = pd.concat([(T[f'detect_{s}'] == 'detected') for s in S3], axis=1).any(axis=1)
st = pd.concat([(T[f'detect_{s}'] == 'detected') & T.struct_ok & T[f'sj_{s}'] for s in S3], axis=1).any(axis=1)
c3 = pd.concat([(T[f'detect_{s}'] == 'detected') & T.struct_ok & T[f'sj_{s}'] & T.tes_ok for s in S3], axis=1).any(axis=1)

steps = [
    ("HC models", pd.Series(True, index=T.index)),
    ("mRNA (protein-coding locus)", mr),
    ("Detected in ≥1 stage", mr & det),
    ("+ structure & splice junctions", mr & st),
    ("+ 3′ end QC", mr & c3),
    ("+ complete CDS", T.full_cds & c3),
]
fun = pd.DataFrame(
    [(n, int(m.sum()), T.loc[m, 'locus'].nunique()) for n, m in steps],
    columns=['step', 'models', 'loci']
)

nucF = T.full_cds & c3 & ~T.chrom.str.contains('API')
fun2 = fun.copy()
fun2.loc[5, ['step', 'models', 'loci']] = [
    "+ complete CDS (nuclear)",
    int(nucF.sum()),
    T.loc[nucF, 'locus'].nunique()
]
fun2['step'] = [
    "High-confidence models",
    "mRNA loci",
    "Detected in ≥1 stage",
    "+ structure & splice junctions",
    "+ 3′ end QC",
    "+ complete CDS"
]

import os
os.makedirs("fig2", exist_ok=True)
fun2.to_csv("fig2/funnel.tsv", sep="\t", index=False)
