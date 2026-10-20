"""Session-derived computation of model_completeness.tsv.

Recovered from the record of the interactive analysis session in which model_completeness.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: v3_model_code_weights.tar.gz, transcript_stage_confirmation.tsv.gz, gene_boundaries_v3_all_genes.tsv, gene_transcript_end_catalogue.tsv, gene_vs_readend_v48.tsv, tss_read_evidence_by_group.tsv
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import LAB_DIR, result  # paths: config/paths.yaml
import numpy as np
import pandas as pd
import sys
import os

# Load transcript stage confirmation data
tsc = pd.read_csv(result('transcript_stage_confirmation.tsv.gz'), sep="\t", low_memory=False)

# Load GFF to get CDS coordinates
import subprocess
sys.path.insert(0, "fig1")

# GFF loading as in annotation/tracks.py (inlined)
LAB = LAB_DIR
GFF = f"{LAB}/PlasmoDB/PlasmoDB-68/PlasmoDB-68_Pfalciparum3D7.gff"

g = pd.read_csv(GFF, sep="\t", comment="#", header=None,
                names="chrom src type start end score strand phase attr".split(),
                dtype={'start': int, 'end': int})
g = g[g.type.isin(['protein_coding_gene', 'ncRNA_gene', 'mRNA', 'CDS',
                   'five_prime_UTR', 'three_prime_UTR', 'exon', 'gene'])].copy()
g['ID'] = g.attr.str.extract(r'ID=([^;]+)')
g['Parent'] = g.attr.str.extract(r'Parent=([^;]+)')

# Get CDS coordinates per gene
cds = g[g.type == 'CDS'].copy()
cds['gene2'] = cds.attr.str.extract(r'Parent=([^;.]+)')[0]
env = cds.groupby('gene2').agg(cds_lo=('start', 'min'), cds_hi=('end', 'max'))
env['cds_lo'] -= 1

S3 = ['ring', 'troph', 'schizont']

T = tsc.join(env, on='locus')
mr = T.unit_class == 'mRNA'
T['full_cds'] = mr & T.cds_lo.notna() & (T.start <= T.cds_lo) & (T.end >= T.cds_hi)

det = pd.concat([(T[f'detect_{s}'] == 'detected') for s in S3], axis=1).any(axis=1)
st = pd.concat([(T[f'detect_{s}'] == 'detected') & T.struct_ok & T[f'sj_{s}'] for s in S3], axis=1).any(axis=1)
c3 = pd.concat([(T[f'detect_{s}'] == 'detected') & T.struct_ok & T[f'sj_{s}'] & T.tes_ok for s in S3], axis=1).any(axis=1)

# Build model completeness for confirmed mRNA models
cm = T[mr & c3 & T.cds_lo.notna() & ~T.chrom.str.contains('API')].copy()
cm['trunc5'] = ~cm.full_cds & np.where(cm.strand == '+', cm.start > cm.cds_lo, cm.end < cm.cds_hi)
cm['cls'] = np.select([cm.full_cds, cm.trunc5], ['full CDS', '5′ truncated'], '3′ truncated')

cmc = cm.cls.value_counts().rename_axis('cls').reset_index(name='n')
cmc['loci'] = cm.locus.nunique()
cmc.to_csv("model_completeness.tsv", sep="\t", index=False)
print(cmc)
