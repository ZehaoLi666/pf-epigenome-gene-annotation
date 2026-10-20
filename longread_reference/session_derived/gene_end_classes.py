"""Session-derived computation of gene_end_classes.tsv.

Recovered from the record of the interactive analysis session in which gene_end_classes.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: gene_boundaries_v3_all_genes.tsv, gene_vs_readend_v48.tsv, siegel2014_pas_per_gene.tsv, gene_transcript_end_catalogue.tsv
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import result  # paths: config/paths.yaml
import pandas as pd
import numpy as np
import os

df = pd.read_csv(result('gene_boundaries_v3_all_genes.tsv'), sep="\t")

Cg = pd.read_csv(result('gene_transcript_end_catalogue.tsv'), sep="\t")

code = lambda s: s.astype(str).str.split(' ').str[0]
Cg['c5'], Cg['c3'] = code(Cg.p5_class), code(Cg.p3_class)

ok = ['3', '4', '5', '7']
fix = ['2a', '2c', '2d', '2e', '6', '8']
ins = ['0', '1']
st = lambda c: np.select([c.isin(ok), c.isin(fix), c.isin(ins)], ['supported', 'needs correction', 'insufficient'], 'other')
Cg['s5'], Cg['s3'] = st(Cg['c5']), st(Cg['c3'])
Cg['endclass'] = np.where((Cg.s5 == 'supported') & (Cg.s3 == 'supported'), 'both supported',
                 np.where((Cg.s5 == 'needs correction') | (Cg.s3 == 'needs correction'), '≥1 end needs correction',
                 np.where((Cg.s5 == 'insufficient') & (Cg.s3 == 'insufficient'), 'both insufficient', 'other')))

grp = df.set_index('gene_id').lr_group
Cg['lr_group'] = Cg.gene.map(grp)

Cg['LR5'] = (Cg.p5_major_end - Cg.cds5).abs()
Cg['LR3'] = (Cg.p3_major_end - Cg.cds3).abs()

sg = pd.read_csv(result('siegel2014_pas_per_gene.tsv'), sep="\t")

M = Cg.merge(sg[['gene_id', 'd_major', 'n_clusters']], left_on='gene', right_on='gene_id', how='left')

os.makedirs("fig1", exist_ok=True)
M.to_csv("fig1/gene_end_classes.tsv", sep="\t", index=False)
print(M.shape)
