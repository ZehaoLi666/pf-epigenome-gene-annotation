"""Session-derived computation of v71_vs_reference_ends.tsv.

Recovered from the record of the interactive analysis session in which v71_vs_reference_ends.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: gene_boundaries_v3_all_genes.tsv, siegel2014_pas_per_gene.tsv, v3_model_code_weights.tar.gz, gene_transcript_end_catalogue.tsv, gene_vs_readend_v48.tsv
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import result  # paths: config/paths.yaml
import pandas as pd
import numpy as np
import os

# Load confirmed gene boundaries
B = pd.read_csv(result('gene_boundaries_v3_all_genes.tsv'), sep="\t")
R = B[B.lr_group == 'confirmed'].copy()

# Extract tar to get ml3_end_positives.tsv
os.makedirs("v3tar", exist_ok=True)
import subprocess
subprocess.run(["tar", "-xzf", result('v3_model_code_weights.tar.gz'),
                "work/data/ml3_end_positives.tsv"], cwd="v3tar")

E = pd.read_csv("v3tar/work/data/ml3_end_positives.tsv", sep="\t")

import sys
sys.path.insert(0, "fig1")
import tracks

LR = tracks.read_bed12(f"{tracks.IGV}/confirmed_transcripts_by_category.bed")
LR['locus'] = LR.name.str.split('|').str[2]
LR['cat'] = LR.name.str.split('|').str[1]
LNC = tracks.read_bed12(f"{tracks.IGV}/lncRNA_catalogue_3datasets.bed")

ref = set(R.gene_id)

Rg = R.set_index('gene_id')

Rg['L5a'] = [(r.atg - r.annot_tss) * (1 if r.strand == '+' else -1) for _, r in Rg.iterrows()]
Rg['L3a'] = [(r.annot_tts - r.stop) * (1 if r.strand == '+' else -1) for _, r in Rg.iterrows()]

T = E[(E.kind == 'terminal') & E.locus.isin(ref)]

def dist_from_cds(g, pos, which):
    r = Rg.loc[g]
    s = 1 if r.strand == '+' else -1
    return ((r.atg - pos) * s) if which == '5' else ((pos - r.stop) * s)

T = T.assign(d=[dist_from_cds(g, p, '5' if t == 'TSS' else '3')
                for g, p, t in zip(T.locus, T.pos, T.end_type)])

D5 = T[T.end_type == 'TSS'].groupby('locus').d.apply(list)
D3 = T[T.end_type == 'TTS'].groupby('locus').d.apply(list)

def cls(L, D, tol):
    if D is None:
        return 'not assessable'
    D = np.array(D)
    if L == 0:
        return 'no annotated UTR' if D.max() > tol else 'agree'
    if np.min(np.abs(L - D)) <= tol:
        return 'agree'
    if L > D.max() + tol:
        return 'over-long'
    if L < D.min() - tol:
        return 'too short'
    return 'between ends'

out = {}
for tol in (100, 200):
    c5 = pd.Series({g: cls(Rg.loc[g, 'L5a'], D5.get(g), tol) for g in Rg.index})
    c3 = pd.Series({g: cls(Rg.loc[g, 'L3a'], D3.get(g), tol) for g in Rg.index})
    out[tol] = (c5, c3)

c5, c3 = out[100]
a5, a3 = out[200]

def excess(g, which):
    r = Rg.loc[g]
    D = np.array((D5 if which == '5' else D3)[g])
    m = D.max()
    if which == '5':
        return (r.annot_tss, r.atg - m) if r.strand == '+' else (r.atg + m, r.annot_tss)
    return (r.stop + m, r.annot_tts) if r.strand == '+' else (r.annot_tts, r.stop - m)

def content(g, which):
    r = Rg.loc[g]
    a, b = sorted(excess(g, which))
    s = r.strand
    lr = LR[(LR.chrom == r.chrom) & (LR.end > a) & (LR.start < b) & (LR.locus != g)]
    ln = LNC[(LNC.chrom == r.chrom) & (LNC.end > a) & (LNC.start < b)]
    if (lr.strand == s).any():
        return 'same-strand long-read unit'
    if (ln.strand == s).any():
        return 'catalogued lncRNA'
    if (lr.strand != s).any() or (ln.strand != s).any():
        return 'antisense unit only'
    return 'no other unit detected'

k5 = pd.Series({g: content(g, '5') for g in c5.index[c5 == 'over-long']})
k3 = pd.Series({g: content(g, '3') for g in c3.index[c3 == 'over-long']})

ex5 = pd.Series({g: abs(np.subtract(*excess(g, '5'))) for g in k5.index})
ex3 = pd.Series({g: abs(np.subtract(*excess(g, '3'))) for g in k3.index})

summ = pd.DataFrame({'c5_100': c5, 'c3_100': c3, 'c5_200': a5, 'c3_200': a3,
                     'excess5_content': k5, 'excess3_content': k3,
                     'excess5_bp': ex5, 'excess3_bp': ex3})
summ['L5_v71'] = Rg.L5a
summ['L3_v71'] = Rg.L3a
summ['ref_TSS_dists'] = D5.reindex(summ.index).apply(
    lambda x: ';'.join(map(str, map(int, x))) if isinstance(x, list) else '')
summ['ref_TTS_dists'] = D3.reindex(summ.index).apply(
    lambda x: ';'.join(map(str, map(int, x))) if isinstance(x, list) else '')
summ.index.name = 'gene_id'
summ.to_csv("v71_vs_reference_ends.tsv", sep="\t")
