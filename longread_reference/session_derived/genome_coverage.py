"""Session-derived computation of genome_coverage.tsv.

Recovered from the record of the interactive analysis session in which genome_coverage.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: gene_end_classes.tsv, gene_boundaries_v3_all_genes.tsv, gene_transcript_end_catalogue.tsv, gene_vs_readend_v48.tsv, transcript_stage_confirmation.tsv.gz, all_cds.bed, all_exon.bed, all_gene.bed, Shaw2022_PacBio.cov1.bed, Shaw2022_PacBio.cov5.bed, Shaw2022_ONT.cov1.bed, Shaw2022_ONT.cov5.bed, Lee2021.cov1.bed, Lee2021.cov5.bed, Hoshizaki2022.cov1.bed, Hoshizaki2022.cov5.bed, Yang2021.cov1.bed, Yang2021.cov5.bed, ALL.cov1.bed, ALL.cov5.bed
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import LAB_DIR, result  # paths: config/paths.yaml
import numpy as np
import pandas as pd

fai = pd.read_csv(f"{LAB_DIR}/IGV_longread_transcript_models/PlasmoDB-68_Pfalciparum3D7_Genome.fasta.fai",
                  sep="\t", header=None, usecols=[0, 1], names=['chrom', 'len'])
fai = fai[~fai.chrom.str.contains('API|M76611|MIT')]

nuc = [f"Pf3D7_{i:02d}_v3" for i in range(1, 15)]
L = {r[0]: r[1] for r in fai.itertuples(index=False) if r[0] in nuc}

def mask(bedfile):
    m = {c: np.zeros(L[c], bool) for c in nuc}
    d = pd.read_csv(bedfile, sep="\t", header=None, usecols=[0, 1, 2])
    for c, s0, e0 in d.itertuples(index=False):
        if c in m:
            m[c][s0:e0] = True
    return m

mC = mask(result('all_cds.bed'))
mE = mask(result('all_exon.bed'))
mG = mask(result('all_gene.bed'))

cat = {c: np.where(mC[c], 0, np.where(mE[c], 1, np.where(mG[c], 2, 3))) for c in nuc}
names = ['CDS', 'UTR / ncRNA exon', 'intron', 'intergenic']
tot = np.array([sum((cat[c] == k).sum() for c in nuc) for k in range(4)])
G = tot.sum()

S = ['Shaw2022_PacBio', 'Shaw2022_ONT', 'Lee2021', 'Hoshizaki2022', 'Yang2021', 'ALL']

cov_files = {
    ('Shaw2022_PacBio', 1): result('Shaw2022_PacBio.cov1.bed'),
    ('Shaw2022_PacBio', 5): result('Shaw2022_PacBio.cov5.bed'),
    ('Shaw2022_ONT', 1): result('Shaw2022_ONT.cov1.bed'),
    ('Shaw2022_ONT', 5): result('Shaw2022_ONT.cov5.bed'),
    ('Lee2021', 1): result('Lee2021.cov1.bed'),
    ('Lee2021', 5): result('Lee2021.cov5.bed'),
    ('Hoshizaki2022', 1): result('Hoshizaki2022.cov1.bed'),
    ('Hoshizaki2022', 5): result('Hoshizaki2022.cov5.bed'),
    ('Yang2021', 1): result('Yang2021.cov1.bed'),
    ('Yang2021', 5): result('Yang2021.cov5.bed'),
    ('ALL', 1): result('ALL.cov1.bed'),
    ('ALL', 5): result('ALL.cov5.bed'),
}

rows = []
for s in S:
    for thr in [1, 5]:
        f = cov_files[(s, thr)]
        m = mask(f)
        cv = np.array([sum((m[c] & (cat[c] == k)).sum() for c in nuc) for k in range(4)])
        rows.append(dict(source=s, min_depth=thr, genome=cv.sum()/G,
                         **{n: cv[k]/tot[k] for k, n in enumerate(names)}))

gcov = pd.DataFrame(rows)
gcov.round(4).to_csv("genome_coverage.tsv", sep="\t", index=False)
