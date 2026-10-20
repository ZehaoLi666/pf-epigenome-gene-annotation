"""Session-derived computation of internal_priming_v71.tsv.

Recovered from the record of the interactive analysis session in which internal_priming_v71.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: gene_boundaries_v3_all_genes.tsv, gene_end_classes.tsv, gene_transcript_end_catalogue.tsv, gene_vs_readend_v48.tsv, transcript_stage_confirmation.tsv.gz
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import LAB_DIR, result  # paths: config/paths.yaml
import numpy as np
import pandas as pd
import glob

# Load gene end classes
Cg = pd.read_csv(result('gene_end_classes.tsv'), sep="\t")
genes = set(Cg.gene)

# Load GFF
g = pd.read_csv(glob.glob(f"{LAB_DIR}/PlasmoDB/PlasmoDB-68/*.gff")[0],
                sep="\t", comment="#", header=None, low_memory=False,
                names=['chrom', 'src', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attr'])

# Get CDS envelope per gene
cds = g[g.type == 'CDS'].copy()
cds['gene'] = cds.attr.str.extract(r'Parent=([^;.]+)')[0]
env = cds.groupby('gene').agg(lo=('start', 'min'), hi=('end', 'max'))
env['lo'] -= 1

# Get mRNA 3' ends from GFF
mr = g[g.type == 'mRNA'].copy()
mr['gene'] = mr.attr.str.extract(r'Parent=([^;]+)')[0]
mr = mr[mr.gene.isin(genes)]
mr['p0'] = np.where(mr.strand == '+', mr.end, mr.start - 1)

# Read genome FASTA
def readfa(p):
    seq = {}; name = None; buf = []
    for line in open(p):
        if line.startswith('>'):
            if name: seq[name] = ''.join(buf)
            name = line[1:].split()[0]; buf = []
        else: buf.append(line.strip().upper())
    seq[name] = ''.join(buf)
    return seq

FA = readfa(f"{LAB_DIR}/IGV_longread_transcript_models/PlasmoDB-68_Pfalciparum3D7_Genome.fasta")
comp = str.maketrans('ACGTN', 'TGCAN')

def winA(ch, p0, strand):
    sq = FA.get(ch)
    if sq is None: return np.nan
    w = sq[p0:p0+20] if strand == '+' else sq[max(0, p0-20):p0].translate(comp)[::-1]
    return w.count('A') / 20 if len(w) == 20 else np.nan

mr['A20'] = [winA(r.chrom, int(r.p0), r.strand) for r in mr.itertuples()]

# Load long-read 3' end table
M = pd.read_csv("cov/master_tes.tsv.gz", sep="\t")
M['A20'] = [winA(r.chrom, int(r.tes), r.strand) for r in M.itertuples()]

# Join with CDS envelope
Mm = M[M.locus.isin(env.index)].join(env, on='locus')
inside = np.where(Mm.strand == '+',
                  (Mm.tes > Mm.lo) & (Mm.tes < Mm.hi),
                  (Mm.tes > Mm.lo) & (Mm.tes < Mm.hi))
down = np.where(Mm.strand == '+',
                Mm.tes >= Mm.hi,
                Mm.tes <= Mm.lo)

# Random genomic sites
fai = pd.read_csv(f"{LAB_DIR}/IGV_longread_transcript_models/PlasmoDB-68_Pfalciparum3D7_Genome.fasta.fai",
                  sep="\t", header=None)
nuc = [f"Pf3D7_{i:02d}_v3" for i in range(1, 15)]
Ls = np.array([len(FA[c]) for c in nuc])
rng = np.random.default_rng(7)
ch = rng.choice(len(nuc), 100000, p=Ls / Ls.sum())
R = pd.DataFrame({
    'chrom': [nuc[i] for i in ch],
    'p0': [rng.integers(100, Ls[i] - 100) for i in ch],
    'strand': rng.choice(['+', '-'], 100000)
})
R['A20'] = [winA(r.chrom, int(r.p0), r.strand) for r in R.itertuples()]

f = lambda x: float((np.asarray(x) >= 0.8).mean())

rates = pd.DataFrame([
    dict(label="Long-read 3′ ends\ninside the CDS", frac=f(Mm.A20[inside]), n=int(inside.sum())),
    dict(label="Long-read 3′ ends\nbeyond the CDS", frac=f(Mm.A20[down]), n=int(down.sum())),
    dict(label="All long-read\n3′ ends", frac=f(M.A20), n=len(M)),
    dict(label="PlasmoDB v71\n3′ ends", frac=f(mr.A20), n=len(mr)),
    dict(label="Random genomic\nsites", frac=f(R.A20.dropna()), n=int(R.A20.notna().sum()))
])

rates.to_csv("cov/internal_priming_v71.tsv", sep="\t", index=False)
