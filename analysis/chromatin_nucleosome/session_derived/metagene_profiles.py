"""Session-derived computation of metagene_profiles.npz.

Recovered from the record of the interactive analysis session in which metagene_profiles.npz was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: gene_boundaries_v3_all_genes.tsv, gene_end_classes.tsv, gene_transcript_end_catalogue.tsv, gene_vs_readend_v48.tsv, transcript_stage_confirmation.tsv.gz
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import LAB_DIR, result  # paths: config/paths.yaml
import numpy as np
import pandas as pd
import pyBigWig
import glob

TR = f"{LAB_DIR}/lncRNA project tracks"
IGV = f"{LAB_DIR}/IGV_longread_transcript_models"

ST = ['Ring', 'Trophozoite', 'Schizont']

daft = {
    '+': f"{IGV}/DAFTseq_pool_plus_CPM_10bp.bw",
    '-': f"{IGV}/DAFTseq_pool_minus_CPM_10bp.bw"
}

RB = pd.read_csv(result('gene_boundaries_v3_all_genes.tsv'), sep="\t")

R = RB[(RB.lr_group == 'confirmed') & RB.lr_TSS.notna() & RB.lr_TTS.notna()].copy()
R = R[(R.lr_TTS - R.lr_TSS).abs() >= 400]

def readfa(p):
    seq = {}; name = None; buf = []
    for line in open(p):
        if line.startswith('>'):
            if name: seq[name] = ''.join(buf)
            name = line[1:].split()[0]; buf = []
        else: buf.append(line.strip().upper())
    seq[name] = ''.join(buf)
    return seq

FA = readfa(f"{IGV}/PlasmoDB-68_Pfalciparum3D7_Genome.fasta")
CL = {c: len(FA[c]) for c in FA}

NU, NB, ND, FLK = 20, 40, 20, 1000

TRK = [
    ('H2A.Z', f"{TR}/histone_marks/H2AZ_{{s}}_minusInput.bw"),
    ('H2B.Z', f"{TR}/histone_marks/H2BZ_{{s}}_minusInput.bw"),
    ('H3K9ac', f"{TR}/histone_marks/H3K9ac_{{s}}_minusInput.bw"),
    ('H3K27ac', f"{TR}/histone_marks/H3K27ac_{{s}}_minusInput.bw"),
    ('H3K4me3', f"{TR}/histone_marks/H3K4me3_{{s}}_minusInput.bw"),
    ('ATAC-seq', f"{TR}/accessibility_nucleosome/ATAC_{{s}}_minusInput.bw"),
    ('MNase-seq', f"{TR}/accessibility_nucleosome/MNase_{{s}}_CPM.bw"),
    ('GRO-seq', f"{TR}/nascent_transcription/GRO_{{s}}_mean_normalized.bw"),
    ('Polysome profiling', f"{TR}/translatome/Polysome_{{s}}_*_CPM_final.bw")
]

def gmean(bw):
    h = bw.header()
    return h['sumData'] / sum(bw.chroms().values())

def prof(bw, rows, strand_sel=None):
    out = np.zeros((len(rows), NU + NB + ND))
    for i, r in enumerate(rows.itertuples()):
        c = r.chrom; L = CL[c]; t, e = int(r.lr_TSS), int(r.lr_TTS)
        if r.strand == '+':
            segs = [(t - FLK, t, NU), (t, e, NB), (e, e + FLK, ND)]
        else:
            segs = [(t, t + FLK, NU), (e, t, NB), (e - FLK, e, ND)]
        vals = []
        for a, b, n in segs:
            a2, b2 = max(0, a), min(L, b)
            v = np.array([x if x is not None else 0.0 for x in bw.stats(c, a2, b2, nBins=n, type='mean')]) if b2 - a2 >= n else np.zeros(n)
            vals.append(v if r.strand == '+' else v[::-1])
        out[i] = np.concatenate(vals)
    return out

PR = {}
for nm, pat in TRK:
    acc = []
    for s in ST:
        p = glob.glob(pat.format(s=s))[0]
        bw = pyBigWig.open(p); m = gmean(bw); acc.append(prof(bw, R) / m); bw.close()
    PR[nm] = np.mean(acc, 0)

pol = sorted(glob.glob(f"{TR}/machinery/PolII_*_minusInput.bw"))
acc = []
for p in pol:
    bw = pyBigWig.open(p); m = gmean(bw); acc.append(prof(bw, R) / m); bw.close()
PR['Pol II Ser2P'] = np.mean(acc, 0)

bwp = pyBigWig.open(daft['+']); bwm = pyBigWig.open(daft['-'])
mp = (gmean(bwp) + gmean(bwm)) / 2
Pp = prof(bwp, R) / mp; Pm = prof(bwm, R) / mp
plus = (R.strand == '+').values[:, None]
PR['DAFT-seq, gene strand'] = np.where(plus, Pp, Pm)
PR['DAFT-seq, opposite strand'] = np.where(plus, Pm, Pp)

np.savez_compressed("metagene_profiles.npz", **{k.replace(' ', '_'): v for k, v in PR.items()})
