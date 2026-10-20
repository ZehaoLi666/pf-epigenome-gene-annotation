"""Strand-specific evidence from unstranded RNA-seq BAMs, 20-bp bins, 14 nuclear chromosomes.
spliced reads: strand from intron motif (GT..AG = +, CT..AC = -); records exon-block coverage and intron-span coverage per strand.
poly(A) reads: soft clip >=6 nt, >=80% A at alignment right end (+ strand TTS at reference_end) or T at left end (- strand TTS at reference_start);
dropped if the genome itself is A-rich (>=7 of 10) immediately downstream (internal priming)."""
from pfconfig import HPCC_BIGDATA, HPCC_LAB  # paths: config/paths.yaml
import sys, numpy as np, pysam
stage = sys.argv[1]; RES = 20
BAM = f"{HPCC_LAB}/SRA_toolkit/RNAseq_7stages/processed_data/{stage}/{stage}_sorted.bam"
FA = f"{HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta"
fa = pysam.FastaFile(FA); bam = pysam.AlignmentFile(BAM)
chroms = [c for c in bam.references if c.startswith("Pf3D7_") and c[6:8].isdigit()]
out = {}; stats = dict(reads=0, spliced=0, motif_plus=0, motif_minus=0, motif_none=0, polyA_plus=0, polyA_minus=0, polyA_internal=0)
for c in chroms:
    L = bam.get_reference_length(c); nb = L // RES; seq = fa.fetch(c).upper()
    ex = np.zeros((2, L + 1), np.float32); it = np.zeros((2, L + 1), np.float32); pa = np.zeros((2, nb + 1), np.float32)
    for r in bam.fetch(c):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 10: continue
        stats["reads"] += 1; cig = r.cigartuples
        if any(op == 3 for op, _ in cig):
            stats["spliced"] += 1; pos = r.reference_start; blocks, introns = [], []
            for op, n in cig:
                if op in (0, 7, 8, 2): blocks.append((pos, pos + n)); pos += n
                elif op == 3: introns.append((pos, pos + n)); pos += n
            votes = 0
            for a, b in introns:
                m = seq[a:a + 2] + seq[b - 2:b]
                votes += 1 if m == "GTAG" else (-1 if m == "CTAC" else 0)
            if any(b - a > 3000 or b - a < 20 for a, b in introns): stats["motif_none"] += 1; continue
            if votes == 0: stats["motif_none"] += 1; continue
            s = 0 if votes > 0 else 1; stats["motif_plus" if s == 0 else "motif_minus"] += 1
            for a, b in blocks: ex[s, a] += 1; ex[s, b] -= 1
            for a, b in introns: it[s, a] += 1; it[s, b] -= 1
        q = r.query_sequence
        if q is None: continue
        if cig[-1][0] == 4 and cig[-1][1] >= 8:
            t = q[-cig[-1][1]:]
            if t.count("A") >= 0.9 * len(t):
                e = r.reference_end; d = seq[e:e + 10]
                if d.count("A") >= 7: stats["polyA_internal"] += 1
                elif e // RES <= nb: pa[0, e // RES] += 1; stats["polyA_plus"] += 1
        if cig[0][0] == 4 and cig[0][1] >= 8:
            t = q[:cig[0][1]]
            if t.count("T") >= 0.9 * len(t):
                st = r.reference_start; d = seq[max(0, st - 10):st]
                if d.count("T") >= 7: stats["polyA_internal"] += 1
                else: pa[1, st // RES] += 1; stats["polyA_minus"] += 1
    ex = np.cumsum(ex, 1)[:, :nb * RES].reshape(2, nb, RES).mean(2); it = np.cumsum(it, 1)[:, :nb * RES].reshape(2, nb, RES).mean(2)
    out[f"{c}__exon"] = ex.astype(np.float16); out[f"{c}__intron"] = it.astype(np.float16); out[f"{c}__polyA"] = pa[:, :nb].astype(np.float16)
    print(stage, c, {k: v for k, v in stats.items()}, flush=True)
np.savez_compressed(f"{HPCC_BIGDATA}/boundary_seg/stranded_v2_{stage}.npz", **out)
print("DONE", stage, stats, flush=True)
