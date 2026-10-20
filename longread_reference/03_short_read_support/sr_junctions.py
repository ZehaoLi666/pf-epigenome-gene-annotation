"""Count spliced short reads per intron; write STAR SJ.out.tab-style files (per stage and pooled)."""
from pfconfig import HPCC_BIGDATA, HPCC_LAB  # paths: config/paths.yaml
import pysam, sys, collections, os
P = f"{HPCC_LAB}/SRA_toolkit/RNAseq_7stages/processed_data"
OUT = f"{HPCC_BIGDATA}/GSE186109_longread/shortread/junctions"; os.makedirs(OUT, exist_ok=True)
GEN = f"{HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta"
fa = pysam.FastaFile(GEN)
MOT = {("GT","AG"):(1,1), ("CT","AC"):(2,2), ("GC","AG"):(1,3), ("CT","GC"):(2,4), ("AT","AC"):(1,5), ("GT","AT"):(2,6)}
stage = sys.argv[1]
bam = pysam.AlignmentFile(f"{P}/{stage}/{stage}_sorted.bam")
uniq = collections.Counter(); multi = collections.Counter(); ovh = collections.Counter()
for r in bam.fetch(until_eof=True):
    if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
    cig = r.cigartuples
    if not cig or not any(op == 3 for op, _ in cig): continue
    nh = r.get_tag("NH") if r.has_tag("NH") else 1
    pos = r.reference_start
    blocks = []  # (op,len,refpos)
    for op, ln in cig:
        blocks.append((op, ln, pos))
        if op in (0, 2, 3, 7, 8): pos += ln
    for i, (op, ln, p0) in enumerate(blocks):
        if op != 3: continue
        left = sum(l for o, l, _ in blocks[:i] if o in (0, 7, 8))
        right = sum(l for o, l, _ in blocks[i+1:] if o in (0, 7, 8))
        # anchor = aligned bases in the adjacent exon block only
        la = 0
        for o, l, _ in reversed(blocks[:i]):
            if o == 3: break
            if o in (0, 7, 8): la += l
        ra = 0
        for o, l, _ in blocks[i+1:]:
            if o == 3: break
            if o in (0, 7, 8): ra += l
        key = (r.reference_name, p0, p0 + ln)
        if nh == 1: uniq[key] += 1
        else: multi[key] += 1
        ovh[key] = max(ovh[key], min(la, ra))
with open(f"{OUT}/{stage}.SJ.out.tab", "w") as fh:
    for key in sorted(set(uniq) | set(multi)):
        ch, s, e = key
        d = fa.fetch(ch, s, s + 2).upper(); a = fa.fetch(ch, e - 2, e).upper()
        st, mo = MOT.get((d, a), (0, 0))
        fh.write(f"{ch}\t{s+1}\t{e}\t{st}\t{mo}\t0\t{uniq[key]}\t{multi[key]}\t{ovh[key]}\n")
print(stage, "junctions", len(set(uniq) | set(multi)), "uniq_reads", sum(uniq.values()))
