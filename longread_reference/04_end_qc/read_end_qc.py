"""Per reference mRNA end: read-level 3' tail and 5' extra-G evidence from long-read BAMs.
ends file: chrom pos strand end_type (0-based genomic position of the end base). Reads assigned when their 5'/3' end lies within +-W bp
of an end on the same strand (alignment strand, reads assumed sense-oriented; orientation checked via tail side).
3': tail = soft clip on the 3' side whose first 15 nt (adjacent to the alignment) have >=8 nt and >=80% A (sense). Also reads ending
    with aligned genomic A-run (last 10 aligned nt >=8 A) and no clip tail = 'templated A end'.
5': extra G = number of consecutive G at the 3' end of the 5' soft clip (sense, adjacent to alignment), capped at 6; also whether a
    5' clip exists (adapter/TSO retained)."""
from pfconfig import HPCC_BIGDATA, HPCC_LAB  # paths: config/paths.yaml
import sys, collections, bisect, numpy as np, pysam
ENDS, OUT = sys.argv[1], sys.argv[2]; W = 30
FA = pysam.FastaFile(f"{HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta")
B = f"{HPCC_BIGDATA}/GSE186109_longread/bam"
BAMS = [(f"{B}/{s}/{t}.sorted.bam", "pacbio", s) for s in ("ring_total", "ring_nuclear", "schizont_total", "schizont_nuclear", "troph_total", "troph_nuclear", "schiztroph_cytoplasmic") for t in ("hq", "lq")]
BAMS += [(f"{B}/ont_mixed_oriented/oriented.sorted.bam", "ont", "ont_mixed")]
import os; BAMS = [b for b in BAMS if os.path.exists(b[0])]
comp = str.maketrans("ACGTN", "TGCAN"); rc = lambda s: s.translate(comp)[::-1]
ends = collections.defaultdict(list)
for i, line in enumerate(open(ENDS)):
    if i == 0: continue
    c, p, s, e = line.rstrip("\n").split("\t"); ends[(c, s, e)].append(int(p))
idx = {k: np.array(sorted(v)) for k, v in ends.items()}
stat = collections.defaultdict(lambda: collections.Counter())
def near(key, pos):
    a = idx.get(key)
    if a is None or not len(a): return None
    j = bisect.bisect_left(a, pos - W)
    best = None
    while j < len(a) and a[j] <= pos + W:
        if best is None or abs(a[j] - pos) < abs(best - pos): best = a[j]
        j += 1
    return best
glob = collections.Counter()
for path, plat, samp in BAMS:
    bam = pysam.AlignmentFile(path)
    for r in bam.fetch():
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.query_sequence is None: continue
        c = r.reference_name; s = "-" if r.is_reverse else "+"; q = r.query_sequence; cig = r.cigartuples
        lc = cig[0][1] if cig[0][0] == 4 else 0; rcl = cig[-1][1] if cig[-1][0] == 4 else 0
        if s == "+":
            p5, p3 = r.reference_start, r.reference_end - 1
            clip5 = q[:lc]; clip3 = q[len(q) - rcl:] if rcl else ""; al3 = q[len(q) - rcl - 10:len(q) - rcl]
        else:
            p5, p3 = r.reference_end - 1, r.reference_start
            clip5 = rc(q[len(q) - rcl:]) if rcl else ""; clip3 = rc(q[:lc]); al3 = rc(q[lc:lc + 10])
        # opposite-side tail (orientation check)
        opp = rc(q[:lc]) if s == "+" else q[len(q) - rcl:] if rcl else ""
        t3 = clip3[:15]; tail = len(t3) >= 8 and t3.count("A") >= 0.8 * len(t3)
        to = opp[:15]; tail_opp = len(to) >= 8 and to.count("T") >= 0.8 * len(to)
        glob[(plat, "tail_sense")] += tail; glob[(plat, "tail_antisense_side")] += tail_opp; glob[(plat, "reads")] += 1
        e3 = near((c, s, "TTS"), p3)
        if e3 is not None:
            k = stat[(c, s, "TTS", int(e3))]; k[f"{plat}_n"] += 1; k[f"{plat}_tail"] += tail
            k[f"{plat}_templatedA_notail"] += (not tail) and al3.count("A") >= 8
        e5 = near((c, s, "TSS"), p5)
        if e5 is not None:
            k = stat[(c, s, "TSS", int(e5))]; k[f"{plat}_n"] += 1; k[f"{plat}_clip5"] += lc > 0 if s == "+" else rcl > 0
            g = len(clip5) - len(clip5.rstrip("G")); g = min(g, 6)
            k[f"{plat}_G1"] += g >= 1; k[f"{plat}_G3"] += g >= 3; k[f"{plat}_Gsum"] += g
    print(path, dict(glob), flush=True)
cols = ["pacbio_n", "pacbio_tail", "pacbio_templatedA_notail", "ont_n", "ont_tail", "ont_templatedA_notail", "pacbio_clip5", "pacbio_G1", "pacbio_G3", "pacbio_Gsum", "ont_clip5", "ont_G1", "ont_G3", "ont_Gsum"]
with open(OUT, "w") as fh:
    fh.write("chrom\tstrand\tend_type\tpos\t" + "\t".join(cols) + "\n")
    for (c, s, e, p), k in sorted(stat.items()):
        fh.write(f"{c}\t{s}\t{e}\t{p}\t" + "\t".join(str(k[x]) for x in cols) + "\n")
with open(OUT + ".global.txt", "w") as fh: fh.write(repr(dict(glob)) + "\n")
print("DONE", len(stat))
