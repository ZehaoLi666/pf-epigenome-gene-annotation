"""Library orientation: fraction of read-1 alignments on the same strand as the overlapping protein-coding gene (v68 GFF)."""
import sys, pysam, collections
gff, bams = sys.argv[1], sys.argv[2:]
G = collections.defaultdict(list)
for l in open(gff):
    if l.startswith('#'): continue
    f = l.split('\t')
    if len(f) > 8 and f[2] == 'protein_coding_gene' and f[0] in ('Pf3D7_01_v3', 'Pf3D7_02_v3', 'Pf3D7_03_v3'): G[f[0]].append((int(f[3]) - 1, int(f[4]), f[6]))
for b in bams:
    bam = pysam.AlignmentFile(b); same = opp = 0; n1 = 0
    for c, gs in G.items():
        for s, e, st in gs:
            for r in bam.fetch(c, s + 200, e - 200) if e - s > 600 else []:
                if r.is_secondary or r.is_supplementary or (r.is_paired and not r.is_read1): continue
                rs = '-' if r.is_reverse else '+'; n1 += 1
                if rs == st: same += 1
                else: opp += 1
            if n1 > 200000: break
    print(b.split('/')[-1], n1, round(same / max(same + opp, 1), 4), flush=True)
