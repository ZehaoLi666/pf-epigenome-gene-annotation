#!/usr/bin/env python3
import sys
from pyfaidx import Fasta

# Usage: python bed_to_fasta.py genome.fa novel_units.sorted.bed out.fa
genome_fa, bed, out_fa = sys.argv[1:4]
fa = Fasta(genome_fa, as_raw=True)  # raw strings for speed

def revcomp(s):
    comp = str.maketrans("ACGTacgtnN", "TGCAtgcanN")
    return s.translate(comp)[::-1]

with open(bed) as fin, open(out_fa, "w") as fout:
    for line in fin:
        if not line.strip() or line.startswith("#"):
            continue
        chrom, start, end, name, score, strand = line.rstrip("\n").split("\t")[:6]
        start = int(start); end = int(end)
        seq = fa[chrom][start:end]
        if strand == "-":
            seq = revcomp(seq)
        fout.write(f">{name}\n")
        # wrap to 60 chars per line
        for i in range(0, len(seq), 60):
            fout.write(seq[i:i+60] + "\n")
