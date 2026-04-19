#!/usr/bin/env python3

import os
import pandas as pd
import shutil   # <-- add this import
import subprocess
import shlex
FEATURE_FILES = {
    #"ATAC-seq": "/rhome/zli529/lab/SRA_toolkit/ATAC-seq_datasets/Toenhake_2018/substract_readcounts/by_antibody/ATAC-seq.txt",
    #"H3R17me2": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Myriam_2023_Chip/coverage/substract_readcounts/by_antibody/H3R17me2.txt",
    #"H2A.Zac": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Myriam_2023_Chip/coverage/substract_readcounts/by_antibody/H2A.Zac.txt",
    #"H3K4me3": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Ashley_2023/coverage/substract_readcounts/by_antibody/H3K4me3.txt",
    #"H3K4me2": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Ashley_2023/coverage/substract_readcounts/by_antibody/H3K4me2.txt",
    #"H3K4me": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Ashley_2023/coverage/substract_readcounts/by_antibody/H3K4me.txt",
    #"H3K18me": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Ashley_2023/coverage/substract_readcounts/by_antibody/H3K18me.txt",
    #"H3K27me": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Ashley_2023/coverage/substract_readcounts/by_antibody/H3K27me.txt",
    #"H3K27ac": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Jingyi_Chip/coverage/substract_readcounts/by_antibody/H3K27ac.txt",
    #"H3K18ac": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Jingyi_Chip/coverage/substract_readcounts/by_antibody/H3K18ac.txt",
    #"H3": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Jingyi_Chip/coverage/substract_readcounts/by_antibody/H3.txt",
    #"H3K4me1": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Jingyi_Chip/coverage/substract_readcounts/by_antibody/H3K4me1.txt",
    #"H2A.Z": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Jingyi_Chip/coverage/substract_readcounts/by_antibody/H2A.Z.txt",
    #"H2B.Z": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Hoeijmakers_2013_Chip/coverage/substract_readcounts/by_antibody/H2B.Z.txt",
    #"MNase": "/rhome/zli529/lab/chip-seq/chipseq/MNase-seq_results/merge_bam/MNseq_sub.txt",
    "H3K9ac": "/rhome/zli529/lab/chip-seq/chipseq/H3K9ac_results/merge_bam/H3K9ac.sub.txt",
    "H3K36me3": "/rhome/zli529/lab/chip-seq/chipseq/H3K36me3_results/H3K36me3.Nor.txt",
}

FAI = "/rhome/zli529/lab/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta.fai"
OUTDIR = "/rhome/zli529/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw"

os.makedirs(OUTDIR, exist_ok=True)
chrom_sizes = os.path.join(OUTDIR, "genome.chrom.sizes")

# Build chrom.sizes
with open(FAI) as fin, open(chrom_sizes, "w") as fout:
    for line in fin:
        toks = line.rstrip("\n").split("\t")
        if len(toks) >= 2:
            fout.write(toks[0] + "\t" + toks[1] + "\n")

def run(cmd):
    print("[cmd]", cmd)
    subprocess.run(cmd, shell=True, check=True)

# Check UCSC tool
if not shutil.which("bedGraphToBigWig"):
    raise SystemExit("ERROR: bedGraphToBigWig not found in PATH. Install: conda install -c bioconda ucsc-bedgraphtobigwig")

for name, txt in FEATURE_FILES.items():
    if not os.path.exists(txt) or os.stat(txt).st_size == 0:
        print(f"[warn] Skipping {name} (missing/empty): {txt}")
        continue
    bg = os.path.join(OUTDIR, f"{name}.bedgraph")
    bgs = os.path.join(OUTDIR, f"{name}.sorted.bedgraph")
    bw = os.path.join(OUTDIR, f"{name}.bw")

    # txt -> bedgraph (skip header if any)
    # awk: col1=chr, col2=pos (1-based), col3=value
    run(f"awk 'BEGIN{{OFS=\"\\t\"}} NR==1{{if($2 !~ /^[0-9]+$/){{next}}}} {{print $1, $2-1, $2, $3}}' {shlex.quote(txt)} > {shlex.quote(bg)}")

    # sort & convert
    run(f"sort -k1,1 -k2,2n {shlex.quote(bg)} -o {shlex.quote(bgs)}")
    run(f"bedGraphToBigWig {shlex.quote(bgs)} {shlex.quote(chrom_sizes)} {shlex.quote(bw)}")

    # cleanup (optional)
    os.remove(bg)
    os.remove(bgs)

    print(f"[✓] Wrote {bw}")

print(f"All done. Output dir: {OUTDIR}")
