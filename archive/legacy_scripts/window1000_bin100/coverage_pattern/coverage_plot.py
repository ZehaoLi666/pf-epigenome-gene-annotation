#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
# --- Input files ---
mRNA_file = "/rhome/zli529/lab/LncRNA_chip_prediction/Final/window1000_bin100/training_features_values/reference_gene_boundaries.csv"
neg_file = "/rhome/zli529/lab/LncRNA_chip_prediction/Final/window1000_bin100/training_features_values/negatives_2k.tsv"

FEATURE_FILES = {
    "ATAC-seq": "/rhome/zli529/lab/SRA_toolkit/ATAC-seq_datasets/Toenhake_2018/substract_readcounts/by_antibody/ATAC-seq.txt"
   # "H3R17me2": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Myriam_2023_Chip/coverage/substract_readcounts/by_antibody/H3R17me2.txt",
   # "H2A.Zac": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Myriam_2023_Chip/coverage/substract_readcounts/by_antibody/H2A.Zac.txt",
   # "H3K4me3": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Ashley_2023/coverage/substract_readcounts/by_antibody/H3K4me3.txt",
   # "H3K4me2": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Ashley_2023/coverage/substract_readcounts/by_antibody/H3K4me2.txt",
   # "H3K4me": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Ashley_2023/coverage/substract_readcounts/by_antibody/H3K4me.txt",
   # "H3K18me": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Ashley_2023/coverage/substract_readcounts/by_antibody/H3K18me.txt",
   # "H3K27me": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Ashley_2023/coverage/substract_readcounts/by_antibody/H3K27me.txt",
   # "H3K27ac": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Jingyi_Chip/coverage/substract_readcounts/by_antibody/H3K27ac.txt",
   # "H3K18ac": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Jingyi_Chip/coverage/substract_readcounts/by_antibody/H3K18ac.txt",
   # "H3": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Jingyi_Chip/coverage/substract_readcounts/by_antibody/H3.txt",
   # "H3K4me1": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Jingyi_Chip/coverage/substract_readcounts/by_antibody/H3K4me1.txt",
   # "H2A.Z": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Jingyi_Chip/coverage/substract_readcounts/by_antibody/H2A.Z.txt",
   # "H2B.Z": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Hoeijmakers_2013_Chip/coverage/substract_readcounts/by_antibody/H2B.Z.txt",
   # "MNase": "/rhome/zli529/lab/chip-seq/chipseq/MNase-seq_results/merge_bam/MNseq_sub.txt",
   # "H3K9ac": "/rhome/zli529/lab/chip-seq/chipseq/H3K9ac_results/merge_bam/H3K9ac.sub.txt",
   # "H3K36me3": "/rhome/zli529/lab/chip-seq/chipseq/H3K36me3_results/H3K36me3.Nor.txt"
}

# Parameters
UP = 1000
DOWN = 1000
BINS_BODY = 1000   # normalized gene body bins

# --- Load gene / neg files ---
mRNA = pd.read_csv(mRNA_file, sep="\t", header=None,
                   names=["chr","start","end","gene_id","strand"])
neg = pd.read_csv(neg_file, sep="\t", header=None,
                  names=["chr","start","end","id","strand"])

# --------------------
# Helper functions
# --------------------
def build_coverage_dict(read_counts_file):
    reads = pd.read_csv(read_counts_file, sep="\t", header=None,
                        names=["chr","site","read_count"])
    coverage_dict = {}
    for chrom, df in reads.groupby("chr"):
        cov = dict(zip(df["site"], df["read_count"]))
        coverage_dict[chrom] = cov
    return coverage_dict

def get_gene_profile(row, coverage_dict):
    chrom = row["chr"]
    if chrom not in coverage_dict:
        return None

    start, end, strand = row["start"], row["end"], row["strand"]
    tss = start if strand == "+" else end
    tts = end if strand == "+" else start

    region_start = tss - UP if strand == "+" else tts - UP
    region_end   = tts + DOWN if strand == "+" else tss + DOWN

    cov = [coverage_dict[chrom].get(pos, 0) for pos in range(region_start, region_end)]
    cov = np.array(cov)

    if strand == "-":
        cov = cov[::-1]

    promoter = cov[:UP]
    body = cov[UP: len(cov)-DOWN]
    downstream = cov[-DOWN:]

    if len(body) > 1:
        body_norm = np.interp(
            np.linspace(0, len(body)-1, BINS_BODY),
            np.arange(len(body)),
            body
        )
    else:
        body_norm = np.zeros(BINS_BODY)

    profile = np.concatenate([promoter, body_norm, downstream])
    return profile

def get_neg_profile(row, coverage_dict):
    chrom = row["chr"]
    if chrom not in coverage_dict:
        return None

    center = (row["start"] + row["end"]) // 2
    region_start = center - UP
    region_end   = center + BINS_BODY + DOWN   # force same total length

    cov = [coverage_dict[chrom].get(pos, 0) for pos in range(region_start, region_end)]
    return np.array(cov)

def compute_average_profile(df, profile_func, coverage_dict):
    profiles = []
    for _, row in df.iterrows():
        prof = profile_func(row, coverage_dict)
        if prof is not None:
            profiles.append(prof)
    return np.vstack(profiles).mean(axis=0)

# --------------------
# Main loop over features
# --------------------
x_axis = np.arange(-UP, BINS_BODY+DOWN)

for feature, filepath in FEATURE_FILES.items():
    if not os.path.exists(filepath):
        print(f"⚠️ Skipping {feature}, file not found")
        continue

    print(f"Processing {feature}...")

    coverage_dict = build_coverage_dict(filepath)

    avg_profile_mRNA = compute_average_profile(mRNA, get_gene_profile, coverage_dict)
    avg_profile_neg  = compute_average_profile(neg, get_neg_profile, coverage_dict)

    # Plot
    plt.figure(figsize=(12,6))
    plt.plot(x_axis, avg_profile_mRNA, color="blue", lw=2, label="genes")
    plt.plot(x_axis, avg_profile_neg, color="orange", lw=2, label="intergenic")

    plt.axvspan(-UP, 0, color="green", alpha=0.2, label="TSS -1000")
    plt.axvspan(BINS_BODY, BINS_BODY+DOWN, color="red", alpha=0.2, label="TTS +1000")

    plt.title(f"{feature} coverage across genes", fontsize=24)
    plt.xlabel("Position (bp / normalized)")
    plt.ylabel("Average coverage")
    plt.xticks([-1000, 0, 1000, 2000], ["-1000", "TSS", "TTS", "+1000"])
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{feature}_gene_vs_neg_coverage.png")
    plt.close()

print("✅ All features processed and plots saved.")
