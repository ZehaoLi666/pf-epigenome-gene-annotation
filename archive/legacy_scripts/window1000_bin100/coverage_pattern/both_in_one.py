#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# --- Input files ---
mRNA_file = "/rhome/zli529/lab/LncRNA_chip_prediction/Final/window1000_bin100/training_features_values/reference_gene_boundaries.csv"
neg_file = "/rhome/zli529/lab/LncRNA_chip_prediction/Final/window1000_bin100/training_features_values/negatives_2k.tsv"
read_counts_file = "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Jingyi_Chip/coverage/substract_readcounts/by_antibody/H2A.Z.txt"

# Parameters
UP = 500
DOWN = 500
BINS_BODY = 1000   # normalized gene body bins

# --- Load data ---
mRNA = pd.read_csv(mRNA_file, sep="\t", header=None,
                   names=["chr","start","end","gene_id","strand"])
neg = pd.read_csv(neg_file, sep="\t", header=None,
                  names=["chr","start","end","id","strand"])
reads = pd.read_csv(read_counts_file, sep="\t", header=None,
                   names=["chr","site","read_count"])

# --- Build coverage dict ---
coverage_dict = {}
for chrom, df in reads.groupby("chr"):
    cov = dict(zip(df["site"], df["read_count"]))
    coverage_dict[chrom] = cov

# --------------------
# mRNA gene coverage
# --------------------
def get_gene_profile(row):
    chrom = row["chr"]
    if chrom not in coverage_dict:
        return None

    start, end, strand = row["start"], row["end"], row["strand"]
    tss = start if strand == "+" else end
    tts = end if strand == "+" else start

    # region from TSS-500 .. TTS+500
    region_start = tss - UP if strand == "+" else tts - UP
    region_end   = tts + DOWN if strand == "+" else tss + DOWN

    cov = [coverage_dict[chrom].get(pos, 0) for pos in range(region_start, region_end)]
    cov = np.array(cov)

    if strand == "-":
        cov = cov[::-1]

    # promoter / body / downstream
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

# --------------------
# Neg window coverage (background, not normalized)
# --------------------
def get_neg_profile(row):
    chrom = row["chr"]
    if chrom not in coverage_dict:
        return None

    # take 2000 bp window around Neg center
    center = (row["start"] + row["end"]) // 2
    region_start = center - UP
    region_end   = center + BINS_BODY + DOWN   # force same length as gene profile

    cov = [coverage_dict[chrom].get(pos, 0) for pos in range(region_start, region_end)]
    return np.array(cov)

# --------------------
# Average profile helper
# --------------------
def compute_average_profile(df, profile_func): 
    profiles = []
    for _, row in df.iterrows():
        prof = profile_func(row)
        if prof is not None:
            profiles.append(prof)
    return np.vstack(profiles).mean(axis=0)

# --- Compute profiles ---
avg_profile_mRNA = compute_average_profile(mRNA, get_gene_profile)
avg_profile_neg = compute_average_profile(neg, get_neg_profile)

# --- Build shared x-axis: -500 .. 0 .. 1000 .. 1500 ---
x_axis = np.arange(-UP, BINS_BODY+DOWN)

# --------------------
# Plot
# --------------------
plt.figure(figsize=(12,6))
plt.plot(x_axis, avg_profile_mRNA, color="blue", lw=2, label="genes")
plt.plot(x_axis, avg_profile_neg, color="orange", lw=2, label="intergenic")

# shaded regions
plt.axvspan(-UP, 0, color="green", alpha=0.2, label="TSS -500")
plt.axvspan(BINS_BODY, BINS_BODY+DOWN, color="red", alpha=0.2, label="TTS +500")

plt.title("H2A.Z coverage across genes", fontsize=18)
plt.xlabel("Position (bp / normalized)")
plt.ylabel("Average coverage")
plt.xticks([-500, 0, 1000, 1500], ["-500", "TSS", "TTS", "+500"], fontsize=14)
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig("H2AZ_gene_vs_neg_coverage.png")
plt.show()
