#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# --- Input: ONLY lncRNAs ---
LNC_FILE = "/rhome/zli529/lab/LncRNA_chip_prediction/Final/lncRNA_model//Marcos2022_lncRNAs_fully_within_v48_genes.csv"

FEATURE_FILES = {
    "H3R17me2": "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Myriam_2023_Chip/coverage/substract_readcounts/by_antibody/H3R17me2.txt",
    "H2A.Zac":  "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Myriam_2023_Chip/coverage/substract_readcounts/by_antibody/H2A.Zac.txt",
    "H3K4me3":  "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Ashley_2023/coverage/substract_readcounts/by_antibody/H3K4me3.txt",
    "H3K4me2":  "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Ashley_2023/coverage/substract_readcounts/by_antibody/H3K4me2.txt",
    "H3K4me":   "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Ashley_2023/coverage/substract_readcounts/by_antibody/H3K4me.txt",
    "H3K18me":  "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Ashley_2023/coverage/substract_readcounts/by_antibody/H3K18me.txt",
    "H3K27me":  "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Ashley_2023/coverage/substract_readcounts/by_antibody/H3K27me.txt",
    "H3K27ac":  "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Jingyi_Chip/coverage/substract_readcounts/by_antibody/H3K27ac.txt",
    "H3K18ac":  "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Jingyi_Chip/coverage/substract_readcounts/by_antibody/H3K18ac.txt",
    "H3":       "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Jingyi_Chip/coverage/substract_readcounts/by_antibody/H3.txt",
    "H3K4me1":  "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Jingyi_Chip/coverage/substract_readcounts/by_antibody/H3K4me1.txt",
    "H2A.Z":    "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Jingyi_Chip/coverage/substract_readcounts/by_antibody/H2A.Z.txt",
    "H2B.Z":    "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Hoeijmakers_2013_Chip/coverage/substract_readcounts/by_antibody/H2B.Z.txt",
    "MNase":    "/rhome/zli529/lab/chip-seq/chipseq/MNase-seq_results/merge_bam/MNseq_sub.txt",
    "H3K9ac":   "/rhome/zli529/lab/chip-seq/chipseq/H3K9ac_results/merge_bam/H3K9ac.sub.txt",
    "H3K36me3": "/rhome/zli529/lab/chip-seq/chipseq/H3K36me3_results/H3K36me3.Nor.txt"
}

# Parameters
UP = 300
DOWN = 200
BINS_BODY = 200   # normalized lncRNA body bins

# --------------------
# Load lncRNA CSV (robust to different column names)
# --------------------
def load_lnc_table(path):
    df = pd.read_csv(path, sep=None, engine="python")
    cols = {c.lower(): c for c in df.columns}

    def pick(*names, required=True):
        for n in names:
            if n in cols:
                return cols[n]
        if required:
            raise ValueError(f"Missing required column among {names} in {path}. Found: {list(df.columns)}")
        return None

    chr_col   = pick("lnc_chr", "chr", "chrom", "seqname", "contig")
    start_col = pick("lnc_start", "start", "begin")
    end_col   = pick("lnc_end", "end", "stop")
    strand_col= pick("lnc_strand", "strand", required=False)
    id_col    = pick("lnc_id", "id", "name", "transcript_id", required=False)

    out = pd.DataFrame({
        "chr":    df[chr_col].astype(str),
        "start":  pd.to_numeric(df[start_col], errors="coerce"),
        "end":    pd.to_numeric(df[end_col], errors="coerce"),
        "strand": df[strand_col].astype(str) if strand_col else ".",
        "lnc_id": df[id_col].astype(str) if id_col else None
    }).dropna(subset=["start","end"]).copy()

    # ensure start <= end
    out[["start","end"]] = np.sort(out[["start","end"]].values, axis=1)
    out["start"] = out["start"].astype(int)
    out["end"]   = out["end"].astype(int)
    return out

lnc = load_lnc_table(LNC_FILE)

# --------------------
# Helper functions
# --------------------
def build_coverage_dict(read_counts_file):
    """Return dict: {chrom -> {position -> value}}"""
    reads = pd.read_csv(read_counts_file, sep="\t", header=None,
                        names=["chr","site","read_count"])
    coverage_dict = {}
    for chrom, df in reads.groupby("chr"):
        coverage_dict[chrom] = dict(zip(df["site"].astype(int), df["read_count"].astype(float)))
    return coverage_dict

def get_lnc_profile(row, coverage_dict):
    """
    Make a strand-aware profile around the lncRNA:
      upstream (UP), normalized body (BINS_BODY bins), downstream (DOWN)
    """
    chrom = row["chr"]
    if chrom not in coverage_dict:
        return None

    s, e = int(row["start"]), int(row["end"])
    strand = row.get("strand", ".")
    # Define lncRNA "5' end" (TSS-like) and "3' end" (TTS-like)
    l5 = s if strand == "+" else e
    l3 = e if strand == "+" else s

    region_start = l5 - UP
    region_end   = l3 + DOWN

    cov = np.array([coverage_dict[chrom].get(pos, 0.0) for pos in range(region_start, region_end)])
    if cov.size == 0:
        return None

    # Align so that body is cov[UP : len(cov)-DOWN]
    promoter   = cov[:UP]
    body_raw   = cov[UP: len(cov)-DOWN] if len(cov) > (UP+DOWN) else np.array([])
    downstream = cov[-DOWN:] if len(cov) >= DOWN else np.array([])

    # Normalize body to fixed bins
    if body_raw.size > 1:
        body_norm = np.interp(
            np.linspace(0, body_raw.size - 1, BINS_BODY),
            np.arange(body_raw.size),
            body_raw
        )
    else:
        body_norm = np.zeros(BINS_BODY, dtype=float)

    prof = np.concatenate([promoter, body_norm, downstream])

    # For minus strand, flip so 5' end is always at x=0
    if strand == "-":
        prof = prof[::-1]
    return prof

def compute_average_profile(df, profile_func, coverage_dict):
    profiles = []
    for _, row in df.iterrows():
        p = profile_func(row, coverage_dict)
        if p is not None:
            profiles.append(p)
    if not profiles:
        return None
    # Ensure equal lengths; pad/truncate if necessary
    L = UP + BINS_BODY + DOWN
    fixed = [pi[:L] if pi.size >= L else np.pad(pi, (0, L - pi.size)) for pi in profiles]
    return np.vstack(fixed).mean(axis=0)

# --------------------
# Main: plot per feature (lncRNAs only)
# --------------------
x_axis = np.arange(-UP, BINS_BODY + DOWN)

for feature, filepath in FEATURE_FILES.items():
    if not os.path.exists(filepath):
        print(f"⚠️ Skipping {feature}, file not found")
        continue

    print(f"Processing {feature}...")
    coverage_dict = build_coverage_dict(filepath)

    avg_profile_lnc = compute_average_profile(lnc, get_lnc_profile, coverage_dict)
    if avg_profile_lnc is None:
        print(f"⚠️ No lncRNA profiles computed for {feature} (check chromosomes/coordinates).")
        continue

    # Plot
    plt.figure(figsize=(12,6))
    plt.plot(x_axis, avg_profile_lnc, lw=2, label="lncRNA (avg)")

    # Shaded regions: upstream of 5' and downstream of 3'
    plt.axvspan(-UP, 0, alpha=0.15, label="5' upstream")
    plt.axvspan(BINS_BODY, BINS_BODY + DOWN, alpha=0.15, label="3' downstream")

    plt.title(f"{feature} coverage across lncRNAs")
    plt.xlabel("Relative position (bp / normalized)")
    plt.ylabel("Average coverage")
    plt.xticks([-UP, 0, BINS_BODY, BINS_BODY + DOWN],
               [f"-{UP}", "lnc 5' start", "lnc 3' end", f"+{DOWN}"])
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{feature}_lncRNA_coverage.png", dpi=200)
    plt.close()

print("✅ Done. Plots saved as <feature>_lncRNA_coverage.png")
