#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
filter_TSS.py
-------------
Constrain candidate TSS to the *nearby intergenic* of the nearby gene with a strand-specific rule:

  • If the nearby gene’s strand is '+': keep only candidates in the *upstream* intergenic
    (i.e., the gap directly before the gene) and with pos < gene.start (BED-style start).

  • If the nearby gene’s strand is '−': keep only candidates in the *downstream intergenic*
    (i.e., the gap directly after the gene) and with pos > gene.end.

Among multiple valid TSS for the same target gene, keep the **farthest** one from the boundary;
ties are broken by **higher score**.

If the TSS table provides a column 'gene_id', we enforce that the selected candidate’s
target gene equals that 'gene_id'.

Inputs:
  - Genes BED (6 cols): chr, start, end, gene_id, score(ignored), strand
  - Chrom sizes (FAI or 2-col TSV: chr\tlength) — enables head/tail gaps at chrom ends
  - TSS table (TSV/CSV): must contain at least
        chr, and one of {pos, mid, (start & end)} to infer position, and score
    Optional columns: strand, gene_id

Outputs:
  - <out-prefix>.tsv    : rich TSV with gap context and chosen target gene
  - <out-prefix>.bed    : 1-bp BED (chr, pos, pos, target_gene, score, strand) 
"""

import sys
import argparse
import numpy as np
import pandas as pd


# -------------------- I/O helpers --------------------

def read_table_auto(path, dtype=None):
    """Auto-detect delimiter (tab vs comma)."""
    with open(path, "r") as fh:
        first = fh.readline()
    sep = "," if first.count(",") > first.count("\t") else "\t"
    return pd.read_csv(path, sep=sep, dtype=dtype, low_memory=False)


def load_genes_bed(path):
    """Genes BED: chr, start, end, gene_id, score, strand -> sorted DataFrame."""
    cols = ["chr","start","end","gene_id","score","strand"]
    g = pd.read_csv(path, sep="\t", header=None, names=cols,
                    usecols=["chr","start","end","gene_id","strand"])
    g = g.sort_values(["chr","start","end"]).reset_index(drop=True)
    return g


def load_chrom_sizes(path):
    """Load chromosome sizes from .fai or 2-col TSV: chr, size."""
    tab = pd.read_csv(path, sep="\t", header=None, usecols=[0,1], names=["chr","size"])
    tab["chr"] = tab["chr"].astype(str)
    tab["size"] = tab["size"].astype(int)
    return tab


# -------------------- intergenic building & indexing --------------------

def build_intergenic_gaps(genes_df, chrom_sizes=None):
    """
    Build intergenic gaps between consecutive genes per chromosome.
    If chrom_sizes provided, add head (0->first.start) and tail (last.end->chrom.size) gaps.

    Returns columns:
      chr, start, end, left_gene_id, left_strand, right_gene_id, right_strand, length
    """
    rows = []
    size_map = dict(zip(chrom_sizes["chr"], chrom_sizes["size"])) if chrom_sizes is not None else {}
    for chrom, g in genes_df.groupby("chr", sort=False):
        g = g.sort_values(["start","end"]).reset_index(drop=True)
        csize = size_map.get(chrom, None)

        # head gap
        if csize is not None and len(g) > 0 and int(g.iloc[0]["start"]) > 0:
            rows.append({
                "chr": chrom, "start": 0, "end": int(g.iloc[0]["start"]),
                "left_gene_id": np.nan, "left_strand": np.nan,
                "right_gene_id": g.iloc[0]["gene_id"], "right_strand": g.iloc[0]["strand"]
            })

        # internal gaps
        for i in range(len(g) - 1):
            L, R = g.iloc[i], g.iloc[i+1]
            s, e = int(L["end"]), int(R["start"])
            if e > s:
                rows.append({
                    "chr": chrom, "start": s, "end": e,
                    "left_gene_id": L["gene_id"], "left_strand": L["strand"],
                    "right_gene_id": R["gene_id"], "right_strand": R["strand"]
                })

        # tail gap
        if csize is not None and len(g) > 0 and int(g.iloc[-1]["end"]) < int(csize):
            last = g.iloc[-1]
            rows.append({
                "chr": chrom, "start": int(last["end"]), "end": int(csize),
                "left_gene_id": last["gene_id"], "left_strand": last["strand"],
                "right_gene_id": np.nan, "right_strand": np.nan
            })

    inter = pd.DataFrame(rows)
    if inter.empty:
        inter["length"] = pd.Series(dtype="int64")
    else:
        inter["length"] = (inter["end"] - inter["start"]).astype(int)
    return inter


def index_intergenic(inter_df):
    """Per-chrom index for fast point-in-interval lookup."""
    idx = {}
    for chrom, df in inter_df.groupby("chr", sort=False):
        df = df.sort_values("start").reset_index(drop=True)
        idx[chrom] = {
            "df": df,
            "starts": df["start"].to_numpy(np.int64),
            "ends": df["end"].to_numpy(np.int64),
        }
    return idx


def assign_gap_point(ig_idx, chrom, pos):
    """Return the gap row for (chrom,pos) or None if outside all gaps."""
    bucket = ig_idx.get(chrom)
    if bucket is None:
        return None
    i = np.searchsorted(bucket["starts"], pos, side="right") - 1
    if i >= 0 and pos < bucket["ends"][i]:
        return bucket["df"].iloc[i]
    return None


# -------------------- core refinement --------------------

def normalize_tss_table(raw: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure required columns exist:
      - chr
      - pos (accept 'mid' or compute from 'start' & 'end')
      - score
      - optional: strand, gene_id
    """
    df = raw.copy()
    cols = set(df.columns)

    if "chr" not in cols and "chrom" in cols:
        df = df.rename(columns={"chrom": "chr"})
        cols = set(df.columns)
    if "pos" not in cols:
        if "mid" in cols:
            df = df.rename(columns={"mid": "pos"})
        elif {"start", "end"} <= cols:
            df["pos"] = ((df["start"].astype(float) + df["end"].astype(float)) / 2.0).astype(int)
        else:
            sys.exit("ERROR: TSS table needs 'pos' or 'mid' or both 'start' and 'end' to infer position.")
    if "score" not in cols:
        if "prob" in cols:
            df = df.rename(columns={"prob": "score"})
        else:
            sys.exit("ERROR: TSS table needs a 'score' column (or 'prob').")

    # Ensure types
    df["chr"] = df["chr"].astype(str)
    df["pos"] = df["pos"].astype(int)
    if "strand" not in df.columns:
        df["strand"] = "+"

    # Keep only relevant columns (but preserve gene_id if present)
    keep = ["chr", "pos", "score", "strand"] + (["gene_id"] if "gene_id" in df.columns else [])
    return df[keep].copy()


def refine_tss_strand_specific(
    genes_df: pd.DataFrame,
    chrom_sizes_df: pd.DataFrame,
    tss_df: pd.DataFrame,
    score_min: float = 0.70,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Keep TSS only in 'nearby intergenic' per strand rule:
      + gene: upstream of gene.start (gap directly before the gene), pos < gene.start
      - gene: downstream of gene.end   (gap directly after  the gene), pos > gene.end

    Among multiple valid TSS for the same target gene, keep the **farthest** one from the boundary;
    ties by higher score. If tss_df has 'gene_id', enforce it matches the selected target gene.

    Returns one TSS per target gene.
    """
    # Maps for gene coords
    gene_start = dict(zip(genes_df["gene_id"], genes_df["start"]))
    gene_end   = dict(zip(genes_df["gene_id"], genes_df["end"]))

    # Build & index intergenic
    inter = build_intergenic_gaps(genes_df, chrom_sizes_df)
    if verbose:
        print(f"[i] Intergenic gaps: {len(inter):,} (mean length {inter['length'].mean():.1f} bp)")
    ig_idx = index_intergenic(inter)

    # Filter TSS by score and assign to gaps
    base_cols = [c for c in ["chr","pos","score","strand","gene_id"] if c in tss_df.columns]
    tss = tss_df.loc[tss_df["score"] >= score_min, base_cols].copy()
    tss["chr"] = tss["chr"].astype(str)
    tss["pos"] = tss["pos"].astype(int)
    if "strand" not in tss.columns:
        tss["strand"] = "+"
    if verbose:
        print(f"[i] TSS above score ≥ {score_min}: {len(tss):,}")

    # Assign each TSS to a gap
    rec = {k: [] for k in ["gap_start","gap_end","left_gene_id","left_strand","right_gene_id","right_strand","gap_len"]}
    keep_mask = []
    for _, r in tss.iterrows():
        g = assign_gap_point(ig_idx, r["chr"], int(r["pos"]))
        if g is None:
            keep_mask.append(False)
            for k in rec: rec[k].append(np.nan)
        else:
            keep_mask.append(True)
            rec["gap_start"].append(int(g["start"]))
            rec["gap_end"].append(int(g["end"]))
            rec["left_gene_id"].append(g["left_gene_id"])
            rec["left_strand"].append(g["left_strand"])
            rec["right_gene_id"].append(g["right_gene_id"])
            rec["right_strand"].append(g["right_strand"])
            rec["gap_len"].append(int(g["length"]))

    tss = tss.reset_index(drop=True)
    for k, v in rec.items():
        tss[k] = v
    tss = tss.loc[keep_mask].copy()
    if verbose:
        print(f"[i] TSS that fall in intergenic gaps: {len(tss):,}")

    # Neighbor gene coordinates (for comparisons)
    tss["right_start"] = tss["right_gene_id"].map(gene_start)  # gene start for right gene
    tss["left_end"]    = tss["left_gene_id"].map(gene_end)     # gene end   for left  gene

    # Decide side per row BEFORE filtering; then filter and reuse the tag (broadcast-safe)
    mask_plus  = (tss["right_gene_id"].notna()) & (tss["right_strand"] == "+") & (tss["pos"] <  tss["right_start"])
    mask_minus = (tss["left_gene_id"].notna())  & (tss["left_strand"]  == "-") & (tss["pos"] >  tss["left_end"])
    tss["keep_side"]   = np.where(mask_plus, "plus", np.where(mask_minus, "minus", "none"))

    tss = tss.loc[tss["keep_side"].isin(["plus","minus"])].copy()
    if verbose:
        print(f"[i] After strand rule (+ upstream / - downstream): {len(tss):,}")

    # Define target gene/strand/boundary per row based on side
    tss["target_gene"]   = np.where(tss["keep_side"] == "plus",  tss["right_gene_id"], tss["left_gene_id"])
    tss["target_strand"] = np.where(tss["keep_side"] == "plus",  tss["right_strand"],  tss["left_strand"])
    # boundary used for distance:
    #   + gene -> boundary = gene.start (right_start)
    #   - gene -> boundary = gene.end   (left_end)
    tss["target_boundary"] = np.where(tss["keep_side"] == "plus", tss["right_start"], tss["left_end"])

    # If provided, enforce TSS.gene_id == target_gene
    if "gene_id" in tss.columns and tss["gene_id"].notna().any():
        before = len(tss)
        tss = tss.loc[tss["gene_id"] == tss["target_gene"]].copy()
        if verbose:
            print(f"[i] Enforcing provided gene_id match: {before} -> {len(tss)}")

    # --- NEW: choose the farthest candidate (primary), then higher score (secondary) ---
    # distance definition (non-negative by construction due to filters above):
    #   + gene: distance = target_boundary - pos  (pos < boundary)
    #   - gene: distance = pos - target_boundary  (pos > boundary)
    dist_plus  = (tss["target_boundary"] - tss["pos"]).where(tss["keep_side"] == "plus", np.nan)
    dist_minus = (tss["pos"] - tss["target_boundary"]).where(tss["keep_side"] == "minus", np.nan)
    tss["distance"] = dist_plus.fillna(dist_minus).astype(int)

    # Sort by: target_gene, distance DESC (farthest first), score DESC
    tss = tss.sort_values(["target_gene", "distance", "score"], ascending=[True, False, False]).copy()

    # One per gene: farthest; tie by score already handled by sort
    final = tss.groupby("target_gene", as_index=False).head(1)

    # Reorder columns
    cols = [
        "chr","pos","strand","score",
        "gap_start","gap_end","gap_len",
        "left_gene_id","left_strand","right_gene_id","right_strand",
        "target_gene","target_strand","target_boundary","distance"
    ]
    return final[cols].copy()


# -------------------- CLI --------------------

def main():
    ap = argparse.ArgumentParser(description="Filter candidate TSS to strand-specific nearby intergenic regions (keep farthest match per gene).")
    ap.add_argument("--genes", required=True, help="Genes BED (6 cols): chr start end gene_id score strand")
    ap.add_argument("--tss",   required=True, help="TSS table (TSV/CSV). Needs chr, and pos|mid|(start&end), and score. Optional: strand, gene_id")
    ap.add_argument("--chrom-sizes", required=True, help="Chrom sizes (.fai or 2-col TSV: chr size)")
    ap.add_argument("--score-min", type=float, default=0.70, help="Minimum TSS score to consider (default: 0.70)")
    ap.add_argument("--out-prefix", required=True, help="Output prefix (writes <prefix>.tsv and <prefix>.bed)")
    args = ap.parse_args()

    # Load inputs
    genes = load_genes_bed(args.genes)
    chrom_sizes = load_chrom_sizes(args.chrom_sizes)
    raw_tss = read_table_auto(args.tss)

    # Normalize TSS columns
    tss = normalize_tss_table(raw_tss)

    # Run refinement
    final = refine_tss_strand_specific(genes, chrom_sizes, tss, score_min=args.score_min, verbose=True)

    # Write TSV
    tsv_path = f"{args.out_prefix}.tsv"
    final.to_csv(tsv_path, sep="\t", index=False)
    print(f"[✓] Wrote: {tsv_path}  (n={len(final):,})")

    # Write 1-bp BED
    bed = final[["chr","pos","pos","target_gene","score","strand"]].copy()
    bed["pos"] = bed["pos"].astype(int)
    bed_path = f"{args.out_prefix}.bed"
    bed.to_csv(bed_path, sep="\t", header=False, index=False)
    print(f"[✓] Wrote: {bed_path}")

if __name__ == "__main__":
    main()
