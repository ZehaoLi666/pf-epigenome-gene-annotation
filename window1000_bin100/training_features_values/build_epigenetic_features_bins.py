
#!/usr/bin/env python3
import argparse
import json
import os
import sys
from typing import Dict, List, Union
import pandas as pd
import numpy as np


def expand_path(path: str) -> str:
    expanded = os.path.expandvars(os.path.expanduser(path))
    if "$" in expanded:
        raise RuntimeError(f"Unresolved environment variable in path: {path}")
    return expanded


def resolve_existing_path(path: str, base_dir: str) -> str:
    expanded = expand_path(path)
    if not os.path.isabs(expanded):
        expanded = os.path.join(base_dir, expanded)
    resolved = os.path.abspath(expanded)
    if not os.path.exists(resolved):
        raise RuntimeError(f"Configured marker file does not exist: {resolved}")
    return resolved

def read_marker_file(path: str) -> pd.DataFrame:
    """
    Read a marker coverage file with columns: chr, pos, value.
    - Whitespace/TSV delimited.
    - Header optional; we force names.
    - 'pos' is 1-based integer coordinate.
    """
    try:
        df = pd.read_csv(
            path, sep=r"\s+", header=None, names=["chr", "pos", "value"],
            dtype={"chr": str, "pos": np.int64, "value": np.float64},
            engine="python", comment="#"
        )
    except Exception as e:
        raise RuntimeError(f"Failed reading marker file {path}: {e}")
    # Ensure sorted for efficient slicing
    df = df.sort_values(["chr", "pos"], kind="mergesort").reset_index(drop=True)
    return df

def merge_replicates(files: List[str]) -> pd.DataFrame:
    """
    Given a list of replicate files for the same marker, return a single tidy DF (chr,pos,value_mean).
    Replicates are merged by outer-join on (chr,pos) and averaged; missing positions treated as 0.
    """
    dfs = []
    for p in files:
        df = read_marker_file(p).rename(columns={"value": os.path.basename(p)})
        dfs.append(df)
    # Outer merge on chr,pos
    merged = dfs[0]
    for df in dfs[1:]:
        merged = merged.merge(df, on=["chr", "pos"], how="outer", sort=True)
    value_cols = [c for c in merged.columns if c not in ("chr", "pos")]
    merged[value_cols] = merged[value_cols].fillna(0.0)
    merged["value"] = merged[value_cols].mean(axis=1)
    merged = merged[["chr", "pos", "value"]].sort_values(["chr", "pos"], kind="mergesort").reset_index(drop=True)
    return merged

def load_markers(marker_config_path: str) -> Dict[str, Dict[str, pd.Series]]:
    """
    Load markers from a JSON config mapping marker names to path(s).
    Returns: dict[marker] -> dict[chr] -> pandas.Series(index=pos, values=coverage)
    """
    marker_config_path = os.path.abspath(expand_path(marker_config_path))
    config_dir = os.path.dirname(marker_config_path)
    with open(marker_config_path, "r") as f:
        cfg = json.load(f)

    markers: Dict[str, Dict[str, pd.Series]] = {}
    for mark, path_or_list in cfg.items():
        if isinstance(path_or_list, list):
            df = merge_replicates([resolve_existing_path(p, config_dir) for p in path_or_list])
        else:
            df = read_marker_file(resolve_existing_path(path_or_list, config_dir))

        # Split by chromosome into Series indexed by position
        chr_groups = {}
        for chrom, sub in df.groupby("chr", sort=False):
            s = pd.Series(sub["value"].to_numpy(), index=sub["pos"].to_numpy(), name=mark)
            chr_groups[chrom] = s
        markers[mark] = chr_groups
    return markers

def compute_bins_for_center(
    s_by_chr: Dict[str, pd.Series],
    chrom: str,
    center: int,
    strand: str,
    window: int,
    bin_size: int
) -> List[float]:
    """
    Compute 20 bins (window/bin_size) mean coverage around a center on a chromosome.
    Missing positions are treated as 0 by dividing sum by bin_size.
    """
    s = s_by_chr.get(chrom, None)
    n_bins = (2 * window) // bin_size
    vals: List[float] = []

    for i in range(n_bins):
        bin_start = center - window + i * bin_size
        bin_end = bin_start + bin_size - 1

        # Clip left edge to 1 to avoid negative coords; right edge left unclipped (missing treated as 0).
        qs = max(1, bin_start)
        qe = bin_end
        if s is None:
            bin_mean = 0.0
        else:
            # sum over available positions
            part = s.loc[qs:qe] if (qs <= qe) else pd.Series(dtype=float)
            bin_sum = float(part.sum()) if not part.empty else 0.0
            # divide by bin_size to treat missing as 0
            bin_mean = bin_sum / float(bin_size)
        vals.append(bin_mean)

    # Reverse for negative strand so upstream->downstream is left->right in returned vector
    if strand == "-":
        vals = vals[::-1]
    return vals

def build_features_for_boundaries(
    boundaries: pd.DataFrame,
    markers: Dict[str, Dict[str, pd.Series]],
    boundary_col: str,
    chr_col: str,
    strand_col: str,
    window: int,
    bin_size: int,
    prefix: str
) -> pd.DataFrame:
    """
    For each row in `boundaries`, compute per-marker binned means around the coordinate in `boundary_col`.
    Returns wide DF: [gene, chr, strand, {mark}_{prefix}_bin01..binN]
    """
    n_bins = (2 * window) // bin_size
    rows = []
    for _, r in boundaries.iterrows():
        chrom = str(r[chr_col])
        center = int(r[boundary_col])
        strand = str(r[strand_col])
        rec = {
            "gene": r.get("gene", r.get("gene_id", r.get("GeneID", ""))),
            "chr": chrom,
            "strand": strand,
            f"{prefix}_coord": center
        }
        for mark, s_by_chr in markers.items():
            vec = compute_bins_for_center(s_by_chr, chrom, center, strand, window, bin_size)
            for j, v in enumerate(vec, 1):
                rec[f"{mark}_{prefix}_bin{j:02d}"] = v
        rows.append(rec)
    df = pd.DataFrame(rows)
    # Order columns: meta, then markers grouped
    meta_cols = ["gene", "chr", "strand", f"{prefix}_coord"]
    mark_cols = [c for c in df.columns if c not in meta_cols]
    df = df[meta_cols + sorted(mark_cols)]
    return df

def main():
    ap = argparse.ArgumentParser(description="Build ±window, binned epigenetic features around TSS/TTS (Path 2).")
    ap.add_argument("--boundaries", required=True, help="CSV with gene boundaries")
    ap.add_argument("--marker-config", required=True, help="JSON: {marker: path or [paths]}")
    ap.add_argument("--out-prefix", required=True, help="Output file prefix (no extension)")
    ap.add_argument("--which", choices=["TSS", "TTS", "BOTH"], default="BOTH", help="Which boundaries to process")
    ap.add_argument("--window", type=int, default=1000, help="Half-window size (bp), default 1000")
    ap.add_argument("--bin-size", type=int, default=100, help="Bin size (bp), default 100")
    # Column mappings
    ap.add_argument("--gene-col", default="gene", help="Gene ID column name")
    ap.add_argument("--chr-col", default="chr", help="Chromosome column name")
    ap.add_argument("--strand-col", default="strand", help="Strand column name (+/-)")
    ap.add_argument("--tss-col", default="tss", help="TSS coordinate column name (1-based)")
    ap.add_argument("--tts-col", default="tts", help="TTS/TES coordinate column name (1-based)")

    args = ap.parse_args()

    # Load boundaries
    try:
        bdf = pd.read_csv(args.boundaries)
    except Exception as e:
        print(f"ERROR: failed to read boundaries CSV {args.boundaries}: {e}", file=sys.stderr)
        sys.exit(1)

    # Normalize required columns; allow alternate casing
    def find_col(df, wanted, default):
        if default in df.columns:
            return default
        cand = None
        for c in df.columns:
            if c.lower() == wanted.lower():
                cand = c; break
        if cand is None:
            raise SystemExit(f"Required column '{default}'/'{wanted}' not found in {args.boundaries}. Found: {list(df.columns)}")
        return cand

    gene_col   = find_col(bdf, args.gene_col, args.gene_col)
    chr_col    = find_col(bdf, args.chr_col, args.chr_col)
    strand_col = find_col(bdf, args.strand_col, args.strand_col)

    if args.which in ("TSS", "BOTH"):
        tss_col = find_col(bdf, args.tss_col, args.tss_col)
    if args.which in ("TTS", "BOTH"):
        tts_col = find_col(bdf, args.tts_col, args.tts_col)

    # Harmonize column names for internal use
    ren = {}
    if gene_col != "gene": ren[gene_col] = "gene"
    if chr_col  != "chr":  ren[chr_col]  = "chr"
    if strand_col != "strand": ren[strand_col] = "strand"
    if args.which in ("TSS", "BOTH") and tss_col != "tss": ren[tss_col] = "tss"
    if args.which in ("TTS", "BOTH") and tts_col != "tts": ren[tts_col] = "tts"
    bdf = bdf.rename(columns=ren)

    # Load markers into dict[marker][chr] -> Series(pos -> coverage)
    markers = load_markers(args.marker_config)

    os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True) if os.path.dirname(args.out_prefix) else None

    if args.which in ("TSS", "BOTH"):
        tss_df = build_features_for_boundaries(
            boundaries=bdf[["gene","chr","strand","tss"]].dropna(),
            markers=markers,
            boundary_col="tss",
            chr_col="chr",
            strand_col="strand",
            window=args.window,
            bin_size=args.bin_size,
            prefix="TSS"
        )
        tss_out = args.out_prefix + "_TSS_features_±{}bp_{}bins.csv".format(args.window, (2*args.window)//args.bin_size)
        tss_df.to_csv(tss_out, index=False)
        print(f"Wrote TSS features to: {tss_out}")

    if args.which in ("TTS", "BOTH"):
        tts_df = build_features_for_boundaries(
            boundaries=bdf[["gene","chr","strand","tts"]].dropna(),
            markers=markers,
            boundary_col="tts",
            chr_col="chr",
            strand_col="strand",
            window=args.window,
            bin_size=args.bin_size,
            prefix="TTS"
        )
        tts_out = args.out_prefix + "_TTS_features_±{}bp_{}bins.csv".format(args.window, (2*args.window)//args.bin_size)
        tts_df.to_csv(tts_out, index=False)
        print(f"Wrote TTS features to: {tts_out}")

if __name__ == "__main__":
    main()
