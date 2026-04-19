#!/usr/bin/env python3
"""
build_bins_for_intergenic_midpoints.py

Assign 20-bin (100 bp each) mean coverage around midpoints (±1000 bp) to intergenic midpoints.
One invocation handles ONE marker (potentially with replicates averaged). Output is a wide CSV:
    chr, mid, ID, strand, {MARK}_bin01,...,{MARK}_bin20

Usage (examples):
  python build_bins_for_intergenic_midpoints.py \
      --midpoints /mnt/data/intergenic_midpoints_100bp.csv \
      --marker-config training_features_values/marker_files.json \
      --marker-index 1 \
      --out-dir prediction/intergenic_bins \
      --window 1000 --bin-size 100

  python build_bins_for_intergenic_midpoints.py \
      --midpoints /mnt/data/intergenic_midpoints_100bp.csv \
      --marker-config training_features_values/marker_files.json \
      --marker-name H3K4me3 \
      --out-dir prediction/intergenic_bins

JSON format (order matters when using --marker-index):
{
  "H3K4me3": ["/path/rep1.txt", "/path/rep2.txt"],
  "H2AZ": "/path/H2AZ.txt",
  ...
}
"""

import argparse
import json
import os
import sys
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
from collections import defaultdict
from tqdm import tqdm


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

# ----------------------------- IO HELPERS -----------------------------

def read_midpoints(path: str) -> pd.DataFrame:
    """Load midpoints CSV and normalize to columns: chr, mid, ID, strand"""
    df = pd.read_csv(path)
    # lower for matching, keep original names too
    orig_cols = list(df.columns)
    low = {c.lower(): c for c in orig_cols}

    # Robustly pick columns
    chr_col = low.get("chr", orig_cols[0])  # assume first if missing
    # prefer the FIRST 'mid' occurrence
    mid_candidates = [c for c in orig_cols if c.lower() == "mid"]
    if not mid_candidates:
        raise SystemExit("No 'mid' column found in midpoints file.")
    mid_col = mid_candidates[0]

    id_col = low.get("id", None)
    if id_col is None:
        # Make one if absent
        df["ID"] = np.arange(1, len(df) + 1)
        id_col = "ID"
    strand_col = low.get("strand", None)
    if strand_col is None:
        raise SystemExit("No 'strand' column found in midpoints file (expected '+'/'-').")

    out = df[[chr_col, mid_col, id_col, strand_col]].copy()
    out.columns = ["chr", "mid", "ID", "strand"]

    # Types
    out["mid"] = out["mid"].astype(int)
    out["ID"] = out["ID"].astype(str)  # keep generic
    out["chr"] = out["chr"].astype(str)
    out["strand"] = out["strand"].astype(str)
    return out

def read_marker_file(path: str) -> pd.DataFrame:
    """Read whitespace/TSV per-base coverage: chr pos value (1-based)."""
    df = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        names=["chr", "pos", "value"],
        dtype={"chr": str, "pos": np.int64, "value": np.float64},
        engine="python",
        comment="#",
    )
    df = df.sort_values(["chr", "pos"], kind="mergesort").reset_index(drop=True)
    return df

def merge_replicates(files: List[str]) -> pd.DataFrame:
    """Outer-join on (chr,pos), fill NA with 0, average across replicates."""
    dfs = []
    for p in files:
        d = read_marker_file(p).rename(columns={"value": os.path.basename(p)})
        dfs.append(d)
    merged = dfs[0]
    for d in dfs[1:]:
        merged = merged.merge(d, on=["chr", "pos"], how="outer", sort=True)
    value_cols = [c for c in merged.columns if c not in ("chr", "pos")]
    merged[value_cols] = merged[value_cols].fillna(0.0)
    merged["value"] = merged[value_cols].mean(axis=1)
    merged = merged[["chr", "pos", "value"]].sort_values(["chr", "pos"], kind="mergesort")
    merged = merged.reset_index(drop=True)
    return merged

def load_marker_series(marker_cfg: dict, marker_name: str) -> Dict[str, pd.Series]:
    """Return dict[chr] -> Series(index=pos, values=value) for the marker."""
    entry = marker_cfg.get(marker_name, None)
    if entry is None:
        raise SystemExit(f"Marker '{marker_name}' not found in config.")
    if isinstance(entry, list):
        df = merge_replicates(entry)
    else:
        df = read_marker_file(entry)
    # Split by chromosome
    per_chr = {}
    for chrom, sub in df.groupby("chr", sort=False):
        s = pd.Series(sub["value"].to_numpy(), index=sub["pos"].to_numpy(), name=marker_name)
        per_chr[chrom] = s
    return per_chr

# ------------------------- PREFIX / BIN MATH --------------------------

def compute_required_max(midpoints: pd.DataFrame, half_window: int) -> Dict[str, int]:
    req = defaultdict(int)
    for r in midpoints.itertuples(index=False):
        chrom = r.chr
        mm = int(r.mid)
        bound = mm + half_window
        if bound > req[chrom]:
            req[chrom] = bound
    return req

def build_prefix_from_series(per_chr_series: Dict[str, pd.Series], req_max_by_chr: Dict[str, int]) -> Dict[str, np.ndarray]:
    """
    Build 1-based prefix sum array per chromosome up to required max.
    prefix[0] = 0, prefix[i] = sum(values[1..i]).
    """
    out = {}
    for chrom, max_needed in req_max_by_chr.items():
        n = int(max_needed)
        arr = np.zeros(n + 1, dtype=np.float64)  # positions 1..n
        s = per_chr_series.get(chrom, None)
        if s is not None:
            # clip to [1..n]
            valid = s.index[(s.index >= 1) & (s.index <= n)]
            if len(valid):
                arr[valid] = s.loc[valid].to_numpy(dtype=np.float64)
        pref = np.zeros(n + 1, dtype=np.float64)
        pref[1:] = np.cumsum(arr[1:])
        out[chrom] = pref
    return out

def range_mean(prefix: np.ndarray, L: int, R: int, bin_size: int) -> float:
    """Mean over [L..R] (inclusive), treating out-of-bound as 0 but dividing by bin_size."""
    if prefix is None:
        return 0.0
    n = len(prefix) - 1
    Lq = max(1, int(L))
    Rq = min(int(R), n)
    if Lq > Rq:
        return 0.0
    s = prefix[Rq] - prefix[Lq - 1]
    return float(s) / float(bin_size)

def compute_bins_for_mid(prefix_by_chr: Dict[str, np.ndarray],
                         chrom: str, mid: int, strand: str,
                         half_window: int, bin_size: int) -> List[float]:
    n_bins = (2 * half_window) // bin_size
    pref = prefix_by_chr.get(chrom, None)
    vals = []
    win_start = mid - half_window
    for b in range(n_bins):
        b_start = win_start + b * bin_size
        b_end = b_start + bin_size - 1
        vals.append(range_mean(pref, b_start, b_end, bin_size))
    if strand == "-":
        vals = vals[::-1]
    return vals

# ------------------------------ MAIN ---------------------------------

def main():
    ap = argparse.ArgumentParser(description="Assign 20-bin means around intergenic midpoints for ONE marker.")
    ap.add_argument("--midpoints", required=True, help="CSV with columns: chr, mid, ID, strand")
    ap.add_argument("--marker-config", required=True, help="JSON: {marker: path or [paths]}")
    ap.add_argument("--marker-index", type=int, default=None, help="1-based index of marker key in JSON (order matters)")
    ap.add_argument("--marker-name", default=None, help="Explicit marker name (overrides --marker-index)")
    ap.add_argument("--out-dir", required=True, help="Output directory")
    ap.add_argument("--window", type=int, default=1000, help="Half-window (bp). Default 1000")
    ap.add_argument("--bin-size", type=int, default=100, help="Bin size (bp). Default 100")
    ap.add_argument("--log1p", action="store_true", help="Apply log1p to bin means")
    args = ap.parse_args()

    # Load midpoints
    mids = read_midpoints(args.midpoints)

    # Load marker config
    cfg_path = os.path.abspath(expand_path(args.marker_config))
    cfg_dir = os.path.dirname(cfg_path)
    with open(cfg_path, "r") as f:
        cfg = json.load(f)
    normalized_cfg = {}
    for mark, entry in cfg.items():
        if isinstance(entry, list):
            normalized_cfg[mark] = [resolve_existing_path(p, cfg_dir) for p in entry]
        else:
            normalized_cfg[mark] = resolve_existing_path(entry, cfg_dir)
    # choose marker
    if args.marker_name:
        marker = args.marker_name
    else:
        if args.marker_index is None:
            raise SystemExit("Provide either --marker-name or --marker-index.")
        # preserve insertion order (Py3.7+)
        keys = list(normalized_cfg.keys())
        if not (1 <= args.marker_index <= len(keys)):
            raise SystemExit(f"--marker-index out of range (1..{len(keys)}).")
        marker = keys[args.marker_index - 1]

    print(f"[INFO] Processing marker: {marker}")

    # Load per-base coverage (merge replicates if needed), split by chr->Series
    per_chr_series = load_marker_series(normalized_cfg, marker)

    # Compute required max per chr and build prefix arrays
    req_max = compute_required_max(mids, args.window)
    prefix = build_prefix_from_series(per_chr_series, req_max)

    # Compute bins
    n_bins = (2 * args.window) // args.bin_size
    if (2 * args.window) % args.bin_size != 0:
        raise SystemExit("Window must be divisible by bin-size.")

    rows = []
    bin_cols = [f"{marker}_bin{j:02d}" for j in range(1, n_bins + 1)]
    out_cols = ["chr", "mid", "ID", "strand"] + bin_cols

    for r in tqdm(mids.itertuples(index=False), total=len(mids), desc=f"Binning {marker}"):
        chrom = r.chr
        mid   = int(r.mid)
        strand = r.strand
        vec = compute_bins_for_mid(prefix, chrom, mid, strand, args.window, args.bin_size)
        if args.log1p:
            vec = list(np.log1p(vec))
        rows.append([chrom, mid, r.ID, strand] + vec)

    os.makedirs(args.out_dir, exist_ok=True)
    out_path = os.path.join(args.out_dir, f"{marker}_intergenic_bins_pm{args.window}bp_{n_bins}bins.csv")
    pd.DataFrame(rows, columns=out_cols).to_csv(out_path, index=False)
    print(f"[OK] Wrote: {out_path}")

if __name__ == "__main__":
    main()
