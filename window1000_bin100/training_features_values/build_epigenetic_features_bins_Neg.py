#!/usr/bin/env python3
import argparse, json, os, sys
from typing import Dict, List
import pandas as pd
import numpy as np

def xp(p: str) -> str:
    expanded = os.path.expandvars(os.path.expanduser(p))
    if "$" in expanded:
        raise RuntimeError(f"Unresolved environment variable in path: {p}")
    return os.path.abspath(expanded)

def resolve_existing_path(path: str, base_dir: str) -> str:
    expanded = os.path.expandvars(os.path.expanduser(path))
    if "$" in expanded:
        raise RuntimeError(f"Unresolved environment variable in path: {path}")
    if not os.path.isabs(expanded):
        expanded = os.path.join(base_dir, expanded)
    resolved = os.path.abspath(expanded)
    if not os.path.exists(resolved):
        raise RuntimeError(f"Configured marker file does not exist: {resolved}")
    return resolved

def read_boundaries(path: str, sep: str, no_header: bool) -> pd.DataFrame:
    """
    Expected columns (headerless by default): chr  tss  tts  gene  strand
    """
    kw = dict(sep=sep, engine="python")
    if no_header:
        kw["header"] = None
        df = pd.read_csv(xp(path), **kw,
                         names=["chr","tss","tts","gene","strand"])
    else:
        df = pd.read_csv(xp(path), **kw)
        # try to normalize column names if a header exists
        ren = {}
        for c in df.columns:
            lc = str(c).lower()
            if lc in ("chrom","chromosome"): ren[c] = "chr"
            elif lc == "tss_coord": ren[c] = "tss"
            elif lc == "tts_coord": ren[c] = "tts"
        if ren: df = df.rename(columns=ren)
        need = ["chr","tss","tts","gene","strand"]
        missing = [c for c in need if c not in df.columns]
        if missing:
            raise SystemExit(f"boundaries missing {missing}. Found {list(df.columns)}")
        df = df[need]
    # types
    df["chr"] = df["chr"].astype(str)
    df["gene"] = df["gene"].astype(str)
    df["strand"] = df["strand"].astype(str)
    df["tss"] = pd.to_numeric(df["tss"], errors="coerce")
    df["tts"] = pd.to_numeric(df["tts"], errors="coerce")
    df = df.dropna(subset=["tss","tts"]).reset_index(drop=True)
    df["tss"] = df["tss"].astype(int)
    df["tts"] = df["tts"].astype(int)
    if df.empty:
        print("No usable rows in boundaries (after parsing).", file=sys.stderr)
    return df

def read_marker_file(path: str) -> pd.DataFrame:
    """
    Coverage file: chr  pos  value   (whitespace-delimited, 1-based pos).
    """
    df = pd.read_csv(xp(path), sep=r"\s+", header=None, comment="#",
                     names=["chr","pos","value"],
                     dtype={"chr":str, "pos":np.int64, "value":np.float64},
                     engine="python")
    return df.sort_values(["chr","pos"], kind="mergesort").reset_index(drop=True)

def merge_reps(files: List[str]) -> pd.DataFrame:
    dfs = []
    for p in files:
        dfi = read_marker_file(p).rename(columns={"value": os.path.basename(xp(p))})
        dfs.append(dfi)
    merged = dfs[0]
    for dfi in dfs[1:]:
        merged = merged.merge(dfi, on=["chr","pos"], how="outer", sort=True)
    val_cols = [c for c in merged.columns if c not in ("chr","pos")]
    merged[val_cols] = merged[val_cols].fillna(0.0)
    merged["value"] = merged[val_cols].mean(axis=1)
    merged = merged[["chr","pos","value"]]
    return merged.sort_values(["chr","pos"], kind="mergesort").reset_index(drop=True)

def load_markers(marker_cfg: str) -> Dict[str, Dict[str, pd.Series]]:
    """
    JSON: { "H3K4me3": "/path/file.txt", "H2A.Z": "/path/file.txt", ... }
    or replicate list: "H3K4me3": ["/path/rep1.txt", "/path/rep2.txt"]
    Returns dict[mark][chrom] -> Series(index=pos, values=value)
    """
    cfg_path = xp(marker_cfg)
    cfg_dir = os.path.dirname(cfg_path)
    cfg = json.load(open(cfg_path))
    out: Dict[str, Dict[str, pd.Series]] = {}
    for mark, src in cfg.items():
        if isinstance(src, list):
            df = merge_reps([resolve_existing_path(p, cfg_dir) for p in src])
        else:
            df = read_marker_file(resolve_existing_path(src, cfg_dir))
        by_chr: Dict[str, pd.Series] = {}
        for chrom, sub in df.groupby("chr", sort=False):
            by_chr[chrom] = pd.Series(sub["value"].to_numpy(),
                                      index=sub["pos"].to_numpy(), name=mark)
        out[mark] = by_chr
    return out

def bin_means(series_by_chr: Dict[str, pd.Series],
              chrom: str, center: int, strand: str,
              window: int, bin_size: int) -> List[float]:
    """
    Compute mean value per-bin across ±window around center.
    Missing positions contribute 0. Reverse bins for '-' strand.
    """
    s = series_by_chr.get(chrom, None)
    n_bins = (2*window)//bin_size
    out: List[float] = []
    for i in range(n_bins):
        b0 = center - window + i*bin_size
        b1 = b0 + bin_size - 1
        qs = max(1, b0); qe = b1
        if s is None or qs > qe:
            m = 0.0
        else:
            part = s.loc[qs:qe]  # inclusive on label index
            m = float(part.sum())/float(bin_size) if not part.empty else 0.0
        out.append(m)
    return out[::-1] if strand == "-" else out

def main():
    ap = argparse.ArgumentParser(description="Assign binned epigenetic values around TSS/TTS (single-pass, no chunking).")
    ap.add_argument("--boundaries", required=True, help="negatives_2k.tsv (chr tss tts gene strand)")
    ap.add_argument("--marker-config", required=True, help="JSON mapping mark -> file or [files]")
    ap.add_argument("--out-csv", required=True, help="Output CSV (combined TSS+TTS features)")
    ap.add_argument("--window", type=int, default=1000)
    ap.add_argument("--bin-size", type=int, default=100)
    ap.add_argument("--sep", default=r"\s+", help="Delimiter for boundaries (default: whitespace)")
    ap.add_argument("--no-header", action="store_true", help="Set if boundaries has no header (default for negatives_2k.tsv)")
    args = ap.parse_args()

    # 1) Read inputs
    bounds = read_boundaries(args.boundaries, args.sep, args.no_header)
    if bounds.empty:
        # write an empty file with headers so pipelines don’t hang
        n_bins = (2*args.window)//args.bin_size
        cols = ["gene","chr","strand","TSS_coord","TTS_coord"]
        for mark in ["MARK"]:  # placeholder to build shape
            pass
        # we cannot know marks yet (need config), so just exit clean with message
        print("No usable rows in boundaries; nothing to do.", file=sys.stderr)
        # produce an empty file with just meta headers (optional)
        pd.DataFrame(columns=["gene","chr","strand","TSS_coord","TTS_coord"]).to_csv(xp(args.out_csv), index=False)
        return

    markers = load_markers(args.marker_config)
    n_bins = (2*args.window)//args.bin_size
    marks = sorted(markers.keys())

    # 2) Build features (combined TSS + TTS)
    rows = []
    for _, r in bounds.iterrows():
        chrom  = r["chr"]; strand = r["strand"]
        tss_c  = int(r["tss"]); tts_c = int(r["tts"])
        rec = {"gene": r["gene"], "chr": chrom, "strand": strand,
               "TSS_coord": tss_c, "TTS_coord": tts_c}
        for mark in marks:
            s_by_chr = markers[mark]
            tss_vec = bin_means(s_by_chr, chrom, tss_c, strand, args.window, args.bin_size)
            tts_vec = bin_means(s_by_chr, chrom, tts_c, strand, args.window, args.bin_size)
            for j, v in enumerate(tss_vec, 1):
                rec[f"{mark}_TSS_bin{j:02d}"] = v
            for j, v in enumerate(tts_vec, 1):
                rec[f"{mark}_TTS_bin{j:02d}"] = v
        rows.append(rec)

    out_df = pd.DataFrame(rows)
    out_df.to_csv(xp(args.out_csv), index=False)
    print(f"Wrote {len(out_df)} rows to {args.out_csv}")

if __name__ == "__main__":
    main()
