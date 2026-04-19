#!/usr/bin/env python3 
import argparse, sys, csv, random
from typing import List, Tuple, Dict, Optional
import pandas as pd
from bisect import bisect_left

def read_chrom_sizes(path: str) -> Dict[str,int]:
    sizes = {}
    with open(path) as fh:
        for line in fh:
            line=line.strip()
            if not line or line.startswith("#"): continue
            chrom, size = line.split()[:2]
            sizes[chrom] = int(size)
    return sizes

def read_tsv(path: str, sep: str, no_header: bool, col_order: str) -> pd.DataFrame:
    kw = dict(sep=sep, engine="python")
    if no_header:
        kw["header"] = None
        df = pd.read_csv(path, **kw)
        cols = [c.strip() for c in col_order.split(",")]
        if len(cols) != df.shape[1]:
            sys.exit(f"--col-order has {len(cols)} names but file has {df.shape[1]} columns")
        df.columns = cols
    else:
        df = pd.read_csv(path, **kw)
    need = ["chr","tss","tts","gene","strand"]
    for c in need:
        if c not in df.columns:
            sys.exit(f"{path}: missing column '{c}'. Found: {list(df.columns)}")
    df["tss"] = pd.to_numeric(df["tss"], errors="coerce").astype("Int64")
    df["tts"] = pd.to_numeric(df["tts"], errors="coerce").astype("Int64")
    df["chr"] = df["chr"].astype(str)
    df["strand"] = df["strand"].astype(str)
    df["gene"] = df["gene"].astype(str)
    return df[need].dropna()

def read_bed(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", header=None, comment="#", engine="python")
    if df.shape[1] < 3: sys.exit(f"{path}: expected 3+ columns")
    df = df.iloc[:, :3].copy()
    df.columns = ["chr","start0","end0"]
    return df

def merge_bed(df: pd.DataFrame) -> pd.DataFrame:
    out = []
    for chrom, sub in df.groupby("chr", sort=False):
        sub = sub.sort_values(["start0","end0"])
        cur_s, cur_e = None, None
        for s,e in zip(sub["start0"], sub["end0"]):
            s = int(s); e = int(e)
            if cur_s is None:
                cur_s, cur_e = s, e
            elif s <= cur_e:
                if e > cur_e: cur_e = e
            else:
                out.append((chrom, cur_s, cur_e))
                cur_s, cur_e = s, e
        if cur_s is not None:
            out.append((chrom, cur_s, cur_e))
    return pd.DataFrame(out, columns=["chr","start0","end0"])

def build_exclude_index(bed: Optional[pd.DataFrame]) -> Dict[str, List[Tuple[int,int]]]:
    idx = {}
    if bed is None or bed.empty: return idx
    for chrom, sub in bed.groupby("chr", sort=False):
        starts = list(map(int, sub["start0"].tolist()))
        ends   = list(map(int, sub["end0"].tolist()))
        idx[chrom] = list(zip(starts, ends))  # assumed sorted from merge_bed
    return idx

def overlaps_exclude(excl_idx: Dict[str,List[Tuple[int,int]]], chrom: str, pos: int) -> bool:
    # BED is 0-based half-open; pos is 1-based center; match to [pos-1,pos)
    if chrom not in excl_idx: return False
    arr = excl_idx[chrom]
    i = bisect_left(arr, (pos-1, -1))
    for j in (i, i-1):
        if 0 <= j < len(arr):
            s0,e0 = arr[j]
            if s0 <= pos-1 < e0:
                return True
    return False

def build_positive_index(df: pd.DataFrame) -> Dict[str, List[int]]:
    idx = {}
    for chrom, sub in df.groupby("chr", sort=False):
        coords = sorted(list(map(int, sub["tss"].tolist())) + list(map(int, sub["tts"].tolist())))
        idx[chrom] = coords
    return idx

def within_keepout(pos_idx: Dict[str,List[int]], chrom: str, pos: int, keepout: int) -> bool:
    arr = pos_idx.get(chrom, [])
    if not arr: return False
    i = bisect_left(arr, pos)
    for j in (i, i-1):
        if 0 <= j < len(arr):
            if abs(arr[j] - pos) <= keepout:
                return True
    return False

def banded_offsets(bands: List[Tuple[int,int]], per_band: int, rng: random.Random) -> List[int]:
    offs = []
    for (a,b) in bands:
        for _ in range(per_band):
            offs.append(rng.randint(a, b))
    return offs

def main():
    ap = argparse.ArgumentParser(description="Generate banded negatives near positives and in quiet intergenic.")
    ap.add_argument("--positives", required=True, help="Positives TSV (chr tss tts gene strand)")
    ap.add_argument("--sep", default=r"\s+", help="Delimiter for positives (default: \\s+)")
    ap.add_argument("--no-header", action="store_true", help="Set if positives has no header")
    ap.add_argument("--col-order", default="chr,tss,tts,gene,strand", help="Order when --no-header")
    ap.add_argument("--chrom-sizes", required=True, help="chrom.sizes")
    ap.add_argument("--quiet-bed", default="", help="quiet intergenic BED (for easy negatives)")
    ap.add_argument("--exclude-bed", default="", help="extra exclude BED (e.g., telomeres)")
    ap.add_argument("--window", type=int, default=1000, help="Half-window to fit")
    ap.add_argument("--keepout", type=int, default=150, help="Avoid placing within +/- keepout of any positive")
    ap.add_argument("--near-bands", default="200-600,600-1500", help="Comma list of bands in bp from boundary")
    ap.add_argument("--per-band-per-side", type=int, default=1, help="# of negs per band per side per positive")
    ap.add_argument("--far-step", type=int, default=200, help="Step for far negatives in quiet BED")
    ap.add_argument("--seed", type=int, default=13)
    ap.add_argument("--out-neg", required=True, help="Output negatives TSV (headerless) chr tss tts gene strand")
    ap.add_argument("--out-meta", required=True, help="Output meta TSV (id, type, band, source, chr, pos)")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    sizes = read_chrom_sizes(args.chrom_sizes)
    pos = read_tsv(args.positives, args.sep, args.no_header, args.col_order)
    pos_idx = build_positive_index(pos)

    excl = None
    if args.exclude_bed:
        excl = merge_bed(read_bed(args.exclude_bed))
    excl_idx = build_exclude_index(excl)

    quiet = None
    if args.quiet_bed:
        quiet = merge_bed(read_bed(args.quiet_bed))

    bands = []
    for tok in args.near_bands.split(","):
        tok = tok.strip()
        if not tok: continue
        a,b = tok.split("-")
        bands.append((int(a), int(b)))

    neg_rows = []
    meta_rows = []

    # 1) Near-boundary (hard) negatives around both TSS and TTS
    for _, r in pos.iterrows():
        chrom = r["chr"]
        for c in (int(r["tss"]), int(r["tts"])):
            for side in (-1, +1):
                for off in banded_offsets(bands, args.per_band_per_side, rng):
                    cand = c + side * off
                    L = sizes.get(chrom, None)
                    if L is None: continue
                    if cand - args.window < 1 or cand + args.window > L:
                        continue
                    if within_keepout(pos_idx, chrom, cand, args.keepout):
                        continue
                    if overlaps_exclude(excl_idx, chrom, cand):
                        continue
                    strand = '+' if rng.random() < 0.5 else '-'
                    neg_id = f"NEG_NEAR_{chrom}_{c}_{cand}"
                    neg_rows.append((chrom, cand, cand, neg_id, strand))
                    meta_rows.append((neg_id, "near", f"{off}", "TSS/TTS", chrom, cand))

    # 2) Far (easy) negatives from quiet intergenic (if provided)
    if quiet is not None and not quiet.empty:
        for chrom, sub in quiet.groupby("chr", sort=False):
            for start0, end0 in zip(sub["start0"], sub["end0"]):
                s = int(start0)+1; e = int(end0)  # to 1-based inclusive
                s2 = s + args.window; e2 = e - args.window
                x = s2
                while x <= e2:
                    if not within_keepout(pos_idx, chrom, x, args.keepout) and not overlaps_exclude(excl_idx, chrom, x):
                        strand = '+' if rng.random() < 0.5 else '-'
                        neg_id = f"NEG_FAR_{chrom}_{x}"
                        neg_rows.append((chrom, x, x, neg_id, strand))
                        meta_rows.append((neg_id, "far", "na", "quiet_intergenic", chrom, x))
                    x += args.far_step

    with open(args.out_neg, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        for row in neg_rows:
            w.writerow(row)

    with open(args.out_meta, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["id","type","band","source","chr","pos"])
        for row in meta_rows:
            w.writerow(row)

    print(f"Wrote {len(neg_rows)} negatives to {args.out_neg}")
    print(f"Wrote meta to {args.out_meta}")

if __name__ == "__main__":
    main()
