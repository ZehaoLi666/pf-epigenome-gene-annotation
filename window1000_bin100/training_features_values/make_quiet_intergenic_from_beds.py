
#!/usr/bin/env python3
import argparse, sys, csv, random
from typing import List, Tuple, Dict
import pandas as pd

# BED is assumed 0-based, half-open [start,end). We'll convert to 1-based inclusive for arithmetic.

def read_bed(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", header=None, comment="#", engine="python")
    if df.shape[1] < 3:
        sys.exit(f"{path}: expected at least 3 BED columns (chr, start, end)")
    df = df.iloc[:, :3]
    df.columns = ["chr", "start0", "end0"]
    # to 1-based inclusive
    df["start"] = df["start0"].astype(int) + 1
    df["end"]   = df["end0"].astype(int)
    df = df[["chr","start","end"]]
    return df

def merge_intervals(df: pd.DataFrame) -> pd.DataFrame:
    out_rows = []
    for chrom, sub in df.groupby("chr", sort=False):
        sub = sub.sort_values(["start","end"])
        cur_s, cur_e = None, None
        for s,e in zip(sub["start"], sub["end"]):
            s = int(s); e = int(e)
            if cur_s is None:
                cur_s, cur_e = s, e
            elif s <= cur_e + 1:
                if e > cur_e: cur_e = e
            else:
                out_rows.append((chrom, cur_s, cur_e))
                cur_s, cur_e = s, e
        if cur_s is not None:
            out_rows.append((chrom, cur_s, cur_e))
    return pd.DataFrame(out_rows, columns=["chr","start","end"])

def subtract_intervals(size: int, merged_exclude: List[Tuple[int,int]]) -> List[Tuple[int,int]]:
    """From [1,size], subtract merged_exclude intervals (1-based inclusive). Return remaining disjoint intervals."""
    keep = []
    cur = 1
    for s,e in merged_exclude:
        s = max(1, s); e = min(size, e)
        if s > e:
            continue
        if cur < s:
            keep.append((cur, s-1))
        cur = max(cur, e+1)
        if cur > size:
            break
    if cur <= size:
        keep.append((cur, size))
    return keep

def shrink_by_margin(intervals: List[Tuple[int,int]], margin: int) -> List[Tuple[int,int]]:
    out = []
    for s,e in intervals:
        s2 = s + margin
        e2 = e - margin
        if e2 >= s2:
            out.append((s2, e2))
    return out

def sample_centers(intervals: List[Tuple[int,int]], window: int, step: int) -> List[int]:
    centers = []
    for s,e in intervals:
        # ensure ±window fits
        s2 = s + window
        e2 = e - window
        x = s2
        while x <= e2:
            centers.append(x)
            x += step
    return centers

def read_chrom_sizes(path: str) -> Dict[str,int]:
    sizes = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            chrom, size = line.split()[:2]
            sizes[chrom] = int(size)
    return sizes

def main():
    ap = argparse.ArgumentParser(description="Quiet intergenic = chromosome - telomeres - merged.bed; sample negatives.")
    ap.add_argument("--chrom-sizes", required=True, help="chrom.sizes (chr<tab>size)")
    ap.add_argument("--telomere-bed", required=True, help="BED of telomere regions (0-based, half-open)")
    ap.add_argument("--merged-bed", required=True, help="BED of all genes/transcripts merged (0-based, half-open)")
    ap.add_argument("--margin", type=int, default=1000, help="Shrink remaining intergenic by this many bp (default 1000)")
    ap.add_argument("--window", type=int, default=1000, help="Half-window for later features (default 1000)")
    ap.add_argument("--step", type=int, default=200, help="Spacing between negative centers (default 200)")
    ap.add_argument("--quiet-bed-out", required=True, help="Output BED (0-based half-open) of quiet intergenic regions")
    ap.add_argument("--neg-tsv-out", required=True, help="Output TSV (headerless) of negatives: chr tss tts gene strand")
    ap.add_argument("--strand-seed", type=int, default=13, help="Seed for random strand assignment")
    args = ap.parse_args()

    sizes = read_chrom_sizes(args.chrom_sizes)

    tel = merge_intervals(read_bed(args.telomere_bed))
    gen = merge_intervals(read_bed(args.merged_bed))

    random.seed(args.strand_seed)

    quiet_rows = []
    neg_rows = []
    for chrom, size in sizes.items():
        tel_c = tel[tel["chr"] == chrom][["start","end"]].values.tolist()
        gen_c = gen[gen["chr"] == chrom][["start","end"]].values.tolist()
        excl = sorted(tel_c + gen_c, key=lambda x: (x[0], x[1]))
        # merge exclusion
        merged = []
        for s,e in excl:
            if not merged or s > merged[-1][1] + 1:
                merged.append([s,e])
            else:
                merged[-1][1] = max(merged[-1][1], e)
        # subtract from chromosome
        keep = subtract_intervals(size, [(s,e) for s,e in merged])
        # shrink by margin
        quiet = shrink_by_margin(keep, args.margin)
        # write quiet bed (convert back to 0-based half-open)
        for s,e in quiet:
            quiet_rows.append((chrom, s-1, e))  # [s-1, e)
        # sample centers
        centers = sample_centers(quiet, args.window, args.step)
        for c in centers:
            strand = '+' if random.random() < 0.5 else '-'
            name = f"NEG_{chrom}_{c}"
            neg_rows.append((chrom, c, c, name, strand))

    qdf = pd.DataFrame(quiet_rows, columns=["chr","start0","end0"])
    qdf.to_csv(args.quiet_bed_out, sep="\t", header=False, index=False)

    with open(args.neg_tsv_out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        for row in neg_rows:
            w.writerow(row)

    print(f"Wrote {len(quiet_rows)} quiet intervals to {args.quiet_bed_out}")
    print(f"Wrote {len(neg_rows)} negatives to {args.neg_tsv_out}")

if __name__ == "__main__":
    main()
