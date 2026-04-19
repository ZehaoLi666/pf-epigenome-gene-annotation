#!/usr/bin/env python3
import argparse, sys, csv, random, os
from typing import List, Tuple, Dict, Optional, Set
import pandas as pd
from bisect import bisect_left

def eprint(*a, **k): print(*a, file=sys.stderr, **k)

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
    df = df[need].dropna()
    return df

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

def overlaps_exclude(excl_idx: Dict[str,List[Tuple[int,int]]], chrom: str, pos_1based: int) -> bool:
    # BED is 0-based half-open; pos is 1-based center; match to [pos-1,pos)
    if chrom not in excl_idx: return False
    arr = excl_idx[chrom]
    i = bisect_left(arr, (pos_1based-1, -1))
    for j in (i, i-1):
        if 0 <= j < len(arr):
            s0,e0 = arr[j]
            if s0 <= pos_1based-1 < e0:
                return True
    return False

def build_positive_index(df: pd.DataFrame) -> Dict[str, List[int]]:
    idx = {}
    for chrom, sub in df.groupby("chr", sort=False):
        coords = sorted(list(map(int, sub["tss"].tolist())) + list(map(int, sub["tts"].tolist())))
        idx[chrom] = coords
    return idx

def within_keepout(pos_idx: Dict[str,List[int]], chrom: str, pos_1based: int, keepout: int) -> bool:
    arr = pos_idx.get(chrom, [])
    if not arr: return False
    i = bisect_left(arr, pos_1based)
    for j in (i, i-1):
        if 0 <= j < len(arr):
            if abs(arr[j] - pos_1based) <= keepout:
                return True
    return False

def banded_offsets(bands: List[Tuple[int,int]], per_band: int, rng: random.Random) -> List[int]:
    offs = []
    for (a,b) in bands:
        if b < a: a, b = b, a
        for _ in range(per_band):
            offs.append(rng.randint(a, b))
    return offs

def stratified_sample(ids: List[str], groups: List[str], k: int, rng: random.Random) -> List[str]:
    df = pd.DataFrame({"id": ids, "grp": groups})
    counts = df["grp"].value_counts()
    if counts.sum() <= k:
        return ids
    alloc = (counts / counts.sum() * k).round().astype(int)
    diff = k - int(alloc.sum())
    if diff != 0:
        order = counts.sort_values(ascending=(diff<0)).index.tolist()
        for g in order:
            if diff == 0: break
            alloc[g] += 1 if diff>0 else -1
            diff += -1 if diff>0 else 1
    chosen = []
    for g, want in alloc.items():
        pool = df.loc[df["grp"]==g, "id"].tolist()
        if want >= len(pool):
            chosen.extend(pool)
        elif want > 0:
            chosen.extend(rng.sample(pool, want))
    if len(chosen) > k:
        chosen = rng.sample(chosen, k)
    elif len(chosen) < k:
        remaining = [i for i in ids if i not in set(chosen)]
        need = k - len(chosen)
        if need > 0 and remaining:
            chosen.extend(rng.sample(remaining, min(need, len(remaining))))
    return chosen

def main():
    ap = argparse.ArgumentParser(description="Generate banded negatives near positives and in quiet intergenic, with optional cap.")
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
    ap.add_argument("--target-neg", type=int, default=0, help="If >0, downsample to this many negatives total")
    ap.add_argument("--near-fraction", type=float, default=0.5, help="Fraction of target negatives to take from 'near' (rest from 'far')")
    ap.add_argument("--near-only", action="store_true", help="Generate only near negatives (overrides near-fraction)")
    ap.add_argument("--far-only", action="store_true", help="Generate only far negatives (overrides near-fraction)")
    ap.add_argument("--out-neg", required=True, help="Output negatives TSV (headerless) chr tss tts gene strand")
    ap.add_argument("--out-meta", required=True, help="Output meta TSV (id, type, band, source, chr, pos)")
    args = ap.parse_args()

    if args.near_only and args.far_only:
        sys.exit("Use at most one of --near-only or --far-only")

    rng = random.Random(args.seed)

    # Load inputs
    if not os.path.exists(args.positives): sys.exit(f"Not found: {args.positives}")
    if not os.path.exists(args.chrom_sizes): sys.exit(f"Not found: {args.chrom_sizes}")
    sizes = read_chrom_sizes(args.chrom_sizes)
    pos = read_tsv(args.positives, args.sep, args.no_header, args.col_order)

    eprint(f"positives: {len(pos)} rows; chroms in positives: {pos['chr'].nunique()}")

    # Chrom name agreement
    pos_chroms = set(pos["chr"].unique().tolist())
    size_chroms = set(sizes.keys())
    missing_in_sizes = sorted(list(pos_chroms - size_chroms))
    if missing_in_sizes:
        sys.exit(f"Chromosomes in positives missing from sizes: {missing_in_sizes}")

    # Optional beds
    def _merge_opt(path: str) -> Optional[pd.DataFrame]:
        if not path: return None
        if not os.path.exists(path):
            eprint(f"WARNING: BED not found: {path} (ignored)")
            return None
        df = read_bed(path)
        if df is None or df.empty:
            eprint(f"WARNING: BED empty: {path} (ignored)")
            return None
        return merge_bed(df)

    excl = _merge_opt(args.exclude_bed)
    excl_idx = build_exclude_index(excl)
    quiet = _merge_opt(args.quiet_bed)
    if args.far_only and (quiet is None or quiet.empty):
        sys.exit("--far-only requested but quiet BED is missing or empty")

    # parse bands
    bands = []
    for tok in args.near_bands.split(","):
        tok = tok.strip()
        if not tok: continue
        a,b = tok.split("-")
        bands.append((int(a), int(b)))
    if not bands and (not args.far_only):
        sys.exit("No valid --near-bands parsed")

    pos_idx = build_positive_index(pos)

    # Generate near candidates
    near_rows: List[Tuple[str,int,int,str,str]] = []
    near_meta: List[Tuple[str,str,str,str,str,int]] = []
    if not args.far_only:
        seen_pos: Set[Tuple[str,int]] = set()
        for _, r in pos.iterrows():
            chrom = r["chr"]
            for c in (int(r["tss"]), int(r["tts"])):
                for side in (-1, +1):
                    for off in banded_offsets(bands, args.per_band_per_side, rng):
                        cand = c + side * off
                        L = sizes[chrom]
                        if cand - args.window < 1 or cand + args.window > L:
                            continue
                        if within_keepout(pos_idx, chrom, cand, args.keepout):
                            continue
                        if overlaps_exclude(excl_idx, chrom, cand):
                            continue
                        if (chrom, cand) in seen_pos:
                            continue
                        seen_pos.add((chrom, cand))
                        strand = '+' if rng.random() < 0.5 else '-'
                        neg_id = f"NEG_NEAR_{chrom}_{c}_{cand}"
                        near_rows.append((chrom, cand, cand, neg_id, strand))
                        # id, type, band, source, chr, pos
                        near_meta.append((neg_id, "near", f"{off}", "TSS/TTS", chrom, cand))

    # Generate far candidates
    far_rows: List[Tuple[str,int,int,str,str]] = []
    far_meta: List[Tuple[str,str,str,str,str,int]] = []
    if (quiet is not None) and (not quiet.empty) and (not args.near_only):
        seen_far: Set[Tuple[str,int]] = set()
        for chrom, sub in quiet.groupby("chr", sort=False):
            L = sizes.get(chrom)
            if L is None:  # chrom not in sizes
                continue
            for start0, end0 in zip(sub["start0"], sub["end0"]):
                s = int(start0)+1; e = int(end0)  # to 1-based inclusive
                s2 = s + args.window; e2 = e - args.window
                x = s2
                while x <= e2:
                    if not within_keepout(pos_idx, chrom, x, args.keepout) and not overlaps_exclude(excl_idx, chrom, x):
                        if (chrom, x) not in seen_far:
                            seen_far.add((chrom, x))
                            strand = '+' if rng.random() < 0.5 else '-'
                            neg_id = f"NEG_FAR_{chrom}_{x}"
                            far_rows.append((chrom, x, x, neg_id, strand))
                            far_meta.append((neg_id, "far", "na", "quiet_intergenic", chrom, x))
                    x += args.far_step

    eprint(f"candidates before cap: near={len(near_rows)}, far={len(far_rows)}")

    # Combine & optionally downsample
    all_rows = near_rows + far_rows
    meta_df = pd.DataFrame(near_meta + far_meta, columns=["id","type","band","source","chr","pos"])

    if args.target_neg and len(all_rows) > args.target_neg:
        near_ids = meta_df.loc[meta_df["type"]=="near","id"].tolist()
        far_ids  = meta_df.loc[meta_df["type"]=="far","id"].tolist()

        if args.near_only:
            want_near, want_far = args.target_neg, 0
        elif args.far_only:
            want_near, want_far = 0, args.target_neg
        else:
            want_near = int(round(args.target_neg * args.near_fraction))
            want_far  = args.target_neg - want_near

        rng.shuffle(near_ids); rng.shuffle(far_ids)

        if want_near > 0 and len(near_ids) > want_near:
            sub = meta_df[meta_df["type"]=="near"][["id","band","chr"]].copy()
            sub["grp"] = sub["band"].astype(str) + "|" + sub["chr"].astype(str)
            chosen_near = stratified_sample(sub["id"].tolist(), sub["grp"].tolist(), want_near, rng)
        else:
            chosen_near = near_ids[:want_near]

        if want_far > 0 and len(far_ids) > want_far:
            sub = meta_df[meta_df["type"]=="far"][["id","chr"]].copy()
            sub["grp"] = sub["chr"].astype(str)
            chosen_far = stratified_sample(sub["id"].tolist(), sub["grp"].tolist(), want_far, rng)
        else:
            chosen_far = far_ids[:want_far]

        chosen: Set[str] = set(chosen_near + chosen_far)
        # Top-up from remaining if needed
        if len(chosen) < args.target_neg:
            remaining = [i for i in meta_df["id"] if i not in chosen]
            need = args.target_neg - len(chosen)
            if need > 0 and remaining:
                chosen.update(rng.sample(remaining, min(need, len(remaining))))

        # Filter and sort
        id_to_row = {row[3]: row for row in all_rows}
        all_rows = [id_to_row[i] for i in chosen if i in id_to_row]
        meta_df = meta_df[meta_df["id"].isin(chosen)].sort_values(["type","chr","pos"]).reset_index(drop=True)

    # Write outputs
    with open(args.out_neg, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        for row in all_rows:
            w.writerow(row)

    meta_df.to_csv(args.out_meta, sep="\t", index=False)

    print(f"Wrote {len(all_rows)} negatives to {args.out_neg}")
    print(f"(near={sum(meta_df['type'].eq('near'))}, far={sum(meta_df['type'].eq('far'))})")

if __name__ == "__main__":
    main()
