#!/usr/bin/env python3
import argparse
import os
import re
from pathlib import Path

import pandas as pd


def collect_files(input_dir: Path, kind: str):
    pattern = re.compile(rf"^features_chunk_(\d+)_{kind}_.*?bins\.csv$")
    matches = []
    for path in input_dir.iterdir():
        if not path.is_file():
            continue
        hit = pattern.match(path.name)
        if hit:
            matches.append((int(hit.group(1)), path))
    matches.sort(key=lambda item: item[0])
    return [path for _, path in matches]


def merge_kind(input_dir: Path, kind: str, out_path: Path):
    files = collect_files(input_dir, kind)
    if not files:
        raise FileNotFoundError(f"No chunk files found for {kind} in {input_dir}")

    columns = None
    parts = []
    for path in files:
        frame = pd.read_csv(path)
        current_columns = list(frame.columns)
        if columns is None:
            columns = current_columns
        elif current_columns != columns:
            raise ValueError(f"Column mismatch in {path}")
        parts.append(frame)

    merged = pd.concat(parts, ignore_index=True)
    merged.to_csv(out_path, index=False)
    print(f"Wrote {out_path} with {len(merged):,} rows from {len(files)} chunk files.")


def main():
    parser = argparse.ArgumentParser(description="Merge chunked positive feature tables.")
    parser.add_argument("--input-dir", default="positive_values", help="Directory containing features_chunk_* CSVs")
    parser.add_argument(
        "--kind",
        choices=["TSS", "TTS", "BOTH"],
        default="BOTH",
        help="Which boundary table(s) to merge"
    )
    parser.add_argument("--out-tss", default="", help="Output CSV for merged TSS features")
    parser.add_argument("--out-tts", default="", help="Output CSV for merged TTS features")
    args = parser.parse_args()

    input_dir = Path(args.input_dir).resolve()
    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")

    out_tss = Path(args.out_tss) if args.out_tss else input_dir / "merged_TSS_500bp_20bins.csv"
    out_tts = Path(args.out_tts) if args.out_tts else input_dir / "merged_TTS_500bp_20bins.csv"

    if args.kind in ("TSS", "BOTH"):
        merge_kind(input_dir, "TSS", out_tss)
    if args.kind in ("TTS", "BOTH"):
        merge_kind(input_dir, "TTS", out_tts)


if __name__ == "__main__":
    main()
