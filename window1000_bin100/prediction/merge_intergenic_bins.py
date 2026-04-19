#!/usr/bin/env python3
import argparse
from functools import reduce
from pathlib import Path

import pandas as pd


KEY_COLS = ["chr", "mid", "ID", "strand"]


def load_one(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    missing = [col for col in KEY_COLS if col not in frame.columns]
    if missing:
        raise ValueError(f"{path.name} is missing key column(s): {missing}")
    for column in frame.columns:
        if column not in KEY_COLS:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame


def main():
    parser = argparse.ArgumentParser(description="Merge per-marker intergenic bin tables into one feature matrix.")
    parser.add_argument("--input-dir", default="intergenic_bins", help="Directory with *_intergenic_bins_*.csv files")
    parser.add_argument("--out-csv", default="", help="Output CSV path")
    args = parser.parse_args()

    input_dir = Path(args.input_dir).resolve()
    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")

    files = sorted(input_dir.glob("*_intergenic_bins_pm*bp_*bins.csv"))
    if not files:
        raise FileNotFoundError(f"No intergenic bin CSVs found in {input_dir}")

    frames = [load_one(path) for path in files]
    merged = reduce(lambda left, right: pd.merge(left, right, on=KEY_COLS, how="outer"), frames)
    for column in merged.columns:
        if column not in KEY_COLS:
            merged[column] = pd.to_numeric(merged[column], errors="coerce").fillna(0.0)

    out_csv = Path(args.out_csv) if args.out_csv else input_dir.parent / "merged_intergenic_features.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv} with shape {merged.shape}")


if __name__ == "__main__":
    main()
