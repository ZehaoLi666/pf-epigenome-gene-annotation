 #!/usr/bin/env python3
import argparse, re, os
import pandas as pd

def main():
    ap = argparse.ArgumentParser(description="Split combined TSS+TTS feature CSV into TSS-only and TTS-only files.")
    ap.add_argument("--in-csv", required=True, help="Input CSV with meta + *_TSS_binXX and *_TTS_binXX columns")
    ap.add_argument("--out-prefix", default="", help="Prefix for outputs (default: input path without .csv)")
    args = ap.parse_args()

    in_csv = args.in_csv
    out_prefix = args.out_prefix or os.path.splitext(in_csv)[0]

    df = pd.read_csv(in_csv)

    # meta columns if present
    meta_candidates = ["gene","chr","strand","TSS_coord","TTS_coord"]
    meta_cols = [c for c in meta_candidates if c in df.columns]

    # feature columns (keep original order)
    tss_cols = [c for c in df.columns if re.search(r"_TSS_bin\d+$", c)]
    tts_cols = [c for c in df.columns if re.search(r"_TTS_bin\d+$", c)]

    if not tss_cols and not tts_cols:
        raise SystemExit("No *_TSS_binXX or *_TTS_binXX columns found.")

    # Build final column orders (drop the other boundary’s coordinate if present)
    tss_keep = [c for c in meta_cols if c != "TTS_coord"] + tss_cols
    tts_keep = [c for c in meta_cols if c != "TSS_coord"] + tts_cols

    # Subset and write
    if tss_cols:
        tss_df = df[tss_keep].copy()
        tss_out = f"{out_prefix}_TSS_only.csv"
        tss_df.to_csv(tss_out, index=False)
        print(f"Wrote {tss_out}  (rows={len(tss_df)}, cols={len(tss_df.columns)})")
    else:
        print("No TSS columns found; skipped TSS file.")

    if tts_cols:
        tts_df = df[tts_keep].copy()
        tts_out = f"{out_prefix}_TTS_only.csv"
        tts_df.to_csv(tts_out, index=False)
        print(f"Wrote {tts_out}  (rows={len(tts_df)}, cols={len(tts_df.columns)})")
    else:
        print("No TTS columns found; skipped TTS file.")

if __name__ == "__main__":
    main()
  