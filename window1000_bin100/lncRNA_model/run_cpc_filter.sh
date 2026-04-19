#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 ]]; then
  echo "Usage: $0 <bed> <genome.fa> <out_prefix> <threads> [--use-cpc1]"
  echo "Example:"
  echo "  $0 novel_units.sorted.bed PlasmoDB-48_Pfalciparum3D7_Genome.fasta out/lncCand 4"
  exit 1
fi

BED=$1
FA=$2
OUT=$3
THREADS=$4
USE_CPC1=${5:-""}

mkdir -p "$(dirname "$OUT")"

# 1) Make FASTA (strand-aware; headers = BED 'name')
FASTA="${OUT}.fa"
bedtools getfasta -s -name -fi "$FA" -bed "$BED" -fo "$FASTA"

# 2) Run CPC2 or CPC1
if [[ -z "$USE_CPC1" ]]; then
  # Try CPC2
  if command -v CPC2.py >/dev/null 2>&1; then
    CPC_OUT="${OUT}.cpc2.txt"
    CPC2.py -i "$FASTA" -o "$CPC_OUT" -t "$THREADS"
    PARSER="cpc2"
  else
    echo "CPC2.py not found. To force CPC1, re-run with --use-cpc1"
    exit 2
  fi
else
  # CPC1 fallback (webserver/old): expects --species=non-vertebrate for Plasmodium
  if command -v CPC.py >/dev/null 2>&1; then
    CPC_OUT="${OUT}.cpc1.txt"
    # Note: adjust --species as appropriate; many use 'non-vertebrate' for Plasmodium
    CPC.py -i "$FASTA" -o "$CPC_OUT" --species=non-vertebrate --threads "$THREADS"
    PARSER="cpc1"
  else
    echo "CPC.py not found."
    exit 2
  fi
fi

# 3) Merge CPC output back to BED and split
python /mnt/data/merge_cpc_to_bed.py "$BED" "$CPC_OUT" "$OUT" "$PARSER"

echo "Done."
