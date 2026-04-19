#!/bin/bash
#SBATCH --job-name=epi_feats_test
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH -p epyc
#SBATCH --array=1-1

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-$(command -v python3 || true)}"

if [[ -z "${PYTHON_BIN}" ]]; then
  echo "python3 was not found on PATH. Set PYTHON_BIN to a working interpreter."
  exit 1
fi

PY_SCRIPT="${PY_SCRIPT:-${SCRIPT_DIR}/build_epigenetic_features_bins_Neg.py}"
BOUNDARIES="${BOUNDARIES:-${SCRIPT_DIR}/neg_near20.tsv}"   # <-- use your 20-row file
MARKER_JSON="${MARKER_JSON:-${SCRIPT_DIR}/marker_files_test.json}"
OUT_DIR="${OUT_DIR:-${SCRIPT_DIR}/negative_values_test}"
WHICH="BOTH"; WINDOW=1000; BIN_SIZE=100

GENE_COL="gene"; CHR_COL="chr"; STRAND_COL="strand"
TSS_COL="TSS_coord"; TTS_COL="TTS_coord"

mkdir -p "$OUT_DIR" logs

ARRAY_MIN=${SLURM_ARRAY_TASK_MIN:-1}
ARRAY_MAX=${SLURM_ARRAY_TASK_MAX:-1}
ARRAY_COUNT=$((ARRAY_MAX - ARRAY_MIN + 1))

TOTAL_LINES=$(wc -l < "$BOUNDARIES")
DATA_LINES=$TOTAL_LINES
FIRST_DATA_LINE=1
(( DATA_LINES <= 0 )) && { echo "No data rows"; exit 1; }

PER_CHUNK=$(( (DATA_LINES + ARRAY_COUNT - 1) / ARRAY_COUNT ))
OFFSET=$((SLURM_ARRAY_TASK_ID - ARRAY_MIN))
START=$(( OFFSET * PER_CHUNK + FIRST_DATA_LINE ))
END=$(( START + PER_CHUNK - 1 ))
LAST=$(( FIRST_DATA_LINE + DATA_LINES - 1 ))
(( START > LAST )) && { echo "Task $SLURM_ARRAY_TASK_ID: nothing to do"; exit 0; }
(( END > LAST )) && END=$LAST

SCRATCH_DIR="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/epi_chunks_${SLURM_JOB_ID}"
mkdir -p "$SCRATCH_DIR"

CHUNK_CSV="$SCRATCH_DIR/boundaries_chunk_${SLURM_ARRAY_TASK_ID}.csv"
{
  echo "gene,chr,strand,TSS_coord,TTS_coord"
  sed -n "${START},${END}p" "$BOUNDARIES" | awk -v OFS=',' 'NF>=5 {print $4,$1,$5,$2,$3}'
} > "$CHUNK_CSV"

echo "[INFO] Chunk rows: $(($(wc -l < "$CHUNK_CSV") - 1))"
echo "[INFO] Chunk preview:"; head -n 3 "$CHUNK_CSV"

OUT_PREFIX="${OUT_DIR}/features_chunk_${SLURM_ARRAY_TASK_ID}"

"${PYTHON_BIN}" "$PY_SCRIPT" \
  --boundaries "$CHUNK_CSV" \
  --marker-config "$MARKER_JSON" \
  --out-prefix "$OUT_PREFIX" \
  --which "$WHICH" \
  --window "$WINDOW" \
  --bin-size "$BIN_SIZE" \
  --gene-col "$GENE_COL" \
  --chr-col "$CHR_COL" \
  --strand-col "$STRAND_COL" \
  --tss-col "$TSS_COL" \
  --tts-col "$TTS_COL"

echo "[INFO] Outputs in: $OUT_DIR"
ls -lh "${OUT_PREFIX}"*
