#!/bin/bash
#SBATCH --job-name=epi_feats
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G  
#SBATCH -p batch
#SBATCH --array=1-17        # <<< adjust number of chunks or override with: sbatch --array=1-20 ...

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-$(command -v python3 || true)}"

if [[ -z "${PYTHON_BIN}" ]]; then
  echo "python3 was not found on PATH. Set PYTHON_BIN to a working interpreter."
  exit 1
fi

PY_SCRIPT="${PY_SCRIPT:-${SCRIPT_DIR}/build_epigenetic_features_bins.py}"
BOUNDARIES="${BOUNDARIES:-${SCRIPT_DIR}/reference_gene_boundaries.csv}"   # headerless TSV: chr tss tts gene strand
DEFAULT_MARKER_JSON="${SCRIPT_DIR}/marker_files.local.json"
if [[ ! -f "${DEFAULT_MARKER_JSON}" ]]; then
  DEFAULT_MARKER_JSON="${SCRIPT_DIR}/marker_files.json"
fi
MARKER_JSON="${MARKER_JSON:-${DEFAULT_MARKER_JSON}}"
OUT_DIR="${OUT_DIR:-${SCRIPT_DIR}/positive_values}"
WHICH="BOTH"; WINDOW=1000; BIN_SIZE=100
GENE_COL="gene"; CHR_COL="chr"; STRAND_COL="strand"; TSS_COL="tss"; TTS_COL="tts"

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

# Convert this chunk from TSV (chr tss tts gene strand) -> CSV with header in the order gene,chr,strand,tss,tts
CHUNK_CSV="$SCRATCH_DIR/boundaries_chunk_${SLURM_ARRAY_TASK_ID}.csv"
{
  echo "gene,chr,strand,tss,tts"
  sed -n "${START},${END}p" "$BOUNDARIES" | \
  awk -v OFS=',' 'NF>=5 {print $4,$1,$5,$2,$3}'
} > "$CHUNK_CSV"

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
