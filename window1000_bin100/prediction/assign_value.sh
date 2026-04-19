#!/bin/bash
#SBATCH --job-name=midbins
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH -p epyc
#SBATCH --array=1-1   # <<< one task per marker; adjust to number of markers in JSON

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-$(command -v python3 || true)}"

if [[ -z "${PYTHON_BIN}" ]]; then
  echo "python3 was not found on PATH. Set PYTHON_BIN to a working interpreter."
  exit 1
fi

# ---- Paths ----
MIDPOINTS="${MIDPOINTS:-${SCRIPT_DIR}/intergenic_midpoints_100bp.csv}"
DEFAULT_MARKER_JSON="${REPO_ROOT}/window1000_bin100/training_features_values/marker_files_ATAC.local.json"
if [[ ! -f "${DEFAULT_MARKER_JSON}" ]]; then
  DEFAULT_MARKER_JSON="${REPO_ROOT}/window1000_bin100/training_features_values/marker_files_ATAC.json"
fi
MARKER_JSON="${MARKER_JSON:-${DEFAULT_MARKER_JSON}}"   # {marker: path or [paths]}
OUT_DIR="${OUT_DIR:-${SCRIPT_DIR}/intergenic_bins}"
PY_SCRIPT="${PY_SCRIPT:-${SCRIPT_DIR}/build_bins_for_intergenic_midpoints.py}"

# ---- Windowing ----
WINDOW=1000      # ±1000 bp
BIN_SIZE=100     # 100-bp bins -> 20 bins

mkdir -p "$(dirname logs/dummy)" "$OUT_DIR"

# If you need a conda env, activate it before submitting the job.

# SLURM_ARRAY_TASK_ID selects the marker index from the JSON (1-based)
"${PYTHON_BIN}" "$PY_SCRIPT" \
  --midpoints "$MIDPOINTS" \
  --marker-config "$MARKER_JSON" \
  --marker-index "${SLURM_ARRAY_TASK_ID}" \
  --out-dir "$OUT_DIR" \
  --window "$WINDOW" \
  --bin-size "$BIN_SIZE" \
