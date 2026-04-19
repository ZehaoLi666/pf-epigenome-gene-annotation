#!/bin/bash
#SBATCH --job-name=epi_bins
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH -p highmem
#SBATCH --array=1-1

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-$(command -v python3 || true)}"

if [[ -z "${PYTHON_BIN}" ]]; then
  echo "python3 was not found on PATH. Set PYTHON_BIN to a working interpreter."
  exit 1
fi

DEFAULT_MARKER_JSON="${SCRIPT_DIR}/marker_files.local.json"
if [[ ! -f "${DEFAULT_MARKER_JSON}" ]]; then
  DEFAULT_MARKER_JSON="${SCRIPT_DIR}/marker_files.json"
fi

"${PYTHON_BIN}" "${SCRIPT_DIR}/build_epigenetic_features_bins_Neg.py" \
  --boundaries "${BOUNDARIES:-${SCRIPT_DIR}/negatives_2k.tsv}" \
  --no-header \
  --sep '\s+' \
  --marker-config "${MARKER_JSON:-${DEFAULT_MARKER_JSON}}" \
  --out-csv "${OUT_CSV:-${SCRIPT_DIR}/ATAC_negatives_2k_features.csv}" \
  --window 1000 \
  --bin-size 100
