#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-$(command -v python3 || true)}"

if [[ -z "${PYTHON_BIN}" ]]; then
  echo "python3 was not found on PATH. Set PYTHON_BIN to a working interpreter."
  exit 1
fi

"${PYTHON_BIN}" "${SCRIPT_DIR}/filter_TSS.py" \
  --genes "${GENES_BED:-${SCRIPT_DIR}/PlasmoDB-48_Pfalciparum3D7_Genes.bed}" \
  --tss "${TSS_INPUT:-${SCRIPT_DIR}/predictions_all.csv}" \
  --chrom-sizes "${CHROM_SIZES:?Set CHROM_SIZES to a genome .fai file}" \
  --score-min "${SCORE_MIN:-0.7}" \
  --out-prefix "${OUT_PREFIX:-${SCRIPT_DIR}/refined_TSS_intergenic}"
