#!/bin/bash -l
#SBATCH --job-name=tss_xgb_cpu
#SBATCH --output=logs/%x_%A.out
#SBATCH --error=logs/%x_%A.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH -p epyc

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-$(command -v python3 || true)}"

if [[ -z "${PYTHON_BIN}" ]]; then
  echo "python3 was not found on PATH. Set PYTHON_BIN to a working interpreter."
  exit 1
fi

mkdir -p "${SCRIPT_DIR}/logs"

POS="${POS:-${REPO_ROOT}/window1000_bin100/training_features_values/positive_values/merged_TSS_500bp_20bins.csv}"
NEG="${NEG:-${REPO_ROOT}/window1000_bin100/training_features_values/negative_values/ATAC_negatives_2k_features_TSS_only.csv}"
OUT="${OUT:-${SCRIPT_DIR}/tss_xgb_out_cpu}"

# Optional sanity check
"${PYTHON_BIN}" - <<'PY'
import sys, sklearn, xgboost
print("Python:", sys.executable)
print("sklearn:", sklearn.__version__)
print("xgboost:", xgboost.__version__)
PY

mkdir -p "$OUT"

"${PYTHON_BIN}" "${SCRIPT_DIR}/model.training.py" \
  --pos "$POS" \
  --neg "$NEG" \
  --out_dir "$OUT" \
  --n_jobs "$SLURM_CPUS_PER_TASK" \
  --device cpu
