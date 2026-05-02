#!/bin/bash
#SBATCH --job-name=secure-skat-1kg
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=08:00:00

set -euo pipefail

# Update this path after copying the repo to Hopper.
REPO_ROOT="${REPO_ROOT:-$HOME/secure-skat}"

# 1KG rare50k dataset/config defaults created in this repo.
DATASET="${DATASET:-.local/datasets/1000g_afr_chr22_rare50k}"
CONFIG_DIR="${CONFIG_DIR:-.local/config_1kg_afr_chr22_rare50k}"
MODE="${MODE:-skat}"
COMPARE_BLOCKS="${COMPARE_BLOCKS:-all}"
COMPARE_SKIP_REFERENCE="${COMPARE_SKIP_REFERENCE:-0}"

# Python environment for preprocessing/export helpers.
CONDA_MODULE="${CONDA_MODULE:-miniconda}"
CONDA_ENV_NAME="${CONDA_ENV_NAME:-skat_env}"

# Optional R module name. Leave empty if Rscript is already on PATH.
R_MODULE="${R_MODULE:-}"

cd "$REPO_ROOT"

module restore skat-mod
module load "$CONDA_MODULE"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV_NAME"

if [[ -n "$R_MODULE" ]]; then
  module load "$R_MODULE"
fi

if ! command -v go >/dev/null 2>&1; then
  echo "go not found in PATH" >&2
  exit 1
fi

if ! command -v plink2 >/dev/null 2>&1; then
  echo "plink2 not found in PATH" >&2
  exit 1
fi

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript not found in PATH" >&2
  exit 1
fi

python - <<'PY'
import numpy
print("numpy ok", numpy.__version__)
PY

RUN_EXAMPLE_LOG="$(mktemp "${TMPDIR:-/tmp}/secure-skat-run.XXXXXX.log")"
trap 'rm -f "$RUN_EXAMPLE_LOG"' EXIT

echo "Submitting secure run from: $REPO_ROOT"
echo "Dataset: $DATASET"
echo "Config directory: $CONFIG_DIR"
echo "Mode: $MODE"

bash run_example.sh \
  --mode "$MODE" \
  --dataset "$DATASET" \
  --config-dir "$CONFIG_DIR" | tee "$RUN_EXAMPLE_LOG"

RUN_ID="$(awk '/^Run ID: / {print $3; exit}' "$RUN_EXAMPLE_LOG")"
RUN_ROOT="$(awk -F': ' '/^Run output directory: / {print $2; exit}' "$RUN_EXAMPLE_LOG")"

if [[ -z "$RUN_ID" || -z "$RUN_ROOT" ]]; then
  echo "Failed to parse run metadata from run_example.sh output" >&2
  exit 1
fi

COMPARE_LOG="${RUN_ROOT}/skat_compare.txt"

COMPARE_CMD=(
  python
  scripts/analysis/skat_compare.py
  compare
  --repo-root .
  --run-id "$RUN_ID"
  --dataset "$DATASET"
  --blocks "$COMPARE_BLOCKS"
)

if [[ "$COMPARE_SKIP_REFERENCE" == "1" ]]; then
  COMPARE_CMD+=(--skip-reference)
fi

echo "Running Python SKAT comparison for run id: $RUN_ID"
"${COMPARE_CMD[@]}" | tee "$COMPARE_LOG"

echo "Secure run root: $RUN_ROOT"
echo "Comparison log: $COMPARE_LOG"
