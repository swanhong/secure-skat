#!/usr/bin/env bash

set -euo pipefail

source scripts/setup/env.sh

config_path="${CONFIG_PATH:-config/aou}"
reference_engine="${REFERENCE_ENGINE:-python}"
TIMEFORMAT='  wall time: %R seconds'

echo "[0/7] Prepare AoU chromosome inputs"
time python3 python/datasets/aou/prepare_aou.py \
  --config "${config_path}"

echo "[1/7] Prepare ancestry-specific secure inputs"
time go run -mod=vendor secure-rvas.go prepare \
  --config "${config_path}" \

echo "[2/7] Generate shared PRG keys"
time go run -mod=vendor secure-rvas.go keygen \
  --config "${config_path}"

echo "[3/7] Run secure Burden/SKAT"
time go run -mod=vendor secure-rvas.go run \
  --config "${config_path}"

echo "[4/7] Run ancestry-specific ${reference_engine} reference"
time python3 python/analysis/run_reference.py \
  --config "${config_path}" \
  --engine "${reference_engine}"

echo "[5/7] Compare secure and reference results"
time python3 python/analysis/compare_secure_to_reference.py \
  --config "${config_path}"

echo "[6/7] Generate scatter and Manhattan plots"
time python3 python/analysis/plot_secure_vs_reference.py \
  --config "${config_path}"

echo "[7/7] Summarize metrics"
time ./scripts/summarize_metrics.sh "${config_path}"

echo "Secure RVAS AoU workflow completed"
