#!/usr/bin/env bash

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/setup/env.sh"

config_path="${CONFIG_PATH:-config/aou}"
reference_engine="${REFERENCE_ENGINE:-python}"

echo "[0/7] Prepare AoU chromosome inputs"
python3 rewrite/testdata/aou/prepare_aou.py \
  --config "${config_path}"

echo "[1/7] Prepare ancestry-specific secure inputs"
go run -mod=vendor secure-rvas.go prepare \
  --config "${config_path}" \
  --clear

echo "[2/7] Generate shared PRG keys"
go run -mod=vendor secure-rvas.go keygen \
  --config "${config_path}"

echo "[3/7] Run secure Burden/SKAT"
go run -mod=vendor secure-rvas.go run \
  --config "${config_path}"

echo "[4/7] Run ancestry-specific ${reference_engine} reference"
python3 rewrite/analysis/run_reference.py \
  --config "${config_path}" \
  --engine "${reference_engine}"

echo "[5/7] Compare secure and reference results"
python3 rewrite/analysis/compare_secure_to_reference.py \
  --config "${config_path}"

echo "[6/7] Generate scatter and Manhattan plots"
python3 rewrite/analysis/plot_secure_vs_reference.py \
  --config "${config_path}"

echo "[7/7] Summarize metrics"
./summarize_metrics.sh "${config_path}"

echo "Secure RVAS AoU workflow completed"
