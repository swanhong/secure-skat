#!/usr/bin/env bash

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/setup/env.sh"

config_path="${CONFIG_PATH:-run.aou.conf}"
reference_engine="${REFERENCE_ENGINE:-python}"

echo "[0/6] Prepare AoU chromosome inputs: ${chromosomes[*]}"
python3 rewrite/testdata/aou/prepare_aou.py \
  --config "${config_path}"

echo "[1/6] Prepare ancestry-specific secure inputs"
go run -mod=vendor secure-rvas.go prepare \
  --config "${config_path}" \
  --clear

echo "[2/6] Run secure Burden/SKAT"
go run -mod=vendor secure-rvas.go run \
  --config "${config_path}"

echo "[3/6] Run ancestry-specific ${reference_engine} reference"
python3 rewrite/analysis/run_reference.py \
  --config "${config_path}" \
  --engine "${reference_engine}"

echo "[4/6] Compare secure and reference results"
python3 rewrite/analysis/compare_secure_to_reference.py \
  --config "${config_path}"

echo "[5/6] Generate scatter and Manhattan plots"
python3 rewrite/analysis/plot_secure_vs_reference.py \
  --config "${config_path}"

echo "[6/6] Summarize metrics"
./summarize_metrics.sh "${config_path}"

echo "Secure RVAS AoU workflow completed"
