#!/usr/bin/env bash

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/setup/env.sh"

config_path="${CONFIG_PATH:-config/1kg}"
reference_engine="${REFERENCE_ENGINE:-python}"
prepare_args=(--config "${config_path}")
if [[ "${CLEAR_RUN_DIR:-0}" == "1" ]]; then
  prepare_args+=(--clear)
fi

echo "[0/7] Generate 1000 Genomes test data"
python3 rewrite/testdata/1kgenome/prepare_1kgenome.py \
  --config "${config_path}" \
  --num-pheno 2

echo "[1/7] Prepare ancestry-specific secure inputs"
go run -mod=vendor secure-rvas.go prepare "${prepare_args[@]}"

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

echo "Secure RVAS 1KG workflow completed"
