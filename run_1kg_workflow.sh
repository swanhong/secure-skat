#!/usr/bin/env bash

set -euo pipefail

config_path="${CONFIG_PATH:-run.1kg.conf}"
prepare_args=(--config "${config_path}")
if [[ "${CLEAR_RUN_DIR:-0}" == "1" ]]; then
  prepare_args+=(--clear)
fi

echo "[0/5] Generate 1000 Genomes test data"
python3 rewrite/testdata/1kgenome/prepare_1kgenome.py \
  --chromosome 21 22 \
  --num-pheno 2

echo "[1/5] Prepare ancestry-specific secure inputs"
go run -mod=vendor secure-rvas.go prepare "${prepare_args[@]}"

echo "[2/5] Run secure Burden/SKAT"
go run -mod=vendor secure-rvas.go run \
  --config "${config_path}"

echo "[3/5] Run ancestry-specific R::SKAT reference"
python3 rewrite/analysis/run_reference.py \
  --config "${config_path}"

echo "[4/5] Compare secure and R results"
python3 rewrite/analysis/compare_results.py \
  --config "${config_path}"

echo "[5/5] Generate scatter and Manhattan plots"
python3 rewrite/analysis/plot_results.py \
  --config "${config_path}"

echo "Secure RVAS 1KG workflow completed"
