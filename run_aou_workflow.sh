#!/usr/bin/env bash

set -euo pipefail

config_path="${CONFIG_PATH:-run.aou.conf}"
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}"

echo "[0/5] Prepare AoU chr21/chr22 local inputs"
python3 rewrite/testdata/aou/prepare_aou.py \
  --chromosome 21 22

echo "[1/5] Prepare ancestry-specific secure inputs"
go run -mod=vendor secure-rvas.go prepare \
  --config "${config_path}" \
  --clear

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

echo "Secure RVAS AoU workflow completed"
