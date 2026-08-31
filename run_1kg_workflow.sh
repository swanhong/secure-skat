#!/usr/bin/env bash

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/setup/env.sh"

config_path="${CONFIG_PATH:-run.1kg.conf}"
reference_engine="${REFERENCE_ENGINE:-python}"
chromosome_values="$(
  python3 -c '
import sys, tomllib
with open(sys.argv[1], "rb") as config_file:
    print(*tomllib.load(config_file)["chromosomes"])
' "${config_path}"
)"
read -r -a chromosomes <<< "${chromosome_values}"
prepare_args=(--config "${config_path}")
if [[ "${CLEAR_RUN_DIR:-0}" == "1" ]]; then
  prepare_args+=(--clear)
fi

echo "[0/5] Generate 1000 Genomes test data"
python3 rewrite/testdata/1kgenome/prepare_1kgenome.py \
  --chromosome "${chromosomes[@]}" \
  --num-pheno 2

echo "[1/5] Prepare ancestry-specific secure inputs"
go run -mod=vendor secure-rvas.go prepare "${prepare_args[@]}"

echo "[2/5] Run secure Burden/SKAT"
go run -mod=vendor secure-rvas.go run \
  --config "${config_path}"

echo "[3/5] Run ancestry-specific ${reference_engine} reference"
python3 rewrite/analysis/run_reference.py \
  --config "${config_path}" \
  --engine "${reference_engine}"

echo "[4/5] Compare secure and reference results"
python3 rewrite/analysis/compare_secure_to_reference.py \
  --config "${config_path}"

echo "[5/5] Generate scatter and Manhattan plots"
python3 rewrite/analysis/plot_secure_vs_reference.py \
  --config "${config_path}"

echo "Secure RVAS 1KG workflow completed"
