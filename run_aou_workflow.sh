#!/usr/bin/env bash

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/setup/env.sh"

config_path="${CONFIG_PATH:-run.aou.conf}"
chromosome_values="$(
  python3 -c '
import sys, tomllib
with open(sys.argv[1], "rb") as config_file:
    print(*tomllib.load(config_file)["chromosomes"])
' "${config_path}"
)"
read -r -a chromosomes <<< "${chromosome_values}"

echo "[0/5] Prepare AoU chromosome inputs: ${chromosomes[*]}"
python3 rewrite/testdata/aou/prepare_aou.py \
  --chromosome "${chromosomes[@]}"

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
