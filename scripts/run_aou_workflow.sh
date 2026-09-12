#!/usr/bin/env bash

set -euo pipefail

source scripts/setup/env.sh

config_path="${CONFIG_PATH:-config/aou}"
reference_engine="${REFERENCE_ENGINE:-python}"
step_labels=()
step_seconds=()

run_step() {
  local label="$1"
  shift
  echo "${label}"
  local started_at="${SECONDS}"
  "$@"
  step_labels+=("${label}")
  step_seconds+=("$((SECONDS - started_at))")
}

run_step "[0/7] Prepare AoU chromosome inputs" \
  python3 python/datasets/aou/prepare_aou.py \
  --config "${config_path}"

run_step "[1/7] Prepare ancestry-specific secure inputs" \
  go run -mod=vendor secure-rvas.go prepare \
  --config "${config_path}"

run_step "[2/7] Generate shared PRG keys" \
  go run -mod=vendor secure-rvas.go keygen \
  --config "${config_path}"

run_step "[3/7] Run secure Burden/SKAT" \
  go run -mod=vendor secure-rvas.go run \
  --config "${config_path}"

run_step "[4/7] Run ancestry-specific ${reference_engine} reference" \
  python3 python/analysis/run_reference.py \
  --config "${config_path}" \
  --engine "${reference_engine}"

run_step "[5/7] Compare secure and reference results" \
  python3 python/analysis/compare_secure_to_reference.py \
  --config "${config_path}"

run_step "[6/7] Generate scatter and Manhattan plots" \
  python3 python/analysis/plot_secure_vs_reference.py \
  --config "${config_path}"

run_step "[7/7] Summarize metrics" \
  ./scripts/summarize_metrics.sh "${config_path}"

echo "Secure RVAS AoU workflow completed"
echo "Workflow timing summary"
for index in "${!step_labels[@]}"; do
  printf "  %-65s %d seconds\n" "${step_labels[index]}" "${step_seconds[index]}"
done
