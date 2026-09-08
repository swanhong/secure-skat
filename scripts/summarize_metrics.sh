#!/usr/bin/env bash

set -euo pipefail

config_path="${1:-config/aou}"
party_id="${2:-1}"
python_bin="${PYTHON_BIN:-python}"

exec "${python_bin}" \
  python/analysis/summarize_metrics.py \
  --config "${config_path}" \
  --party-id "${party_id}"
