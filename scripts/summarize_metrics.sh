#!/usr/bin/env bash

set -euo pipefail

config_path="${1:-config/aou}"
party_id="${2:-1}"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
python_bin="${PYTHON_BIN:-python}"

exec "${python_bin}" \
  "${script_dir}/rewrite/analysis/summarize_metrics.py" \
  --config "${config_path}" \
  --party-id "${party_id}"
