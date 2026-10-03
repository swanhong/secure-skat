#!/usr/bin/env bash

set -euo pipefail

cd "$(dirname "$0")/.."
source scripts/setup/env.sh

config_a="${CONFIG_A:-config/1kg-A}"
config_b="${CONFIG_B:-config/1kg-B}"
config_0="${CONFIG_0:-config/1kg-cp0}"
reference_engine="${REFERENCE_ENGINE:-python}"
log_dir="${LOG_DIR:-output/1kg-logs}"

echo "[0/8] Generate 1000 Genomes test data"
python3 python/datasets/onekg/prepare_1kgenome.py --config "$config_a" --num-pheno 2

echo "[1/8] Prepare A and its public variant lists"
go run -mod=vendor secure-rvas.go prepare --config "$config_a" --party 1

echo "[2/8] Prepare B using A's public variant lists"
go run -mod=vendor secure-rvas.go prepare --config "$config_b" --party 2

echo "[3/8] Generate shared PRG keys"
go run -mod=vendor secure-rvas.go keygen \
  --config-party0 "$config_0" --config-party1 "$config_a" --config-party2 "$config_b"

echo "[4/8] Run secure Burden/SKAT; logs: $log_dir"
mkdir -p "$log_dir"
go run -mod=vendor secure-rvas.go party --config "$config_0" --party 0 > "$log_dir/party0.log" 2>&1 &
pid_0=$!
go run -mod=vendor secure-rvas.go party --config "$config_a" --party 1 > "$log_dir/party1.log" 2>&1 &
pid_a=$!
go run -mod=vendor secure-rvas.go party --config "$config_b" --party 2 > "$log_dir/party2.log" 2>&1 &
pid_b=$!
wait "$pid_0"
wait "$pid_a"
wait "$pid_b"

echo "[5/8] Run ancestry-specific $reference_engine reference"
python3 python/analysis/run_reference.py --config "$config_a" --config-b "$config_b" --engine "$reference_engine"

echo "[6/8] Compare secure and reference results"
python3 python/analysis/compare_secure_to_reference.py --config "$config_a"

echo "[7/8] Generate scatter and Manhattan plots"
python3 python/analysis/plot_secure_vs_reference.py --config "$config_a"

echo "[8/8] Summarize metrics"
"${PYTHON_BIN:-python}" python/analysis/summarize_metrics.py --config "$config_a" --config-b "$config_b"

echo "Secure RVAS 1KG workflow completed"
