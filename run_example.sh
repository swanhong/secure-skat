#!/bin/bash

set -euo pipefail

if [[ -n "${GOFLAGS:-}" ]]; then
  export GOFLAGS="${GOFLAGS} -mod=vendor"
else
  export GOFLAGS="-mod=vendor"
fi

unset GOROOT
export GOCACHE="$(pwd)/.local/go-build-cache"
mkdir -p "${GOCACHE}"

NUM_MAIN_PARTY=2
MODE="gwas"
SKATO_RHO="0.5"
DATASET="example_data"
CONFIG_DIR="config"
RUN_ROOT=""
RUN_NAME=""
RUN_ID=""
ORIGINAL_ARGS=("$@")

usage() {
  cat <<'EOF'
Usage: bash run_example.sh [options]

Options:
  --mode <gwas|skat|burden|skato>
  --dataset <dataset-root>
  --config-dir <config-dir>
  --skato-rho <0..1>
  --help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode)
      MODE="${2:?missing value for --mode}"
      shift 2
      ;;
    --dataset)
      DATASET="${2:?missing value for --dataset}"
      shift 2
      ;;
    --config-dir)
      CONFIG_DIR="${2:?missing value for --config-dir}"
      shift 2
      ;;
    --skato-rho)
      SKATO_RHO="${2:?missing value for --skato-rho}"
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

timestamp="$(date '+%y%m%d_%H%M%S')"
RUN_ID="$(od -An -N2 -tx1 /dev/urandom | tr -d ' \n' | cut -c1-4)"
RUN_NAME="output_${timestamp}_${RUN_ID}"
RUN_ROOT="out/${RUN_NAME}"
mkdir -p "$RUN_ROOT"
mkdir -p "$RUN_ROOT/cache"

metadata_file="${RUN_ROOT}/run_metadata.txt"
{
  echo "run_name=${RUN_NAME}"
  echo "run_id=${RUN_ID}"
  echo "started_at=$(date '+%Y-%m-%d %H:%M:%S %Z')"
  echo "cwd=$(pwd)"
  echo "command=bash run_example.sh ${ORIGINAL_ARGS[*]}"
  echo "mode=${MODE}"
  echo "dataset=${DATASET}"
  echo "config_dir=${CONFIG_DIR}"
  echo "skato_rho=${SKATO_RHO}"
  echo "run_root=${RUN_ROOT}"
  echo "pid_count=$((NUM_MAIN_PARTY + 1))"
  echo "git_branch=$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)"
  echo "git_commit=$(git rev-parse HEAD 2>/dev/null || echo unknown)"
} > "$metadata_file"

echo "Run output directory: ${RUN_ROOT}"
echo "Run ID: ${RUN_ID}"
echo "Dataset: ${DATASET}"
echo "Config directory: ${CONFIG_DIR}"
echo "Metadata written to: ${metadata_file}"

filter_terminal_output() {
  awk '
    /Running rare-variant protocol in mode:/ ||
    /Starting GWAS protocol/ ||
    /Finished QC/ ||
    /SKAT Step [1-4]\/4:/ ||
    /SKAT Progress:/ ||
    /Finished rare-variant statistic computation/ ||
    /Output collectively decrypted and saved to:/ ||
    /MatMult: block [0-9]+ \/ [0-9]+ elapsed time/ ||
    /panic:/ ||
    /fatal:/ ||
    /exit status/ ||
    /Error:/ ||
    /error:/ ||
    /listen tcp/ ||
    /Connection failed/ ||
    /unsupported rare-variant mode/ {
      print
      fflush()
    }
  '
}

run_party() {
  local pid="$1"
  local log_file="${RUN_ROOT}/stdout_party${pid}.txt"

  echo "Running PID=${pid} mode=${MODE}"
  (
    PID="$pid" \
    SFGWAS_MODE="$MODE" \
    SFGWAS_DATASET="$DATASET" \
    SFGWAS_CONFIG_PATH="$CONFIG_DIR" \
    SFGWAS_RUN_ROOT="$RUN_ROOT" \
    SKATO_RHO="$SKATO_RHO" \
    go run -mod=vendor sfgwas.go 2>&1
  ) | awk -v pid="$pid" '{ print "[PID=" pid "] " $0; fflush() }' | tee "$log_file" | filter_terminal_output
}

status=0
pids=()
for (( i = 0; i <= NUM_MAIN_PARTY; i++ ))
do
  run_party "$i" &
  pids+=("$!")
done

for pid in "${pids[@]}"; do
  if ! wait "$pid"; then
    status=1
  fi
done

{
  echo "finished_at=$(date '+%Y-%m-%d %H:%M:%S %Z')"
  echo "exit_status=${status}"
} >> "$metadata_file"

echo "Run finished with status ${status}"
echo "Run ID: ${RUN_ID}"
echo "Dataset: ${DATASET}"
echo "Config directory: ${CONFIG_DIR}"
echo "Outputs are under: ${RUN_ROOT}"

exit "$status"
