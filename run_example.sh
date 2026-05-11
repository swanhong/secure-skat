#!/bin/bash

# Fail early on script errors, undefined variables, and failed commands inside pipes.
set -euo pipefail

# Force vendored Go dependencies so local runs do not depend on the caller's module settings.
if [[ -n "${GOFLAGS:-}" ]]; then
  export GOFLAGS="${GOFLAGS} -mod=vendor"
else
  export GOFLAGS="-mod=vendor"
fi

# Normalize the Go toolchain environment and keep caches inside the repo.
# This makes repeated local runs faster and avoids surprises from stale global Go state.
unset GOROOT
export GOCACHE="$(pwd)/.local/go-build-cache"
export GOMODCACHE="$(pwd)/.local/go-mod-cache"
mkdir -p "${GOCACHE}"
mkdir -p "${GOMODCACHE}"

# Default knobs for the example run. These are the only user-facing settings this wrapper exposes.
# The actual protocol code reads them through environment variables set in run_party().
NUM_MAIN_PARTY=2
MODE="skat"
SKATO_RHO="0.5"
DATASET="example_data"
CONFIG_DIR="config"
RUN_ROOT=""
RUN_BASE="${SFGWAS_RUN_BASE:-out}"
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
  --run-base <output-parent-dir>
  --skato-rho <0..1>
  --help
EOF
}

# Parse a small CLI surface area and leave the rest of the configuration in TOML/env vars.
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
    --run-base)
      RUN_BASE="${2:?missing value for --run-base}"
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

# Create a unique output root for this invocation so logs/results do not overwrite prior runs.
# As a result, the name of root directory will be set as output_YYMMDD_HHMMSS_RAND
#   - timestamp is in the local timezone
#   - RAND is a 4-digit random string, in order to enable easy recall
#   - e.g.) output_251231_235959_a1b2 (later, we can call this result via the RUN_ID "a1b2")
run_started_epoch="$(date '+%s')"
run_started_at="$(date '+%Y-%m-%d %H:%M:%S %Z')"
run_started_at_iso="$(date '+%Y-%m-%dT%H:%M:%S%z')"
timestamp="$(date '+%y%m%d_%H%M%S')"
RUN_ID="$(od -An -N2 -tx1 /dev/urandom | tr -d ' \n' | cut -c1-4)"
RUN_NAME="output_${timestamp}_${RUN_ID}"
RUN_ROOT="${RUN_BASE%/}/${RUN_NAME}"
mkdir -p "$RUN_ROOT"
mkdir -p "$RUN_ROOT/cache"

# Record run metadata for later debugging/comparison scripts.
# This is not required for the protocol itself, but it helps track what was executed.
metadata_file="${RUN_ROOT}/run_metadata.txt"
{
  echo "run_name=${RUN_NAME}"
  echo "run_id=${RUN_ID}"
  echo "started_at=${run_started_at}"
  echo "started_at_iso=${run_started_at_iso}"
  echo "started_epoch=${run_started_epoch}"
  echo "cwd=$(pwd)"
  echo "command=bash run_example.sh ${ORIGINAL_ARGS[*]}"
  echo "mode=${MODE}"
  echo "dataset=${DATASET}"
  echo "config_dir=${CONFIG_DIR}"
  echo "skato_rho=${SKATO_RHO}"
  echo "run_base=${RUN_BASE}"
  echo "run_root=${RUN_ROOT}"
  echo "pid_count=$((NUM_MAIN_PARTY + 1))"
  echo "git_branch=$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)"
  echo "git_commit=$(git rev-parse HEAD 2>/dev/null || echo unknown)"
} > "$metadata_file"

# Print key run settings for easy reference in the terminal and logs.
echo "Run output directory: ${RUN_ROOT}"
echo "Run ID: ${RUN_ID}"
echo "Dataset: ${DATASET}"
echo "Config directory: ${CONFIG_DIR}"
echo "Metadata written to: ${metadata_file}"

# Keep the terminal readable by showing only high-signal progress/error lines.
# Full unfiltered logs are still written to stdout_party*.txt by tee in run_party().
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

# ==========================================
# ======= Main execution starts here =======
# ==========================================

# Launch one local process per party.
# sfgwas.go uses PID plus the SFGWAS_* env vars to pick config files, dataset paths, mode, and output dirs.
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
  # Prefix each log line with the party ID, save the full stream to a file,
  # and show only the filtered subset on the terminal.
  ) | awk -v pid="$pid" '{ print "[PID=" pid "] " $0; fflush() }' | tee "$log_file" | filter_terminal_output
}

# Spawn the auxiliary party (PID=0) and the data parties (PID=1..NUM_MAIN_PARTY) in parallel.
status=0
pids=()
for (( i = 0; i <= NUM_MAIN_PARTY; i++ ))
do
  run_party "$i" &
  pids+=("$!")
done

# Wait for every background process and return nonzero if any one of them failed.
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then
    status=1
  fi
done

# Append final status so other scripts can inspect the run after completion.
run_finished_epoch="$(date '+%s')"
run_finished_at="$(date '+%Y-%m-%d %H:%M:%S %Z')"
run_finished_at_iso="$(date '+%Y-%m-%dT%H:%M:%S%z')"
{
  echo "finished_at=${run_finished_at}"
  echo "finished_at_iso=${run_finished_at_iso}"
  echo "finished_epoch=${run_finished_epoch}"
  echo "exit_status=${status}"
} >> "$metadata_file"

echo "Run finished with status ${status}"
echo "Run ID: ${RUN_ID}"
echo "Dataset: ${DATASET}"
echo "Config directory: ${CONFIG_DIR}"
echo "Outputs are under: ${RUN_ROOT}"

timing_summary_file="${RUN_ROOT}/timing_summary.txt"
if python3 scripts/analysis/summarize_run_timing.py "$RUN_ROOT" | tee "$timing_summary_file"; then
  echo "Timing summary saved to: ${timing_summary_file}"
else
  echo "Timing summary unavailable; inspect ${RUN_ROOT}/stdout_party*.txt directly." >&2
fi

analysis_dir="${RUN_ROOT}/analysis"
case "$MODE" in
  skat|burden|skato)
    if [[ "$status" -eq 0 ]]; then
      if python3 scripts/analysis/skat_compare.py compare \
        --run-root "$RUN_ROOT" \
        --dataset "$DATASET" \
        --skip-reference; then
        echo "Analysis outputs are under: ${analysis_dir}"
      else
        echo "Plain-vs-secure analysis unavailable; secure run outputs are still saved under ${RUN_ROOT}." >&2
      fi
    else
      echo "Plain-vs-secure analysis skipped because the secure run failed." >&2
    fi
    ;;
  *)
    echo "Plain-vs-secure SKAT analysis skipped for mode=${MODE}."
    ;;
esac

if [[ -d "$analysis_dir" ]]; then
  echo "Analysis directory: ${analysis_dir}"
fi

exit "$status"
