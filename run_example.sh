#!/bin/bash

set -euo pipefail

NUM_MAIN_PARTY=2
MODE="gwas"
OUTPUT_SUFFIX=""
SKATO_RHO="0.5"
LOG_PREFIX="stdout"

usage() {
  cat <<'EOF'
Usage: bash run_example.sh [options]

Options:
  --mode <gwas|skat|burden|skato>
  --output-suffix <suffix>
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
    --output-suffix)
      OUTPUT_SUFFIX="${2:?missing value for --output-suffix}"
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

if [[ -n "$OUTPUT_SUFFIX" ]]; then
  LOG_PREFIX="${LOG_PREFIX}_${OUTPUT_SUFFIX}"
fi

for (( i = 0; i <= NUM_MAIN_PARTY; i++ ))
do
  echo "Running PID=$i mode=$MODE"
  CMD="PID=$i SFGWAS_MODE=$MODE SFGWAS_OUTPUT_SUFFIX=$OUTPUT_SUFFIX SKATO_RHO=$SKATO_RHO go run sfgwas.go > ${LOG_PREFIX}_party${i}.txt 2>&1"
  if [[ $i -eq $NUM_MAIN_PARTY ]]; then
    eval "$CMD"
  else
    eval "$CMD" &
  fi
done
