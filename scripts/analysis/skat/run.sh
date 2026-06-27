#!/bin/bash
# Secure-vs-plain SKAT comparison driver.
#
#   run.sh [--secure-only|--plain-only] [DATASET]
#     DATASET: example (default) | toy | medium
#
#   default: secure (3-party) + plain (R::SKAT) + per-block compare (+ PNG).
#
# Notes: toy must be generated once (python3 scripts/make_toy_secure.py).
# medium secure is impractical at v2 (huge + stale config paths); use --plain-only.
set -e

MODE=both
case "$1" in
  --secure-only) MODE=secure; shift ;;
  --plain-only)  MODE=plain;  shift ;;
esac
DATASET="${1:-example}"

case "$DATASET" in
  example|example_data) CONFIG="config";                         PDATA="example_data";                                  FMT="pgen" ;;
  toy)                  CONFIG="test_data/toy/secure/config";    PDATA="test_data/toy/secure";                          FMT="blocks" ;;
  medium)               CONFIG="test_data/medium/config_medium"; PDATA="test_data/medium/1000g_chr22_anchor50kb_merged"; FMT="pgen" ;;
  *) echo "unknown dataset: $DATASET (example|toy|medium)" >&2; exit 1 ;;
esac

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$ROOT"
NUM_MAIN_PARTY=2
SECURE_DIR="out/party1"
PLAIN_DIR="out/plain"

run_secure() {
  echo "== secure ($DATASET, config=$CONFIG, SFGWAS_DEBUG=1 for per-block dump) =="
  go test -c ./gwas/ -o skat_test.test
  pkill -f skat_test.test 2>/dev/null || true
  rm -f test_stdout_party*.txt out/party*/qBlock_block*.txt
  for ((i = 0; i <= NUM_MAIN_PARTY; i++)); do
    PID=$i SFGWAS_DEBUG=1 SFGWAS_CONFIG_PATH="$CONFIG" \
      ./skat_test.test -test.v -test.run TestSecureSKATEndToEnd > test_stdout_party${i}.txt 2>&1 &
    if [ $i -lt $NUM_MAIN_PARTY ]; then sleep 2; else wait; fi
  done
  grep -m1 -E "PASS|FAIL" test_stdout_party1.txt || true
}

run_plain() {
  # pass secure_dir only when comparing (both mode), so --plain-only skips the diff
  local sdir=""; [ "$MODE" = both ] && sdir="$SECURE_DIR"
  echo "== plain R::SKAT ($DATASET, $FMT) =="
  Rscript scripts/analysis/skat/skat_compare.R "$ROOT" "$PDATA" "$FMT" "$PLAIN_DIR" "$sdir"
}

if [ "$MODE" != plain ];  then run_secure; fi
if [ "$MODE" != secure ]; then run_plain;  fi
