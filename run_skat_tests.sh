#!/bin/bash

set -euo pipefail

if [[ -n "${GOFLAGS:-}" ]]; then
  export GOFLAGS="${GOFLAGS} -mod=vendor"
else
  export GOFLAGS="-mod=vendor"
fi

unset GOROOT
export GOCACHE="$(pwd)/.local/go-build-cache"
export GOMODCACHE="$(pwd)/.local/go-mod-cache"
mkdir -p "${GOCACHE}"
mkdir -p "${GOMODCACHE}"

TEST_NAME="TestSecureSKATEndToEnd"
RUN_SUFFIX=""
NUM_MAIN_PARTY=2
LOG_PREFIX="test_stdout"

usage() {
  cat <<'EOF'
Usage: bash run_skat_tests.sh [options]

Options:
  --test-name <go-test-regex>
  --output-suffix <suffix>
  --help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --test-name)
      TEST_NAME="${2:?missing value for --test-name}"
      shift 2
      ;;
    --output-suffix)
      RUN_SUFFIX="${2:?missing value for --output-suffix}"
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

if [[ -n "$RUN_SUFFIX" ]]; then
  LOG_PREFIX="${LOG_PREFIX}_${RUN_SUFFIX}_party"
else
  LOG_PREFIX="${LOG_PREFIX}_party"
fi

pkill -f "skat_test.test" 2>/dev/null || true

echo "Compiling test binary..."
go test -mod=vendor -c ./gwas/ -o skat_test.test || { echo "Compilation failed"; exit 1; }

rm -f ${LOG_PREFIX}*.txt

for (( i = 0; i <= NUM_MAIN_PARTY; i++ ))
do
  echo "Running PID=$i"
  CMD="PID=$i TEST_RUN_SUFFIX=$RUN_SUFFIX ./skat_test.test -test.v -test.run $TEST_NAME > ${LOG_PREFIX}${i}.txt 2>&1"
  if [[ $i -eq $NUM_MAIN_PARTY ]]; then
    eval "$CMD"
  else
    eval "$CMD" &
    sleep 2
  fi
done

wait
if [[ -n "$RUN_SUFFIX" ]]; then
  echo "Tests completed. Logs: ${LOG_PREFIX}*.txt"
  echo "Outputs: out/party*_${RUN_SUFFIX}, cache/party*_${RUN_SUFFIX}"
else
  echo "Tests completed. Check ${LOG_PREFIX}*.txt for outputs."
fi
