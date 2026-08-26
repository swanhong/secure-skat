#!/usr/bin/env bash
set -euo pipefail

: "${CODE_BUNDLE:?CODE_BUNDLE is required}"
: "${CHROMOSOME_RESULTS:?CHROMOSOME_RESULTS is required}"
: "${FINAL_RESULTS:?FINAL_RESULTS is required}"
: "${chromosomes:?chromosomes is required}"

WORK=/mnt/data/secure-skat-aggregate
REPO=$WORK/repo
mkdir -p "$REPO" "$FINAL_RESULTS" "$WORK/matplotlib"
tar -xzf "$CODE_BUNDLE" -C "$REPO"
IFS=, read -r -a CHROMOSOME_NUMBER_LIST <<< "$chromosomes"
PYTHONPATH="$REPO" MPLCONFIGDIR="$WORK/matplotlib" \
  python3 -m rewrite.workbench aggregate \
  --input-root "$CHROMOSOME_RESULTS" \
  --output-dir "$FINAL_RESULTS" \
  --chromosomes "${CHROMOSOME_NUMBER_LIST[@]}"

touch "$FINAL_RESULTS/_SUCCESS"
