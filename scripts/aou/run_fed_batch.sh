#!/usr/bin/env bash
set -euo pipefail

echo "=== [0/3] prepare Batch workspace ==="
: "${CODE_BUNDLE:?missing CODE_BUNDLE}"
: "${PGEN_PGEN:?missing PGEN_PGEN}"
: "${PGEN_PVAR:?missing PGEN_PVAR}"
: "${PGEN_PSAM:?missing PGEN_PSAM}"
: "${ANCESTRY_TSV:?missing ANCESTRY_TSV}"
: "${PHENO_CSV:?missing PHENO_CSV}"
: "${RESULTS:?missing RESULTS}"

WORK_ROOT=${BATCH_WORK_ROOT:-/mnt/data/secure-skat-work}
REPO=$WORK_ROOT/repo
RUN=$WORK_ROOT/run
mkdir -p "$REPO" "$RUN" "$RESULTS"
tar -xzf "$CODE_BUNDLE" -C "$REPO"

set -a
. "$REPO/fed_aou.conf"
set +a
PLINK2=$(command -v plink2)
python3 -c 'import numpy'
export PLINK2

echo "=== [1/3] connect localized inputs ==="
FED_PGEN=${PGEN_PGEN%.pgen}
[ "$PGEN_PVAR" = "$FED_PGEN.pvar" ] && [ "$PGEN_PSAM" = "$FED_PGEN.psam" ] \
  || { echo "error: localized PGEN files do not share one prefix" >&2; exit 1; }

FED_ANCESTRY=$ANCESTRY_TSV
FED_PHENO=$PHENO_CSV
FED_OUT=$RUN
FED_KEYS=$REPO/example_data/keys
SKIP_BUILD=1
export FED_PGEN FED_ANCESTRY FED_PHENO FED_OUT FED_KEYS SKIP_BUILD

echo "=== [2/3] run secure-SKAT ==="
(cd "$REPO" && bash run_fed.sh)

echo "=== [3/3] collect aggregate results ==="
while IFS= read -r -d '' src; do
  rel=${src#"$RUN"/}
  mkdir -p "$RESULTS/$(dirname "$rel")"
  cp "$src" "$RESULTS/$rel"
done < <(
  find "$RUN" -type f \( \
    -name fed_results.csv -o \
    -name timing_steps.csv -o \
    -name communication_summary.csv -o \
    -name manifest.json -o \
    -name 'manhattan_*.png' -o \
    -name 'scatter_*.png' -o \
    -name skat_fed_out.txt -o \
    -name skat_fed_burden_p_out.txt -o \
    -name skat_fed_skat_p_out.txt \
  \) -print0
)
cp "$REPO/fed_aou.conf" "$RESULTS/fed_aou.conf"
