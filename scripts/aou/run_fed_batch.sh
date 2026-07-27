#!/usr/bin/env bash
set -euo pipefail

echo "=== [0/3] prepare Batch workspace ==="
: "${CODE_BUNDLE_GCS:?missing CODE_BUNDLE_GCS}"
: "${GENO_GCS:?missing GENO_GCS}"
: "${COV_GCS:?missing COV_GCS}"
: "${PHENO_GCS:?missing PHENO_GCS}"
: "${RESULTS:?missing RESULTS}"
: "${GOOGLE_CLOUD_PROJECT:?missing GOOGLE_CLOUD_PROJECT}"

WORK_ROOT=${BATCH_WORK_ROOT:-/mnt/data/secure-skat-work}
REPO=$WORK_ROOT/repo
RUN=$WORK_ROOT/run
DL=$WORK_ROOT/inputs
mkdir -p "$REPO" "$RUN" "$RESULTS" "$DL"

echo "=== [1/3] download inputs (requester-pays bucket) ==="
# dsub --input localization was failing before this script ran (job died in a
# pre-script runnable, no logs). The controlled bucket is requester-pays; pull
# inputs here with gsutil -u so the download runs in this image, under our control.
gsutil -u "$GOOGLE_CLOUD_PROJECT" -m cp \
  "$CODE_BUNDLE_GCS" \
  "$GENO_GCS".pgen "$GENO_GCS".pvar "$GENO_GCS".psam \
  "$COV_GCS" "$PHENO_GCS" \
  "$DL/"

CODE_BUNDLE=$DL/$(basename "$CODE_BUNDLE_GCS")
tar -xzf "$CODE_BUNDLE" -C "$REPO"

set -a
. "$REPO/fed_aou.conf"
set +a
PLINK2=$(command -v plink2)
python3 -c 'import numpy'
export PLINK2

echo "=== [2/3] run secure-SKAT ==="
export FED_PGEN=$DL/$(basename "$GENO_GCS")
export FED_ANCESTRY=$DL/$(basename "$COV_GCS")
export FED_PHENO=$DL/$(basename "$PHENO_GCS")
export FED_OUT=$RUN
export FED_KEYS=$REPO/example_data/keys
export SKIP_BUILD=1
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
