#!/usr/bin/env bash
set -euo pipefail

DIAG=${BATCH_DIAG:-0}
WORK=/mnt/data/secure-skat
REPO=$WORK/repo
INPUT=$WORK/input
RUN=$WORK/run
mkdir -p "$REPO" "$INPUT" "$RUN" "$RESULTS"

finish_diag() {
  local rc=$?
  trap - EXIT
  cp "$RUN"/prep.log "$RUN"/party*.log "$RESULTS"/ 2>/dev/null || true
  echo "exit=$rc"
  exit 0
}

if [ "$DIAG" = 1 ]; then
  exec > "$RESULTS/diag_report.txt" 2>&1
  trap finish_diag EXIT
fi

gsutil -u "$GOOGLE_CLOUD_PROJECT" cp "$CODE_BUNDLE_GCS" "$WORK/code.tar.gz"
tar -xzf "$WORK/code.tar.gz" -C "$REPO"

set -a
. "$REPO/fed_aou.conf"
set +a

gsutil -u "$GOOGLE_CLOUD_PROJECT" -m cp \
  "$GENO_GCS".pgen "$GENO_GCS".pvar "$GENO_GCS".psam \
  "$COV_GCS" "$PHENO_GCS" "$INPUT/"

# plink2 ships in the bundle (the runner image has no plink2).
chmod +x "$REPO/plink2"
export PLINK2="$REPO/plink2"
export FED_PGEN="$INPUT/$(basename "$GENO_GCS")"
export FED_ANCESTRY="$INPUT/$(basename "$COV_GCS")"
export FED_PHENO="$INPUT/$(basename "$PHENO_GCS")"
export FED_OUT=$RUN
export FED_KEYS=$REPO/example_data/keys
export SKIP_BUILD=1
echo ">>> setup done (plink2=$PLINK2) — launching run_fed.sh"

if [ "$DIAG" = 1 ]; then
  export FED_SPLIT_ANCESTRY=0 FED_ANCESTRY_GROUP=EUR
  export FED_NGENES=1 FED_NSUB=1000 FED_PROBES=1
  timeout 480 bash "$REPO/run_fed.sh"
else
  bash "$REPO/run_fed.sh"
fi

# Only aggregate outputs leave the Batch VM.
cd "$RUN"
while IFS= read -r -d '' file; do
  rel=${file#./}
  mkdir -p "$RESULTS/$(dirname "$rel")"
  cp "$file" "$RESULTS/$rel"
done < <(
  find . -type f \( \
    -name fed_results.csv -o \
    -name timing_steps.csv -o \
    -name communication_summary.csv -o \
    -name manifest.json -o \
    -name '*.png' -o \
    -name 'skat_fed*_out.txt' \
  \) -print0
)
cp "$REPO/fed_aou.conf" "$RESULTS/fed_aou.conf"
