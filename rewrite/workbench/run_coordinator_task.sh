#!/usr/bin/env bash
set -euo pipefail

: "${CODE_BUNDLE:?CODE_BUNDLE is required}"
: "${R_ENV_ARCHIVE:?R_ENV_ARCHIVE is required}"
: "${PHENOTYPE:?PHENOTYPE is required}"
: "${COVARIATE:?COVARIATE is required}"
: "${ANNOTATIONS:?ANNOTATIONS is required}"
: "${COORDINATOR_RESULTS:?COORDINATOR_RESULTS is required}"
: "${project:?project is required}"
: "${genotype_prefix:?genotype_prefix is required}"
: "${chromosomes:?chromosomes is required}"

WORK=${WORK:-/mnt/data/secure-skat-coordinator}
REPO=$WORK/shared/repo
R_ENV=$WORK/shared/r-env
CHROMOSOME_RESULTS=$COORDINATOR_RESULTS/chromosomes
FINAL_RESULTS=$COORDINATOR_RESULTS/final
mkdir -p "$REPO" "$R_ENV" "$CHROMOSOME_RESULTS" "$FINAL_RESULTS"

LOG=${LOG:-$COORDINATOR_RESULTS/coordinator.log}
if [ "$LOG" != "-" ]; then
  exec > >(tee -a "$LOG") 2>&1
fi

finish() {
  local rc=$?
  trap - EXIT
  if [ "$rc" -ne 0 ] && [ -n "${RESULT_GCS:-}" ]; then
    gcloud storage rsync --recursive "$COORDINATOR_RESULTS" "$RESULT_GCS" \
      --billing-project "$project" 2>/dev/null || true
  fi
  exit "$rc"
}
trap finish EXIT

echo ">>> coordinator: shared runtime"
tar -xzf "$CODE_BUNDLE" -C "$REPO"
if [ ! -x "$REPO/plink2" ]; then
  echo "plink2 binary not found: $REPO/plink2" >&2
  exit 1
fi
tar -xzf "$R_ENV_ARCHIVE" -C "$R_ENV"
"$R_ENV/bin/conda-unpack"
"$R_ENV/bin/Rscript" -e \
  'stopifnot(requireNamespace("SKAT", quietly = TRUE))'

IFS=, read -r -a CHROMOSOME_NUMBERS <<< "$chromosomes"
for chromosome_number in "${CHROMOSOME_NUMBERS[@]}"; do
  chromosome=chr$chromosome_number
  chromosome_input=$WORK/$chromosome/input
  chromosome_work=$WORK/$chromosome/work
  chromosome_output=$CHROMOSOME_RESULTS/$chromosome
  chromosome_genotype=$genotype_prefix.$chromosome
  annotation=$ANNOTATIONS/${chromosome}_annotation.tsv
  pgen=$chromosome_input/$(basename "$chromosome_genotype.pgen")
  pvar=$chromosome_input/$(basename "$chromosome_genotype.pvar")
  psam=$chromosome_input/$(basename "$chromosome_genotype.psam")
  mkdir -p "$chromosome_input" "$chromosome_output"

  if [ ! -s "$annotation" ]; then
    echo "annotation not found: $annotation" >&2
    exit 1
  fi

  echo ">>> $chromosome: genotype localization"
  gcloud storage cp --billing-project "$project" \
    "$chromosome_genotype.pgen" \
    "$chromosome_genotype.pvar" \
    "$chromosome_genotype.psam" \
    "$chromosome_input/"

  PGEN=$pgen \
  PVAR=$pvar \
  PSAM=$psam \
  ANNOTATION=$annotation \
  CHROMOSOME=$chromosome \
  CHROMOSOME_RESULTS=$chromosome_output \
  RESULT_GCS=${RESULT_GCS:+$RESULT_GCS/chromosomes/$chromosome} \
  WORK=$chromosome_work \
  REPO=$REPO \
  R_ENV=$R_ENV \
    bash "$REPO/rewrite/workbench/run_chromosome_task.sh"

  rm -rf "$chromosome_input" "$chromosome_work"
done

echo ">>> coordinator: aggregate results"
mkdir -p "$WORK/matplotlib"
PYTHONPATH="$REPO" MPLCONFIGDIR="$WORK/matplotlib" \
  python3 -m rewrite.workbench aggregate \
  --input-root "$CHROMOSOME_RESULTS" \
  --output-dir "$FINAL_RESULTS" \
  --chromosomes "${CHROMOSOME_NUMBERS[@]}"

touch "$FINAL_RESULTS/_SUCCESS" "$COORDINATOR_RESULTS/_SUCCESS"
echo ">>> coordinator: complete"
