#!/usr/bin/env bash
set -euo pipefail

: "${CHROMOSOME:?CHROMOSOME is required}"
: "${PGEN:?PGEN is required}"
: "${PVAR:?PVAR is required}"
: "${PSAM:?PSAM is required}"
: "${ANNOTATION:?ANNOTATION is required}"
: "${PHENOTYPE:?PHENOTYPE is required}"
: "${COVARIATE:?COVARIATE is required}"
: "${CODE_BUNDLE:?CODE_BUNDLE is required}"
: "${R_ENV_ARCHIVE:?R_ENV_ARCHIVE is required}"
: "${CHROMOSOME_RESULTS:?CHROMOSOME_RESULTS is required}"
: "${project:?project is required}"

WORK=${WORK:-/mnt/data/secure-skat}
REPO=${REPO:-$WORK/repo}
INPUT=${INPUT:-$WORK/input}
PREPROCESSED=${PREPROCESSED:-$WORK/preprocessed}
RUN=${RUN:-$WORK/run}
R_ENV=${R_ENV:-$WORK/r-env}
PORT_BASE=${PORT_BASE:-18000}
LANE_INDEX=${LANE_INDEX:-0}
mkdir -p "$REPO" "$INPUT" "$RUN" "$CHROMOSOME_RESULTS"

LOG=$CHROMOSOME_RESULTS/run.log
exec > >(tee -a "$LOG") 2>&1

WORKFLOW_TIMING=$CHROMOSOME_RESULTS/workflow_timing.csv
clock_ns() {
  python3 -c 'import time; print(time.perf_counter_ns())'
}
record_workflow_timing() {
  local phase=$1
  local started_ns=$2
  local status=${3:-success}
  PYTHONPATH="$REPO" python3 -m rewrite.workbench timing-event \
    --output "$WORKFLOW_TIMING" \
    --started-ns "$started_ns" \
    --component workflow \
    --scope chromosome \
    --chromosome "$CHROMOSOME" \
    --lane "$LANE_INDEX" \
    --phase "$phase" \
    --parent-phase chromosome_total \
    --status "$status"
}
chromosome_started_ns=${CHROMOSOME_STARTED_NS:-$(clock_ns)}

PARTY_PIDS=()
finish() {
  local rc=$?
  trap - EXIT
  for process in "${PARTY_PIDS[@]:-}"; do
    kill "$process" 2>/dev/null || true
  done
  for process in "${PARTY_PIDS[@]:-}"; do
    wait "$process" 2>/dev/null || true
  done
  cp "$RUN"/party*.log "$CHROMOSOME_RESULTS/" 2>/dev/null || true
  if [ "$rc" -ne 0 ] && [ -n "${RESULT_GCS:-}" ]; then
    gcloud storage rsync --recursive "$CHROMOSOME_RESULTS" "$RESULT_GCS" \
      --billing-project "$project" 2>/dev/null || true
  fi
  exit "$rc"
}
trap finish EXIT

if [ ! -x "$REPO/secure-skat" ]; then
  tar -xzf "$CODE_BUNDLE" -C "$REPO"
fi
if [ ! -x "$REPO/plink2" ]; then
  echo "plink2 binary not found: $REPO/plink2" >&2
  exit 1
fi

echo ">>> $CHROMOSOME: R environment"
r_environment_started_ns=$(clock_ns)
if [ ! -x "$R_ENV/bin/Rscript" ]; then
  mkdir -p "$R_ENV"
  tar -xzf "$R_ENV_ARCHIVE" -C "$R_ENV"
  "$R_ENV/bin/conda-unpack"
fi
"$R_ENV/bin/Rscript" -e \
  'stopifnot(requireNamespace("SKAT", quietly = TRUE))'
record_workflow_timing r_environment "$r_environment_started_ns"

ln -s "$PGEN" "$INPUT/genotype.pgen"
ln -s "$PVAR" "$INPUT/genotype.pvar"
ln -s "$PSAM" "$INPUT/genotype.psam"
ln -s "$ANNOTATION" "$INPUT/annotation.csv"
ln -s "$PHENOTYPE" "$INPUT/phenotype.csv"
ln -s "$COVARIATE" "$INPUT/covariate.csv"

IFS=, read -r -a PHENOTYPE_COLUMN_LIST <<< "${phenotype_columns:?}"
IFS=, read -r -a COVARIATE_COLUMN_LIST <<< "${covariate_columns:?}"
IFS=, read -r -a MASKS <<< "${mask:?}"
genes_path=$genes
if [ "$genes_path" != "all" ] && [[ "$genes_path" != /* ]]; then
  genes_path=$REPO/$genes_path
fi

PREPARE=(
  python3 -m rewrite.preprocessing prepare
  --pgen-prefix "$INPUT/genotype"
  --annotation "$ANNOTATION"
  --phenotype "$PHENOTYPE"
  --covariate "$COVARIATE"
  --phenotype-id-column "${phenotype_id_column:?}"
  --covariate-id-column "${covariate_id_column:?}"
  --phenotype-columns "${PHENOTYPE_COLUMN_LIST[@]}"
  --covariate-columns "${COVARIATE_COLUMN_LIST[@]}"
  --chromosome "$CHROMOSOME"
  --genes "$genes_path"
  --samples-per-cohort "${samples_per_cohort:?}"
  --sample-seed "${sample_seed:?}"
  --role-seed "${role_seed:?}"
  --out "$PREPROCESSED"
  --timing-output "$CHROMOSOME_RESULTS/preprocessing_timing.csv"
)
if [ -n "${covariate_array_column:-}" ]; then
  PREPARE+=(--covariate-array-column "$covariate_array_column")
fi
for mask in "${MASKS[@]}"; do
  PREPARE+=(--mask "$mask")
done

echo ">>> $CHROMOSOME: preprocessing"
preprocessing_started_ns=$(clock_ns)
if ! PATH="$REPO:$PATH" PYTHONPATH="$REPO" "${PREPARE[@]}"; then
  record_workflow_timing preprocessing "$preprocessing_started_ns" failure
  cp "$PREPROCESSED/validation_errors.csv" "$CHROMOSOME_RESULTS/" 2>/dev/null || true
  exit 1
fi
record_workflow_timing preprocessing "$preprocessing_started_ns"

SHARED_KEYS=$RUN/shared-keys
shared_keys_started_ns=$(clock_ns)
if ! PYTHONPATH="$REPO" python3 -m rewrite.workbench keys --output-dir "$SHARED_KEYS"; then
  record_workflow_timing shared_key_generation "$shared_keys_started_ns" failure
  exit 1
fi
record_workflow_timing shared_key_generation "$shared_keys_started_ns"

echo ">>> $CHROMOSOME: secure protocol"
secure_started_ns=$(clock_ns)
for party in 0 1 2; do
  "$REPO/secure-skat" run \
    --party "$party" \
    --lane "$LANE_INDEX" \
    --input "$PREPROCESSED" \
    --output "$RUN/secure_results.csv" \
    --timing-output "$CHROMOSOME_RESULTS/secure_timing_party${party}.csv" \
    --port-base "$PORT_BASE" \
    --shared-keys "$SHARED_KEYS" \
    --ckks "${ckks:?}" \
    --data-bits "${data_bits:?}" \
    --frac-bits "${frac_bits:?}" \
    --probes "${probes:?}" \
    --seed "${seed:?}" \
    > "$RUN/party${party}.log" 2>&1 &
  PARTY_PIDS+=("$!")
done

remaining_parties=${#PARTY_PIDS[@]}
while [ "$remaining_parties" -gt 0 ]; do
  if wait -n "${PARTY_PIDS[@]}"; then
    remaining_parties=$((remaining_parties - 1))
  else
    secure_rc=$?
    record_workflow_timing secure_protocol "$secure_started_ns" failure
    exit "$secure_rc"
  fi
done
PARTY_PIDS=()
record_workflow_timing secure_protocol "$secure_started_ns"

echo ">>> $CHROMOSOME: R::SKAT Burden and Davies"
r_started_ns=$(clock_ns)
if ! "$R_ENV/bin/Rscript" "$REPO/rewrite/workbench/run_reference_skat.R" \
  "$PREPROCESSED" "$RUN/r_results.csv" \
  "$CHROMOSOME_RESULTS/r_timing.csv"; then
  record_workflow_timing r_reference "$r_started_ns" failure
  exit 1
fi
record_workflow_timing r_reference "$r_started_ns"

echo ">>> $CHROMOSOME: join secure and R results"
join_started_ns=$(clock_ns)
if ! PYTHONPATH="$REPO" python3 -m rewrite.workbench join \
  --secure "$RUN/secure_results.csv" \
  --r-results "$RUN/r_results.csv" \
  --output "$CHROMOSOME_RESULTS/gene_results.csv" \
  --summary "$CHROMOSOME_RESULTS/error_summary.csv"; then
  record_workflow_timing join_results "$join_started_ns" failure
  exit 1
fi
record_workflow_timing join_results "$join_started_ns"

cp "$PREPROCESSED/metadata.json" "$CHROMOSOME_RESULTS/metadata.json"
cp "$RUN/secure_results.csv" "$CHROMOSOME_RESULTS/secure_results.csv"
cp "$RUN/r_results.csv" "$CHROMOSOME_RESULTS/r_results.csv"
record_workflow_timing chromosome_total "$chromosome_started_ns"
PYTHONPATH="$REPO" python3 -m rewrite.workbench timing-combine \
  --output "$CHROMOSOME_RESULTS/timing_events.csv" \
  --inputs \
    "$WORKFLOW_TIMING" \
    "$CHROMOSOME_RESULTS/preprocessing_timing.csv" \
    "$CHROMOSOME_RESULTS/secure_timing_party0.csv" \
    "$CHROMOSOME_RESULTS/secure_timing_party1.csv" \
    "$CHROMOSOME_RESULTS/secure_timing_party2.csv" \
    "$CHROMOSOME_RESULTS/r_timing.csv"
touch "$CHROMOSOME_RESULTS/_SUCCESS"
echo ">>> $CHROMOSOME: complete"
