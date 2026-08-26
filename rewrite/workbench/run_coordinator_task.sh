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
: "${max_parallel_chromosomes:?max_parallel_chromosomes is required}"

if [[ ! "$max_parallel_chromosomes" =~ ^[1-9][0-9]*$ ]]; then
  echo "max_parallel_chromosomes must be a positive integer" >&2
  exit 1
fi

WORK=${WORK:-/mnt/data/secure-skat-coordinator}
REPO=$WORK/shared/repo
R_ENV=$WORK/shared/r-env
CHROMOSOME_RESULTS=$COORDINATOR_RESULTS/chromosomes
FINAL_RESULTS=$COORDINATOR_RESULTS/final
mkdir -p "$REPO" "$R_ENV" "$CHROMOSOME_RESULTS" "$FINAL_RESULTS"

COORDINATOR_TIMING=$COORDINATOR_RESULTS/coordinator_timing.csv
clock_ns() {
  python3 -c 'import time; print(time.perf_counter_ns())'
}
record_coordinator_timing() {
  local phase=$1
  local started_ns=$2
  local status=${3:-success}
  PYTHONPATH="$REPO" python3 -m rewrite.workbench timing-event \
    --output "$COORDINATOR_TIMING" \
    --started-ns "$started_ns" \
    --component coordinator \
    --scope run \
    --phase "$phase" \
    --status "$status"
}
record_chromosome_timing() {
  local output=$1
  local chromosome=$2
  local lane_index=$3
  local phase=$4
  local started_ns=$5
  local status=${6:-success}
  PYTHONPATH="$REPO" python3 -m rewrite.workbench timing-event \
    --output "$output" \
    --started-ns "$started_ns" \
    --component workflow \
    --scope chromosome \
    --chromosome "$chromosome" \
    --lane "$lane_index" \
    --phase "$phase" \
    --parent-phase chromosome_total \
    --status "$status"
}
coordinator_started_ns=$(clock_ns)

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
shared_runtime_started_ns=$(clock_ns)
tar -xzf "$CODE_BUNDLE" -C "$REPO"
if [ ! -x "$REPO/plink2" ]; then
  echo "plink2 binary not found: $REPO/plink2" >&2
  exit 1
fi
tar -xzf "$R_ENV_ARCHIVE" -C "$R_ENV"
"$R_ENV/bin/conda-unpack"
"$R_ENV/bin/Rscript" -e \
  'stopifnot(requireNamespace("SKAT", quietly = TRUE))'
record_coordinator_timing shared_runtime "$shared_runtime_started_ns"

IFS=, read -r -a CHROMOSOME_NUMBERS <<< "$chromosomes"
run_chromosome() {
  local chromosome_number=$1
  local port_base=$2
  local lane_index=$3
  local chromosome=chr$chromosome_number
  local chromosome_input=$WORK/$chromosome/input
  local chromosome_work=$WORK/$chromosome/work
  local chromosome_output=$CHROMOSOME_RESULTS/$chromosome
  local gcloud_config=$chromosome_work/gcloud
  local chromosome_genotype=$genotype_prefix.$chromosome
  local annotation=$ANNOTATIONS/${chromosome}_annotation.tsv
  local pgen=$chromosome_input/$(basename "$chromosome_genotype.pgen")
  local pvar=$chromosome_input/$(basename "$chromosome_genotype.pvar")
  local psam=$chromosome_input/$(basename "$chromosome_genotype.psam")
  mkdir -p "$chromosome_input" "$chromosome_output" "$gcloud_config"
  export CLOUDSDK_CONFIG=$gcloud_config

  if [ ! -s "$annotation" ]; then
    echo "annotation not found: $annotation" >&2
    exit 1
  fi

  echo ">>> $chromosome: genotype localization"
  local localization_started_ns
  localization_started_ns=$(clock_ns)
  if ! gcloud storage cp --billing-project "$project" \
    "$chromosome_genotype.pgen" \
    "$chromosome_genotype.pvar" \
    "$chromosome_genotype.psam" \
    "$chromosome_input/"; then
    record_chromosome_timing \
      "$chromosome_output/workflow_timing.csv" \
      "$chromosome" "$lane_index" \
      genotype_localization "$localization_started_ns" failure
    return 1
  fi
  record_chromosome_timing \
    "$chromosome_output/workflow_timing.csv" \
    "$chromosome" "$lane_index" \
    genotype_localization "$localization_started_ns"

  PGEN=$pgen \
  PVAR=$pvar \
  PSAM=$psam \
  ANNOTATION=$annotation \
  CHROMOSOME=$chromosome \
  CHROMOSOME_RESULTS=$chromosome_output \
  RESULT_GCS=${RESULT_GCS:+$RESULT_GCS/chromosomes/$chromosome} \
  PORT_BASE=$port_base \
  LANE_INDEX=$lane_index \
  CHROMOSOME_STARTED_NS=$localization_started_ns \
  WORK=$chromosome_work \
  REPO=$REPO \
  R_ENV=$R_ENV \
    bash "$REPO/rewrite/workbench/run_chromosome_task.sh"

  rm -rf "$chromosome_input" "$chromosome_work"
}

run_chromosome_lane() {
  local lane_index=$1
  local port_base=$((18000 + 10 * lane_index))
  local chromosome_index

  chromosome_index=$lane_index
  while [ "$chromosome_index" -lt "${#CHROMOSOME_NUMBERS[@]}" ]; do
    run_chromosome \
      "${CHROMOSOME_NUMBERS[$chromosome_index]}" "$port_base" "$lane_index"
    chromosome_index=$((chromosome_index + max_parallel_chromosomes))
  done
}

lane_count=$max_parallel_chromosomes
if [ "$lane_count" -gt "${#CHROMOSOME_NUMBERS[@]}" ]; then
  lane_count=${#CHROMOSOME_NUMBERS[@]}
fi

CHROMOSOME_PIDS=()
for ((lane_index = 0; lane_index < lane_count; lane_index++)); do
  run_chromosome_lane "$lane_index" &
  CHROMOSOME_PIDS+=("$!")
done

chromosome_failure=0
for lane_index in "${!CHROMOSOME_PIDS[@]}"; do
  if ! wait "${CHROMOSOME_PIDS[$lane_index]}"; then
    echo "chromosome lane failed: $lane_index" >&2
    chromosome_failure=1
  fi
done
if [ "$chromosome_failure" -ne 0 ]; then
  exit 1
fi

echo ">>> coordinator: aggregate results"
aggregate_started_ns=$(clock_ns)
mkdir -p "$WORK/matplotlib"
if ! PYTHONPATH="$REPO" MPLCONFIGDIR="$WORK/matplotlib" \
  python3 -m rewrite.workbench aggregate \
  --input-root "$CHROMOSOME_RESULTS" \
  --output-dir "$FINAL_RESULTS" \
  --chromosomes "${CHROMOSOME_NUMBERS[@]}"; then
  record_coordinator_timing aggregate_results "$aggregate_started_ns" failure
  exit 1
fi
record_coordinator_timing aggregate_results "$aggregate_started_ns"
record_coordinator_timing run_total "$coordinator_started_ns"

PYTHONPATH="$REPO" python3 -m rewrite.workbench timing-summary \
  --input-root "$CHROMOSOME_RESULTS" \
  --output-dir "$FINAL_RESULTS" \
  --coordinator-timing "$COORDINATOR_TIMING" \
  --chromosomes "${CHROMOSOME_NUMBERS[@]}"

touch "$FINAL_RESULTS/_SUCCESS" "$COORDINATOR_RESULTS/_SUCCESS"
echo ">>> coordinator: complete"
