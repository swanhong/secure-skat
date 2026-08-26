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

WORK=/mnt/data/secure-skat
REPO=$WORK/repo
INPUT=$WORK/input
PREPROCESSED=$WORK/preprocessed
RUN=$WORK/run
R_ENV=$WORK/r-env
mkdir -p "$REPO" "$INPUT" "$RUN" "$CHROMOSOME_RESULTS"

LOG=$CHROMOSOME_RESULTS/run.log
exec > >(tee -a "$LOG") 2>&1

PARTY_PIDS=()
finish() {
  local rc=$?
  trap - EXIT
  for process in "${PARTY_PIDS[@]:-}"; do
    kill "$process" 2>/dev/null || true
  done
  cp "$RUN"/party*.log "$CHROMOSOME_RESULTS/" 2>/dev/null || true
  if [ "$rc" -ne 0 ] && [ -n "${RESULT_GCS:-}" ]; then
    gcloud storage rsync --recursive "$CHROMOSOME_RESULTS" "$RESULT_GCS" \
      --billing-project "$project" 2>/dev/null || true
  fi
  exit "$rc"
}
trap finish EXIT

tar -xzf "$CODE_BUNDLE" -C "$REPO"
ln -s "$PGEN" "$INPUT/genotype.pgen"
ln -s "$PVAR" "$INPUT/genotype.pvar"
ln -s "$PSAM" "$INPUT/genotype.psam"

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
)
for mask in "${MASKS[@]}"; do
  PREPARE+=(--mask "$mask")
done

echo ">>> $CHROMOSOME: preprocessing"
if ! PYTHONPATH="$REPO" "${PREPARE[@]}"; then
  cp "$PREPROCESSED/validation_errors.csv" "$CHROMOSOME_RESULTS/" 2>/dev/null || true
  exit 1
fi

SHARED_KEYS=$RUN/shared-keys
PYTHONPATH="$REPO" python3 -m rewrite.workbench keys --output-dir "$SHARED_KEYS"

echo ">>> $CHROMOSOME: secure protocol"
for party in 0 1 2; do
  "$REPO/secure-skat" run \
    --party "$party" \
    --input "$PREPROCESSED" \
    --output "$RUN/secure_results.csv" \
    --port-base 18000 \
    --shared-keys "$SHARED_KEYS" \
    --ckks "${ckks:?}" \
    --data-bits "${data_bits:?}" \
    --frac-bits "${frac_bits:?}" \
    --probes "${probes:?}" \
    --seed "${seed:?}" \
    > "$RUN/party${party}.log" 2>&1 &
  PARTY_PIDS+=("$!")
done

SECURE_RC=0
for process in "${PARTY_PIDS[@]}"; do
  if ! wait "$process"; then
    SECURE_RC=1
  fi
done
PARTY_PIDS=()
if [ "$SECURE_RC" -ne 0 ]; then
  exit "$SECURE_RC"
fi

echo ">>> $CHROMOSOME: R::SKAT Burden and Davies"
mkdir -p "$R_ENV"
tar -xzf "$R_ENV_ARCHIVE" -C "$R_ENV"
"$R_ENV/bin/conda-unpack"
"$R_ENV/bin/Rscript" "$REPO/rewrite/workbench/run_reference_skat.R" \
  "$PREPROCESSED" "$RUN/r_results.csv"

echo ">>> $CHROMOSOME: join secure and R results"
PYTHONPATH="$REPO" python3 -m rewrite.workbench join \
  --secure "$RUN/secure_results.csv" \
  --r-results "$RUN/r_results.csv" \
  --output "$CHROMOSOME_RESULTS/gene_results.csv" \
  --summary "$CHROMOSOME_RESULTS/error_summary.csv"

cp "$PREPROCESSED/metadata.json" "$CHROMOSOME_RESULTS/metadata.json"
cp "$RUN/secure_results.csv" "$CHROMOSOME_RESULTS/secure_results.csv"
cp "$RUN/r_results.csv" "$CHROMOSOME_RESULTS/r_results.csv"
touch "$CHROMOSOME_RESULTS/_SUCCESS"
echo ">>> $CHROMOSOME: complete"
