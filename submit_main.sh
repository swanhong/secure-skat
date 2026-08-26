#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")" && pwd)

usage() {
  echo "usage: bash submit_main.sh run.conf"
}

case ${1:-} in
  --help|-h)
    usage
    exit 0
    ;;
  "")
    usage >&2
    exit 2
    ;;
esac
if [ "$#" -ne 1 ]; then
  usage >&2
  exit 2
fi

config_file=$1
if [ ! -f "$config_file" ]; then
  echo "config file not found: $config_file" >&2
  exit 1
fi

# Required settings. run.conf must provide these values.
annotation_dir=
data_bits=

# Optional analysis and protocol defaults.
phenotype_columns=LDLC_final_mgdl_6sd_masked,HDLC_mgdl_6sd_masked,TotChol_corrected_mvp_explicit_duration_mgdl_6sd_masked,ln_Trig_6sd_masked,nonHDL_corrected_mvp_explicit_duration_mgdl_6sd_masked
phenotype_id_column=person_id
covariate_columns=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16
covariate_id_column=research_id
covariate_array_column=pca_features
mask=annotation=pLoF
genes=all
samples_per_cohort=all
sample_seed=42
role_seed=42
ckks=PN14QP436S45
frac_bits=30
probes=50
seed=42
chromosomes=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22

# Optional Workbench and storage defaults.
project=${GOOGLE_CLOUD_PROJECT:-}
service_account=${PET_SA_EMAIL:-}
genotype_prefix=gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/exome/pgen/exome
phenotype=gs://gwas-data-wgs-wb-jaunty-blueberry-8679/pheno/v9_final_lipid_med_corrected_short_read_tot.csv
covariate=gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv
output_root=
r_env_archive=
region=us-west1
image=gcr.io/deeplearning-platform-release/base-cpu@sha256:394be02c6b020a39837e0719c546d56b829994ce091fde7970c232b7f16a6640
timeout=24h
min_cores=16
min_ram=64
boot_disk_size=50
disk_size=500
plink2_bin=

# Load the trusted local run configuration after defaults are set.
# shellcheck source=/dev/null
. "$config_file"

: "${annotation_dir:?run.conf must set annotation_dir}"
: "${covariate_columns:?covariate_columns must not be empty}"
: "${data_bits:?run.conf must set data_bits}"

IFS=, read -r -a CHROMOSOME_NUMBERS <<< "$chromosomes"
if [ "${#CHROMOSOME_NUMBERS[@]}" -eq 0 ]; then
  echo "chromosomes must contain at least one chromosome" >&2
  exit 1
fi
for chromosome_number in "${CHROMOSOME_NUMBERS[@]}"; do
  if [[ ! "$chromosome_number" =~ ^[0-9]+$ ]] ||
     [ "$chromosome_number" -lt 1 ] ||
     [ "$chromosome_number" -gt 22 ]; then
    echo "invalid chromosome number: $chromosome_number" >&2
    exit 1
  fi
done

# Validate required commands and settings.
for command_name in go gcloud dsub tar; do
  command -v "$command_name" >/dev/null || {
    echo "required command not found: $command_name" >&2
    exit 1
  }
done

if [ -z "$project" ]; then
  project=$(gcloud config get-value project 2>/dev/null || true)
fi
if [ -z "$plink2_bin" ]; then
  plink2_bin=$(command -v plink2 || true)
fi
if [ -z "$plink2_bin" ] && [ -n "${HOME:-}" ] && [ -x "$HOME/plink2" ]; then
  plink2_bin=$HOME/plink2
fi

: "${project:?set project in run.conf or configure a gcloud default project}"
: "${service_account:?set service_account in run.conf outside Workbench}"
: "${plink2_bin:?set plink2_bin in run.conf, add plink2 to PATH, or install it at ~/plink2}"

if [ -z "$output_root" ]; then
  output_root=gs://dataproc-staging-${project}/secure-skat
fi
if [ -z "$r_env_archive" ]; then
  r_env_archive=$output_root/runtime/skat-r44-skat225.tar.gz
fi

if ! gcloud storage objects describe \
  "$r_env_archive" \
  --billing-project "$project" >/dev/null 2>&1; then
  echo "R environment archive not found: $r_env_archive" >&2
  exit 1
fi

run_id=$(date -u +%Y%m%d-%H%M%S)
run_gcs=$output_root/runs/$run_id
code_gcs=$run_gcs/code.tar.gz
temporary=$(mktemp -d /tmp/secure-skat-submit.XXXXXX)
trap 'rm -rf "$temporary"' EXIT

mkdir -p "$temporary/repo/rewrite"
GOCACHE="$temporary/go-cache" CGO_ENABLED=0 GOOS=linux GOARCH=amd64 \
  go -C "$root" build -trimpath -mod=vendor \
  -o "$temporary/repo/secure-skat" ./rewrite/cmd/secure-skat

cp "$root/rewrite/__init__.py" "$temporary/repo/rewrite/"
cp -a "$root/rewrite/preprocessing" "$temporary/repo/rewrite/"
cp -a "$root/rewrite/workbench" "$temporary/repo/rewrite/"
cp -a "$root/rewrite/analysis" "$temporary/repo/rewrite/"
cp "$plink2_bin" "$temporary/repo/plink2"
chmod +x "$temporary/repo/secure-skat" "$temporary/repo/plink2"

genes_argument=all
if [ "$genes" != "all" ] && [ "$genes" != "ALL" ]; then
  cp "$genes" "$temporary/repo/genes.txt"
  genes_argument=genes.txt
fi

tar --exclude='__pycache__' --exclude='*.pyc' \
  -C "$temporary/repo" -czf "$temporary/code.tar.gz" .
gcloud storage cp "$temporary/code.tar.gz" "$code_gcs"
gcloud storage cp "$config_file" "$run_gcs/run.conf"

annotations_gcs=$run_gcs/inputs/annotations
for chromosome_number in "${CHROMOSOME_NUMBERS[@]}"; do
  chromosome=chr$chromosome_number
  annotation="$annotation_dir/${chromosome}_annotation.tsv"
  [ -f "$annotation" ] || {
    echo "missing annotation: $annotation" >&2
    exit 1
  }
  gcloud storage cp \
    "$annotation" \
    "$annotations_gcs/${chromosome}_annotation.tsv"
done

common_dsub=(
  --provider google-batch
  --project "$project"
  --location "$region"
  --service-account "$service_account"
  --network "projects/${project}/global/networks/network"
  --subnetwork "projects/${project}/regions/${region}/subnetworks/subnetwork"
  --use-private-address
  --user-project "$project"
  --image "$image"
  --timeout "$timeout"
  --min-cores "$min_cores"
  --min-ram "$min_ram"
  --boot-disk-size "$boot_disk_size"
  --disk-size "$disk_size"
)

echo ">>> submitting coordinator for ${#CHROMOSOME_NUMBERS[@]} chromosomes"
dsub "${common_dsub[@]}" \
  --name secure-skat-coordinator \
  --logging "$run_gcs/logs/coordinator" \
  --input CODE_BUNDLE="$code_gcs" \
          R_ENV_ARCHIVE="$r_env_archive" \
          PHENOTYPE="$phenotype" \
          COVARIATE="$covariate" \
  --input-recursive ANNOTATIONS="$annotations_gcs" \
  --output-recursive COORDINATOR_RESULTS="$run_gcs/results" \
  --env project="$project" \
        RESULT_GCS="$run_gcs/results" \
        genotype_prefix="$genotype_prefix" \
        chromosomes="$chromosomes" \
        phenotype_columns="$phenotype_columns" \
        covariate_columns="$covariate_columns" \
        phenotype_id_column="$phenotype_id_column" \
        covariate_id_column="$covariate_id_column" \
        covariate_array_column="$covariate_array_column" \
        mask="$mask" \
        genes="$genes_argument" \
        samples_per_cohort="$samples_per_cohort" \
        sample_seed="$sample_seed" \
        role_seed="$role_seed" \
        ckks="$ckks" \
        data_bits="$data_bits" \
        frac_bits="$frac_bits" \
        probes="$probes" \
        seed="$seed" \
  --script "$root/rewrite/workbench/run_coordinator_task.sh"

echo "run config:         $run_gcs/run.conf"
echo "chromosome results: $run_gcs/results/chromosomes"
echo "final results:      $run_gcs/results/final"
echo "logs:               $run_gcs/logs"
