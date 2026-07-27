#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")" && pwd)
CONFIG=$ROOT/fed_aou.conf
MODE=${1:-}
[ -z "$MODE" ] || [ "$MODE" = "--dry-run" ] \
  || { echo "usage: bash submit_fed_dsub.sh [--dry-run]" >&2; exit 2; }

resolve_gcs_root() {
  local root=$1
  local resource_id=$2

  if [ -z "$root" ]; then
    root=$(wb resource resolve --id "$resource_id" --format=TEXT \
      | awk '/^gs:\/\// { value=$0 } END { print value }')
  fi
  root=${root%/}
  [[ "$root" = gs://* ]] \
    || { echo "error: could not resolve GCS resource: $resource_id" >&2; exit 1; }
  printf '%s\n' "$root"
}

echo "=== [0/3] load parameters ==="
[ -f "$CONFIG" ] || { echo "error: config not found: $CONFIG" >&2; exit 1; }
: "${GOOGLE_CLOUD_PROJECT:?GOOGLE_CLOUD_PROJECT is not set}"
: "${PET_SA_EMAIL:?PET_SA_EMAIL is not set}"
# shellcheck disable=SC1090
. "$CONFIG"

if [ -z "${BATCH_IMAGE:-}" ]; then
  IMAGE_DIGEST=$(gcloud container images list-tags "$BATCH_IMAGE_REPOSITORY" \
    --limit=1 --sort-by='~timestamp' --format='value(digest)')
  : "${IMAGE_DIGEST:?could not resolve a Terra AoU image}"
  BATCH_IMAGE=$BATCH_IMAGE_REPOSITORY@$IMAGE_DIGEST
fi
printf '  config: %s\n  image:  %s\n' "$CONFIG" "$BATCH_IMAGE"

echo "=== [1/3] resolve GCS paths ==="
AOU_DATA_ROOT=$(resolve_gcs_root "${AOU_DATA_GCS_ROOT:-}" "$AOU_DATA_RESOURCE_ID")
PHENO_ROOT=$(resolve_gcs_root "${PHENO_GCS_ROOT:-}" "$PHENO_RESOURCE_ID")
GENO_GCS=$AOU_DATA_ROOT/${GENO_PGEN_PREFIX_RELATIVE#/}.$FED_CHR
COV_GCS=$AOU_DATA_ROOT/${COV_ANCESTRY_TSV_RELATIVE#/}
PHENO_GCS=$PHENO_ROOT/${PHENO_CSV_RELATIVE#/}

RUN_TAG=$(date -u +%Y%m%d-%H%M%S)
RUN_GCS=${BATCH_BUCKET%/}/${BATCH_RUN_PREFIX#/}/$RUN_TAG
BUNDLE_GCS=${BATCH_BUCKET%/}/${BATCH_CODE_PREFIX#/}/secure-skat-$RUN_TAG.tar.gz

printf '  genotype: %s.{pgen,pvar,psam}\n' "$GENO_GCS"
printf '  ancestry: %s\n  phenotype: %s\n' "$COV_GCS" "$PHENO_GCS"
printf '  results:  %s/results\n' "$RUN_GCS"

echo "=== [2/3] build and upload code ==="
if [ "$MODE" = "--dry-run" ]; then
  echo "  skipped"
else
  TMP_ROOT=$(mktemp -d)
  trap 'rm -rf "$TMP_ROOT"' EXIT
  STAGE=$TMP_ROOT/repo
  mkdir -p "$STAGE"

  CGO_ENABLED=0 GOOS=linux GOARCH=amd64 \
    go -C "$ROOT" build -trimpath -mod=vendor -o "$TMP_ROOT/sfgwas" .

  while IFS= read -r -d '' path; do
    mkdir -p "$STAGE/$(dirname "$path")"
    cp -a "$ROOT/$path" "$STAGE/$path"
  done < <(git -C "$ROOT" ls-files -z -- run_fed.sh scripts example_data/keys)
  sed "s|^BATCH_IMAGE=.*|BATCH_IMAGE=$BATCH_IMAGE|" "$CONFIG" > "$STAGE/fed_aou.conf"
  cp "$TMP_ROOT/sfgwas" "$STAGE/sfgwas"

  tar -C "$STAGE" -czf "$TMP_ROOT/secure-skat.tar.gz" .
  gcloud storage cp "$TMP_ROOT/secure-skat.tar.gz" "$BUNDLE_GCS" \
    --billing-project="$GOOGLE_CLOUD_PROJECT"
fi

echo "=== [3/3] submit dsub ==="
dsub_args=(
  --provider google-batch
  --project "$GOOGLE_CLOUD_PROJECT"
  --location "$BATCH_REGION"
  --service-account "$PET_SA_EMAIL"
  --network "projects/${GOOGLE_CLOUD_PROJECT}/global/networks/${BATCH_NETWORK}"
  --subnetwork "projects/${GOOGLE_CLOUD_PROJECT}/regions/${BATCH_REGION}/subnetworks/${BATCH_SUBNETWORK}"
  --use-private-address
  --user-project "$GOOGLE_CLOUD_PROJECT"
  --name "$BATCH_JOB_NAME"
  --timeout "${BATCH_TIMEOUT:-24h}"
  --logging "$RUN_GCS/logs"
  --image "$BATCH_IMAGE"
  --min-cores "$BATCH_MIN_CORES"
  --min-ram "$BATCH_MIN_RAM"
  --boot-disk-size "$BATCH_BOOT_DISK_SIZE"
  --disk-size "$BATCH_DISK_SIZE"
  --output-recursive "RESULTS=$RUN_GCS/results"
  --env
    "GOOGLE_CLOUD_PROJECT=$GOOGLE_CLOUD_PROJECT"
    "BATCH_WORK_ROOT=$BATCH_WORK_ROOT"
    "CODE_BUNDLE_GCS=$BUNDLE_GCS"
    "GENO_GCS=$GENO_GCS"
    "COV_GCS=$COV_GCS"
    "PHENO_GCS=$PHENO_GCS"
  --script "$ROOT/${BATCH_SCRIPT:-scripts/aou/run_fed_batch.sh}"
)

if [ "$MODE" = "--dry-run" ]; then
  printf '  '
  printf '%q ' dsub "${dsub_args[@]}"
  printf '\n'
else
  dsub "${dsub_args[@]}"
fi

echo "results:   $RUN_GCS/results"
echo "logs:      $RUN_GCS/logs"
echo "monitor: dstat --provider google-batch --project $GOOGLE_CLOUD_PROJECT --location $BATCH_REGION --status '*' --age 1d --summary"
