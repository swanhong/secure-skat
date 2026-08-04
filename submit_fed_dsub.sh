#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")" && pwd)
: "${GOOGLE_CLOUD_PROJECT:?GOOGLE_CLOUD_PROJECT is not set}"
: "${PET_SA_EMAIL:?PET_SA_EMAIL is not set}"
. "$ROOT/fed_aou.conf"

case ${1:-} in
  "") DIAG=0 ;;
  --diag) DIAG=1 ;;
  *) echo "usage: bash submit_fed_dsub.sh [--diag]" >&2; exit 2 ;;
esac

RUN_GCS=$BATCH_ROOT/runs/$(date -u +%Y%m%d-%H%M%S)
CODE_GCS=$RUN_GCS/code.tar.gz
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP/repo/example_data"

CGO_ENABLED=0 GOOS=linux GOARCH=amd64 \
  go -C "$ROOT" build -trimpath -mod=vendor -o "$TMP/repo/sfgwas" .
cp -a "$ROOT/run_fed.sh" "$ROOT/fed_aou.conf" "$ROOT/scripts" "$TMP/repo/"
cp -a "$ROOT/example_data/keys" "$TMP/repo/example_data/"
cp -a "${PLINK2:-$HOME/plink2}" "$TMP/repo/plink2"   # runner image has no plink2
tar --exclude='__pycache__' --exclude='*.pyc' \
  -C "$TMP/repo" -czf "$TMP/code.tar.gz" .
gcloud storage cp "$TMP/code.tar.gz" "$CODE_GCS"

# Per-run inputs kept out of the code bundle: the annotation table is VAT-derived (Controlled Tier)
# and can be hundreds of MB, so it goes to the same in-project bucket as everything else here.
EXTRA_ENV=()
if [ -n "${FED_CHRS:-}" ]; then
  : "${FED_ANNOT_DIR:?FED_ANNOT_DIR is required with FED_CHRS}"
  IFS=, read -ra CHRS <<< "$FED_CHRS"
  ANNOT_GCS_LIST=
  for chr in "${CHRS[@]}"; do
    annot="$FED_ANNOT_DIR/${chr}_annotation.tsv"
    [ -f "$annot" ] || { echo "missing annotation: $annot" >&2; exit 1; }
    annot_gcs="$RUN_GCS/${chr}_annotation.tsv"
    gcloud storage cp "$annot" "$annot_gcs"
    ANNOT_GCS_LIST="${ANNOT_GCS_LIST:+$ANNOT_GCS_LIST,}$annot_gcs"
  done
  EXTRA_ENV+=("ANNOT_GCS_LIST=$ANNOT_GCS_LIST")
elif [ -n "${FED_ANNOT:-}" ]; then
  gcloud storage cp "$FED_ANNOT" "$RUN_GCS/$(basename "$FED_ANNOT")"
  EXTRA_ENV+=("ANNOT_GCS=$RUN_GCS/$(basename "$FED_ANNOT")")
fi
if [ -n "${FED_GENES:-}" ] && [ -f "$FED_GENES" ]; then
  gcloud storage cp "$FED_GENES" "$RUN_GCS/$(basename "$FED_GENES")"
  EXTRA_ENV+=("GENES_GCS=$RUN_GCS/$(basename "$FED_GENES")")
elif [ -n "${FED_GENES:-}" ]; then
  EXTRA_ENV+=("FED_GENES=$FED_GENES")   # "ALL", or a path that already exists on the VM
fi

# Forward every FED_* the caller set so fed_aou.conf's defaults can be overridden per run.
for v in $(compgen -v FED_); do
  case $v in FED_ANNOT|FED_ANNOT_DIR|FED_GENES) continue ;; esac
  [ -n "${!v}" ] && EXTRA_ENV+=("$v=${!v}")
done

dsub \
  --provider google-batch \
  --project "$GOOGLE_CLOUD_PROJECT" --location "$BATCH_REGION" \
  --service-account "$PET_SA_EMAIL" \
  --network "projects/${GOOGLE_CLOUD_PROJECT}/global/networks/network" \
  --subnetwork "projects/${GOOGLE_CLOUD_PROJECT}/regions/${BATCH_REGION}/subnetworks/subnetwork" \
  --use-private-address --user-project "$GOOGLE_CLOUD_PROJECT" \
  --name "${BATCH_NAME:-secure-skat-fed}" --timeout "${BATCH_TIMEOUT_OVERRIDE:-$BATCH_TIMEOUT}" \
  --logging "$RUN_GCS/logs" --image "${BATCH_IMAGE_OVERRIDE:-$BATCH_IMAGE}" \
  --min-cores "$BATCH_MIN_CORES" --min-ram "$BATCH_MIN_RAM" \
  --boot-disk-size "$BATCH_BOOT_DISK_SIZE" --disk-size "$BATCH_DISK_SIZE" \
  --output-recursive "RESULTS=$RUN_GCS/results" \
  --env "GOOGLE_CLOUD_PROJECT=$GOOGLE_CLOUD_PROJECT" \
        "CODE_BUNDLE_GCS=$CODE_GCS" "BATCH_DIAG=$DIAG" \
        "LIVE_LOG_GCS=$RUN_GCS/live.log" \
        ${EXTRA_ENV[@]+"${EXTRA_ENV[@]}"} \
  --script "$ROOT/scripts/aou/run_fed_batch.sh"

echo "results: $RUN_GCS/results"
echo "logs:    $RUN_GCS/logs"
echo "live:    $RUN_GCS/live.log"
