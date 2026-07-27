#!/usr/bin/env bash
set -euo pipefail
: "${GOOGLE_CLOUD_PROJECT:?missing GOOGLE_CLOUD_PROJECT}"
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
. "$ROOT/fed_aou.conf"

PLINK2=${PLINK2:-$HOME/plink2}
[ -x "$PLINK2" ] || { echo "missing executable plink2: $PLINK2" >&2; exit 1; }

CTX=$(mktemp -d)
trap 'rm -rf "$CTX"' EXIT
cp "$ROOT/docker/Dockerfile" "$CTX/"
cp "$PLINK2" "$CTX/plink2"

# GCR (gcr.io) auto-creates its backing bucket on first push, so no repo-create
# permission is needed (the PET SA is denied artifactregistry.repositories.create).
gcloud builds submit "$CTX" --tag "$BATCH_IMAGE" \
  --project "$GOOGLE_CLOUD_PROJECT" --region "$BATCH_REGION" \
  --gcs-source-staging-dir "$BATCH_ROOT/cloudbuild/source" \
  --gcs-log-dir "$BATCH_ROOT/cloudbuild/logs"
