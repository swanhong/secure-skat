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

# Artifact Registry doesn't auto-create the repo on push (unlike gcr.io); make it.
gcloud artifacts repositories create "$(echo "$BATCH_IMAGE" | cut -d/ -f3)" \
  --repository-format=docker --location="$BATCH_REGION" \
  --project "$GOOGLE_CLOUD_PROJECT" 2>/dev/null || true

gcloud builds submit "$CTX" --tag "$BATCH_IMAGE" --project "$GOOGLE_CLOUD_PROJECT"
