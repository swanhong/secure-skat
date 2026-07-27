#!/usr/bin/env bash
# Diagnostic drop-in for run_fed_batch.sh. Downloads inputs in-script (gsutil -u)
# and runs the pipeline at tiny scale. Writes a full trace to $RESULTS/diag_report.txt
# with a plain group redirect (NO process substitution, which was crashing earlier
# versions before anything got written) and ALWAYS exits 0 so dsub delocalizes
# $RESULTS every time. Every blocking op is wrapped in `timeout` so a hang can't
# stop the upload. Read: gs://.../results/diag_report.txt
RESULTS=${RESULTS:-/tmp/out}
mkdir -p "$RESULTS" 2>/dev/null
LOG=$RESULTS/diag_report.txt

{
  echo "=== runnable-3 RAN $(date -u +%FT%TZ) ==="
  echo "shell=$0 bash=${BASH_VERSION:-none}"
  command -v gsutil || echo "NO gsutil"
  command -v gcloud || echo "NO gcloud"
  echo "host: nproc=$(nproc 2>/dev/null) mem=$(free -g 2>/dev/null | awk '/Mem:/{print $2}')GB"
  echo "GOOGLE_CLOUD_PROJECT=$GOOGLE_CLOUD_PROJECT"
  echo "CODE_BUNDLE_GCS=$CODE_BUNDLE_GCS"
  echo "GENO_GCS=$GENO_GCS  COV_GCS=$COV_GCS  PHENO_GCS=$PHENO_GCS"
  echo "RESULTS=$RESULTS  BATCH_WORK_ROOT=${BATCH_WORK_ROOT:-unset}"

  WORK_ROOT=${BATCH_WORK_ROOT:-/mnt/data/secure-skat-work}
  REPO=$WORK_ROOT/repo; RUN=$WORK_ROOT/run; DL=$WORK_ROOT/inputs
  mkdir -p "$REPO" "$RUN" "$DL"

  echo "=== download inputs (timeout 300, gsutil -u) ==="
  timeout 300 gsutil -u "$GOOGLE_CLOUD_PROJECT" -m cp \
    "$CODE_BUNDLE_GCS" \
    "$GENO_GCS".pgen "$GENO_GCS".pvar "$GENO_GCS".psam \
    "$COV_GCS" "$PHENO_GCS" \
    "$DL/"
  echo "=== gsutil rc=$? (124 = timed out) ==="
  ls -la "$DL"

  CODE_BUNDLE=$DL/$(basename "$CODE_BUNDLE_GCS")
  tar -xzf "$CODE_BUNDLE" -C "$REPO" && echo "extract OK" || echo "extract rc=$?"
  . "$REPO/fed_aou.conf" 2>/dev/null
  command -v plink2 || echo "NO plink2"
  python3 -c 'import numpy;print("numpy",numpy.__version__)' 2>&1
  [ -x "$REPO/sfgwas" ] && echo "sfgwas present" || echo "sfgwas MISSING"

  echo "=== tiny pipeline (split=0 EUR, NGENES=1, NSUB=1000, timeout 480) ==="
  export PLINK2=$(command -v plink2)
  export FED_PGEN=$DL/$(basename "$GENO_GCS")
  export FED_ANCESTRY=$DL/$(basename "$COV_GCS")
  export FED_PHENO=$DL/$(basename "$PHENO_GCS")
  export FED_OUT=$RUN FED_KEYS=$REPO/example_data/keys SKIP_BUILD=1
  export FED_SPLIT_ANCESTRY=0 FED_ANCESTRY_GROUP=EUR
  export FED_NGENES=1 FED_NSUB=1000 FED_PROBES=0
  ( cd "$REPO" && timeout 480 bash run_fed.sh )
  echo "=== run_fed.sh rc=$? (124 = 8-min timeout = HUNG) ==="

  cp "$RUN"/party*.log "$RUN"/prep.log "$RESULTS"/ 2>/dev/null
  echo "=== DIAG END $(date -u +%FT%TZ) ==="
} >> "$LOG" 2>&1

exit 0
