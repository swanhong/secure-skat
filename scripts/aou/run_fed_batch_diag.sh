#!/usr/bin/env bash
# Diagnostic drop-in for run_fed_batch.sh. Downloads inputs in-script (requester-
# pays, via gsutil -u) instead of dsub --input, then runs the pipeline at tiny
# scale. NEVER fails (set +e; exit 0) so dsub delocalizes $RESULTS every time;
# full trace -> $RESULTS/diag_report.txt. Read: gs://.../results/diag_report.txt
set +e
RESULTS=${RESULTS:-/tmp/out}
mkdir -p "$RESULTS"
exec > >(tee -a "$RESULTS/diag_report.txt") 2>&1
set -x

echo "=== DIAG START $(date -u +%FT%TZ) ==="
command -v gsutil || echo "NO gsutil"
command -v gcloud || echo "NO gcloud"
echo "host: nproc=$(nproc) mem=$(free -g 2>/dev/null | awk '/Mem:/{print $2}')GB"
for v in GOOGLE_CLOUD_PROJECT CODE_BUNDLE_GCS GENO_GCS COV_GCS PHENO_GCS RESULTS BATCH_WORK_ROOT; do
  echo "env $v=${!v:-<unset>}"
done

WORK_ROOT=${BATCH_WORK_ROOT:-/mnt/data/secure-skat-work}
REPO=$WORK_ROOT/repo; RUN=$WORK_ROOT/run; DL=$WORK_ROOT/inputs
mkdir -p "$REPO" "$RUN" "$DL"

echo "=== download inputs (gsutil -u, requester-pays) ==="
gsutil -u "$GOOGLE_CLOUD_PROJECT" -m cp \
  "$CODE_BUNDLE_GCS" \
  "$GENO_GCS".pgen "$GENO_GCS".pvar "$GENO_GCS".psam \
  "$COV_GCS" "$PHENO_GCS" \
  "$DL/"
echo "=== gsutil rc=$? ==="
ls -la "$DL"

CODE_BUNDLE=$DL/$(basename "$CODE_BUNDLE_GCS")
tar -xzf "$CODE_BUNDLE" -C "$REPO" && echo "extract OK" || echo "extract FAILED rc=$?"
set -a; . "$REPO/fed_aou.conf"; set +a
command -v plink2 && echo "plink2 OK" || echo "NO plink2"
python3 -c 'import numpy;print("numpy",numpy.__version__)' || echo "NO numpy"
[ -x "$REPO/sfgwas" ] && echo "sfgwas present" || echo "sfgwas MISSING"
export PLINK2=$(command -v plink2)

echo "=== tiny pipeline (split=0 EUR, NGENES=1, NSUB=1000), 8-min cap ==="
export FED_PGEN=$DL/$(basename "$GENO_GCS")
export FED_ANCESTRY=$DL/$(basename "$COV_GCS")
export FED_PHENO=$DL/$(basename "$PHENO_GCS")
export FED_OUT=$RUN FED_KEYS=$REPO/example_data/keys SKIP_BUILD=1
export FED_SPLIT_ANCESTRY=0 FED_ANCESTRY_GROUP=EUR
export FED_NGENES=1 FED_NSUB=1000 FED_PROBES=0
( cd "$REPO" && timeout 480 bash run_fed.sh )
echo "=== run_fed.sh rc=$? (124 = 8-min timeout = HUNG) ==="

cp "$RUN"/party*.log "$RUN"/prep.log "$RESULTS"/ 2>/dev/null || true
cp -r "$RUN"/out "$RESULTS"/out 2>/dev/null || true
echo "=== DIAG END $(date -u +%FT%TZ) ==="
exit 0
