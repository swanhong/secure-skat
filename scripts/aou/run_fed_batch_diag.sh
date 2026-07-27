#!/usr/bin/env bash
# Diagnostic drop-in for run_fed_batch.sh (pass as dsub --script).
# On AoU we can't read failure logs: dsub uploads its GCS logs only on SUCCESS,
# and the PET service account is denied Cloud Logging. So this script:
#   1. NEVER fails (set +e; always exit 0) -> dsub delocalizes $RESULTS every time.
#   2. Traces everything (set -x) into $RESULTS/diag_report.txt via the $RESULTS
#      mount, which dsub uploads with NO gcloud needed inside the container.
#   3. Bounds the inner pipeline with `timeout` so a hang can't block the upload.
# Read the result afterwards from:  gs://.../results/diag_report.txt
set +e

RESULTS=${RESULTS:-/tmp/out}
mkdir -p "$RESULTS"
REPORT=$RESULTS/diag_report.txt
exec > >(tee -a "$REPORT") 2>&1   # everything (incl. set -x) -> report + stdout
set -x

echo "=== DIAG START $(date -u +%FT%TZ) ==="

# --- can this image even talk to GCS? (answers the gcloud-missing question) ---
echo "--- tooling ---"
command -v gcloud && gcloud --version 2>&1 | head -1 || echo "NO gcloud"
command -v gsutil && gsutil version 2>&1 | head -1 || echo "NO gsutil"
command -v timeout || echo "NO timeout"
echo "host: nproc=$(nproc) mem=$(free -g 2>/dev/null | awk '/Mem:/{print $2}')GB disk=$(df -h /mnt/data 2>/dev/null | tail -1)"

echo "--- env ---"
for v in CODE_BUNDLE PGEN_PGEN PGEN_PVAR PGEN_PSAM ANCESTRY_TSV PHENO_CSV RESULTS BATCH_WORK_ROOT; do
  echo "$v=${!v:-<unset>}"
done

echo "--- inputs on disk (did localization actually place them?) ---"
for f in "${PGEN_PGEN:-}" "${PGEN_PVAR:-}" "${PGEN_PSAM:-}" "${ANCESTRY_TSV:-}" "${PHENO_CSV:-}" "${CODE_BUNDLE:-}"; do
  [ -n "$f" ] || continue
  if [ -e "$f" ]; then echo "OK      $(du -h "$f" 2>/dev/null | cut -f1)  $f"
  else                 echo "MISSING $f"; fi
done

WORK_ROOT=${BATCH_WORK_ROOT:-/mnt/data/secure-skat-work}
REPO=$WORK_ROOT/repo; RUN=$WORK_ROOT/run
mkdir -p "$REPO" "$RUN"
echo "--- extract bundle ---"
tar -xzf "$CODE_BUNDLE" -C "$REPO" && echo "extract OK" || echo "extract FAILED rc=$?"
ls -la "$REPO" | head -30

echo "--- conf + deps ---"
set -a; . "$REPO/fed_aou.conf"; set +a
PLINK2=$(command -v plink2) && echo "plink2 $PLINK2" || echo "NO plink2"
python3 -c 'import numpy;print("numpy",numpy.__version__)' || echo "NO numpy"
[ -x "$REPO/sfgwas" ] && echo "sfgwas present" || echo "sfgwas MISSING/NOT-EXEC"
"$REPO/sfgwas" 2>&1 | head -3   # does the Go binary even load on this VM?
export PLINK2

echo "=== tiny pipeline (split=0 EUR, NGENES=1, NSUB=1000), capped at 8 min ==="
export FED_PGEN=${PGEN_PGEN%.pgen}
export FED_ANCESTRY=$ANCESTRY_TSV FED_PHENO=$PHENO_CSV FED_OUT=$RUN
export FED_KEYS=$REPO/example_data/keys SKIP_BUILD=1
export FED_SPLIT_ANCESTRY=0 FED_ANCESTRY_GROUP=EUR
export FED_NGENES=1 FED_NSUB=1000 FED_PROBES=0

( cd "$REPO" && timeout 480 bash run_fed.sh )
echo "=== run_fed.sh rc=$? (124 = hit the 8-min timeout = HUNG) ==="

echo "--- collect party logs ---"
cp "$RUN"/party*.log "$RUN"/prep.log "$RESULTS"/ 2>/dev/null || true
cp -r "$RUN"/out "$RESULTS"/out 2>/dev/null || true

echo "=== DIAG END $(date -u +%FT%TZ) ==="
exit 0   # force success so dsub delocalizes $RESULTS
