#!/usr/bin/env bash
set -euo pipefail

DIAG=${BATCH_DIAG:-0}
WORK=/mnt/data/secure-skat
REPO=$WORK/repo
INPUT=$WORK/input
RUN=$WORK/run
mkdir -p "$REPO" "$INPUT" "$RUN" "$RESULTS"

LOG=$RESULTS/$([ "$DIAG" = 1 ] && echo diag_report.txt || echo run_log.txt)
exec > "$LOG" 2>&1

# dsub uploads $RESULTS only at the end (and only on success), so live-stream the
# log to GCS every 20s. Watch: watch -n20 "gcloud storage cat $LIVE_LOG_GCS | tail -40"
if [ -n "${LIVE_LOG_GCS:-}" ]; then
  ( while :; do sleep 20; gsutil -q cp "$LOG" "$LIVE_LOG_GCS" 2>/dev/null || true; done ) &
  STREAM_PID=$!
fi

finish() {
  local rc=$?
  trap - EXIT
  [ -n "${STREAM_PID:-}" ] && kill "$STREAM_PID" 2>/dev/null || true
  cp "$RUN"/prep.log "$RUN"/party*.log "$RESULTS"/ 2>/dev/null || true
  echo "exit=$rc"
  [ -n "${LIVE_LOG_GCS:-}" ] && gsutil -q cp "$LOG" "$LIVE_LOG_GCS" 2>/dev/null || true
  # diag always succeeds so $RESULTS delocalizes; real run keeps its real exit code
  # (its log already lives in $LIVE_LOG_GCS even on failure).
  [ "$DIAG" = 1 ] && exit 0 || exit "$rc"
}
trap finish EXIT

gsutil -u "$GOOGLE_CLOUD_PROJECT" cp "$CODE_BUNDLE_GCS" "$WORK/code.tar.gz"
tar -xzf "$WORK/code.tar.gz" -C "$REPO"

set -a
. "$REPO/fed_aou.conf"
set +a

gsutil -u "$GOOGLE_CLOUD_PROJECT" -m cp \
  "$GENO_GCS".pgen "$GENO_GCS".pvar "$GENO_GCS".psam \
  "$COV_GCS" "$PHENO_GCS" "$INPUT/"

# Per-run inputs uploaded by submit_fed_dsub.sh.
# gcloud storage, not gsutil: a >150 MB upload lands as a COMPOSITE object, and gsutil refuses to
# download those unless crcmod's C extension is installed -- which it is not on the runner image.
if [ -n "${ANNOT_GCS:-}" ]; then
  gcloud storage cp --billing-project "$GOOGLE_CLOUD_PROJECT" "$ANNOT_GCS" "$INPUT/annotation.tsv"
  export FED_ANNOT="$INPUT/annotation.tsv"
fi
if [ -n "${GENES_GCS:-}" ]; then
  gcloud storage cp --billing-project "$GOOGLE_CLOUD_PROJECT" "$GENES_GCS" "$INPUT/genes.txt"
  export FED_GENES="$INPUT/genes.txt"
fi

# plink2 ships in the bundle (the runner image has no plink2).
chmod +x "$REPO/plink2"
export PLINK2="$REPO/plink2"
export FED_PGEN="$INPUT/$(basename "$GENO_GCS")"
export FED_ANCESTRY="$INPUT/$(basename "$COV_GCS")"
export FED_PHENO="$INPUT/$(basename "$PHENO_GCS")"
export FED_OUT=$RUN
export FED_KEYS=$REPO/example_data/keys
export SKIP_BUILD=1
echo ">>> setup done (plink2=$PLINK2) — launching run_fed.sh"

if [ "$DIAG" = 1 ]; then
  export FED_SPLIT_ANCESTRY=0 FED_ANCESTRY_GROUP=EUR
  export FED_NGENES=1 FED_NSUB=1000 FED_PROBES=1   # fast smoke; for full-n validation run normally with FED_NGENES=1
  timeout 480 bash "$REPO/run_fed.sh"
else
  bash "$REPO/run_fed.sh"
fi

# Aggregate result files leave the Batch VM here; diagnostic logs (prep/party) also
# exit via finish(), and live.log streams during the run — all intended output.
cd "$RUN"
while IFS= read -r -d '' file; do
  rel=${file#./}
  mkdir -p "$RESULTS/$(dirname "$rel")"
  cp "$file" "$RESULTS/$rel"
done < <(
  find . -type f \( \
    -name fed_results.csv -o \
    -name timing_steps.csv -o \
    -name communication_summary.csv -o \
    -name manifest.json -o \
    -name '*.png' -o \
    -name 'skat_fed*_out.txt' \
  \) -print0
)
cp "$REPO/fed_aou.conf" "$RESULTS/fed_aou.conf"
