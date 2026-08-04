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

if [ "${FED_PROBES:-0}" -gt 0 ]; then
  R_ENV=$WORK/skat-r
  mkdir -p "$R_ENV"
  gcloud storage cp --billing-project "$GOOGLE_CLOUD_PROJECT" \
    "$FED_R_ENV_GCS" "$WORK/skat-r.tar.gz"
  tar -xzf "$WORK/skat-r.tar.gz" -C "$R_ENV"
  "$R_ENV/bin/conda-unpack"
  export PATH="$R_ENV/bin:$PATH"
fi

gsutil -u "$GOOGLE_CLOUD_PROJECT" -m cp "$COV_GCS" "$PHENO_GCS" "$INPUT/"

# Per-run inputs uploaded by submit_fed_dsub.sh.
# gcloud storage, not gsutil: a >150 MB upload lands as a COMPOSITE object, and gsutil refuses to
# download those unless crcmod's C extension is installed -- which it is not on the runner image.
if [ -n "${GENES_GCS:-}" ]; then
  gcloud storage cp --billing-project "$GOOGLE_CLOUD_PROJECT" "$GENES_GCS" "$INPUT/genes.txt"
  export FED_GENES="$INPUT/genes.txt"
fi

# plink2 ships in the bundle (the runner image has no plink2).
chmod +x "$REPO/plink2"
export PLINK2="$REPO/plink2"
export FED_ANCESTRY="$INPUT/$(basename "$COV_GCS")"
export FED_PHENO="$INPUT/$(basename "$PHENO_GCS")"
export FED_KEYS=$REPO/example_data/keys
export SKIP_BUILD=1

run_chromosome() {
  local chr=$1 annot_gcs=$2 chr_input=$INPUT/$1
  local chr_geno_gcs=${GENO_GCS%.*}.$chr
  mkdir -p "$chr_input"
  gsutil -u "$GOOGLE_CLOUD_PROJECT" -m cp \
    "$chr_geno_gcs".pgen "$chr_geno_gcs".pvar "$chr_geno_gcs".psam "$chr_input/"
  if [ -n "$annot_gcs" ]; then
    gcloud storage cp --billing-project "$GOOGLE_CLOUD_PROJECT" \
      "$annot_gcs" "$chr_input/annotation.tsv"
    export FED_ANNOT="$chr_input/annotation.tsv"
  fi
  export FED_CHR=$chr
  export FED_PGEN="$chr_input/$(basename "$chr_geno_gcs")"
  export FED_OUT=$RUN/$chr
  echo ">>> launching $chr (out=$FED_OUT)"
  bash "$REPO/run_fed.sh"
}

if [ -n "${FED_CHRS:-}" ]; then
  IFS=, read -ra CHRS <<< "$FED_CHRS"
  IFS=, read -ra ANNOTS <<< "${ANNOT_GCS_LIST:-}"
  [ "${#CHRS[@]}" -eq "${#ANNOTS[@]}" ] || {
    echo "FED_CHRS and uploaded annotations do not match" >&2
    exit 1
  }
  export FED_SKIP_PLOT=1
  CSV_FILES=()
  for i in "${!CHRS[@]}"; do
    run_chromosome "${CHRS[$i]}" "${ANNOTS[$i]}"
    CSV_FILES+=("$RUN/${CHRS[$i]}/fed_results.csv")
  done
  unset FED_SKIP_PLOT
  awk -F, 'BEGIN { OFS=","; gene=-1 } NR==1 { print; next } FNR==1 { next } { $1=++gene; print }' \
    "${CSV_FILES[@]}" > "$RUN/fed_results.csv"
  python3 "$REPO/scripts/analysis/fed_plot.py" "$RUN/fed_results.csv" "$RUN"
  echo "=== MULTI-CHROMOSOME COST ==="
  printf '  %-6s %10s %10s %10s %10s %10s %12s\n' chromosome prep build secure compare total network_MiB
  for chr in "${CHRS[@]}"; do
    network_bytes=$(awk -F, '$2=="protocol_total" && $7=="network_total" { print $8; exit }' \
      "$RUN/$chr/communication_summary.csv")
    awk -F, -v chr="$chr" -v bytes="${network_bytes:-0}" '
      $7=="step" && $1=="prep"    { prep=$3 }
      $7=="step" && $1=="build"   { build=$3 }
      $7=="step" && $1=="secure"  { secure=$3 }
      $7=="step" && $1=="compare" { compare=$3 }
      $7=="step" && $1=="total"   { total=$3 }
      END { printf "  %-6s %9.3fs %9.3fs %9.3fs %9.3fs %9.3fs %12.3f\n",
                   chr, prep/1000, build/1000, secure/1000, compare/1000, total/1000, bytes/1048576 }
    ' "$RUN/$chr/timing_steps.csv"
  done
else
  gsutil -u "$GOOGLE_CLOUD_PROJECT" -m cp \
    "$GENO_GCS".pgen "$GENO_GCS".pvar "$GENO_GCS".psam "$INPUT/"
  if [ -n "${ANNOT_GCS:-}" ]; then
    gcloud storage cp --billing-project "$GOOGLE_CLOUD_PROJECT" \
      "$ANNOT_GCS" "$INPUT/annotation.tsv"
    export FED_ANNOT="$INPUT/annotation.tsv"
  fi
  export FED_PGEN="$INPUT/$(basename "$GENO_GCS")"
  export FED_OUT=$RUN
  echo ">>> setup done (plink2=$PLINK2) — launching run_fed.sh"
  if [ "$DIAG" = 1 ]; then
    export FED_SPLIT_ANCESTRY=0 FED_ANCESTRY_GROUP=EUR
    export FED_NGENES=1 FED_NSUB=1000 FED_PROBES=1 FED_DATABITS=70
    timeout 480 bash "$REPO/run_fed.sh"
  else
    bash "$REPO/run_fed.sh"
  fi
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
    -name parameter.txt -o \
    -name '*.png' -o \
    -name 'skat_fed*_out.txt' \
  \) -print0
)
cp "$REPO/fed_aou.conf" "$RESULTS/fed_aou.conf"
