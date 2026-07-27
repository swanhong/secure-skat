#!/usr/bin/env bash
# Diagnostic drop-in for run_fed_batch.sh (pass as dsub --script).
# dsub uploads its own logs ONLY on success, so a Batch max-runtime kill (exit
# 50005) leaves nothing to read. This script pushes a growing log to
# $HEARTBEAT_GCS every ~15s, so even a killed job leaves a breadcrumb:
#   - "BEACON 0" present in GCS  => localization finished, the script ran; read later beacons.
#   - "BEACON 0" absent          => the hour was spent BEFORE the script (input download / setup).
# It also probes env + inputs, then runs the REAL pipeline at tiny scale so a
# healthy run finishes in minutes; if it still hangs, it is our code/network.
# No set -e: every probe must run and beacon its own status.

HB_LOCAL=/tmp/diag.log
: > "$HB_LOCAL"
push()   { [ -n "${HEARTBEAT_GCS:-}" ] && gcloud storage cp "$HB_LOCAL" "$HEARTBEAT_GCS/diag.log" >/dev/null 2>&1 || true; }
beacon() { printf '%s | %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$HB_LOCAL"; push; }
( while :; do sleep 15; push; done ) & HB_PID=$!
trap 'rc=$?; beacon "EXIT (code=$rc)"; kill $HB_PID 2>/dev/null; push' EXIT

beacon "BEACON 0: script started => localization runnables completed"
beacon "host: nproc=$(nproc) mem=$(free -g 2>/dev/null | awk '/Mem:/{print $2}')GB"

echo "--- env ---" >> "$HB_LOCAL"
for v in CODE_BUNDLE PGEN_PGEN PGEN_PVAR PGEN_PSAM ANCESTRY_TSV PHENO_CSV RESULTS BATCH_WORK_ROOT HEARTBEAT_GCS; do
  beacon "env $v=${!v:-<unset>}"
done

# Are the localized inputs actually present on disk, and how big?
for f in "${PGEN_PGEN:-}" "${PGEN_PVAR:-}" "${PGEN_PSAM:-}" "${ANCESTRY_TSV:-}" "${PHENO_CSV:-}" "${CODE_BUNDLE:-}"; do
  [ -n "$f" ] || continue
  if [ -e "$f" ]; then beacon "input OK  $(du -h "$f" 2>/dev/null | cut -f1)  $f"
  else                 beacon "input MISSING  $f"; fi
done

WORK_ROOT=${BATCH_WORK_ROOT:-/mnt/data/secure-skat-work}
REPO=$WORK_ROOT/repo; RUN=$WORK_ROOT/run
mkdir -p "$REPO" "$RUN" "$RESULTS"
beacon "extracting bundle -> $REPO"
if tar -xzf "$CODE_BUNDLE" -C "$REPO"; then beacon "bundle extracted OK"; else beacon "bundle extract FAILED"; fi

set -a; . "$REPO/fed_aou.conf"; set +a
if PLINK2=$(command -v plink2); then beacon "plink2 OK $PLINK2"; else beacon "plink2 MISSING"; fi
if python3 -c 'import numpy,sys;print(numpy.__version__)' >/tmp/np 2>&1; then beacon "numpy OK $(cat /tmp/np)"; else beacon "numpy FAIL $(cat /tmp/np)"; fi
[ -x "$REPO/sfgwas" ] && beacon "sfgwas binary present" || beacon "sfgwas binary MISSING/NOT-EXEC"
export PLINK2

# --- tiny real run: shrink so a healthy pipeline is minutes, not an hour ---
FED_PGEN=${PGEN_PGEN%.pgen}
export FED_PGEN
export FED_ANCESTRY=$ANCESTRY_TSV FED_PHENO=$PHENO_CSV FED_OUT=$RUN
export FED_KEYS=$REPO/example_data/keys SKIP_BUILD=1
export FED_SPLIT_ANCESTRY=0 FED_ANCESTRY_GROUP=EUR
export FED_NGENES=1 FED_NSUB=1000 FED_PROBES=0

beacon "BEACON 1: starting tiny run_fed.sh (split=0 EUR, NGENES=1, NSUB=1000)"
( cd "$REPO" && bash run_fed.sh ) 2>&1 | tee -a "$HB_LOCAL"
rc=${PIPESTATUS[0]}
beacon "BEACON 2: run_fed.sh exited rc=$rc"

cp -r "$RUN"/. "$RESULTS"/ 2>/dev/null || true
beacon "results copied to \$RESULTS; diag complete"
