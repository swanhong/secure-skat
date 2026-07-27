#!/bin/bash
# End-to-end federated-private SKAT (#4 benchmark): prep -> build -> secure 3-party -> compare.
#
#   PLINK2=$HOME/plink2 bash run_fed.sh                                  # full run → $HOME/runs/out<YYMMDDHHMMSS>
#   FED_SPLIT_ANCESTRY=1 FED_OUT=$HOME/runs/by_ancestry bash run_fed.sh  # run EUR/AFR/AMR separately
#   PLINK2=$HOME/plink2 FED_NSUB=38000 FED_DATABITS=100 bash run_fed.sh  # scale n (needs more fixed-point range)
#   FED_PREP_SRC=$HOME/fed_prep_out SKIP_PREP=1 SKIP_BUILD=1 bash run_fed.sh   # reuse a prior prep, archive results in a fresh dir
#   FED_OUT=$HOME/runs/out260712202430 SKIP_PREP=1 SKIP_BUILD=1 bash run_fed.sh # re-run in an existing archived dir
#
# Each run lands in its own timestamped dir ($HOME/runs/out<YYMMDDHHMMSS>) unless FED_OUT is set, so past
# results are kept: fed_in (blocks/config) + fed_out (secure out/, logs, fed_results.csv) live together.
# FED_* / PLINK2 env vars are BAKED into <rundir>/config by fed_prep (data_bits, ports, ckks, dims), so
# FED_DATABITS etc. only take effect when fed_prep runs; with SKIP_PREP=1 the existing config is reused as-is
# (FED_PREP_SRC seeds it into the fresh dir). SKIP_BUILD=1 skips go build. FED_CSV=1 makes step [4/4] also dump
# fed_results.csv (per-gene positions + secure/plain p-values) for scripts/analysis/fed_plot.py.
# FED_PROBES=N (N>0) turns ON the SKAT p-value (exact traces when m_pub<=N, else Hutchinson → WH z;
# 0 = Q+burden p only);
set -eo pipefail

REPO=$(cd "$(dirname "$0")" && pwd)
if [ "${FED_SPLIT_ANCESTRY:-0}" = "1" ]; then
  : "${FED_OUT:=$HOME/runs/out$(date +%y%m%d%H%M%S)}"
  split_skip_build=${SKIP_BUILD:-}
  echo "=== ancestry split: EUR, AFR, AMR -> $FED_OUT ==="
  for ancestry in EUR AFR AMR; do
    echo "=== ancestry: $ancestry ==="
    FED_SPLIT_ANCESTRY=0 \
    FED_ANCESTRY_GROUP=$ancestry \
    FED_OUT="$FED_OUT/$ancestry" \
    FED_PREP_SRC="${FED_PREP_SRC:+$FED_PREP_SRC/$ancestry}" \
    SKIP_BUILD="$split_skip_build" \
      bash "$REPO/run_fed.sh"
    split_skip_build=1
  done
  echo "=== ancestry split complete: $FED_OUT ==="
  exit 0
fi

# Each run lives in its own timestamped dir so past results are kept: FED_OUT unset → $HOME/runs/out<YYMMDDHHMMSS>
# holds fed_in (prep blocks/config) + fed_out (secure results, compare csv). Set FED_OUT to reuse/target a dir.
# FED_PREP_SRC=<prior run dir>: with SKIP_PREP, seed this run's fed_in from there (symlink blocks, copy+repath
# config) so a slow prep is reused while results still land in the fresh archived dir.
: "${FED_OUT:=$HOME/runs/out$(date +%y%m%d%H%M%S)}"
export FED_OUT # fed_prep.py reads it too
OUT=$FED_OUT
CFG=$OUT/config
PREP=$REPO/scripts/preprocessing
T_PREP_MS=0 T_BUILD_MS=0 T_SECURE_MS=0 T_COMPARE_MS=0
S_PREP=not_run S_BUILD=not_run S_SECURE=not_run S_COMPARE=not_run
CURRENT_STEP="" STEP_START_MS=0 TIMING_WRITTEN=0
PARTY_PIDS=()
PARTY_LOG_TAIL_PID=""
mkdir -p "$OUT"
rm -f "$OUT/prep.log" "$OUT/compare.log" \
  "$OUT/party0.log" "$OUT/party1.log" "$OUT/party2.log" \
  "$OUT/communication_summary.csv"

now_ms() {
  python3 -c 'import time; print(time.monotonic_ns() // 1000000)' 2>/dev/null \
    || echo "$(($(date +%s) * 1000))"
}

start_step() {
  CURRENT_STEP=$1
  STEP_START_MS=$(now_ms)
  case "$CURRENT_STEP" in
    prep) S_PREP=running ;;
    build) S_BUILD=running ;;
    secure) S_SECURE=running ;;
    compare) S_COMPARE=running ;;
  esac
}

finish_step() {
  local status=${2:-done}
  local elapsed=$(( $(now_ms) - STEP_START_MS ))
  case "$1" in
    prep) T_PREP_MS=$elapsed; S_PREP=$status ;;
    build) T_BUILD_MS=$elapsed; S_BUILD=$status ;;
    secure) T_SECURE_MS=$elapsed; S_SECURE=$status ;;
    compare) T_COMPARE_MS=$elapsed; S_COMPARE=$status ;;
  esac
  CURRENT_STEP=""
  STEP_START_MS=0
}

duration_s() {
  awk -v ms="$1" 'BEGIN { printf "%.3fs", ms / 1000 }'
}

write_timing_summary() {
  local exit_code=${1:-0}
  local total_ms=$((T_PREP_MS + T_BUILD_MS + T_SECURE_MS + T_COMPARE_MS))
  echo "=== TIMING (steps; status is explicit for skipped/failed phases) ==="
  printf '  [1] prep     %-8s %10s\n' "$S_PREP" "$(duration_s "$T_PREP_MS")"
  printf '  [2] build    %-8s %10s\n' "$S_BUILD" "$(duration_s "$T_BUILD_MS")"
  printf '  [3] secure   %-8s %10s\n' "$S_SECURE" "$(duration_s "$T_SECURE_MS")"
  printf '  [4] compare  %-8s %10s\n' "$S_COMPARE" "$(duration_s "$T_COMPARE_MS")"
  printf '  --------------------------------\n  total               %10s\n' "$(duration_s "$total_ms")"
  {
    echo "step,status,milliseconds"
    printf 'prep,%s,%d\n' "$S_PREP" "$T_PREP_MS"
    printf 'build,%s,%d\n' "$S_BUILD" "$T_BUILD_MS"
    printf 'secure,%s,%d\n' "$S_SECURE" "$T_SECURE_MS"
    printf 'compare,%s,%d\n' "$S_COMPARE" "$T_COMPARE_MS"
    printf 'total,exit_%d,%d\n' "$exit_code" "$total_ms"
  } > "$OUT/timing_steps.csv"

  # Merge machine-readable substeps emitted by prep, each secure party, and
  # comparison.  A parser failure must not hide the original run status.
  local timing_args=("$OUT/timing_steps.csv" --output "$OUT/timing_steps.csv")
  [ -f "$OUT/prep.log" ] && timing_args+=(--prep-log "$OUT/prep.log")
  [ -f "$OUT/party0.log" ] && timing_args+=(--party-log "$OUT/party0.log")
  [ -f "$OUT/party1.log" ] && timing_args+=(--party-log "$OUT/party1.log")
  [ -f "$OUT/party2.log" ] && timing_args+=(--party-log "$OUT/party2.log")
  [ -f "$OUT/compare.log" ] && timing_args+=(--compare-log "$OUT/compare.log")
  python3 "$REPO/scripts/analysis/timing_summary.py" "${timing_args[@]}" \
    || echo "warning: failed to append detailed timing rows" >&2
  echo "  timing CSV -> $OUT/timing_steps.csv"
  TIMING_WRITTEN=1
}

cleanup_party_processes() {
  local child_pid
  if [ -n "$PARTY_LOG_TAIL_PID" ]; then
    kill "$PARTY_LOG_TAIL_PID" 2>/dev/null || true
    kill -9 "$PARTY_LOG_TAIL_PID" 2>/dev/null || true
    wait "$PARTY_LOG_TAIL_PID" 2>/dev/null || true
    PARTY_LOG_TAIL_PID=""
  fi
  if [ "${#PARTY_PIDS[@]}" -eq 0 ]; then
    return
  fi

  # A failed/interrupting party can leave its peers blocked in network I/O.
  # Do not wait in an EXIT/signal path: a stuck child must not prevent the
  # timing summary or shell exit. SIGKILL makes the cleanup bounded.
  for child_pid in "${PARTY_PIDS[@]}"; do
    kill "$child_pid" 2>/dev/null || true
  done
  for child_pid in "${PARTY_PIDS[@]}"; do
    kill -9 "$child_pid" 2>/dev/null || true
  done
  PARTY_PIDS=()
}

on_exit() {
  local exit_code=$1
  trap - EXIT INT TERM
  cleanup_party_processes
  if [ "$TIMING_WRITTEN" -eq 0 ]; then
    if [ -n "$CURRENT_STEP" ]; then
      finish_step "$CURRENT_STEP" failed
    fi
    write_timing_summary "$exit_code"
  fi
  exit "$exit_code"
}
trap 'on_exit $?' EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

if [ -n "$SKIP_PREP" ] && [ -n "$FED_PREP_SRC" ] && [ ! -f "$CFG/configGlobal.toml" ]; then
  echo "=== seeding fed_in from $FED_PREP_SRC -> $OUT (symlink blocks, repath config) ==="
  ln -sfn "$FED_PREP_SRC/A" "$OUT/A"
  ln -sfn "$FED_PREP_SRC/B" "$OUT/B"
  cp -r "$FED_PREP_SRC/config" "$CFG"
  cp "$FED_PREP_SRC"/*.txt "$FED_PREP_SRC"/*.json "$OUT"/ 2>/dev/null || true
  find "$CFG" -name '*.toml' -exec sed -i.bak "s|$FED_PREP_SRC|$OUT|g" {} + && rm -f "$CFG"/*.bak
fi

if [ -n "$SKIP_PREP" ] && [ -n "$FED_ANCESTRY_GROUP" ]; then
  [ -f "$OUT/manifest.json" ] || { echo "error: missing $OUT/manifest.json" >&2; exit 1; }
  PREP_ANCESTRY=$(python3 -c 'import json, sys; print(json.load(open(sys.argv[1])).get("ancestry_group", ""))' "$OUT/manifest.json")
  REQUESTED_ANCESTRY=$(printf '%s' "$FED_ANCESTRY_GROUP" | tr '[:upper:]' '[:lower:]')
  [ "$PREP_ANCESTRY" = "$REQUESTED_ANCESTRY" ] || {
    echo "error: requested ancestry $FED_ANCESTRY_GROUP, but prep contains ${PREP_ANCESTRY:-unknown}" >&2
    exit 1
  }
fi

echo "=== run dir: $OUT ==="
echo "=== run knobs (env; blank = fed_prep default) ==="
printf '  PLINK2=%s\n  FED_SPLIT_ANCESTRY=%s FED_ANCESTRY_GROUP=%s FED_CHR=%s FED_NSUB=%s FED_NGENES=%s FED_NPCS=%s\n  FED_CKKS=%s FED_DATABITS=%s FED_FRACBITS=%s FED_PHENO_COL=%s\n   FED_PROBES=%s\n  SKIP_PREP=%s SKIP_BUILD=%s FED_CSV=%s FED_PREP_SRC=%s\n' \
  "${PLINK2:-plink2}" "${FED_SPLIT_ANCESTRY:-0}" "${FED_ANCESTRY_GROUP:-}" "${FED_CHR:-}" "${FED_NSUB:-}" "${FED_NGENES:-}" "${FED_NPCS:-}" \
  "${FED_CKKS:-}" "${FED_DATABITS:-}" "${FED_FRACBITS:-}" "${FED_PHENO_COL:-}" \
  "${FED_PROBES:-}" \
  "${SKIP_PREP:-}" "${SKIP_BUILD:-}" "${FED_CSV:-}" "${FED_PREP_SRC:-}"

if [ -z "$SKIP_PREP" ]; then
  echo "=== [1/4] fed_prep (blocks + cov + pheno + config) ==="
  start_step prep
  python3 "$PREP/fed_prep.py" 2>&1 | tee "$OUT/prep.log"
  finish_step prep
  echo "  [1/4] fed_prep done: $(duration_s "$T_PREP_MS")"
else
  S_PREP=skipped
  echo "=== [1/4] fed_prep skipped (SKIP_PREP is set) ==="
fi

if [ -z "$SKIP_BUILD" ]; then
  echo "=== [2/4] build sfgwas ==="
  start_step build
  ( cd "$REPO" && { go build -mod=vendor -o sfgwas || go build -mod=mod -o sfgwas; } )
  finish_step build
  echo "  [2/4] build done: $(duration_s "$T_BUILD_MS")"
else
  S_BUILD=skipped
  echo "=== [2/4] build skipped (SKIP_BUILD is set) ==="
fi

if [ -f "$CFG/configGlobal.toml" ]; then
  echo "=== resolved params (n=num_inds, m=num_snps, c=num_covs; from $CFG) ==="
  grep -E '^(num_inds|num_snps|num_covs|geno_num_blocks|ckks_params|mpc_data_bits|mpc_frac_bits|rotkey_pow2only|private_pid)\b' "$CFG/configGlobal.toml" | sed 's/^/  /'
fi

echo "=== [3/4] secure skat_fed (3 parties) ==="
pkill -9 -f "$REPO/sfgwas" 2>/dev/null || true
sleep 1
# A reused FED_OUT can contain mode-exclusive files from an older raw-Q or
# SKAT-p run. Clear generated results before launching so comparison/plotting
# can only consume outputs from this run.
for party_dir in "$OUT"/out/party*; do
  [ -d "$party_dir" ] || continue
  rm -f "$party_dir/skat_fed_out.txt" \
    "$party_dir/skat_fed_burden_p_out.txt" \
    "$party_dir/skat_fed_skat_p_out.txt"
done
rm -f "$OUT/communication_summary.csv" "$OUT/fed_results.csv" \
  "$OUT/manhattan_burden.png" "$OUT/scatter_burden.png" \
  "$OUT/manhattan_skat.png" "$OUT/scatter_skat.png"
start_step secure
party_ids=(0 1 2)
for i in "${party_ids[@]}"; do
  PID=$i SFGWAS_MODE=skat_fed SFGWAS_CONFIG_PATH="$CFG" "$REPO/sfgwas" > "$OUT/party$i.log" 2>&1 &
  PARTY_PIDS+=("$!")
done
# Keep the long-running gene progress/ETA visible while retaining a complete
# party-2 log. This portable follower is tracked and terminated on every exit.
tail -n +1 -f "$OUT/party2.log" &
PARTY_LOG_TAIL_PID=$!

# macOS ships Bash 3.2 without `wait -n`. Poll all three children so a failure
# from any party is noticed promptly instead of waiting forever on a peer that
# is blocked in network I/O.
party_done=(0 0 0)
remaining=${#PARTY_PIDS[@]}
while [ "$remaining" -gt 0 ]; do
  observed_exit=0
  for idx in 0 1 2; do
    if [ "${party_done[$idx]}" -eq 1 ]; then
      continue
    fi
    child_pid=${PARTY_PIDS[$idx]}
    if kill -0 "$child_pid" 2>/dev/null; then
      continue
    fi
    set +e
    wait "$child_pid"
    child_status=$?
    set -e
    party_done[$idx]=1
    remaining=$((remaining - 1))
    observed_exit=1
    if [ "$child_status" -ne 0 ]; then
      cat "$OUT/party2.log"
      echo "secure party ${party_ids[$idx]} failed: status=$child_status" >&2
      exit "$child_status"
    fi
  done
  if [ "$remaining" -gt 0 ] && [ "$observed_exit" -eq 0 ]; then
    sleep 0.1
  fi
done
PARTY_PIDS=()
# Stop at protocol completion; log rendering and communication validation are
# post-processing and intentionally excluded from the secure runtime.
finish_step secure
sleep 0.1 # allow tail to print the final buffered lines; excluded from secure time
kill "$PARTY_LOG_TAIL_PID" 2>/dev/null || true
set +e
wait "$PARTY_LOG_TAIL_PID" 2>/dev/null
set -e
PARTY_LOG_TAIL_PID=""
echo "  [3/4] secure done: $(duration_s "$T_SECURE_MS")"

python3 "$REPO/scripts/analysis/communication_summary.py" \
  "$OUT/party0.log" "$OUT/party1.log" "$OUT/party2.log" \
  --scope skat_fed_total --expected-parties 0 1 2 \
  --output "$OUT/communication_summary.csv"

echo "=== [4/4] fed_compare (secure vs plaintext) ==="
start_step compare
python3 "$PREP/fed_compare.py" 2>&1 | tee "$OUT/compare.log"
finish_step compare
echo "  [4/4] compare done: $(duration_s "$T_COMPARE_MS")"

echo "=== [plot] fed_plot (Manhattan + secure-vs-plain scatter) ==="
if [ -f "$OUT/fed_results.csv" ]; then
  python3 "$REPO/scripts/analysis/fed_plot.py" "$OUT/fed_results.csv" "$OUT" \
    && echo "  [plot] PNGs -> $OUT" || echo "  [plot] fed_plot failed (skipped)"
else
  echo "  [plot] no fed_results.csv (run with FED_CSV=1) -> skipped"
fi

write_timing_summary 0
