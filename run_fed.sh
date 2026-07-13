#!/bin/bash
# End-to-end federated-private SKAT (#4 benchmark): prep -> build -> secure 3-party -> compare.
#
#   PLINK2=$HOME/plink2 bash run_fed.sh                                  # full run → $HOME/runs/out<YYMMDDHHMMSS>
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
# FED_PROBES=N (N>0) turns ON the SKAT p-value (Hutchinson probes → WH pivot z; 0 = Q+burden p only);
set -eo pipefail

REPO=$(cd "$(dirname "$0")" && pwd)
# Each run lives in its own timestamped dir so past results are kept: FED_OUT unset → $HOME/runs/out<YYMMDDHHMMSS>
# holds fed_in (prep blocks/config) + fed_out (secure results, compare csv). Set FED_OUT to reuse/target a dir.
# FED_PREP_SRC=<prior run dir>: with SKIP_PREP, seed this run's fed_in from there (symlink blocks, copy+repath
# config) so a slow prep is reused while results still land in the fresh archived dir.
: "${FED_OUT:=$HOME/runs/out$(date +%y%m%d%H%M%S)}"
export FED_OUT # fed_prep.py reads it too
OUT=$FED_OUT
CFG=$OUT/config
PREP=$REPO/scripts/preprocessing
T_PREP=0 T_BUILD=0 T_SECURE=0 T_COMPARE=0
mkdir -p "$OUT"

if [ -n "$SKIP_PREP" ] && [ -n "$FED_PREP_SRC" ] && [ ! -f "$CFG/configGlobal.toml" ]; then
  echo "=== seeding fed_in from $FED_PREP_SRC -> $OUT (symlink blocks, repath config) ==="
  ln -sfn "$FED_PREP_SRC/A" "$OUT/A"
  ln -sfn "$FED_PREP_SRC/B" "$OUT/B"
  cp -r "$FED_PREP_SRC/config" "$CFG"
  cp "$FED_PREP_SRC"/*.txt "$FED_PREP_SRC"/*.json "$OUT"/ 2>/dev/null || true
  find "$CFG" -name '*.toml' -exec sed -i.bak "s|$FED_PREP_SRC|$OUT|g" {} + && rm -f "$CFG"/*.bak
fi

echo "=== run dir: $OUT ==="
echo "=== run knobs (env; blank = fed_prep default) ==="
printf '  PLINK2=%s\n  FED_CHR=%s FED_NSUB=%s FED_NGENES=%s FED_NPCS=%s\n  FED_CKKS=%s FED_DATABITS=%s FED_FRACBITS=%s FED_PHENO_COL=%s\n   FED_PROBES=%s\n  SKIP_PREP=%s SKIP_BUILD=%s FED_CSV=%s FED_PREP_SRC=%s\n' \
  "${PLINK2:-plink2}" "${FED_CHR:-}" "${FED_NSUB:-}" "${FED_NGENES:-}" "${FED_NPCS:-}" \
  "${FED_CKKS:-}" "${FED_DATABITS:-}" "${FED_FRACBITS:-}" "${FED_PHENO_COL:-}" \
  "${FED_PROBES:-}" \
  "${SKIP_PREP:-}" "${SKIP_BUILD:-}" "${FED_CSV:-}" "${FED_PREP_SRC:-}"

if [ -z "$SKIP_PREP" ]; then
  echo "=== [1/4] fed_prep (blocks + cov + pheno + config) ==="
  s=$SECONDS; python3 "$PREP/fed_prep.py"; T_PREP=$((SECONDS-s))
  echo "  [1/4] fed_prep done: ${T_PREP}s"
fi

if [ -z "$SKIP_BUILD" ]; then
  echo "=== [2/4] build sfgwas ==="
  s=$SECONDS; ( cd "$REPO" && { go build -mod=vendor -o sfgwas || go build -mod=mod -o sfgwas; } ); T_BUILD=$((SECONDS-s))
  echo "  [2/4] build done: ${T_BUILD}s"
fi

if [ -f "$CFG/configGlobal.toml" ]; then
  echo "=== resolved params (n=num_inds, m=num_snps, c=num_covs; from $CFG) ==="
  grep -E '^(num_inds|num_snps|num_covs|geno_num_blocks|ckks_params|mpc_data_bits|mpc_frac_bits|rotkey_pow2only|private_pid)\b' "$CFG/configGlobal.toml" | sed 's/^/  /'
fi

echo "=== [3/4] secure skat_fed (3 parties) ==="
s=$SECONDS
pkill -9 -f "$REPO/sfgwas" 2>/dev/null || true
sleep 1
for i in 0 1; do
  PID=$i SFGWAS_MODE=skat_fed SFGWAS_CONFIG_PATH="$CFG" "$REPO/sfgwas" > "$OUT/party$i.log" 2>&1 &
done
PID=2 SFGWAS_MODE=skat_fed SFGWAS_CONFIG_PATH="$CFG" "$REPO/sfgwas" 2>&1 | tee "$OUT/party2.log"
wait   # let party0/1 finish writing their logs
T_SECURE=$((SECONDS-s))
echo "  [3/4] secure done: ${T_SECURE}s"

echo "=== [4/4] fed_compare (secure vs plaintext) ==="
s=$SECONDS; python3 "$PREP/fed_compare.py"; T_COMPARE=$((SECONDS-s))
echo "  [4/4] compare done: ${T_COMPARE}s"

echo "=== TIMING (steps) ==="
printf '  [1] prep    %5ds\n  [2] build   %5ds\n  [3] secure  %5ds\n  [4] compare %5ds\n  ---------------------\n  total     %7ds\n' \
  "$T_PREP" "$T_BUILD" "$T_SECURE" "$T_COMPARE" "$((T_PREP+T_BUILD+T_SECURE+T_COMPARE))"
