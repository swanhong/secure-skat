#!/bin/bash
# End-to-end federated-private SKAT (#4 benchmark): prep -> build -> secure 3-party -> compare.
#
#   PLINK2=$HOME/plink2 bash run_fed.sh                                  # full run (default N_SUB etc.)
#   PLINK2=$HOME/plink2 FED_NSUB=38000 FED_DATABITS=100 bash run_fed.sh  # scale n (needs more fixed-point range)
#   SKIP_PREP=1 SKIP_BUILD=1 bash run_fed.sh                             # reuse existing blocks+config, rerun secure
#
# FED_* / PLINK2 env vars are BAKED into ~/fed_prep_out/config by fed_prep (data_bits, ports, ckks,
# dims), so FED_DATABITS etc. only take effect when fed_prep runs; with SKIP_PREP=1 the existing
# config is reused as-is. SKIP_BUILD=1 skips go build.
set -eo pipefail

REPO=$(cd "$(dirname "$0")" && pwd)
OUT=${FED_OUT:-$HOME/fed_prep_out}
CFG=$OUT/config
PREP=$REPO/scripts/preprocessing
T_PREP=0 T_BUILD=0 T_SECURE=0 T_COMPARE=0

if [ -z "$SKIP_PREP" ]; then
  echo "=== [1/4] fed_prep (blocks + cov + pheno + config) ==="
  s=$SECONDS; python3 "$PREP/fed_prep.py"; T_PREP=$((SECONDS-s))
fi

if [ -z "$SKIP_BUILD" ]; then
  echo "=== [2/4] build sfgwas ==="
  s=$SECONDS; ( cd "$REPO" && { go build -mod=vendor -o sfgwas || go build -mod=mod -o sfgwas; } ); T_BUILD=$((SECONDS-s))
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

echo "=== [4/4] fed_compare (secure vs plaintext) ==="
s=$SECONDS; python3 "$PREP/fed_compare.py"; T_COMPARE=$((SECONDS-s))

echo "=== TIMING (steps) ==="
printf '  [1] prep    %5ds\n  [2] build   %5ds\n  [3] secure  %5ds\n  [4] compare %5ds\n  ---------------------\n  total     %7ds\n' \
  "$T_PREP" "$T_BUILD" "$T_SECURE" "$T_COMPARE" "$((T_PREP+T_BUILD+T_SECURE+T_COMPARE))"
