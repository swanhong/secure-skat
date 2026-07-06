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

if [ -z "$SKIP_PREP" ]; then
  echo "=== [1/4] fed_prep (blocks + cov + pheno + config) ==="
  python3 "$PREP/fed_prep.py"
fi

if [ -z "$SKIP_BUILD" ]; then
  echo "=== [2/4] build sfgwas ==="
  ( cd "$REPO" && { go build -mod=vendor -o sfgwas || go build -mod=mod -o sfgwas; } )
fi

echo "=== [3/4] secure skat_fed (3 parties) ==="
pkill -9 -f "$REPO/sfgwas" 2>/dev/null || true
sleep 1
for i in 0 1; do
  PID=$i SFGWAS_MODE=skat_fed SFGWAS_CONFIG_PATH="$CFG" "$REPO/sfgwas" > "$OUT/party$i.log" 2>&1 &
done
PID=2 SFGWAS_MODE=skat_fed SFGWAS_CONFIG_PATH="$CFG" "$REPO/sfgwas" 2>&1 | tee "$OUT/party2.log"
wait   # let party0/1 finish writing their logs

echo "=== [4/4] fed_compare (secure vs plaintext) ==="
python3 "$PREP/fed_compare.py"
