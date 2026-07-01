#!/bin/bash
# End-to-end federated-private SKAT (#4 benchmark): prep -> build -> secure 3-party -> compare.
#
#   PLINK2=$HOME/plink2 bash run_fed.sh                 # full run (default N_SUB etc.)
#   PLINK2=$HOME/plink2 FED_NSUB=38000 bash run_fed.sh  # scale n
#   FED_DATABITS=100 SKIP_PREP=1 bash run_fed.sh        # reuse blocks, rerun secure w/ more fixed-point range
#
# All FED_* / PLINK2 env vars pass through to fed_prep (data + config are generated to match).
# SKIP_PREP=1 skips fed_prep (reuse existing ~/fed_prep_out); SKIP_BUILD=1 skips go build.
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
