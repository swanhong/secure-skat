#!/bin/bash
# Federated-SKAT diagnostic sweep: run the secure protocol for RAW and STANDARDIZED covariates,
# and for each print the full fed_compare diagnostics (per-gene Q, PART A/B split, secure β̂ vs
# plaintext, and plaintext-Q-built-on-secure-β̂). One command → both variants side by side.
#
#   PLINK2=$HOME/plink2 bash run_fed_diag.sh            # prep (if missing) + build + both variants
#   VARIANTS="raw std" SKIP_PREP=1 bash run_fed_diag.sh # reuse blocks, run chosen variants
#   VARIANTS=raw ...                                    # just one variant
set -eo pipefail

REPO=$(cd "$(dirname "$0")" && pwd)
OUT=${FED_OUT:-$HOME/fed_prep_out}
CFG=$OUT/config
PREP=$REPO/scripts/preprocessing
VARIANTS=${VARIANTS:-"raw std"}

# 1. prep (blocks + RAW cov + pheno + config) unless skipped; snapshot the raw cov right after
if [ -z "$SKIP_PREP" ]; then
  echo "=== fed_prep (blocks + raw cov + pheno + config) ==="
  python3 "$PREP/fed_prep.py"
  cp "$OUT/A/cov.txt" "$OUT/A/cov_raw.txt"; cp "$OUT/B/cov.txt" "$OUT/B/cov_raw.txt"  # fresh raw snapshot
fi

# 2. build
if [ -z "$SKIP_BUILD" ]; then
  echo "=== build sfgwas ==="
  ( cd "$REPO" && { go build -mod=vendor -o sfgwas || go build -mod=mod -o sfgwas; } )
fi

# 3. the variants rewrite cov.txt from the raw snapshot — it must exist
if [ ! -f "$OUT/A/cov_raw.txt" ]; then
  echo "ERROR: $OUT/A/cov_raw.txt missing (current cov.txt may be standardized). Run once WITHOUT SKIP_PREP." >&2
  exit 1
fi

set_cov() {  # $1 = raw|std ; rewrite {A,B}/cov.txt from the raw snapshot
  if [ "$1" = std ]; then
    python3 - "$OUT" <<'PY'
import sys, numpy as np
O = sys.argv[1]
A = np.loadtxt(f"{O}/A/cov_raw.txt"); B = np.loadtxt(f"{O}/B/cov_raw.txt")
allc = np.vstack([A, B]); m, s = allc.mean(0), allc.std(0); s[s == 0] = 1.0  # global (A+B) mean/std
np.savetxt(f"{O}/A/cov.txt", (A - m) / s, delimiter="\t")
np.savetxt(f"{O}/B/cov.txt", (B - m) / s, delimiter="\t")
PY
  else
    cp "$OUT/A/cov_raw.txt" "$OUT/A/cov.txt"; cp "$OUT/B/cov_raw.txt" "$OUT/B/cov.txt"
  fi
}

run_secure() {  # 3-party skat_fed; logs to party*.log
  pkill -9 -f "$REPO/sfgwas" 2>/dev/null || true; sleep 1
  for i in 0 1; do
    PID=$i SFGWAS_MODE=skat_fed SFGWAS_CONFIG_PATH="$CFG" "$REPO/sfgwas" > "$OUT/party$i.log" 2>&1 &
  done
  PID=2 SFGWAS_MODE=skat_fed SFGWAS_CONFIG_PATH="$CFG" "$REPO/sfgwas" > "$OUT/party2.log" 2>&1
  wait
}

for v in $VARIANTS; do
  echo ""
  echo "############################## VARIANT: $v covariates ##############################"
  set_cov "$v"
  s=$SECONDS; run_secure; echo "  [$v] secure done: $((SECONDS-s))s"
  python3 "$PREP/fed_compare.py" | tee "$OUT/diag_$v.txt"
done

set_cov raw  # leave cov.txt in the known (raw) state
echo ""
echo "=== SWEEP SUMMARY (max rel per variant) ==="
for v in $VARIANTS; do
  printf '  %-4s : %s\n' "$v" "$(grep 'max rel' "$OUT/diag_$v.txt" | tail -1)"
done
echo "  full logs: $OUT/diag_{$(echo $VARIANTS | tr ' ' ,)}.txt"
