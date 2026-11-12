#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
source scripts/setup/env.sh

export CONFIG_A="${CONFIG_A:-config/1kg-A}"
export CONFIG_B="${CONFIG_B:-config/1kg-B}"
export CONFIG_0="${CONFIG_0:-config/1kg-cp0}"
export LOG_DIR="${LOG_DIR:-output/1kg-logs}"

case "${JOINT_STAT_ONLY:-0}" in
  0) bash scripts/run_1kg_workflow.sh ;;
  1) echo "Using existing prepared inputs and keys for joint-stat" ;;
  *) echo "JOINT_STAT_ONLY must be 0 or 1" >&2; exit 1 ;;
esac

# Run real binary processes so a failed party can stop its peers reliably.
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/1kg-joint-stat.XXXXXX")
pids=""
cleanup() {
  for pid in $pids; do kill "$pid" 2>/dev/null || true; done
  for pid in $pids; do wait "$pid" 2>/dev/null || true; done
  rm -f "$build_dir/secure-rvas"
  rmdir "$build_dir"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

go build -mod=vendor -o "$build_dir/secure-rvas" secure-rvas.go
stat_logs="$LOG_DIR/joint-stat"
mkdir -p "$stat_logs"
echo "Run figure statistics after the main workflow; logs: $stat_logs"
configs=("$CONFIG_0" "$CONFIG_A" "$CONFIG_B")
for party in 0 1 2; do
  "$build_dir/secure-rvas" joint-stat --config "${configs[$party]}" --party "$party" "$@" \
    > "$stat_logs/party$party.log" 2>&1 &
  pids="$pids $!"
done

# Poll all children: waiting for party 0 first can hang if another party fails.
while [[ -n "$pids" ]]; do
  running=""
  for pid in $pids; do
    if kill -0 "$pid" 2>/dev/null; then
      running="$running $pid"
    elif ! wait "$pid"; then
      echo "joint-stat failed; inspect $stat_logs/party{0,1,2}.log" >&2
      exit 1
    fi
  done
  pids="$running"
  if [[ -n "$pids" ]]; then sleep 1; fi
done

python3 - "$CONFIG_A" "$CONFIG_B" <<'PY'
from pathlib import Path
import sys
import tomllib

for directory in sys.argv[1:]:
    with (Path(directory) / "configPrepare.toml").open("rb") as stream:
        run_dir = Path(tomllib.load(stream)["run_dir"])
    print(f"Figure table: {run_dir / 'resource_stats/resource_stats.tsv'}")
PY
echo "1KG analysis and joint-stat workflow completed"
