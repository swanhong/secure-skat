#!/usr/bin/env bash

set -euo pipefail

input_path="${1:-output/secure-rvas-aou}"
party_id="${2:-1}"

if [[ -d "${input_path}/metrics" ]]; then
    metrics_dir="${input_path}/metrics"
else
    metrics_dir="${input_path}"
fi

python3 - "${metrics_dir}" "${party_id}" <<'PY'
import csv
import sys
from collections import defaultdict
from pathlib import Path


STAGES = (
    ("Total", "chromosome_total"),
    ("Weights", "compute_weights"),
    ("Packed", "packed_statistics"),
    ("GtG1", "first_gtg_action"),
    ("GtG2", "second_gtg_action"),
    ("Trace", "private_trace_correction"),
    ("Finalize", "finalize"),
    ("Release", "release"),
)


def format_duration(seconds: float) -> str:
    if seconds >= 3600:
        hours, remainder = divmod(seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        return f"{int(hours)}h{int(minutes):02d}m{seconds:06.3f}s"
    if seconds >= 60:
        minutes, seconds = divmod(seconds, 60)
        return f"{int(minutes)}m{seconds:06.3f}s"
    return f"{seconds:.3f}s"


def print_table(headers: list[str], rows: list[list[str]]) -> None:
    widths = [
        max(len(header), *(len(row[index]) for row in rows))
        for index, header in enumerate(headers)
    ]
    formats = [f"{{:<{widths[0]}}}", f"{{:>{widths[1]}}}"] + [
        f"{{:>{width}}}" for width in widths[2:]
    ]

    def print_row(row: list[str]) -> None:
        print("  ".join(template.format(value) for template, value in zip(formats, row)))

    print_row(headers)
    print("  ".join("-" * width for width in widths))
    for row in rows:
        print_row(row)


metrics_dir = Path(sys.argv[1])
party_id = sys.argv[2]
metric_files = sorted(metrics_dir.rglob(f"metrics_party{party_id}.csv"))
if not metric_files:
    raise SystemExit(
        f"no metrics_party{party_id}.csv files found under {metrics_dir}"
    )

durations: dict[tuple[str, int], dict[str, float]] = defaultdict(
    lambda: defaultdict(float)
)
for metric_path in metric_files:
    with metric_path.open(newline="") as input_file:
        reader = csv.DictReader(input_file)
        required = {"stage", "chromosome", "duration_seconds"}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise SystemExit(
                f"{metric_path} is missing columns: {', '.join(sorted(missing))}"
            )

        for row in reader:
            if not row["chromosome"]:
                continue
            ancestry = row.get("ancestry") or metric_path.parent.name
            key = (ancestry, int(row["chromosome"]))
            durations[key][row["stage"]] += float(row["duration_seconds"])

chromosomes = sorted(
    key for key, values in durations.items()
    if "chromosome_total" in values
)
if not chromosomes:
    raise SystemExit(f"no chromosome_total rows found under {metrics_dir}")

rows = []
ancestry_totals: dict[str, dict[str, float]] = defaultdict(
    lambda: defaultdict(float)
)
all_totals: dict[str, float] = defaultdict(float)

for ancestry, chromosome in chromosomes:
    stage_values = durations[(ancestry, chromosome)]
    rows.append([
        ancestry,
        str(chromosome),
        *(format_duration(stage_values[stage]) for _, stage in STAGES),
    ])
    for _, stage in STAGES:
        ancestry_totals[ancestry][stage] += stage_values[stage]
        all_totals[stage] += stage_values[stage]

for ancestry in sorted(ancestry_totals):
    rows.append([
        ancestry,
        "TOTAL",
        *(
            format_duration(ancestry_totals[ancestry][stage])
            for _, stage in STAGES
        ),
    ])

if len(ancestry_totals) > 1:
    rows.append([
        "ALL",
        "TOTAL",
        *(format_duration(all_totals[stage]) for _, stage in STAGES),
    ])

print(f"Metrics timing summary (party{party_id})")
print(f"Source: {metrics_dir}")
print()
print_table(
    ["Ancestry", "Chr", *(label for label, _ in STAGES)],
    rows,
)
print()
print("Stage timings are measured independently and can overlap during parallel execution.")
PY
