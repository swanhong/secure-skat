#!/usr/bin/env bash

set -euo pipefail

config_path="${1:-run.aou.conf}"
party_id="${2:-1}"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

python3 - "${config_path}" "${party_id}" "${script_dir}" <<'PY'
import csv
import sys
import tomllib
from collections import defaultdict
from pathlib import Path


config_path = Path(sys.argv[1])
party_id = sys.argv[2]
repository_root = Path(sys.argv[3])
sys.path.insert(0, str(repository_root / "rewrite" / "analysis"))

from compare_secure_to_reference import R2_COMPARISONS, read_rows, r_squared


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

COMMUNICATION_FIELDS = (
    "sent_bytes",
    "received_bytes",
)

SETUP_COMMUNICATION_STAGES = {
    "sample_count_exchange",
    "collective_setup",
    "null_model",
}


def format_duration(seconds: float) -> str:
    if seconds >= 3600:
        hours, remainder = divmod(seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        return f"{int(hours)}h{int(minutes):02d}m{seconds:06.3f}s"
    if seconds >= 60:
        minutes, seconds = divmod(seconds, 60)
        return f"{int(minutes)}m{seconds:06.3f}s"
    return f"{seconds:.3f}s"


def format_bytes(byte_count: int) -> str:
    value = float(byte_count)
    units = ("B", "KiB", "MiB", "GiB", "TiB")
    for unit in units:
        if value < 1024 or unit == units[-1]:
            if unit == "B":
                return f"{int(value)} B"
            return f"{value:.2f} {unit}"
        value /= 1024
    raise AssertionError("unreachable")


def format_r_squared(score: float | None) -> str:
    if score is None:
        return "NA"
    return f"{score:.6f}"


def print_table(
    headers: list[str],
    rows: list[list[str]],
    left_columns: set[int],
) -> None:
    widths = [
        max(len(header), *(len(row[index]) for row in rows))
        for index, header in enumerate(headers)
    ]

    def print_row(row: list[str]) -> None:
        values = []
        for index, (value, width) in enumerate(zip(row, widths)):
            if index in left_columns:
                values.append(f"{value:<{width}}")
            else:
                values.append(f"{value:>{width}}")
        print("  ".join(values))

    print_row(headers)
    print("  ".join("─" * width for width in widths))
    for row in rows:
        print_row(row)


def add_communication(
    destination: dict[str, int],
    source: dict[str, int],
) -> None:
    for field in COMMUNICATION_FIELDS:
        destination[field] += source[field]


def communication_row(
    ancestry: str,
    scope: str,
    values: dict[str, int],
) -> list[str]:
    sent_bytes = values["sent_bytes"]
    received_bytes = values["received_bytes"]
    return [
        ancestry,
        scope,
        format_bytes(sent_bytes),
        format_bytes(received_bytes),
        format_bytes(sent_bytes + received_bytes),
    ]


def eligible_rows(
    rows: list[dict[str, str]],
    convergence_column: str | None,
) -> list[dict[str, str]]:
    if convergence_column is None:
        return rows
    return [row for row in rows if row[convergence_column] == "1"]


def r_squared_summary(
    run_dir: Path,
    ancestries: list[str],
) -> tuple[list[list[str]], list[str]]:
    summary_rows = []
    missing_ancestries = []

    for ancestry in ancestries:
        comparison_path = (
            run_dir / "comparison" / ancestry / "all_comparison.csv"
        )
        if not comparison_path.is_file():
            missing_ancestries.append(ancestry)
            continue

        comparison_rows = read_rows(comparison_path)
        groups: dict[
            tuple[int, str, int],
            list[dict[str, str]],
        ] = defaultdict(list)
        for row in comparison_rows:
            key = (
                int(row["phenotype_index"]),
                row["phenotype_name"],
                int(row["chromosome"]),
            )
            groups[key].append(row)

        for (
            label,
            secure_column,
            reference_column,
            convergence_column,
        ) in R2_COMPARISONS:
            pooled_rows = eligible_rows(
                comparison_rows,
                convergence_column,
            )
            pooled_count, pooled_score = r_squared(
                pooled_rows,
                secure_column,
                reference_column,
            )

            worst = None
            for group_key, group_rows in sorted(groups.items()):
                count, score = r_squared(
                    eligible_rows(group_rows, convergence_column),
                    secure_column,
                    reference_column,
                )
                if score is None:
                    continue
                candidate = (score, *group_key, count)
                if worst is None or candidate < worst:
                    worst = candidate

            if worst is None:
                worst_score = None
                worst_phenotype = "NA"
                worst_chromosome = "NA"
                worst_count = "0"
            else:
                (
                    worst_score,
                    phenotype_index,
                    phenotype_name,
                    chromosome,
                    count,
                ) = worst
                worst_phenotype = f"{phenotype_index}:{phenotype_name}"
                worst_chromosome = str(chromosome)
                worst_count = str(count)

            summary_rows.append([
                ancestry,
                label,
                str(pooled_count),
                format_r_squared(pooled_score),
                format_r_squared(worst_score),
                worst_phenotype,
                worst_chromosome,
                worst_count,
            ])

    return summary_rows, missing_ancestries


try:
    with config_path.open("rb") as config_file:
        config = tomllib.load(config_file)
except OSError as error:
    raise SystemExit(f"cannot read config {config_path}: {error}") from error
except tomllib.TOMLDecodeError as error:
    raise SystemExit(f"cannot parse config {config_path}: {error}") from error

configured_run_dir = str(config.get("run_dir", "")).strip()
if not configured_run_dir:
    raise SystemExit(f"run_dir is missing from {config_path}")

run_dir = Path(configured_run_dir)
metrics_dir = run_dir / "metrics"
metric_files = sorted(metrics_dir.rglob(f"metrics_party{party_id}.csv"))
if not metric_files:
    raise SystemExit(
        f"no metrics_party{party_id}.csv files found under {metrics_dir}"
    )

durations: dict[tuple[str, int], dict[str, float]] = defaultdict(
    lambda: defaultdict(float)
)
communications: dict[tuple[str, int], dict[str, int]] = defaultdict(
    lambda: defaultdict(int)
)

for metric_path in metric_files:
    with metric_path.open(newline="", encoding="utf-8") as input_file:
        reader = csv.DictReader(input_file)
        required = {
            "stage",
            "chromosome",
            "duration_seconds",
            *COMMUNICATION_FIELDS,
        }
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise SystemExit(
                f"{metric_path} is missing columns: "
                f"{', '.join(sorted(missing))}"
            )

        for row in reader:
            ancestry = row.get("ancestry") or metric_path.parent.name
            stage = row["stage"]
            chromosome_text = row["chromosome"]
            chromosome = int(chromosome_text) if chromosome_text else 0

            if chromosome:
                durations[(ancestry, chromosome)][stage] += float(
                    row["duration_seconds"]
                )

            include_communication = (
                chromosome == 0 and stage in SETUP_COMMUNICATION_STAGES
            ) or (
                chromosome != 0 and stage == "chromosome_total"
            )
            if not include_communication:
                continue

            values = communications[(ancestry, chromosome)]
            for field in COMMUNICATION_FIELDS:
                values[field] += int(row[field] or 0)

chromosomes = sorted(
    key
    for key, values in durations.items()
    if "chromosome_total" in values
)
if not chromosomes:
    raise SystemExit(f"no chromosome_total rows found under {metrics_dir}")

timing_rows = []
ancestry_timing_totals: dict[str, dict[str, float]] = defaultdict(
    lambda: defaultdict(float)
)
all_timing_totals: dict[str, float] = defaultdict(float)

for ancestry, chromosome in chromosomes:
    stage_values = durations[(ancestry, chromosome)]
    timing_rows.append([
        ancestry,
        str(chromosome),
        *(format_duration(stage_values[stage]) for _, stage in STAGES),
    ])
    for _, stage in STAGES:
        ancestry_timing_totals[ancestry][stage] += stage_values[stage]
        all_timing_totals[stage] += stage_values[stage]

for ancestry in sorted(ancestry_timing_totals):
    timing_rows.append([
        ancestry,
        "TOTAL",
        *(
            format_duration(ancestry_timing_totals[ancestry][stage])
            for _, stage in STAGES
        ),
    ])

if len(ancestry_timing_totals) > 1:
    timing_rows.append([
        "ALL",
        "TOTAL",
        *(format_duration(all_timing_totals[stage]) for _, stage in STAGES),
    ])

communication_rows = []
communication_ancestries = sorted(
    {ancestry for ancestry, _ in communications}
)
all_communication: dict[str, int] = defaultdict(int)

for ancestry in communication_ancestries:
    ancestry_total: dict[str, int] = defaultdict(int)
    scopes = sorted(
        scope
        for current_ancestry, scope in communications
        if current_ancestry == ancestry
    )
    for scope in scopes:
        values = communications[(ancestry, scope)]
        scope_label = "Setup" if scope == 0 else f"chr{scope}"
        communication_rows.append(
            communication_row(ancestry, scope_label, values)
        )
        add_communication(ancestry_total, values)
        add_communication(all_communication, values)

    communication_rows.append(
        communication_row(ancestry, "TOTAL", ancestry_total)
    )

if len(communication_ancestries) > 1:
    communication_rows.append(
        communication_row("ALL", "TOTAL", all_communication)
    )

configured_ancestries = [
    str(ancestry).strip().upper()
    for ancestry in config.get("ancestries", [])
]
if not configured_ancestries:
    configured_ancestries = sorted(ancestry_timing_totals)

r2_rows, missing_r2_ancestries = r_squared_summary(
    run_dir,
    configured_ancestries,
)

print("Secure RVAS metrics summary")
print(f"Config: {config_path}")
print(f"Run directory: {run_dir}")

print(f"\nTiming summary (party{party_id})")
print_table(
    ["Ancestry", "Chr", *(label for label, _ in STAGES)],
    timing_rows,
    {0},
)
print(
    "Stage timings are measured independently and can overlap during "
    "parallel execution."
)

print(f"\nCommunication summary (party{party_id})")
print_table(
    [
        "Ancestry",
        "Scope",
        "Sent",
        "Received",
        "Total I/O",
    ],
    communication_rows,
    {0, 1},
)
print(
    "Totals use non-overlapping setup scopes and chromosome_total rows; "
    "Total I/O is sent + received for this party."
)

print("\nR^2 summary on -log10(p)")
if r2_rows:
    print_table(
        [
            "Ancestry",
            "Comparison",
            "Pooled n",
            "Pooled R^2",
            "Worst R^2",
            "Worst phenotype",
            "Chr",
            "n",
        ],
        r2_rows,
        {0, 1, 5},
    )
else:
    print("Unavailable: no comparison CSV files were found.")

if missing_r2_ancestries:
    print(
        "Unavailable ancestries: "
        + ", ".join(missing_r2_ancestries)
        + "; run compare_secure_to_reference.py first."
    )
print(
    "SKAT WH vs Davies includes approximation/reference differences and "
    "uses only converged Davies rows."
)
PY
