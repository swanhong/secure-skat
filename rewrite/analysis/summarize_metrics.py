#!/usr/bin/env python3

import argparse
import tomllib
from pathlib import Path

import pandas as pd

from compare_secure_to_reference import R2_COMPARISONS, comparison_r_squared


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

REQUIRED_METRIC_COLUMNS = {
    "stage",
    "chromosome",
    "duration_seconds",
    *COMMUNICATION_FIELDS,
}

def print_table(
    table: pd.DataFrame,
    left_columns: set[str],
) -> None:
    headers = [str(column) for column in table.columns]
    rows = [
        ["NA" if pd.isna(value) else str(value) for value in row]
        for row in table.itertuples(index=False, name=None)
    ]
    widths = [
        max([len(header), *(len(row[index]) for row in rows)])
        for index, header in enumerate(headers)
    ]

    def print_row(row: list[str]) -> None:
        values = []
        for header, value, width in zip(headers, row, widths):
            if header in left_columns:
                values.append(f"{value:<{width}}")
            else:
                values.append(f"{value:>{width}}")
        print("  ".join(values))

    print_row(headers)
    print("  ".join("─" * width for width in widths))
    for row in rows:
        print_row(row)

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

def read_metrics(metric_paths: list[Path]) -> pd.DataFrame:
    frames = []

    for path in metric_paths:
        frame = pd.read_csv(path)

        missing_columns = REQUIRED_METRIC_COLUMNS - set(frame.columns)
        if missing_columns:
            missing = ", ".join(sorted(missing_columns))
            raise ValueError(f"{path} is missing columns: {missing}")

        if "ancestry" not in frame.columns:
            frame["ancestry"] = path.parent.name
        else:
            ancestry = frame["ancestry"].astype("string").str.strip()
            frame["ancestry"] = ancestry.mask(
                ancestry.isna() | ancestry.eq(""),
                path.parent.name,
            )

        frame["chromosome"] = (
            pd.to_numeric(frame["chromosome"])
            .fillna(0)
            .astype(int)
        )
        frame["duration_seconds"] = pd.to_numeric(
            frame["duration_seconds"]
        )

        for column in COMMUNICATION_FIELDS:
            frame[column] = (
                pd.to_numeric(frame[column])
                .fillna(0)
                .astype(int)
            )

        frames.append(frame)

    return pd.concat(frames, ignore_index=True)

def make_time_table(metrics: pd.DataFrame) -> pd.DataFrame:
    stage_names = [stage for _, stage in STAGES]
    columns = ["Ancestry", "Chr", *(label for label, _ in STAGES)]

    completed = metrics.loc[
        (metrics["chromosome"] != 0)
        & (metrics["stage"] == "chromosome_total"),
        ["ancestry", "chromosome"],
    ].drop_duplicates()
    if completed.empty:
        raise ValueError("no chromosome_total rows found")

    completed = completed.sort_values(["ancestry", "chromosome"])
    completed_index = pd.MultiIndex.from_frame(completed)

    durations = (
        metrics.loc[metrics["chromosome"] != 0]
        .groupby(["ancestry", "chromosome", "stage"])[
            "duration_seconds"
        ]
        .sum()
        .unstack(fill_value=0)
        .reindex(completed_index)
        .reindex(columns=stage_names, fill_value=0)
    )

    def time_row(
        ancestry: str,
        chromosome: str,
        values: pd.Series,
    ) -> list[str]:
        return [
            ancestry,
            chromosome,
            *(format_duration(float(values[stage])) for stage in stage_names),
        ]

    rows = [
        time_row(ancestry, str(chromosome), values)
        for (ancestry, chromosome), values in durations.iterrows()
    ]

    ancestry_totals = durations.groupby(level="ancestry", sort=True).sum()
    rows.extend(
        time_row(ancestry, "TOTAL", values)
        for ancestry, values in ancestry_totals.iterrows()
    )

    if len(ancestry_totals) > 1:
        rows.append(time_row("ALL", "TOTAL", durations.sum()))

    return pd.DataFrame(rows, columns=columns)


def make_comm_table(metrics: pd.DataFrame) -> pd.DataFrame:
    columns = ["Ancestry", "Scope", "Sent", "Received", "Total I/O"]
    setup_rows = (
        (metrics["chromosome"] == 0)
        & metrics["stage"].isin(SETUP_COMMUNICATION_STAGES)
    )
    chromosome_rows = (
        (metrics["chromosome"] != 0)
        & (metrics["stage"] == "chromosome_total")
    )
    selected = metrics.loc[setup_rows | chromosome_rows]
    if selected.empty:
        return pd.DataFrame(columns=columns)

    communication = selected.groupby(
        ["ancestry", "chromosome"],
        sort=True,
    )[list(COMMUNICATION_FIELDS)].sum()

    def comm_row(
        ancestry: str,
        scope: str,
        values: pd.Series,
    ) -> list[str]:
        sent = int(values["sent_bytes"])
        received = int(values["received_bytes"])
        return [
            ancestry,
            scope,
            format_bytes(sent),
            format_bytes(received),
            format_bytes(sent + received),
        ]

    rows = []
    ancestries = sorted(
        communication.index.get_level_values("ancestry").unique()
    )
    for ancestry in ancestries:
        ancestry_values = communication.loc[ancestry].sort_index()
        for chromosome, values in ancestry_values.iterrows():
            scope = "Setup" if chromosome == 0 else f"chr{chromosome}"
            rows.append(comm_row(ancestry, scope, values))
        rows.append(comm_row(ancestry, "TOTAL", ancestry_values.sum()))

    if len(ancestries) > 1:
        rows.append(comm_row("ALL", "TOTAL", communication.sum()))

    return pd.DataFrame(rows, columns=columns)


def read_accuracy(
    comparison_specs: list[str],
) -> pd.DataFrame:
    frames = []

    required_columns = {
        "phenotype_index",
        "phenotype_name",
        "chromosome",
    }
    for _, secure, reference, convergence in R2_COMPARISONS:
        required_columns.update((secure, reference))
        if convergence is not None:
            required_columns.add(convergence)

    for spec in comparison_specs:
        ancestry, separator, path_text = spec.partition("=")
        ancestry = ancestry.strip()
        path_text = path_text.strip()
        if not separator or not ancestry or not path_text:
            raise ValueError(
                f"invalid comparison file {spec!r}; expected ANCESTRY=PATH"
            )

        path = Path(path_text)
        if not path.is_file():
            continue

        frame = pd.read_csv(path, dtype=str)
        missing_columns = required_columns - set(frame.columns)
        if missing_columns:
            missing = ", ".join(sorted(missing_columns))
            raise ValueError(f"{path} is missing columns: {missing}")

        frame["ancestry"] = ancestry
        frames.append(frame)

    if not frames:
        return pd.DataFrame()

    return pd.concat(frames, ignore_index=True)


def make_accuracy_table(comparisons: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "Ancestry",
        "Comparison",
        "n total",
        "n failed",
        "R^2 (all, R::SKAT)",
        "R^2 (non-failed)",
        "Worst R^2 (all)",
        "Worst phenotype",
        "Chr",
        "Worst n",
        "Worst failed",
    ]
    if comparisons.empty:
        return pd.DataFrame(columns=columns)

    rows = []
    for ancestry in comparisons["ancestry"].drop_duplicates():
        ancestry_rows = comparisons.loc[
            comparisons["ancestry"] == ancestry
        ]
        groups = [
            (
                (
                    int(group_key[0]),
                    str(group_key[1]),
                    int(group_key[2]),
                ),
                group,
            )
            for group_key, group in ancestry_rows.groupby(
                ["phenotype_index", "phenotype_name", "chromosome"],
                sort=True,
            )
        ]

        for label, secure, reference, convergence in R2_COMPARISONS:
            (
                pooled_count,
                pooled_failed,
                pooled_score,
                pooled_non_failed_score,
            ) = comparison_r_squared(
                ancestry_rows.to_dict("records"),
                secure,
                reference,
                convergence,
            )

            worst = None
            for group_key, group in groups:
                count, failed, score, _ = comparison_r_squared(
                    group.to_dict("records"),
                    secure,
                    reference,
                    convergence,
                )
                if score is None:
                    continue

                candidate = (score, *group_key, count, failed)
                if worst is None or candidate < worst:
                    worst = candidate

            if worst is None:
                worst_score = None
                worst_phenotype = "NA"
                worst_chromosome = "NA"
                worst_count = 0
                worst_failed = 0
            else:
                (
                    worst_score,
                    phenotype_index,
                    phenotype_name,
                    chromosome,
                    worst_count,
                    worst_failed,
                ) = worst
                worst_phenotype = f"{phenotype_index}:{phenotype_name}"
                worst_chromosome = str(chromosome)

            rows.append([
                ancestry,
                label,
                str(pooled_count),
                str(pooled_failed),
                format_r_squared(pooled_score),
                format_r_squared(pooled_non_failed_score),
                format_r_squared(worst_score),
                worst_phenotype,
                worst_chromosome,
                str(worst_count),
                str(worst_failed),
            ])

    return pd.DataFrame(rows, columns=columns)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--party-id", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    try:
        with args.config.open("rb") as config_file:
            config = tomllib.load(config_file)

        configured_run_dir = str(config.get("run_dir", "")).strip()
        if not configured_run_dir:
            raise ValueError(f"run_dir is missing from {args.config}")

        run_dir = Path(configured_run_dir)
        metrics_dir = run_dir / "metrics"
        metric_paths = sorted(
            metrics_dir.rglob(
                f"metrics_party{args.party_id}.csv"
            )
        )
        if not metric_paths:
            raise ValueError(
                f"no metrics_party{args.party_id}.csv files "
                f"found under {metrics_dir}"
            )

        ancestries = [
            str(ancestry).strip().upper()
            for ancestry in config.get("ancestries", [])
        ]
        comparison_specs = [
            (
                f"{ancestry}="
                f"{run_dir}/comparison/{ancestry}/all_comparison.csv"
            )
            for ancestry in ancestries
        ]

        metrics = read_metrics(metric_paths)
        time_table = make_time_table(metrics)
        comm_table = make_comm_table(metrics)

        accuracy = read_accuracy(comparison_specs)
        accuracy_table = make_accuracy_table(accuracy)
    except (OSError, ValueError, pd.errors.ParserError) as error:
        raise SystemExit(str(error)) from error

    print("Secure RVAS metrics summary")
    print(f"Config: {args.config}")
    print(f"Run directory: {run_dir}")

    print(f"\nTiming summary (party{args.party_id})")
    print_table(time_table, {"Ancestry"})

    print(f"\nCommunication summary (party{args.party_id})")
    print_table(comm_table, {"Ancestry", "Scope"})

    print("\nR^2 summary on -log10(p)")
    if accuracy_table.empty:
        print("Unavailable: no comparison CSV files were found.")
    else:
        print_table(
            accuracy_table,
            {"Ancestry", "Comparison", "Worst phenotype"},
        )

if __name__ == "__main__":
    main()
