#!/usr/bin/env python3

import argparse
import csv
import math
import tomllib
from pathlib import Path


KEY_COLUMNS = [
    "chromosome",
    "gene_index",
    "gene_id",
    "phenotype_index",
    "phenotype_name",
]

OUTPUT_COLUMNS = [
    *KEY_COLUMNS,
    "secure_burden_p",
    "r_burden_p",
    "burden_abs_error",
    "secure_skat_wh_p",
    "r_skat_liu_p",
    "skat_liu_abs_error",
    "r_skat_davies_p",
    "r_skat_davies_converged",
    "skat_davies_abs_error",
]

R2_COMPARISONS = [
    (
        "Burden",
        "secure_burden_p",
        "r_burden_p",
        None,
    ),
    (
        "SKAT WH vs Liu",
        "secure_skat_wh_p",
        "r_skat_liu_p",
        None,
    ),
    (
        "SKAT WH vs Davies",
        "secure_skat_wh_p",
        "r_skat_davies_p",
        "r_skat_davies_converged",
    ),
]

def read_rows(
    path: Path,
    delimiter: str = ",",
) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as input_file:
        return list(csv.DictReader(input_file, delimiter=delimiter))


def row_key(row: dict[str, str]) -> tuple[str, ...]:
    return tuple(row[column] for column in KEY_COLUMNS)


def index_rows(
    rows: list[dict[str, str]],
    label: str,
) -> dict[tuple[str, ...], dict[str, str]]:
    indexed: dict[tuple[str, ...], dict[str, str]] = {}

    for row in rows:
        key = row_key(row)
        if key in indexed:
            raise ValueError(f"duplicate {label} result key: {key}")
        indexed[key] = row

    return indexed


def absolute_error(left: str, right: str) -> str:
    try:
        left_value = float(left)
        right_value = float(right)
    except ValueError:
        return ""

    if not math.isfinite(left_value) or not math.isfinite(right_value):
        return ""

    return format(abs(left_value - right_value), ".17g")

def r_squared(
        rows: list[dict[str, str]],
        secure_column: str,
        reference_column: str,
) -> tuple[int, float | None]:
    secure_values = []
    reference_values = []

    for row in rows:
        try:
            secure_value = float(row[secure_column])
            reference_value = float(row[reference_column])
        except ValueError:
            continue

        if math.isfinite(secure_value) and math.isfinite(reference_value):
            secure_values.append(secure_value)
            reference_values.append(reference_value)

    if len(reference_values) < 2:
        return len(reference_values), None

    reference_mean = sum(reference_values) / len(reference_values)
    total_sum_squares = sum(
        (value - reference_mean) ** 2 for value in reference_values
    )
    if total_sum_squares == 0:
        return len(reference_values), None

    residual_sum_squares = sum(
        (secure - reference) ** 2
        for secure, reference in zip(secure_values, reference_values)
    )
    return (len(reference_values), 1 - residual_sum_squares / total_sum_squares)

def print_r_squared_by_pheno_and_chr(
        rows: list[dict[str, str]],
) -> None:
    grouped_rows = {}

    for row in rows:
        key = (
            int(row["phenotype_index"]),
            row["phenotype_name"],
            int(row["chromosome"]),
        )
        grouped_rows.setdefault(key, []).append(row)

    print("\n=== R^2 by phenotype and chromosome ===")
    print(
        f"  {'phenotype':<20} {'chr':>3} "
        f"{'comparison':<20} {'n':>5} {'R^2':>12}"
    )

    for (_, phenotype_name, chromosome), group in sorted(
        grouped_rows.items()
    ):
        for (
            label,
            secure_column,
            reference_column,
            convergence_column,
        ) in R2_COMPARISONS:
            eligible_rows = group
            if convergence_column is not None:
                eligible_rows = [
                    row
                    for row in group
                    if row[convergence_column] == "1"
                ]

            count, score = r_squared(
                eligible_rows,
                secure_column,
                reference_column,
            )
            score_text = "NA" if score is None else f"{score:.6f}"

            print(
                f"  {phenotype_name:<20} {chromosome:>3} "
                f"{label:<20} {count:>5} {score_text:>12}"
            )

def compare_results(config_path: Path) -> None:
    with config_path.open("rb") as config_file:
        config = tomllib.load(config_file)

    run_dir = Path(config["run_dir"])
    comparison_dir = run_dir / "comparison"
    comparison_dir.mkdir(parents=True, exist_ok=True)

    success_path = comparison_dir / "_SUCCESS"
    success_path.unlink(missing_ok=True)

    secure_rows = read_rows(
        run_dir / "secure" / "all_secure_results.tsv",
        delimiter="\t",
    )
    reference_rows = read_rows(
        run_dir / "reference" / "all_r_results.csv"
    )

    secure_by_key = index_rows(secure_rows, "secure")
    reference_by_key = index_rows(reference_rows, "R")

    if set(secure_by_key) != set(reference_by_key):
        raise ValueError("secure and R result keys do not match")

    comparison_rows: list[dict[str, str]] = []

    for secure in secure_rows:
        reference = reference_by_key[row_key(secure)]

        davies_error = ""
        if reference["r_skat_davies_converged"] == "1":
            davies_error = absolute_error(
                secure["secure_skat_wh_p"],
                reference["r_skat_davies_p"],
            )

        comparison_rows.append(
            {
                **{
                    column: secure[column]
                    for column in KEY_COLUMNS
                },
                "secure_burden_p": secure["secure_burden_p"],
                "r_burden_p": reference["r_burden_p"],
                "burden_abs_error": absolute_error(
                    secure["secure_burden_p"],
                    reference["r_burden_p"],
                ),
                "secure_skat_wh_p": secure["secure_skat_wh_p"],
                "r_skat_liu_p": reference["r_skat_liu_p"],
                "skat_liu_abs_error": absolute_error(
                    secure["secure_skat_wh_p"],
                    reference["r_skat_liu_p"],
                ),
                "r_skat_davies_p":
                    reference["r_skat_davies_p"],
                "r_skat_davies_converged":
                    reference["r_skat_davies_converged"],
                "skat_davies_abs_error": davies_error,
            }
        )

    output_path = comparison_dir / "all_comparison.csv"
    with output_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as output_file:
        writer = csv.DictWriter(
            output_file,
            fieldnames=OUTPUT_COLUMNS,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(comparison_rows)

    success_path.touch()
    print("Wrote secure/plain comparison to", output_path)
    print_r_squared_by_pheno_and_chr(comparison_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        default="run.1kg.conf",
        help="path to the run configuration",
    )
    args = parser.parse_args()

    compare_results(Path(args.config))


if __name__ == "__main__":
    main()