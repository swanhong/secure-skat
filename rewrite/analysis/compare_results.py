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