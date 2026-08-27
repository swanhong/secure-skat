from __future__ import annotations

import csv
import math
from collections import defaultdict
from collections.abc import Sequence
from pathlib import Path


RESULT_FIELDS = [
    "chromosome",
    "gene_index",
    "gene_id",
    "gene_symbol",
    "gene_order",
    "phenotype_index",
    "phenotype_name",
    "secure_burden_p",
    "r_burden_p",
    "burden_abs_error",
    "burden_rel_error",
    "secure_skat_wh_p",
    "r_skat_davies_p",
    "r_skat_davies_converged",
    "skat_abs_error",
    "skat_rel_error",
]

SUMMARY_FIELDS = [
    "scope",
    "chromosome",
    "phenotype_index",
    "phenotype_name",
    "method",
    "valid_count",
    "excluded_count",
    "max_abs_error",
    "max_abs_error_gene_id",
    "max_rel_error",
    "max_rel_error_gene_id",
]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as file:
        return list(csv.DictReader(file))


def write_rows(
    path: Path,
    rows: list[dict[str, object]],
    fieldnames: list[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def result_key(row: dict[str, str]) -> tuple[str, str, int]:
    return row["chromosome"], row["gene_id"], int(row["phenotype_index"])


def index_unique(
    rows: list[dict[str, str]],
    label: str,
) -> dict[tuple[str, str, int], dict[str, str]]:
    indexed = {}
    for row in rows:
        key = result_key(row)
        if key in indexed:
            raise ValueError(f"duplicate {label} result for {key}")
        indexed[key] = row
    return indexed


def parse_converged(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "t", "yes"}


def format_number(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return "NA"
    return format(value, ".17g")


def errors(secure: float, reference: float) -> tuple[str, str]:
    absolute = abs(secure - reference)
    relative = None if reference == 0 else absolute / abs(reference)
    return format_number(absolute), format_number(relative)


def join_results(
    secure_path: Path,
    r_path: Path,
    output_path: Path,
) -> list[dict[str, str]]:
    secure_rows = read_rows(secure_path)
    r_rows = read_rows(r_path)
    secure_by_key = index_unique(secure_rows, "secure")
    r_by_key = index_unique(r_rows, "R")

    if secure_by_key.keys() != r_by_key.keys():
        secure_only = sorted(secure_by_key.keys() - r_by_key.keys())
        r_only = sorted(r_by_key.keys() - secure_by_key.keys())
        raise ValueError(f"secure/R key mismatch: secure-only={secure_only}, R-only={r_only}")

    output = []
    for secure in secure_rows:
        reference = r_by_key[result_key(secure)]
        secure_burden = float(secure["secure_burden_p"])
        r_burden = float(reference["r_burden_p"])
        burden_abs, burden_rel = errors(secure_burden, r_burden)

        converged = parse_converged(reference["r_skat_davies_converged"])
        if converged:
            secure_skat = float(secure["secure_skat_wh_p"])
            r_skat = float(reference["r_skat_davies_p"])
            skat_abs, skat_rel = errors(secure_skat, r_skat)
        else:
            skat_abs, skat_rel = "NA", "NA"

        output.append({
            "chromosome": secure["chromosome"],
            "gene_index": secure["gene_index"],
            "gene_id": secure["gene_id"],
            "gene_symbol": secure["gene_symbol"],
            "gene_order": secure["gene_order"],
            "phenotype_index": secure["phenotype_index"],
            "phenotype_name": secure["phenotype_name"],
            "secure_burden_p": secure["secure_burden_p"],
            "r_burden_p": reference["r_burden_p"],
            "burden_abs_error": burden_abs,
            "burden_rel_error": burden_rel,
            "secure_skat_wh_p": secure["secure_skat_wh_p"],
            "r_skat_davies_p": reference["r_skat_davies_p"],
            "r_skat_davies_converged": "true" if converged else "false",
            "skat_abs_error": skat_abs,
            "skat_rel_error": skat_rel,
        })

    write_rows(output_path, output, RESULT_FIELDS)
    return output


def valid_number(value: str) -> float | None:
    if value in {"", "NA", "NaN", "nan"}:
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def summarize_group(
    rows: list[dict[str, str]],
    scope: str,
    chromosome: str,
    method: str,
) -> dict[str, object]:
    if method == "burden":
        abs_column = "burden_abs_error"
        rel_column = "burden_rel_error"
        eligible = rows
    else:
        abs_column = "skat_abs_error"
        rel_column = "skat_rel_error"
        eligible = [row for row in rows if parse_converged(row["r_skat_davies_converged"])]

    absolute_values = []
    relative_values = []
    for row in eligible:
        absolute = valid_number(row[abs_column])
        relative = valid_number(row[rel_column])
        if absolute is not None:
            absolute_values.append((row, absolute))
        if relative is not None:
            relative_values.append((row, relative))

    max_absolute = max(absolute_values, key=lambda item: item[1], default=None)
    max_relative = max(relative_values, key=lambda item: item[1], default=None)
    first = rows[0]
    return {
        "scope": scope,
        "chromosome": chromosome,
        "phenotype_index": first["phenotype_index"],
        "phenotype_name": first["phenotype_name"],
        "method": method,
        "valid_count": len(absolute_values),
        "excluded_count": len(rows) - len(absolute_values),
        "max_abs_error": format_number(max_absolute[1] if max_absolute else None),
        "max_abs_error_gene_id": max_absolute[0]["gene_id"] if max_absolute else "",
        "max_rel_error": format_number(max_relative[1] if max_relative else None),
        "max_rel_error_gene_id": max_relative[0]["gene_id"] if max_relative else "",
    }


def build_error_summary(
    rows: list[dict[str, str]],
    include_all: bool,
) -> list[dict[str, object]]:
    groups = defaultdict(list)
    for row in rows:
        groups[(row["chromosome"], row["phenotype_index"])].append(row)

    summary = []
    for chromosome, phenotype in sorted(
        groups,
        key=lambda key: (chromosome_number(key[0]), int(key[1])),
    ):
        grouped = groups[(chromosome, phenotype)]
        for method in ("burden", "skat_davies"):
            summary.append(summarize_group(grouped, "chromosome", chromosome, method))

    if not include_all:
        return summary

    by_phenotype = defaultdict(list)
    for row in rows:
        by_phenotype[row["phenotype_index"]].append(row)
    for phenotype in sorted(by_phenotype, key=int):
        grouped = by_phenotype[phenotype]
        for method in ("burden", "skat_davies"):
            summary.append(summarize_group(grouped, "all", "all", method))
    return summary


def chromosome_number(value: str) -> int:
    value = value.lower()
    return int(value[3:] if value.startswith("chr") else value)


def aggregate_results(
    input_root: Path,
    output_dir: Path,
    chromosomes: Sequence[int] = tuple(range(1, 23)),
) -> list[dict[str, str]]:
    rows = []
    for chromosome in chromosomes:
        directory = input_root / f"chr{chromosome}"
        if not (directory / "_SUCCESS").exists():
            raise ValueError(f"missing {directory / '_SUCCESS'}")
        rows.extend(read_rows(directory / "gene_results.csv"))

    rows.sort(key=lambda row: (
        chromosome_number(row["chromosome"]),
        int(row["gene_order"]),
        int(row["phenotype_index"]),
    ))

    global_index = {}
    for row in rows:
        key = row["chromosome"], row["gene_id"]
        if key not in global_index:
            global_index[key] = len(global_index)
        row["global_gene_index"] = str(global_index[key])

    output_fields = ["global_gene_index", *RESULT_FIELDS]
    write_rows(output_dir / "all_gene_results.csv", rows, output_fields)
    write_rows(
        output_dir / "error_summary.csv",
        build_error_summary(rows, include_all=True),
        SUMMARY_FIELDS,
    )
    return rows
