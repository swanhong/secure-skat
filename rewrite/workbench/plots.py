from __future__ import annotations

import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .results import chromosome_number, parse_converged, read_rows


PLOT_FLOOR = 1e-300


def safe_name(value: str) -> str:
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")
    return name or "phenotype"


def minus_log10(value: str) -> float | None:
    try:
        number = float(value)
    except ValueError:
        return None
    if not math.isfinite(number) or number < 0 or number > 1:
        return None
    return -math.log10(max(number, PLOT_FLOOR))


def scatter_plot(
    rows: list[dict[str, str]],
    secure_column: str,
    reference_column: str,
    title: str,
    output_path: Path,
) -> None:
    points = []
    for row in rows:
        x = minus_log10(row[reference_column])
        y = minus_log10(row[secure_column])
        if x is not None and y is not None:
            points.append((x, y))
    x_values = [point[0] for point in points]
    y_values = [point[1] for point in points]
    high = max(x_values + y_values, default=1.0) * 1.05 + 0.1
    figure, axis = plt.subplots(figsize=(5, 5))
    axis.scatter(x_values, y_values, s=20, color="#333333")
    axis.plot([0, high], [0, high], "--", color="red", linewidth=1)
    axis.set_xlim(0, high)
    axis.set_ylim(0, high)
    axis.set_xlabel("-log10(R p)")
    axis.set_ylabel("-log10(secure p)")
    axis.set_title(f"{title} (n={len(points)})")
    figure.tight_layout()
    figure.savefig(output_path, dpi=150)
    plt.close(figure)


def write_scatter_plots(rows: list[dict[str, str]], output_dir: Path) -> None:
    groups = defaultdict(list)
    for row in rows:
        groups[(row["chromosome"], row["phenotype_index"])].append(row)

    for (chromosome, _), grouped in groups.items():
        phenotype = safe_name(grouped[0]["phenotype_name"])
        chromosome_label = f"chr{chromosome_number(chromosome)}"
        scatter_plot(
            grouped,
            "secure_burden_p",
            "r_burden_p",
            f"Burden: {chromosome_label} {phenotype}",
            output_dir / f"scatter_burden_{chromosome_label}_{phenotype}.png",
        )
        converged = [
            row for row in grouped
            if parse_converged(row["r_skat_davies_converged"])
        ]
        scatter_plot(
            converged,
            "secure_skat_wh_p",
            "r_skat_davies_p",
            f"SKAT-WH vs Davies: {chromosome_label} {phenotype}",
            output_dir / f"scatter_skat_davies_{chromosome_label}_{phenotype}.png",
        )


def manhattan_plot(
    rows: list[dict[str, str]],
    p_column: str,
    title: str,
    output_path: Path,
    test_count: int | None = None,
) -> None:
    points = []
    for row in rows:
        value = minus_log10(row[p_column])
        if value is not None:
            points.append((
                int(row["global_gene_index"]),
                chromosome_number(row["chromosome"]),
                value,
            ))
    x_values = [point[0] for point in points]
    colors = ["#1f77b4" if point[1] % 2 == 0 else "#ff7f0e" for point in points]
    y_values = [point[2] for point in points]
    chromosomes = sorted({point[1] for point in points})
    centers = [
        sum(point[0] for point in points if point[1] == chromosome)
        / sum(1 for point in points if point[1] == chromosome)
        for chromosome in chromosomes
    ]

    figure, axis = plt.subplots(figsize=(12, 4))
    axis.scatter(x_values, y_values, c=colors, s=14)
    correction_count = test_count if test_count is not None else len(rows)
    if correction_count:
        axis.axhline(
            -math.log10(0.05 / correction_count),
            color="red",
            linestyle="--",
            linewidth=1,
        )
    axis.set_xticks(centers)
    axis.set_xticklabels([f"chr{chromosome}" for chromosome in chromosomes])
    axis.set_ylabel("-log10(p)")
    axis.set_title(title)
    figure.tight_layout()
    figure.savefig(output_path, dpi=150)
    plt.close(figure)


def write_manhattan_plots(rows: list[dict[str, str]], output_dir: Path) -> None:
    by_phenotype = defaultdict(list)
    for row in rows:
        by_phenotype[row["phenotype_index"]].append(row)

    for grouped in by_phenotype.values():
        phenotype = safe_name(grouped[0]["phenotype_name"])
        manhattan_plot(
            grouped,
            "secure_burden_p",
            f"Secure Burden: {phenotype}",
            output_dir / f"manhattan_secure_burden_{phenotype}.png",
        )
        manhattan_plot(
            grouped,
            "secure_skat_wh_p",
            f"Secure SKAT-WH: {phenotype}",
            output_dir / f"manhattan_secure_skat_wh_{phenotype}.png",
        )
        manhattan_plot(
            grouped,
            "r_burden_p",
            f"R Burden: {phenotype}",
            output_dir / f"manhattan_r_burden_{phenotype}.png",
        )
        converged = [
            row for row in grouped
            if parse_converged(row["r_skat_davies_converged"])
        ]
        manhattan_plot(
            converged,
            "r_skat_davies_p",
            f"R SKAT-Davies: {phenotype}",
            output_dir / f"manhattan_r_skat_davies_{phenotype}.png",
            test_count=len(grouped),
        )


def write_all_plots(results_path: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = read_rows(results_path)
    write_scatter_plots(rows, output_dir)
    write_manhattan_plots(rows, output_dir)
