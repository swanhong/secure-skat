#!/usr/bin/env python3

import argparse
import math
import tomllib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from compare_secure_to_reference import read_rows, r_squared


SCATTER_COMPARISONS = [
    (
        "burden",
        "Burden",
        "secure_burden_p",
        "r_burden_p",
    ),
    (
        "skat_liu",
        "SKAT WH vs Reference Liu",
        "secure_skat_wh_p",
        "r_skat_liu_p",
    ),
]
MANHATTAN_COMPARISONS = [
    (
        "burden",
        "Burden",
        "secure_burden_p",
        "Secure Burden",
        "r_burden_p",
        "Reference Burden",
    ),
    (
        "skat_liu",
        "SKAT",
        "secure_skat_wh_p",
        "Secure SKAT WH",
        "r_skat_liu_p",
        "Reference SKAT Liu",
    ),
]

def negative_log10(p_value: float) -> float:
    return -math.log10(max(p_value, 1e-300))


def negative_log10_p_value_pairs(
    rows: list[dict[str, str]],
    secure_column: str,
    reference_column: str,
) -> tuple[list[float], list[float]]:
    secure_values = []
    reference_values = []

    for row in rows:
        try:
            secure_value = float(row[secure_column])
            reference_value = float(row[reference_column])
        except ValueError:
            continue

        if (
            math.isfinite(secure_value)
            and math.isfinite(reference_value)
            and secure_value >= 0
            and reference_value >= 0
        ):
            secure_values.append(negative_log10(secure_value))
            reference_values.append(negative_log10(reference_value))

    return secure_values, reference_values


def write_scatter_plot(
    rows: list[dict[str, str]],
    secure_column: str,
    reference_column: str,
    title: str,
    output_path: Path,
) -> None:
    secure_values, reference_values = negative_log10_p_value_pairs(
        rows,
        secure_column,
        reference_column,
    )
    count, score = r_squared(
        rows,
        secure_column,
        reference_column,
    )
    score_text = "NA" if score is None else f"{score:.6f}"

    plot_values = [0.0, *secure_values, *reference_values]
    lower = min(plot_values)
    upper = max(plot_values)
    padding = max(0.02 * (upper - lower), 0.02)
    limits = [lower - padding, upper + padding]

    figure, axis = plt.subplots(figsize=(5, 5))
    axis.scatter(
        reference_values,
        secure_values,
        color="#2f5597",
        s=24,
    )
    axis.plot(
        limits,
        limits,
        color="#c00000",
        linestyle="--",
        linewidth=1,
    )
    axis.set_xlim(limits)
    axis.set_ylim(limits)
    axis.set_xlabel("Reference -log10(p)")
    axis.set_ylabel("Secure -log10(p)")
    axis.set_title(
        f"{title}\nn={count}, $R^2$ (-log10 p)={score_text}"
    )
    axis.grid(alpha=0.2)

    figure.tight_layout()
    figure.savefig(output_path, dpi=150)
    plt.close(figure)


def write_manhattan_plot(
    rows: list[dict[str, str]],
    secure_column: str,
    secure_label: str,
    reference_column: str,
    reference_label: str,
    title: str,
    output_path: Path,
) -> None:
    positioned_rows = [
        (
            int(row["chromosome"]),
            int(row["gene_index"]),
            row,
        )
        for row in rows
    ]

    positioned_rows.sort(
        key=lambda item: (item[0], item[1])
    )

    chromosomes = sorted({
        chromosome
        for chromosome, _, _ in positioned_rows
    })
    chromosome_indices = {
        chromosome: [
            index
            for index, (row_chromosome, _, _) in enumerate(
                positioned_rows
            )
            if row_chromosome == chromosome
        ]
        for chromosome in chromosomes
    }
    chromosome_colors = {
        chromosome: (
            "#2f5597" if index % 2 == 0 else "#70ad47"
        )
        for index, chromosome in enumerate(chromosomes)
    }

    figure, axes = plt.subplots(
        2,
        1,
        figsize=(11, 7),
        sharex=True,
        sharey=True,
    )

    panels = [
        (secure_column, secure_label),
        (reference_column, reference_label),
    ]
    bonferroni = negative_log10(
        0.05 / len(positioned_rows)
    )

    for axis, (p_value_column, panel_label) in zip(
        axes,
        panels,
    ):
        x_values = []
        y_values = []
        colors = []

        for index, (chromosome, _, row) in enumerate(
            positioned_rows
        ):
            try:
                p_value = float(row[p_value_column])
            except ValueError:
                continue

            if not math.isfinite(p_value) or p_value < 0:
                continue

            x_values.append(index)
            y_values.append(negative_log10(p_value))
            colors.append(chromosome_colors[chromosome])

        axis.scatter(
            x_values,
            y_values,
            c=colors,
            s=22,
        )
        axis.axhline(
            bonferroni,
            color="#c00000",
            linestyle="--",
            linewidth=1,
            label=f"Bonferroni 0.05/{len(positioned_rows)}",
        )
        axis.set_ylabel("-log10(p)")
        axis.set_title(panel_label)
        axis.grid(axis="y", alpha=0.2)

        for chromosome in chromosomes[:-1]:
            boundary = max(chromosome_indices[chromosome]) + 0.5
            axis.axvline(
                boundary,
                color="#bfbfbf",
                linewidth=0.7,
            )

    axes[0].legend(loc="upper right")
    axes[-1].set_xticks([
        sum(chromosome_indices[chromosome])
        / len(chromosome_indices[chromosome])
        for chromosome in chromosomes
    ])
    axes[-1].set_xticklabels([
        f"chr{chromosome}"
        for chromosome in chromosomes
    ])
    axes[-1].set_xlabel("Gene order")

    figure.suptitle(title)
    figure.tight_layout()
    figure.savefig(output_path, dpi=150)
    plt.close(figure)


def plot_ancestry(run_dir: Path, ancestry: str) -> None:
    comparison_dir = run_dir / "comparison" / ancestry
    plots_dir = comparison_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    success_path = plots_dir / "_SUCCESS"
    success_path.unlink(missing_ok=True)

    rows = read_rows(comparison_dir / "all_comparison.csv")

    chromosome_groups = {}
    phenotype_groups = {}

    for row in rows:
        chromosome_key = (
            int(row["phenotype_index"]),
            row["phenotype_name"],
            int(row["chromosome"]),
        )
        chromosome_groups.setdefault(
            chromosome_key,
            [],
        ).append(row)

        phenotype_key = (
            int(row["phenotype_index"]),
            row["phenotype_name"],
        )
        phenotype_groups.setdefault(
            phenotype_key,
            [],
        ).append(row)

    written = []

    for (
        phenotype_index,
        phenotype_name,
        chromosome,
    ), group in sorted(chromosome_groups.items()):
        for (
            tag,
            label,
            secure_column,
            reference_column,
        ) in SCATTER_COMPARISONS:
            output_path = plots_dir / (
                f"scatter_{tag}_pheno{phenotype_index}_"
                f"chr{chromosome}.png"
            )
            write_scatter_plot(
                rows=group,
                secure_column=secure_column,
                reference_column=reference_column,
                title=(
                    f"{label}: {phenotype_name}, "
                    f"chromosome {chromosome}"
                ),
                output_path=output_path,
            )
            written.append(output_path)

    for (
        phenotype_index,
        phenotype_name,
    ), group in sorted(phenotype_groups.items()):
        for (
            tag,
            label,
            secure_column,
            secure_label,
            reference_column,
            reference_label,
        ) in MANHATTAN_COMPARISONS:
            output_path = plots_dir / (
                f"manhattan_{tag}_pheno{phenotype_index}.png"
            )
            write_manhattan_plot(
                rows=group,
                secure_column=secure_column,
                secure_label=secure_label,
                reference_column=reference_column,
                reference_label=reference_label,
                title=f"{label}: {phenotype_name}",
                output_path=output_path,
            )
            written.append(output_path)

    success_path.touch()
    print(f"Wrote {len(written)} plots to {plots_dir}")


def plot_secure_vs_reference(config_path: Path) -> None:
    with config_path.open("rb") as config_file:
        config = tomllib.load(config_file)

    run_dir = Path(config["run_dir"])
    success_path = run_dir / "comparison" / "_PLOTS_SUCCESS"
    success_path.unlink(missing_ok=True)

    for ancestry in config["ancestries"]:
        plot_ancestry(run_dir, ancestry)

    success_path.touch()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        default="run.1kg.conf",
        help="path to the run configuration",
    )
    args = parser.parse_args()

    plot_secure_vs_reference(Path(args.config))


if __name__ == "__main__":
    main()
