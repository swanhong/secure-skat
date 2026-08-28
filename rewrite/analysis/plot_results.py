#!/usr/bin/env python3

import argparse
import math
import tomllib
import statistics
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from compare_results import read_rows, r_squared


SCATTER_COMPARISONS = [
    (
        "burden",
        "Burden",
        "secure_burden_p",
        "r_burden_p",
    ),
    (
        "skat_liu",
        "SKAT WH vs R Liu",
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
        "R Burden",
    ),
    (
        "skat_liu",
        "SKAT",
        "secure_skat_wh_p",
        "Secure SKAT WH",
        "r_skat_liu_p",
        "R SKAT Liu",
    ),
]

def finite_p_value_pairs(
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
        ):
            secure_values.append(secure_value)
            reference_values.append(reference_value)

    return secure_values, reference_values


def write_scatter_plot(
    rows: list[dict[str, str]],
    secure_column: str,
    reference_column: str,
    title: str,
    output_path: Path,
) -> None:
    secure_values, reference_values = finite_p_value_pairs(
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

    figure, axis = plt.subplots(figsize=(5, 5))
    axis.scatter(
        reference_values,
        secure_values,
        color="#2f5597",
        s=24,
    )
    axis.plot(
        [0, 1],
        [0, 1],
        color="#c00000",
        linestyle="--",
        linewidth=1,
    )
    axis.set_xlim(-0.02, 1.02)
    axis.set_ylim(-0.02, 1.02)
    axis.set_xlabel("R p-value")
    axis.set_ylabel("Secure p-value")
    axis.set_title(f"{title}\nn={count}, $R^2$={score_text}")
    axis.grid(alpha=0.2)

    figure.tight_layout()
    figure.savefig(output_path, dpi=150)
    plt.close(figure)

def read_gene_positions(
    run_dir: Path,
    chromosomes: list[int],
) -> dict[tuple[str, str, str], float]:
    gene_positions = {}

    for chromosome in chromosomes:
        chromosome_dir = run_dir / "prepared" / f"chr{chromosome}"

        gene_ids = (
            chromosome_dir / "genes.txt"
        ).read_text(encoding="utf-8").splitlines()
        block_sizes = [
            int(value)
            for value in (
                chromosome_dir / "block_sizes.txt"
            ).read_text(encoding="utf-8").splitlines()
        ]
        position_rows = [
            line.split()
            for line in (
                chromosome_dir / "pos.txt"
            ).read_text(encoding="utf-8").splitlines()
        ]

        if len(gene_ids) != len(block_sizes):
            raise ValueError(
                f"chr{chromosome}: genes and block sizes do not match"
            )
        if sum(block_sizes) != len(position_rows):
            raise ValueError(
                f"chr{chromosome}: block sizes and positions do not match"
            )

        offset = 0
        for gene_index, (gene_id, block_size) in enumerate(
            zip(gene_ids, block_sizes)
        ):
            block_positions = position_rows[
                offset:offset + block_size
            ]
            offset += block_size

            if not block_positions:
                raise ValueError(
                    f"chr{chromosome}: gene {gene_id} has no position"
                )
            if any(
                int(row[0]) != chromosome
                for row in block_positions
            ):
                raise ValueError(
                    f"chr{chromosome}: position chromosome mismatch"
                )

            gene_positions[
                (str(chromosome), str(gene_index), gene_id)
            ] = float(
                statistics.median(
                    int(row[1])
                    for row in block_positions
                )
            )

    return gene_positions


def negative_log10(p_value: float) -> float:
    return -math.log10(max(p_value, 1e-300))


def write_manhattan_plot(
    rows: list[dict[str, str]],
    gene_positions: dict[tuple[str, str, str], float],
    secure_column: str,
    secure_label: str,
    reference_column: str,
    reference_label: str,
    title: str,
    output_path: Path,
) -> None:
    positioned_rows = []

    for row in rows:
        key = (
            row["chromosome"],
            row["gene_index"],
            row["gene_id"],
        )
        if key not in gene_positions:
            raise ValueError(
                f"gene position not found: {key}"
            )

        positioned_rows.append(
            (
                int(row["chromosome"]),
                gene_positions[key],
                int(row["gene_index"]),
                row,
            )
        )

    positioned_rows.sort(
        key=lambda item: (item[0], item[1], item[2])
    )

    chromosomes = sorted({
        chromosome
        for chromosome, _, _, _ in positioned_rows
    })
    chromosome_indices = {
        chromosome: [
            index
            for index, (row_chromosome, _, _, _) in enumerate(
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

        for index, (chromosome, _, _, row) in enumerate(
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

def plot_results(config_path: Path) -> None:
    with config_path.open("rb") as config_file:
        config = tomllib.load(config_file)

    run_dir = Path(config["run_dir"])
    comparison_dir = run_dir / "comparison"
    plots_dir = comparison_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    success_path = plots_dir / "_SUCCESS"
    success_path.unlink(missing_ok=True)

    rows = read_rows(comparison_dir / "all_comparison.csv")
    chromosomes = sorted({
        int(row["chromosome"])
        for row in rows
    })
    gene_positions = read_gene_positions(
        run_dir,
        chromosomes,
    )

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
                gene_positions=gene_positions,
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        default="run.1kg.conf",
        help="path to the run configuration",
    )
    args = parser.parse_args()

    plot_results(Path(args.config))


if __name__ == "__main__":
    main()