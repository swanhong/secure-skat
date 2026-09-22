"""Plot three count distributions from gene_summary.csv."""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", nargs="?", type=Path, default=Path("metadata/gene_summary.csv"))
    args = parser.parse_args()
    with args.csv.open(newline="") as file:
        rows = list(csv.DictReader(file))
    if not rows:
        raise ValueError(f"No genes in {args.csv}")

    figure, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True, layout="constrained")
    for axis, column, color in zip(
        axes, ("N_variants", "nz_samples", "total_MAC"), ("#4477AA", "#228833", "#AA3377")
    ):
        values = np.array([int(row[column]) for row in rows])
        median = np.median(values)
        axis.hist(np.log10(values + 1), bins=40, color=color, alpha=0.85)
        axis.axvline(np.log10(median + 1), color="#333333", linestyle="--", linewidth=1)
        axis.set(title=column, xlabel=f"log10({column} + 1)")
        axis.text(0.97, 0.97, f"Median: {median:,.1f}\nZero: {np.mean(values == 0):.1%}",
                  transform=axis.transAxes, ha="right", va="top")
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", alpha=0.15)
    axes[0].set_ylabel("Number of genes")
    figure.suptitle(f"Gene summary ({len(rows):,} genes)")
    for extension in (".png", ".pdf"):
        output = args.csv.with_suffix(extension)
        figure.savefig(output, dpi=200)
        print(f"Saved {output}")
    plt.close(figure)


if __name__ == "__main__":
    main()
