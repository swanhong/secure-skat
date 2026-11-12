#!/usr/bin/env python3
"""Plot panels A/B/C from one site's resource_stats.tsv; no other inputs required."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np

COLORS = {"MVP": "#2878A2", "AoU": "#E58B37", "joint": "#B34255", "local_any": "#717782", "shared": "#459B82"}
SERIES = {"A": {"MVP", "AoU", "local_any", "joint"}, "B": {"MVP", "AoU", "joint"}, "C": {"MVP_only", "AoU_only", "shared"}}


def load_tables(path, chromosomes=None):
    """Keep ancestry/MAF strata separate and sum matching chromosome aggregates."""
    groups = defaultdict(lambda: defaultdict(dict))
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"ancestry", "maf_threshold", "annotation", "chromosome", "gene_set", "panel", "series", "lower", "upper", "count", "status"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("input is missing resource_stats.tsv columns")
        for row in reader:
            chromosome = int(row["chromosome"])
            if chromosomes and chromosome not in chromosomes:
                continue
            if row["status"] != "released" or int(row["count"]) < 0:
                raise ValueError("only released, nonnegative counts can be plotted")
            group = (row["ancestry"], float(row["maf_threshold"]), row["annotation"], row["gene_set"])
            cell = (row["panel"], row["series"], row["lower"], row["upper"])
            table = groups[group][chromosome]
            if cell in table:
                raise ValueError(f"duplicate cell: {group}, chr{chromosome}, {cell}")
            table[cell] = int(row["count"])
    if not groups:
        raise ValueError("no rows selected")
    for group, by_chr in sorted(groups.items()):
        if chromosomes and set(by_chr) != set(chromosomes):
            raise ValueError(f"missing requested chromosomes in {group}")
        cells = set(next(iter(by_chr.values())))
        if any(set(table) != cells for table in by_chr.values()):
            raise ValueError(f"chromosomes have different bins/series in {group}")
        for panel, series in SERIES.items():
            if {s for p, s, lo, hi in cells if p == panel} != series:
                raise ValueError(f"missing or unexpected series in panel {panel}: {group}")
        totals = {cell: sum(table[cell] for table in by_chr.values()) for cell in cells}
        yield group, sorted(by_chr), totals


def plot_panels(group, chromosomes, table, labels, title_prefix):
    ancestry, maf, annotation, gene_set = group
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), gridspec_kw={"width_ratios": [1, 1.3, 1]})
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_axisbelow(True)
        ax.grid(axis="y", color="#E8EBEE", linewidth=0.7)
        ax.tick_params(labelsize=9)

    # A: exact counts at the released thresholds, without fitting a smooth curve.
    ax = axes[0]
    thresholds = sorted({int(lo) for p, s, lo, hi in table if p == "A"})
    styles = {"MVP": ("o", "-", 7), "AoU": ("s", "--", 5), "local_any": ("x", ":", 8), "joint": ("D", "-", 4)}
    for series in ("MVP", "AoU", "local_any", "joint"):
        values = [table["A", series, str(t), "NA"] for t in thresholds]
        marker, line, size = styles[series]
        ax.plot(thresholds, values, color=COLORS[series], marker=marker, linestyle=line,
                markersize=size, markerfacecolor="white", linewidth=1.6, label=labels[series])
    ax.set_xscale("log")
    ax.set_xticks(thresholds, [str(t) for t in thresholds], rotation=40, ha="right")
    ax.minorticks_off()
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_ylim(0, max(1, max(v for (p, s, lo, hi), v in table.items() if p == "A")) * 1.25)
    ax.set_xlabel("Carrier threshold per gene (log scale)")
    ax.set_ylabel("Genes with carriers ≥ threshold")
    ax.set_title("A   Carrier coverage", loc="left", fontweight="bold", pad=25)
    ax.text(0, 1.02, "Evaluated at the released thresholds", transform=ax.transAxes, fontsize=9, color="#555555")
    ax.legend(frameon=False, fontsize=8, loc="upper right")

    # B: all three histograms use the same eligible genes and existing bin edges.
    ax = axes[1]
    bins = sorted({(int(lo), hi) for p, s, lo, hi in table if p == "B"})
    heights = {s: [table["B", s, str(lo), hi] for lo, hi in bins] for s in ("MVP", "AoU", "joint")}
    eligible = {sum(values) for values in heights.values()}
    if len(eligible) != 1:
        raise ValueError("panel B series do not cover the same number of genes")
    eligible = eligible.pop()
    positions = np.arange(len(bins))
    for i, (series, values) in enumerate(heights.items()):
        ax.bar(positions + (i - 1) * 0.25, values, width=0.24, color=COLORS[series], label=labels[series])
    bin_labels = [f"≥{lo:,}" if hi == "Inf" else f"{lo:,}–{int(hi)-1:,}" for lo, hi in bins]
    ax.set_xticks(positions, bin_labels, rotation=45, ha="right")
    ax.set_ylim(0, max(1, max(max(v) for v in heights.values())) * 1.3)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("Gene cMAC (predefined bins)")
    ax.set_ylabel("Gene count")
    ax.set_title("B   cMAC distribution", loc="left", fontweight="bold", pad=25)
    ax.text(0, 1.02, f"Carrier ≥20 in both sites · {eligible:,} eligible gene{'s' if eligible != 1 else ''}", transform=ax.transAxes, fontsize=9, color="#555555")
    if eligible == 0:
        ax.text(0.5, 0.5, "No genes meet the\nboth-site carrier threshold", transform=ax.transAxes,
                ha="center", va="center", fontsize=10, color="#555555")
    ax.legend(frameon=False, fontsize=8, loc="upper right")

    # C: stack disjoint overlap components; shared variants appear once in joint.
    ax = axes[2]
    shared = table["C", "shared", "NA", "NA"]
    only_a = table["C", "MVP_only", "NA", "NA"]
    only_b = table["C", "AoU_only", "NA", "NA"]
    totals = np.array([shared + only_a, shared + only_b, shared + only_a + only_b])
    components = [("Shared", [shared] * 3, COLORS["shared"]),
                  (f"{labels['MVP']} only", [only_a, 0, only_a], COLORS["MVP"]),
                  (f"{labels['AoU']} only", [0, only_b, only_b], COLORS["AoU"])]
    left = np.zeros(3)
    for name, values, color in components:
        ax.barh([2, 1, 0], values, left=left, height=0.55, color=color, label=name)
        for i, value in enumerate(values):
            if value > 0 and value >= max(totals.max(), 1) * 0.09:
                ax.text(left[i] + value / 2, 2-i, f"{value:,}", ha="center", va="center", color="white", fontsize=9)
        left += values
    for i, value in enumerate(totals):
        ax.text(value + max(totals.max(), 1) * 0.025, 2-i, f"{value:,}", va="center", fontsize=10)
    ax.set_yticks([2, 1, 0], [labels["MVP"], labels["AoU"], labels["joint"]])
    ax.set_xlim(0, max(totals.max(), 1) * 1.25)
    ax.set_ylim(-0.6, 3.8)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", color="#E8EBEE", linewidth=0.7)
    ax.set_xlabel("Distinct observed variants")
    ax.set_title("C   Variant coverage", loc="left", fontweight="bold", pad=25)
    ax.text(0, 1.02, "Shared variants counted once in joint", transform=ax.transAxes, fontsize=9, color="#555555")
    ax.legend(frameon=False, fontsize=8, loc="upper left", ncol=1)

    prefix = f"{title_prefix}  |  " if title_prefix else ""
    fig.suptitle(f"{prefix}{ancestry}  ·  pooled MAF ≤ {maf:.3g}  ·  chr{', '.join(map(str, chromosomes))}",
                 fontsize=13, fontweight="bold", y=0.98)
    fig.text(0.06, 0.015, f"Input annotation label: {annotation}  |  Gene set: {gene_set}  |  Chromosomes summed; ancestry and MAF strata kept separate",
             fontsize=8, color="#555555")
    fig.tight_layout(rect=(0, 0.055, 1, 0.91), w_pad=2.2)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="one site's resource_stats.tsv")
    parser.add_argument("--out-dir", type=Path, help="default: plots/ beside the input TSV")
    parser.add_argument("--chromosomes", nargs="+", type=int, help="sum only these chromosomes")
    parser.add_argument("--label-a", default="MVP")
    parser.add_argument("--label-b", default="AoU")
    parser.add_argument("--title-prefix", default="")
    args = parser.parse_args()
    out = args.out_dir or args.input.parent / "plots"
    out.mkdir(parents=True, exist_ok=True)
    labels = {"MVP": args.label_a, "AoU": args.label_b, "local_any": "Either site", "joint": "Joint"}
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "pdf.fonttype": 42})
    for group, chromosomes, table in load_tables(args.input, args.chromosomes):
        fig = plot_panels(group, chromosomes, table, labels, args.title_prefix)
        ancestry, maf, annotation, gene_set = group
        stem = f"joint_stats_{ancestry}_maf{maf:g}_{annotation}_{gene_set}"
        for extension in ("png", "pdf"):
            path = out / f"{stem}.{extension}"
            fig.savefig(path, dpi=220, facecolor="white")
            print(path)
        plt.close(fig)


if __name__ == "__main__":
    main()
