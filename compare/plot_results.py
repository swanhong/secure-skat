#!/usr/bin/env python3
"""Plot merged R v9 / All-by-All v8 results without changing input files."""

import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from common import index_rows, read_rows, valid_p


def paired_points(rows, qc, test):
    points, missing, zero = [], 0, 0
    for row in rows:
        x, y = valid_p(row[f"axa_{test}"]), valid_p(row[f"R_{test}"])
        if x is None or y is None:
            missing += 1
            continue
        if x == 0 or y == 0:
            zero += 1
            continue
        status = qc.get((row["gene_id"], row["axa_phenotype_id"]), {})
        flagged = status.get(f"R_{test}_status") != "ok"
        if test == "skat":
            flagged |= status.get("R_skat_converged") != "1"
        points.append((-math.log10(x), -math.log10(y),
                       row["gene_symbol"] or row["gene_id"], flagged))
    return points, missing, zero


def scatter(ax, rows, qc, test, name, threshold, label_top):
    points, missing, zero = paired_points(rows, qc, test)
    line = -math.log10(threshold) if threshold is not None else 0
    limit = max([1, line] + [max(x, y) for x, y, _, _ in points]) * 1.18
    ax.plot([0, limit], [0, limit], color="0.6", lw=1, zorder=1)
    if threshold is not None:
        ax.axvline(line, color="0.4", ls="--", lw=0.8)
        ax.axhline(line, color="0.4", ls="--", lw=0.8)
    for flagged, color, marker in ((False, "#2375AB", "o"), (True, "#C75D16", "x")):
        selected = [(x, y) for x, y, _, bad in points if bad == flagged]
        if selected:
            xs, ys = zip(*selected)
            ax.scatter(xs, ys, c=color, marker=marker, s=35, alpha=0.75, zorder=3)
    labels = sorted(points, key=lambda p: max(p[:2]), reverse=True)[:label_top]
    placed = []
    gap = min(0.065, 0.75 / max(1, len(labels))) * limit
    for rank, (x, y, symbol, _) in enumerate(sorted(labels, key=lambda p: p[1])):
        # Spread nearby labels vertically; a connector keeps their point unambiguous.
        label_y = min(y + limit * 0.025, limit * 0.94 - gap * (len(labels) - rank - 1))
        for old_x, old_y in placed:
            if abs(x - old_x) < limit * 0.3 and abs(label_y - old_y) < gap:
                label_y = old_y + gap
        placed.append((x, label_y))
        right = x > limit * 0.7
        ax.annotate(symbol, (x, y), (x + limit * (-0.025 if right else 0.025), label_y),
                    fontsize=8, ha="right" if right else "left",
                    arrowprops={"arrowstyle": "-", "color": "0.65", "lw": 0.5})
    ax.set(xlim=(0, limit), ylim=(0, limit), xlabel="All-by-All v8  -log10(p)",
           ylabel="R v9  -log10(p)")
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(f"{name}  |  pairs={len(points)}/{len(rows)}", fontsize=11)
    flagged = sum(p[3] for p in points)
    ax.text(0.02, 0.98, f"Missing/invalid: {missing}  |  zero p: {zero}  |  flagged: {flagged}",
            transform=ax.transAxes, va="top", fontsize=7,
            bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"})
    ax.grid(alpha=0.15)
    print(f"{test}/{name}: pairs={len(points)}, missing/invalid={missing}, zero={zero}, flagged={flagged}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r-dir", type=Path, required=True, help="Directory containing comparison.csv")
    parser.add_argument("--threshold", type=float, help="Optional common p-value cutoff; no default")
    parser.add_argument("--label-top", type=int, default=7, help="Labels per panel, ranked by stronger signal")
    args = parser.parse_args()
    if args.threshold is not None and not 0 < args.threshold <= 1:
        parser.error("--threshold must be in (0, 1]")
    if args.label_top < 0:
        parser.error("--label-top must be nonnegative")
    directory = args.r_dir.expanduser().resolve()
    config = json.loads((directory / "run_info.json").read_text())["settings"]
    rows = read_rows(directory / "comparison.csv")
    key = lambda r: (r["gene_id"], r["axa_phenotype_id"])
    index_rows(rows, key)
    qc = index_rows(read_rows(directory / "comparison_qc.csv"), key)
    phenotypes = config["phenotypes"]
    for test in ("burden", "skat"):
        nrows = math.ceil((len(phenotypes) + 1) / 3)
        fig, axes = plt.subplots(nrows, 3, figsize=(15, 5 * nrows), squeeze=False,
                                 constrained_layout=True)
        panels = list(axes.flat)
        for ax, pheno in zip(panels, phenotypes):
            selected = [r for r in rows if r["axa_phenotype_id"] == pheno["axa_id"]]
            scatter(ax, selected, qc, test, pheno["name"], args.threshold, args.label_top)
        for ax in panels[len(phenotypes):]:
            ax.set_axis_off()
        legend = panels[len(phenotypes)]
        legend.legend(handles=[
            Line2D([], [], color="#2375AB", marker="o", ls="", label="R status OK"),
            Line2D([], [], color="#C75D16", marker="x", ls="", label="R flagged / QC unknown"),
            Line2D([], [], color="0.6", label="Equal p-values (y = x)"),
        ], loc="upper left", frameon=False)
        note = "No significance cutoff applied."
        if args.threshold is not None:
            note = f"Common cutoff: p = {args.threshold:g}\nNot an All-by-All official cutoff."
        legend.text(0.03, 0.57, note + "\n\nMissing/invalid pairs are omitted.\n"
                    "Remaining pairs with p=0 are omitted\n(log10 is undefined); counts are shown.\n"
                    "Flagged points are retained.\nAxes use the same scale within each panel.\n\n"
                    "Cross-release comparison; selected genes\nare not a genome-wide validation.",
                    transform=legend.transAxes, va="top", fontsize=10, linespacing=1.6)
        fig.suptitle(f"{test.upper()} | {config['ancestry']} | {config['annotation']} | "
                     f"max MAF {config['max_maf']:g}\nR v9 vs All-by-All v8", fontsize=17)
        path = directory / f"{test}_scatter.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        print(f"Saved: {path}")


if __name__ == "__main__":
    main()
