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
    # Match the main plot's R²: R predictions against the All-by-All reference.
    reference = [x for x, _, _, _ in points]
    score = None
    if len(reference) >= 2:
        mean = sum(reference) / len(reference)
        total = sum((x - mean) ** 2 for x in reference)
        if total:
            score = 1 - sum((y - x) ** 2 for x, y, _, _ in points) / total
    score_text = "NA" if score is None else f"{score:.6f}"
    line = -math.log10(threshold) if threshold is not None else 0
    limit = max([1, line] + [max(x, y) for x, y, _, _ in points]) * 1.18
    lower = -0.02 * limit
    ax.plot([lower, limit], [lower, limit], color="#c00000", ls="--", lw=1, zorder=1)
    if threshold is not None:
        ax.axvline(line, color="#a6a6a6", ls=":", lw=0.8)
        ax.axhline(line, color="#a6a6a6", ls=":", lw=0.8)
    for flagged, color, marker in ((False, "#2f5597", "o"), (True, "#777777", "x")):
        selected = [(x, y) for x, y, _, bad in points if bad == flagged]
        if selected:
            xs, ys = zip(*selected)
            ax.scatter(xs, ys, c=color, marker=marker, s=24, zorder=3)
    significant = [p for p in points if threshold is not None and max(p[:2]) > line]
    labels = sorted(significant, key=lambda p: max(p[:2]), reverse=True)[:label_top]
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
                    fontsize=9, ha="right" if right else "left",
                    arrowprops={"arrowstyle": "-", "color": "#a6a6a6", "lw": 0.6})
    ax.set(xlim=(lower, limit), ylim=(lower, limit), xlabel=r"All-by-All v8  $-\log_{10}(p)$",
           ylabel=r"R v9  $-\log_{10}(p)$")
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(f"{name}\nn={len(points)}, $R^2$={score_text}", fontsize=12)
    flagged = sum(p[3] for p in points)
    ax.grid(alpha=0.2)
    print(f"{test}/{name}: pairs={len(points)}, R2={score_text}, labels={len(labels)}, "
          f"missing/invalid={missing}, zero={zero}, flagged={flagged}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r-dir", type=Path, required=True, help="Directory containing comparison.csv")
    parser.add_argument("--threshold", type=float, help="Label genes with either p below this cutoff; no default")
    parser.add_argument("--label-top", type=int, help="Optional maximum number of significant labels per panel")
    args = parser.parse_args()
    if args.threshold is not None and not 0 < args.threshold <= 1:
        parser.error("--threshold must be in (0, 1]")
    if args.label_top is not None and args.label_top < 0:
        parser.error("--label-top must be nonnegative")
    directory = args.r_dir.expanduser().resolve()
    config = json.loads((directory / "run_info.json").read_text())["settings"]
    rows = read_rows(directory / "comparison.csv")
    key = lambda r: (r["gene_id"], r["axa_phenotype_id"])
    index_rows(rows, key)
    qc = index_rows(read_rows(directory / "comparison_qc.csv"), key)
    phenotypes = config["phenotypes"]
    for test in ("burden", "skat"):
        nrows = math.ceil(len(phenotypes) / 3)
        fig, axes = plt.subplots(nrows, 3, figsize=(15, 5 * nrows), squeeze=False,
                                 constrained_layout=True)
        panels = list(axes.flat)
        for ax, pheno in zip(panels, phenotypes):
            selected = [r for r in rows if r["axa_phenotype_id"] == pheno["axa_id"]]
            scatter(ax, selected, qc, test, pheno["name"], args.threshold, args.label_top)
        for ax in panels[len(phenotypes):]:
            ax.set_axis_off()
        handles = [Line2D([], [], color="#2f5597", marker="o", ls="", label="R status OK"),
                   Line2D([], [], color="#777777", marker="x", ls="", label="R flagged / QC unknown"),
                   Line2D([], [], color="#c00000", ls="--", label="Equal p-values")]
        if args.threshold is not None:
            handles.append(Line2D([], [], color="#a6a6a6", ls=":", label=f"p = {args.threshold:.3g}"))
        fig.legend(handles=handles, loc="lower right", bbox_to_anchor=(0.98, 0.14),
                   frameon=False, fontsize=10)
        title = "Burden" if test == "burden" else "SKAT"
        fig.suptitle(f"{title} · {config['ancestry']} · {config['annotation']} · "
                     f"MAF {config['max_maf']:.1%}", fontsize=16)
        path = directory / f"{test}_scatter.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        print(f"Saved: {path}")


if __name__ == "__main__":
    main()
