"""Plot MAC categories and, when variable, nonmissing sample counts."""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", nargs="?", type=Path, default=Path("metadata/variant_summary.csv"))
    args = parser.parse_args()
    with args.csv.open(newline="") as file:
        variants = {r["variant_id"]: (int(r["MAC"]), int(r["N"])) for r in csv.DictReader(file)}
    if not variants:
        raise ValueError(f"No variants in {args.csv}")
    mac, n = np.array(list(variants.values())).T
    panels = 1 + int(np.ptp(n) > 0)
    figure, axes = plt.subplots(1, panels, figsize=(6 * panels, 4), squeeze=False, layout="constrained")
    axes = axes[0]
    counts, _ = np.histogram(mac, bins=[0, 1, 2, 3, 6, 11, np.inf])
    bars = axes[0].bar(("0", "1", "2", "3–5", "6–10", ">10"), counts, color="#4477AA")
    axes[0].bar_label(bars, labels=[f"{v:,}\n({v / len(mac):.1%})" for v in counts], padding=3)
    axes[0].set(title="Minor allele count", xlabel="MAC", ylim=(0, counts.max() * 1.25))
    if panels == 2:
        axes[1].hist(n, bins=min(40, int(np.ptp(n)) + 1), color="#228833", alpha=0.85)
        axes[1].set(title="Available samples per variant", xlabel="N (nonmissing samples)")
    else:
        axes[0].set_title(f"Minor allele count (N = {n[0]:,} for all variants)")
    for axis in axes:
        axis.set_ylabel("Number of unique variants")
        axis.ticklabel_format(axis="y", style="plain", useOffset=False)
        axis.spines[["top", "right"]].set_visible(False)
        axis.set_axisbelow(True)
        axis.grid(axis="y", alpha=0.15)
    figure.suptitle(f"Variant summary ({len(mac):,} unique variants)")
    for extension in (".png", ".pdf"):
        output = args.csv.with_suffix(extension)
        figure.savefig(output, dpi=200)
        print(f"Saved {output}")
    plt.close(figure)


if __name__ == "__main__":
    main()
