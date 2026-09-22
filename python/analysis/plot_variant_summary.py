"""Plot variant counts and percentages by MAC category."""

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
        variants = {r["variant_id"]: int(r["MAC"]) for r in csv.DictReader(file)}
    if not variants:
        raise ValueError(f"No variants in {args.csv}")
    mac = list(variants.values())
    figure, axis = plt.subplots(figsize=(6, 4), layout="constrained")
    counts, _ = np.histogram(mac, bins=[0, 1, 2, 3, 6, 11, np.inf])
    bars = axis.bar(("0", "1", "2", "3–5", "6–10", ">10"), counts, color="#4477AA")
    axis.bar_label(bars, labels=[f"{v:,}\n({v / len(mac):.1%})" for v in counts], padding=3)
    axis.set(title="Minor allele count", xlabel="MAC", ylabel="Number of unique variants",
             ylim=(0, counts.max() * 1.25))
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
