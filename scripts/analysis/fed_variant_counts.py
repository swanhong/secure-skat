#!/usr/bin/env python3
"""Count filtered variants per protein-coding gene across chromosomes.

Uses the same annotation and PVAR assignment functions as fed_prep.py, but does
not extract genotypes or run the secure protocol.
"""

import argparse
import csv
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PREP_DIR = Path(__file__).resolve().parents[1] / "preprocessing"
sys.path.insert(0, str(PREP_DIR))
from fed_prep import (  # noqa: E402
    FRAC_PUBONLY,
    FRAC_SHARED,
    load_annotation,
    scan_pvar_into_genes,
)

DEFAULT_PVAR_DIR = Path(
    "~/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/"
    "snpindel/exome/pgen"
).expanduser()
DEFAULT_GENCODE = PREP_DIR / "gencode_v44_pc_genes.bed"


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chroms", nargs="+", default=[f"chr{i}" for i in range(1, 23)])
    parser.add_argument("--annotation-dir", default="~/fed_prep_out")
    parser.add_argument("--pvar-dir", default=str(DEFAULT_PVAR_DIR))
    parser.add_argument("--gencode", default=str(DEFAULT_GENCODE))
    parser.add_argument("--output-dir", default="~/fed_variant_counts")
    parser.add_argument("--mask", default="pLoF")
    parser.add_argument("--max-maf", type=float, default=0.001)
    parser.add_argument("--af", choices=("gnomad", "gvs"), default="gvs")
    parser.add_argument("--probes", type=int, default=50)
    return parser.parse_args()


def chromosome_name(value):
    value = value.strip().lower()
    return value if value.startswith("chr") else f"chr{value}"


def load_genes(path, chrom):
    genes = []
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == chrom:
                genes.append((fields[3], int(fields[1]), int(fields[2])))
    genes.sort(key=lambda row: row[1])
    return {name: (start, end) for name, start, end in genes}


def quantile(values, fraction):
    if not values:
        return 0
    ordered = sorted(values)
    return ordered[round(fraction * (len(ordered) - 1))]


def write_csv(path, fieldnames, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def histogram(values, title, output, zero_count=0):
    positive = np.asarray([value for value in values if value > 0], dtype=float)
    fig, ax = plt.subplots(figsize=(8, 5))
    if positive.size:
        upper = max(2.0, float(positive.max()))
        bins = np.unique(np.geomspace(1, upper + 1, 32).astype(int))
        bins = np.unique(np.append(bins, int(upper) + 2)) - 0.5
        ax.hist(positive, bins=bins, color="#35618f", edgecolor="white")
        ax.set_xscale("log")
    ax.set_xlabel("Filtered variants per nonempty gene (log scale)")
    ax.set_ylabel("Genes")
    ax.set_title(f"{title} (empty genes={zero_count:,})")
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)


def chromosome_histograms(rows, chroms, output):
    columns = 4
    rows_n = math.ceil(len(chroms) / columns)
    fig, axes = plt.subplots(rows_n, columns, figsize=(16, 3.2 * rows_n), squeeze=False)
    for ax, chrom in zip(axes.flat, chroms):
        values = [row["filtered_variants"] for row in rows if row["chromosome"] == chrom]
        positive = np.asarray([value for value in values if value > 0], dtype=float)
        if positive.size:
            upper = max(2.0, float(positive.max()))
            bins = np.unique(np.geomspace(1, upper + 1, 22).astype(int))
            bins = np.unique(np.append(bins, int(upper) + 2)) - 0.5
            ax.hist(positive, bins=bins, color="#35618f", edgecolor="white")
            ax.set_xscale("log")
        ax.set_title(f"{chrom}: {len(positive)}/{len(values)} nonempty")
        ax.set_xlabel("variants")
        ax.set_ylabel("genes")
    for ax in axes.flat[len(chroms):]:
        ax.axis("off")
    fig.suptitle("Filtered variants per gene by chromosome", y=1.0)
    fig.tight_layout()
    fig.savefig(output, dpi=160, bbox_inches="tight")
    plt.close(fig)


def main():
    args = arguments()
    chroms = [chromosome_name(chrom) for chrom in args.chroms]
    annotation_dir = Path(args.annotation_dir).expanduser()
    pvar_dir = Path(args.pvar_dir).expanduser()
    gencode = Path(args.gencode).expanduser()
    output_dir = Path(args.output_dir).expanduser()

    inputs = []
    for chrom in chroms:
        inputs.extend((annotation_dir / f"{chrom}_annotation.tsv", pvar_dir / f"exome.{chrom}.pvar"))
    missing = [path for path in [gencode, *inputs] if not path.is_file()]
    if missing:
        raise SystemExit("missing required files:\n  " + "\n  ".join(map(str, missing)))

    output_dir.mkdir(parents=True, exist_ok=True)
    gene_rows = []
    summary_rows = []
    for chrom in chroms:
        print(f"\n=== {chrom} ===", flush=True)
        genes = load_genes(gencode, chrom)
        annotation = load_annotation(
            annotation_dir / f"{chrom}_annotation.tsv", args.mask, args.max_maf, args.af
        )
        names, keys = scan_pvar_into_genes(pvar_dir / f"exome.{chrom}.pvar", genes, annotation)
        counts = dict(zip(names, map(len, keys)))

        for symbol, (start, end) in genes.items():
            total = counts.get(symbol, 0)
            public = round(FRAC_SHARED * total) + round(FRAC_PUBONLY * total)
            private = total - public
            mode = "excluded" if total == 0 else ("exact" if public <= args.probes else "hutchinson")
            gene_rows.append({
                "chromosome": chrom,
                "gene_symbol": symbol,
                "start": start,
                "end": end,
                "filtered_variants": total,
                "public_variants": public,
                "private_variants": private,
                "trace_mode": mode,
            })

        positive = list(counts.values())
        public_counts = [round(FRAC_SHARED * count) + round(FRAC_PUBONLY * count) for count in positive]
        summary_rows.append({
            "chromosome": chrom,
            "gencode_genes": len(genes),
            "genes_with_variants": len(positive),
            "genes_without_variants": len(genes) - len(positive),
            "filtered_gene_variant_slots": sum(positive),
            "exact_genes": sum(count <= args.probes for count in public_counts),
            "hutchinson_genes": sum(count > args.probes for count in public_counts),
            "min_variants": min(positive, default=0),
            "p25_variants": quantile(positive, 0.25),
            "median_variants": quantile(positive, 0.5),
            "p75_variants": quantile(positive, 0.75),
            "max_variants": max(positive, default=0),
            "mean_variants": f"{np.mean(positive):.2f}" if positive else "0.00",
        })

    gene_fields = list(gene_rows[0]) if gene_rows else []
    summary_fields = list(summary_rows[0]) if summary_rows else []
    write_csv(output_dir / "gene_variant_counts.csv", gene_fields, gene_rows)
    write_csv(output_dir / "chromosome_summary.csv", summary_fields, summary_rows)

    print("\nchromosome  genes  nonempty  filtered slots  exact  hutchinson")
    for row in summary_rows:
        print(
            f"{row['chromosome']:>10}  {row['gencode_genes']:>5}  "
            f"{row['genes_with_variants']:>8}  {row['filtered_gene_variant_slots']:>14}  "
            f"{row['exact_genes']:>5}  {row['hutchinson_genes']:>11}"
        )

    values = [row["filtered_variants"] for row in gene_rows]
    histogram(
        values,
        "Filtered variants per gene, chr1-22",
        output_dir / "variant_histogram_all.png",
        sum(value == 0 for value in values),
    )
    chromosome_histograms(gene_rows, chroms, output_dir / "variant_histogram_by_chromosome.png")

    print(f"\nwrote {output_dir / 'chromosome_summary.csv'}")
    print(f"wrote {output_dir / 'gene_variant_counts.csv'}")
    print(f"wrote {output_dir / 'variant_histogram_all.png'}")
    print(f"wrote {output_dir / 'variant_histogram_by_chromosome.png'}")


if __name__ == "__main__":
    main()
