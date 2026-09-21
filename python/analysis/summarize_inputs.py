"""Summarize Step 0 genotypes without A/B or phenotype/covariate selection.

Run from the repository root:
    python3 -m python.analysis.summarize_inputs --config config/aou
"""

import argparse
import csv
from pathlib import Path
from tempfile import TemporaryDirectory
import tomllib

import numpy as np
import pgenlib

from python.preprocessing.input import load_inputs
from python.preprocessing.pipeline import select_gene_variants


def write_summaries(config: dict, output_dir: Path) -> None:
    chromosomes = config["chromosomes"]
    if not chromosomes or len(set(chromosomes)) != len(chromosomes):
        raise ValueError("chromosomes must be nonempty and contain no duplicates")
    if any(chromosome not in range(1, 23) for chromosome in chromosomes):
        raise ValueError("summaries support diploid autosomes 1-22 only")
    mask = dict([part.strip() for part in item.split("=", 1)] for item in config.get("masks", []))
    sample_ids = None
    variant_stats = {}
    genes = set()

    with (
        (output_dir / "variant_summary.csv").open("w", newline="") as variant_file,
        (output_dir / "gene_summary.csv").open("w", newline="") as gene_file,
    ):
        variant_rows = csv.writer(variant_file)
        variant_rows.writerow(("variant_id", "gene", "N", "MAC", "MAF"))
        gene_rows = csv.writer(gene_file)
        gene_rows.writerow(("chromosome", "gene", "N_variants", "nz_samples", "total_MAC", "N_variants_MAC_gt_0"))

        for chromosome in chromosomes:
            print(f"Summarizing chromosome {chromosome}", flush=True)
            inputs = load_inputs(
                pgen_prefix=Path(config["genotype"].format(chromosome=chromosome)),
                gene_panel_path=Path(config["gene_panel"].format(chromosome=chromosome)),
                annotation_path=Path(config["annotation"].format(chromosome=chromosome)),
            )
            if sample_ids is None:
                sample_ids = inputs.psam_ids
                if not sample_ids or len(set(sample_ids)) != len(sample_ids):
                    raise ValueError("PSAM must contain nonempty, unique sample IDs")
            elif inputs.psam_ids != sample_ids:
                raise ValueError("ordered PSAM sample IDs differ across chromosomes")

            groups = select_gene_variants(
                gene_panel=inputs.gene_panel, variants=inputs.variants,
                annotations=inputs.annotations, annotation_columns=inputs.annotation_columns,
                chromosome=str(chromosome), gene_selection="all",
                mask=mask, max_maf=config.get("max_maf"),
            )
            wanted = {variant.key for group in groups for variant in group.variants}
            prefix = str(inputs.pgen_prefix).encode()
            with pgenlib.PvarReader(prefix + b".pvar") as pvar:
                variant_indices = {}
                for index in range(pvar.get_variant_ct()):
                    key = pvar.get_variant_id(index).decode()
                    if key not in wanted or pvar.get_allele_ct(index) != 2:
                        continue
                    if key in variant_indices:
                        raise ValueError(f"duplicate variant ID: {key}")
                    if pvar.get_allele_code(index, 1).decode() != key.split(":", 3)[-1]:
                        raise ValueError(f"{key}: counted allele is not the key's ALT")
                    variant_indices[key] = index

                with pgenlib.PgenReader(
                    prefix + b".pgen", pvar=pvar, raw_sample_ct=len(sample_ids),
                ) as reader:
                    dosages = np.empty(len(sample_ids), dtype=np.float32)
                    for group in groups:
                        keys = [v.key for v in group.variants if v.key in variant_indices]
                        if not keys:
                            continue
                        gene = group.gene.gene_id
                        if gene in genes:
                            raise ValueError(f"duplicate gene ID: {gene}")
                        genes.add(gene)
                        nonzero_any = np.zeros(len(sample_ids), dtype=bool)
                        total_mac = 0
                        for key in keys:
                            reader.read_dosages(variant_indices[key], dosages)
                            observed = dosages != -9
                            nonzero_any |= dosages > 0
                            if key not in variant_stats:
                                if not np.isin(dosages, (-9, 0, 1, 2)).all():
                                    raise ValueError(f"{key}: fractional dosages are unsupported")
                                n = int(np.count_nonzero(observed))
                                alt_count = int(dosages.sum(where=observed, dtype=np.int64))
                                mac = min(alt_count, 2 * n - alt_count)
                                variant_stats[key] = (n, mac, mac / (2 * n) if n else "")
                            n, mac, maf = variant_stats[key]
                            variant_rows.writerow((key, gene, n, mac, maf))
                            total_mac += mac
                        gene_rows.writerow((
                            chromosome, gene, len(keys), int(np.count_nonzero(nonzero_any)), total_mac,
                            sum(variant_stats[key][1] > 0 for key in keys),
                        ))

    with (output_dir / "dataset_summary.csv").open("w", newline="") as dataset_file:
        csv.writer(dataset_file).writerows((
            ("N_samples_total", "N_variants_total", "N_genes_total", "MAC_total"),
            (len(sample_ids), len(variant_stats), len(genes),
             sum(mac for _, mac, _ in variant_stats.values())),
        ))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("config/aou"))
    parser.add_argument("--output-dir", type=Path, default=Path("metadata"))
    args = parser.parse_args()
    try:
        config = {}
        for name in ("configGlobal.toml", "configPrepare.toml"):
            with (args.config / name).open("rb") as config_file:
                config.update(tomllib.load(config_file))
        args.output_dir.mkdir(parents=True, exist_ok=True)
        # Publish only after every chromosome has been summarized successfully.
        with TemporaryDirectory(prefix=".summary-", dir=args.output_dir) as temporary:
            directory = Path(temporary)
            write_summaries(config, directory)
            for name in ("variant_summary.csv", "gene_summary.csv", "dataset_summary.csv"):
                (directory / name).replace(args.output_dir / name)
    except (OSError, ValueError) as error:
        parser.exit(1, f"Summary failed: {error}\n")
    print(f"Wrote three summary CSVs to {args.output_dir}")


if __name__ == "__main__":
    main()
