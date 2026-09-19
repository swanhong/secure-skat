#!/usr/bin/env python3
"""Run on the Hail cluster. Export both tests for all five phenotypes."""

import argparse
import subprocess

from common import settings


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", help="TSV path, including gs://; default: rvas-shared")
    args = parser.parse_args()
    import hail as hl

    config = settings()
    output = args.output
    if output is None:
        bucket = subprocess.check_output(
            ["wb", "resource", "resolve", "--id=rvas-shared"], text=True
        ).strip().rstrip("/")
        if not bucket.startswith("gs://"):
            raise ValueError(f"Expected a GCS bucket URI, got {bucket!r}")
        output = f"{bucket}/allxall-comparison/axa_test/allxall_v8.tsv"

    hl.init()
    try:
        tables = []
        for phenotype in config["phenotypes"]:
            uri = (f"{config['allxall_root']}/{config['ancestry']}/"
                   f"phenotype_{phenotype['axa_id']}_gene_results.ht")
            print(f"Reading {phenotype['name']}: {uri}", flush=True)
            ht = hl.read_table(uri)
            ht = ht.filter((ht.annotation == config["annotation"]) &
                           (ht.max_MAF == config["max_maf"]))
            if ht.take(1) == []:
                raise ValueError(f"No matching rows: {uri}")
            ht = ht.key_by().select(
                "gene_id", "gene_symbol", "annotation", "max_MAF",
                phenotype=phenotype["name"], axa_phenotype_id=phenotype["axa_id"],
                ancestry=config["ancestry"], release="v8", source_uri=uri,
                axa_burden=ht.Pvalue_Burden, axa_skat=ht.Pvalue_SKAT,
            ).select_globals()
            tables.append(ht)
        tables[0].union(*tables[1:]).export(output)
        print(f"Saved: {output}")
    finally:
        hl.stop()


if __name__ == "__main__":
    main()
