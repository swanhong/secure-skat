#!/usr/bin/env python3

import argparse
import csv
import subprocess
import tomllib
from pathlib import Path


R_COLUMNS = [
    "gene_index",
    "gene_id",
    "phenotype_index",
    "burden_p",
    "skat_davies_p",
    "skat_davies_converged",
    "skat_liu_p",
]

OUTPUT_COLUMNS = [
    "chromosome",
    "gene_index",
    "gene_id",
    "phenotype_index",
    "phenotype_name",
    "r_burden_p",
    "r_skat_davies_p",
    "r_skat_davies_converged",
    "r_skat_liu_p",
]


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(
            output,
            fieldnames=OUTPUT_COLUMNS,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def run_ancestry_reference(
    run_dir: Path,
    r_script: Path,
    ancestry: str,
    chromosomes: list[int],
    phenotype_names: list[str],
) -> None:
    reference_dir = run_dir / "reference" / ancestry
    reference_dir.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict[str, str]] = []

    for chromosome in chromosomes:
        print(f"Running R::SKAT for {ancestry} chromosome {chromosome}")

        completed = subprocess.run(
            [
                "Rscript",
                str(r_script),
                str(
                    run_dir
                    / "prepared"
                    / ancestry
                    / f"chr{chromosome}"
                ),
            ],
            check=True,
            stdout=subprocess.PIPE,
            text=True,
        )

        reader = csv.DictReader(
            completed.stdout.splitlines(),
            delimiter="\t",
        )
        if reader.fieldnames != R_COLUMNS:
            raise ValueError(
                f"unexpected R output columns: {reader.fieldnames}"
            )

        chromosome_rows: list[dict[str, str]] = []
        for row in reader:
            phenotype_index = int(row["phenotype_index"])
            if not 0 <= phenotype_index < len(phenotype_names):
                raise ValueError(
                    f"invalid phenotype index {phenotype_index}"
                )

            chromosome_rows.append(
                {
                    "chromosome": str(chromosome),
                    "gene_index": row["gene_index"],
                    "gene_id": row["gene_id"],
                    "phenotype_index": row["phenotype_index"],
                    "phenotype_name": phenotype_names[phenotype_index],
                    "r_burden_p": row["burden_p"],
                    "r_skat_davies_p": row["skat_davies_p"],
                    "r_skat_davies_converged":
                        row["skat_davies_converged"],
                    "r_skat_liu_p": row["skat_liu_p"],
                }
            )

        write_csv(
            reference_dir / f"chr{chromosome}.csv",
            chromosome_rows,
        )
        all_rows.extend(chromosome_rows)

    output_path = reference_dir / "all_r_results.csv"
    write_csv(output_path, all_rows)
    print("Wrote R reference results to", output_path)


def run_reference(config_path: Path) -> None:
    with config_path.open("rb") as config_file:
        config = tomllib.load(config_file)

    run_dir = Path(config["run_dir"])
    reference_root = run_dir / "reference"
    reference_root.mkdir(parents=True, exist_ok=True)

    success_path = reference_root / "_SUCCESS"
    success_path.unlink(missing_ok=True)

    r_script = Path(__file__).resolve().with_name("main.R")
    for ancestry in config["ancestries"]:
        run_ancestry_reference(
            run_dir,
            r_script,
            ancestry,
            config["chromosomes"],
            config["phenotype_columns"],
        )

    success_path.touch()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        default="run.1kg.conf",
        help="path to the run configuration",
    )
    args = parser.parse_args()

    run_reference(Path(args.config))


if __name__ == "__main__":
    main()
