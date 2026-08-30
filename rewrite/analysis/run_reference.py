#!/usr/bin/env python3

import argparse
import csv
import os
import subprocess
import tomllib
from concurrent.futures import ThreadPoolExecutor
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

TARGET_BLAS_THREADS_PER_PROCESS = 8

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


def run_chromosome_reference(
    run_dir: Path,
    r_script: Path,
    ancestry: str,
    chromosome: int,
    phenotype_names: list[str],
    environment: dict[str, str],
) -> list[dict[str, str]]:
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
        env=environment,
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

    reference_dir = run_dir / "reference" / ancestry
    reference_dir.mkdir(parents=True, exist_ok=True)
    write_csv(
        reference_dir / f"chr{chromosome}.csv",
        chromosome_rows,
    )
    return chromosome_rows


def run_reference(config_path: Path) -> None:
    with config_path.open("rb") as config_file:
        config = tomllib.load(config_file)

    run_dir = Path(config["run_dir"])
    reference_root = run_dir / "reference"
    reference_root.mkdir(parents=True, exist_ok=True)

    success_path = reference_root / "_SUCCESS"
    success_path.unlink(missing_ok=True)

    r_script = Path(__file__).resolve().with_name("main.R")
    tasks = [
        (ancestry, chromosome)
        for ancestry in config["ancestries"]
        for chromosome in config["chromosomes"]
    ]

    cpu_count = os.cpu_count() or 1
    automatic_workers = max(
        1,
        cpu_count // TARGET_BLAS_THREADS_PER_PROCESS,
    )
    requested_workers = int(
        os.environ.get("R_REFERENCE_WORKERS", automatic_workers)
    )
    if requested_workers < 1:
        raise ValueError("R_REFERENCE_WORKERS must be positive")

    worker_count = min(requested_workers, len(tasks), cpu_count)
    blas_threads = max(1, cpu_count // worker_count)
    r_environment = os.environ.copy()
    r_environment["OPENBLAS_NUM_THREADS"] = str(blas_threads)

    print(
        f"Running {len(tasks)} R::SKAT tasks with "
        f"{worker_count} workers and {blas_threads} BLAS threads per worker"
    )

    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        futures = [
            executor.submit(
                run_chromosome_reference,
                run_dir,
                r_script,
                ancestry,
                chromosome,
                config["phenotype_columns"],
                r_environment,
            )
            for ancestry, chromosome in tasks
        ]
        task_rows = [future.result() for future in futures]

    rows_by_ancestry = {
        ancestry: [] for ancestry in config["ancestries"]
    }
    for (ancestry, _), chromosome_rows in zip(tasks, task_rows):
        rows_by_ancestry[ancestry].extend(chromosome_rows)

    for ancestry, ancestry_rows in rows_by_ancestry.items():
        output_path = reference_root / ancestry / "all_r_results.csv"
        write_csv(output_path, ancestry_rows)
        print("Wrote R reference results to", output_path)

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
