"""Run on the AoU host from ~/secure-skat; export only complete result CSVs.

Default: print a plan. --run: reuse validated fixed-local chromosome CSVs,
compute missing tasks, assemble all chromosomes, and archive aggregate CSVs.
"""
import argparse
import csv
import math
import os
import shutil
import sys
import tarfile
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "python/analysis"))
from config import analysis_groups, load_party_config
from run_reference import OUTPUT_COLUMNS, run_chromosome_reference, write_csv


def checked_rows(path, genes, chrom, phenotypes):
    if not path.is_file():
        return None
    try:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames != OUTPUT_COLUMNS:
                return None
            rows = list(reader)
        if len(rows) != len(genes) * len(phenotypes):
            return None
        seen = set()
        for row in rows:
            gi, pi = int(row["gene_index"]), int(row["phenotype_index"])
            if not (0 <= gi < len(genes) and 0 <= pi < len(phenotypes)):
                return None
            if (int(row["chromosome"]) != chrom or row["gene_id"] != genes[gi]
                    or row["phenotype_name"] != phenotypes[pi] or (gi, pi) in seen):
                return None
            seen.add((gi, pi))
            for col in ("r_burden_p", "r_skat_davies_p", "r_skat_liu_p"):
                value = float(row[col])
                if not math.isfinite(value) or not 0 <= value <= 1:
                    return None
        return sorted(rows, key=lambda row: (int(row["gene_index"]), int(row["phenotype_index"])))
    except (ValueError, KeyError, TypeError, csv.Error):
        return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()
    config = load_party_config(REPO / "config/aou-local-fixed-1111", party_id=2)
    run = (REPO / config["run_dir"]).resolve()
    assert run == REPO / "output/aou-local-fixed-1111", run
    assert config["chromosomes"] == list(range(1, 23)), "Expected all 22 chromosomes"
    assert set(config["ancestries"]) == {"AFR", "EUR", "AMR", "ALL"}
    assert set(config["global_max_maf"]) == {0.01, 0.001}
    phenotypes = config["phenotype_columns"]
    destination = run / "local-resumed22"
    sources = [destination, run / "local", run / "local-partial-20261117-201935"]
    tasks, reused, genes_by_task = [], {}, {}
    groups = analysis_groups(config)
    for group in groups:
        ancestry, _, label = group
        for chrom in config["chromosomes"]:
            prepared = run / "prepared/local" / ancestry / f"chr{chrom}"
            assert (prepared / "config.json").is_file(), f"Missing prepared input: {prepared}"
            genes = (prepared / "genes.txt").read_text().splitlines()
            assert genes and len(set(genes)) == len(genes), prepared
            task = (group, chrom)
            genes_by_task[task] = genes
            tasks.append(task)
            for source in sources:
                path = source / label / f"chr{chrom}.csv"
                if checked_rows(path, genes, chrom, phenotypes) is not None:
                    reused[task] = path
                    break
                if path.exists():
                    print(f"Incomplete/incompatible CSV; preserved, not reused: {path}", flush=True)
    missing = [task for task in tasks if task not in reused]
    print(f"Total {len(tasks)}; reusable {len(reused)}; to compute {len(missing)}", flush=True)
    for group in groups:
        print(group[2], "missing chromosomes:", [c for g, c in missing if g == group], flush=True)
    if not args.run:
        print("Plan only. Add --run to execute.")
        return
    if not reused:
        raise SystemExit("No reusable fixed-local CSVs found. Locate the existing partial snapshot before running; no analysis started.")
    workers = min(max(1, int(os.environ.get("REFERENCE_WORKERS", "2"))), os.cpu_count() or 1)
    env = os.environ.copy()
    for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        env[key] = "1"
    destination.mkdir(parents=True, exist_ok=True)
    script = REPO / "python/analysis/python_skat/compute_reference.py"
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")

    def publish(task, source):
        group, chrom = task
        target = destination / group[2] / f"chr{chrom}.csv"
        target.parent.mkdir(parents=True, exist_ok=True)
        if source == target:
            return
        if target.exists():
            target.rename(target.with_name(target.name + f".before-{stamp}"))
        temp = target.with_suffix(".csv.tmp")
        shutil.copyfile(source, temp)
        os.replace(temp, target)

    for task, path in reused.items():
        publish(task, path)

    def compute(task):
        group, chrom = task
        with tempfile.TemporaryDirectory(prefix=".task-", dir=destination) as folder:
            root = Path(folder)
            run_chromosome_reference(run, run, root, "Python", sys.executable,
                                     script, group, chrom, phenotypes, env, local=True)
            path = root / group[2] / f"chr{chrom}.csv"
            assert checked_rows(path, genes_by_task[task], chrom, phenotypes) is not None, path
            publish(task, path)

    print(f"Computing with {workers} workers, 1 BLAS thread each", flush=True)
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(compute, task) for task in missing]
        for count, future in enumerate(as_completed(futures), 1):
            future.result()
            print(f"New tasks completed: {count}/{len(missing)}", flush=True)
    aggregates = []
    for group in groups:
        rows = []
        for chrom in config["chromosomes"]:
            path = destination / group[2] / f"chr{chrom}.csv"
            part = checked_rows(path, genes_by_task[(group, chrom)], chrom, phenotypes)
            assert part is not None, path
            rows.extend(part)
        aggregate = destination / group[2] / "all_local_results.csv"
        temp = aggregate.with_suffix(".csv.tmp")
        write_csv(temp, rows)
        os.replace(temp, aggregate)
        aggregates.append(aggregate)
    archive = run / f"aou-local-full-csv-{datetime.now():%Y%m%d-%H%M%S}.tar.gz"
    temp_archive = archive.with_suffix(".gz.tmp")
    with tarfile.open(temp_archive, "w:gz") as tar:
        for path in aggregates:
            tar.add(path, arcname=str(path.relative_to(destination)))
    with tarfile.open(temp_archive) as tar:
        assert len(tar.getmembers()) == 8 and all(m.isfile() for m in tar.getmembers())
    os.replace(temp_archive, archive)
    print(f"ALL_LOCAL_CSV_COMPLETE\nARCHIVE={archive}", flush=True)


if __name__ == "__main__":
    main()
