#!/usr/bin/env python3
"""Run R Burden and SKAT on existing v9 prepared inputs (pilot by default)."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

from common import DEFAULT_OUTPUT, HERE, gene_id, index_rows, read_rows, settings, write_rows


def loader_path(repo):
    for relative in ("rewrite/analysis/r_skat/load_preprocessed.R",
                     "python/analysis/r_skat/load_preprocessed.R"):
        path = repo / relative
        if path.is_file():
            return path
    raise FileNotFoundError(f"Cannot find the existing R prepared-input loader under {repo}")


def plan_inputs(run_dir, config, mode):
    root = run_dir / "prepared" / config["ancestry"]
    targets = {g["gene_id"]: g for g in config["pilot_genes"]}
    needed = set(range(1, 23)) if mode == "all" else {g["chromosome"] for g in targets.values()}
    found_chromosomes, found_genes, plans = set(), set(), []
    for directory in sorted(root.glob("chr*"), key=lambda p: int(p.name[3:])):
        chromosome = int(directory.name[3:])
        if chromosome not in needed:
            continue
        if not (directory / "pos.txt").is_file():
            raise ValueError(f"Incomplete prepared inputs: {directory}")
        saved = json.loads((directory / "config.json").read_text())
        expected = {"chromosome": chromosome, "ancestry_group": config["ancestry"],
                    "max_maf": config["max_maf"], "shared_rate": 1.0,
                    "mask": {"annotation": [config["annotation"]]}}
        for key, value in expected.items():
            if saved.get(key) != value:
                raise ValueError(f"{directory}: expected {key}={value!r}, got {saved.get(key)!r}")
        columns = saved["phenotype_columns"]
        if len(columns) != len(set(columns)):
            raise ValueError(f"Duplicate phenotype columns: {directory}")
        phenotypes = [{**p, "index": columns.index(p["column"]) + 1}
                      for p in config["phenotypes"]]
        original_genes = (directory / "genes.txt").read_text().splitlines()
        normalized = [gene_id(g) for g in original_genes]
        if len(normalized) != len(set(normalized)) or "" in normalized:
            raise ValueError(f"Duplicate/empty normalized gene IDs: {directory}")
        genes = []
        for index, (original, stable) in enumerate(zip(original_genes, normalized), 1):
            if mode == "pilot" and stable not in targets:
                continue
            if stable in targets and targets[stable]["chromosome"] != chromosome:
                raise ValueError(f"Unexpected chromosome for {stable}")
            if stable in found_genes:
                raise ValueError(f"Gene occurs in multiple chromosomes: {stable}")
            found_genes.add(stable)
            genes.append({"index": index, "prepared_gene_id": original, "gene_id": stable,
                          "gene_symbol": targets.get(stable, {}).get("symbol", "")})
        found_chromosomes.add(chromosome)
        if not genes:
            raise ValueError(f"No selected genes: {directory}")
        plans.append({"chromosome": chromosome, "prepared": str(directory.resolve()),
                      "saved_config": saved, "genes": genes, "phenotypes": phenotypes})
    if found_chromosomes != needed:
        raise ValueError(f"Missing prepared chromosomes: {sorted(needed - found_chromosomes)}")
    if mode == "pilot" and found_genes != targets.keys():
        raise ValueError(f"Missing pilot genes: {sorted(targets.keys() - found_genes)}")
    return plans


def validate_rows(rows, plan):
    expected = {(g["gene_id"], p["axa_id"]) for g in plan["genes"] for p in plan["phenotypes"]}
    actual = index_rows(rows, lambda row: (row["gene_id"], row["axa_phenotype_id"]))
    if actual.keys() != expected:
        raise ValueError(f"Incomplete or unexpected R results for chr{plan['chromosome']}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.home() / "secure-skat")
    parser.add_argument("--run-dir", type=Path,
                        help="Original run directory containing prepared/EUR; otherwise use old run_info.json")
    parser.add_argument("--plain-directory", type=Path, default=Path.home() / "allxall-comparison")
    parser.add_argument("--mode", choices=("pilot", "all"), default="pilot")
    parser.add_argument("--out", type=Path, help="Default: ~/allxall-comparison/axa_test/<mode>")
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    config = settings()
    repo = args.repo.expanduser().resolve()
    run_dir = args.run_dir
    if run_dir is None:
        run_dir = Path(json.loads((args.plain_directory.expanduser() / "run_info.json").read_text())["run_dir"])
    run_dir = run_dir.expanduser()
    if not run_dir.is_absolute():
        run_dir = repo / run_dir
    loader = loader_path(repo)
    plans = plan_inputs(run_dir, config, args.mode)
    out = (args.out or DEFAULT_OUTPUT / args.mode).expanduser().resolve()
    row_count = sum(len(p["genes"]) * len(p["phenotypes"]) for p in plans)
    print(f"Mode: {args.mode}; chromosomes: {len(plans)}; gene-phenotype rows: {row_count}")
    print(f"Output: {out}")
    print("Prepared covariates are preserved; this does not harmonize with All-by-All.")
    if args.dry_run:
        for p in plans:
            s = p["saved_config"]
            print(f"chr{p['chromosome']}: {len(p['genes'])} genes; num_cov={s.get('num_cov')}; "
                  f"covariate={s.get('covariate')}")
        return

    versions = subprocess.check_output(
        [args.rscript, "-e", 'cat(R.version.string, "\\n", as.character(packageVersion("SKAT")))'],
        text=True,
    ).strip()
    manifest = {"settings": config, "mode": args.mode, "run_dir": str(run_dir.resolve()),
                "R": versions, "plans": plans,
                "code_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in (loader, HERE / "run_genes.R", HERE / "run_r.py", HERE / "common.py")}}
    out.mkdir(parents=True, exist_ok=True)
    manifest_path = out / "run_info.json"
    if manifest_path.exists() and json.loads(manifest_path.read_text()) != manifest:
        raise ValueError("Inputs/code changed. Choose a new --out directory to preserve the previous run.")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    (out / "_SUCCESS").unlink(missing_ok=True)
    environment = {**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1",
                   "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1"}
    combined = []
    for plan in plans:
        chromosome = str(plan["chromosome"])
        folder = out / f"chr{chromosome}"
        folder.mkdir(exist_ok=True)
        result = folder / "r_results.csv"
        if result.exists():
            print(f"Reusing completed chr{chromosome}", flush=True)
            rows = read_rows(result)
            validate_rows(rows, plan)
        else:
            for name in ("genes", "phenotypes"):
                write_rows(folder / f"{name}.tsv", plan[name], delimiter="\t")
            temporary = folder / "r_results.csv.tmp"
            subprocess.run([args.rscript, str(HERE / "run_genes.R"), str(loader),
                            plan["prepared"], str(folder / "genes.tsv"), str(folder / "phenotypes.tsv"), chromosome,
                            str(temporary)], check=True, env=environment)
            rows = read_rows(temporary)
            validate_rows(rows, plan)
            temporary.replace(result)
        combined.extend(rows)
    write_rows(out / "r_v9.csv", combined)
    (out / "_SUCCESS").touch()
    flagged = sum(r["R_burden_status"] != "ok" or r["R_skat_status"] != "ok" for r in combined)
    print(f"Saved {len(combined)} rows: {out / 'r_v9.csv'}; diagnostic rows: {flagged}")


if __name__ == "__main__":
    main()
