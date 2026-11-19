#!/usr/bin/env python3
"""Prepare AFR chr20/22 local inputs, analyze at MAF 1%, and package result tables."""
import argparse
import csv
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import tomllib

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "python/analysis"))
from config import load_party_config


def set_value(text, name, value):
    pattern = rf"(?m)^{name}\s*=\s*" + (r"\[[^\]]*\]" if isinstance(value, list) else r'"[^"\n]*"')
    text, count = re.subn(pattern, lambda _: f"{name} = {json.dumps(value)}", text)
    if count != 1:
        raise ValueError(f"Expected one {name} setting; config was not rewritten")
    return text


def check_scope(path, delimiter, phenotypes):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter=delimiter)
        seen = set()
        for row in reader:
            key = (int(row["chromosome"]), row["gene_id"], row["phenotype_name"])
            index = int(row["phenotype_index"])
            if (key in seen or key[0] not in (20, 22) or not key[1]
                    or not 0 <= index < len(phenotypes) or phenotypes[index] != key[2]):
                raise ValueError("Result scope, phenotype mapping, or unique keys are inconsistent")
            seen.add(key)
    genes = {(chrom, gene) for chrom, gene, _ in seen}
    if {chrom for chrom, _, _ in seen} != {20, 22} or len(seen) != len(genes) * len(phenotypes):
        raise ValueError("Results must cover chr20 and chr22 with every configured phenotype")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("config/aou"))
    parser.add_argument("--secure-run", type=Path, default=Path("output/aou-joint-keyfix-01"))
    args = parser.parse_args()
    os.chdir(REPO)
    os.umask(0o077)
    config = load_party_config(args.config, party_id=2)
    secure = args.secure_run / "secure/AFR/maf0.01/all_secure_results.tsv"
    check_scope(secure, "\t", config["phenotype_columns"])
    if config["samples_per_cohort"] != 0 or config.get("debug_split_num_parties", 0):
        raise ValueError("This pilot requires all eligible samples and real-site inputs")

    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    (REPO / "output").mkdir(exist_ok=True)
    run = Path(tempfile.mkdtemp(prefix=f"aou-local-afr20-22-{stamp}-", dir=REPO / "output"))
    cfg = run / "input-config"
    cfg.mkdir()
    global_text = (args.config / "configGlobal.toml").read_text()
    for name, value in (("chromosomes", [20, 22]), ("ancestries", ["AFR"]), ("global_max_maf", [0.01])):
        global_text = set_value(global_text, name, value)
    (cfg / "configGlobal.toml").write_text(global_text)
    prepare_text = set_value((args.config / "configPrepare.toml").read_text(), "run_dir", str(run))
    (cfg / "configPrepare.toml").write_text(prepare_text)
    local_config = args.config / "configLocal.toml"
    if not local_config.is_file():
        local_config = args.config / "configLocal.Party2.toml"
    shutil.copyfile(local_config, cfg / "configLocal.toml")
    effective = load_party_config(cfg, party_id=2)
    for name, expected in (("chromosomes", [20, 22]), ("ancestries", ["AFR"]),
                           ("global_max_maf", [0.01]), ("run_dir", str(run))):
        if effective[name] != expected:
            raise ValueError(f"Conflicting {name} in the copied configuration")

    bundle = run / "download"
    bundle.mkdir()
    shutil.copyfile(secure, bundle / "all_secure_results.tsv")
    fields = ("chromosomes", "ancestries", "global_max_maf", "phenotype_columns",
              "covariate_columns", "local_max_maf", "samples_per_cohort",
              "ckks", "mpc_num_threads", "data_bits", "fractional_bits", "probes")
    # No source paths, sample records, variant lists, network settings, or auth keys.
    metadata = {
        "local_settings": {key: effective[key] for key in fields if key in effective},
        "local_gene_selection_mode": effective["gene_selection"]["mode"],
        "secure_result_scope": {"ancestry": "AFR", "chromosomes": [20, 22], "global_max_maf": 0.01},
        "secure_scope_source": "ancestry/MAF directory and validated chromosome/phenotype columns",
        "local_maf_rule": "site ALL cohort at global_max_maf; local_max_maf is joint-only",
        "skat_methods": {"secure": "Wilson-Hilferty", "local": "Davies with recorded fallback/flags"},
        "source_sha256": {
            name: hashlib.sha256((REPO / name).read_bytes()).hexdigest()
            for name in ("python/preprocessing/plink.py", "python/preprocessing/prepare.py",
                         "python/analysis/run_reference.py")
        },
    }
    snapshot = args.secure_run / "config/configGlobal.toml"
    if snapshot.is_file():
        recorded = tomllib.loads(snapshot.read_text())
        metadata["secure_saved_config_not_runtime_verified"] = {
            key: recorded[key] for key in fields if key in recorded
        }

    print(f"Local run: {run}", flush=True)
    subprocess.run(["go", "run", "-mod=vendor", "secure-rvas.go", "prepare",
                    "--config", str(cfg), "--party", "2", "--mode", "local"], check=True)
    env = os.environ.copy()
    env.setdefault("REFERENCE_WORKERS", "1")
    env.setdefault("REFERENCE_BLAS_THREADS", "1")
    subprocess.run([sys.executable, "-u", "python/analysis/run_reference.py",
                    "--config", str(cfg), "--mode", "local"], env=env, check=True)
    local = run / "local/AFR/maf0.01/all_local_results.csv"
    check_scope(local, ",", effective["phenotype_columns"])
    shutil.copyfile(local, bundle / "all_local_results.csv")

    with (bundle / "gene_labels.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=("chromosome", "gene_id", "gene_symbol"),
                                delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for chromosome in (20, 22):
            panel = Path(effective["gene_panel"].replace("{chromosome}", str(chromosome)))
            with panel.open(newline="") as source:
                writer.writerows(csv.DictReader(source, delimiter="\t"))
    names = ("all_secure_results.tsv", "all_local_results.csv", "gene_labels.tsv")
    metadata["result_sha256"] = {
        name: hashlib.sha256((bundle / name).read_bytes()).hexdigest() for name in names
    }
    (bundle / "manifest.json").write_text(json.dumps(metadata, indent=2) + "\n")
    archive = run / "afr-chr20-22-maf0.01-results.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        for name in (*names, "manifest.json"):
            tar.add(bundle / name, arcname=name, recursive=False)
    print(f"LOCAL_PILOT_COMPLETED\nDownload: {archive}", flush=True)


if __name__ == "__main__":
    main()
