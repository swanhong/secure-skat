#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path

from plink import run_vcf_to_pgen


def default_1000g_vcf_name(chromosome: int) -> str:
    return f"ALL.chr{chromosome}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"


def default_1000g_panel_name() -> str:
    return "integrated_call_samples_v3.20130502.ALL.panel"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare reusable AoU-style 1000 Genomes source inputs: "
            "PGEN trio plus split phenotype/covariate CSV files."
        )
    )
    parser.add_argument(
        "--source-dir",
        required=True,
        help="Reusable 1000G source root. Defaults read from raw/ and outputs to pgen/ and pheno/ under this root.",
    )
    parser.add_argument("--chromosome", type=int, default=22, help="Chromosome to convert (default: 22).")
    parser.add_argument(
        "--vcf",
        help=(
            "Input 1000G VCF. If omitted, uses raw/ALL.chr<chromosome>...vcf.gz under --source-dir "
            "when present."
        ),
    )
    parser.add_argument(
        "--panel-file",
        help="1000G integrated_call_samples panel file. If omitted, uses raw/integrated_call_samples...panel.",
    )
    parser.add_argument(
        "--sample-file",
        help="Optional sample list/keep file for pheno/cov generation. Uses all panel samples if omitted.",
    )
    parser.add_argument(
        "--pgen-prefix-out",
        help="Output PGEN prefix without extension (default: <source-dir>/pgen/1000g.chr<chromosome>).",
    )
    parser.add_argument(
        "--pheno-out",
        help="Output phenotype CSV path (default: <source-dir>/pheno/phenotype_data.csv).",
    )
    parser.add_argument(
        "--cov-out",
        help="Output covariate CSV path (default: <source-dir>/pheno/covariate_data.csv).",
    )
    parser.add_argument("--seed", default="270318_pheno", help="Deterministic phenotype seed.")
    parser.add_argument(
        "--cov-mode",
        choices=("superpop", "afr-pop"),
        default="superpop",
        help="Covariates to emit. superpop is for all 1000G; afr-pop matches the old AFR subset.",
    )
    parser.add_argument("--plink2", default="plink2", help="Path to plink2 executable.")
    parser.add_argument(
        "--new-id-max-allele-len",
        type=int,
        default=1000,
        help="PLINK2 --new-id-max-allele-len value used with --set-all-var-ids.",
    )
    parser.add_argument(
        "--vcf-no-double-id",
        action="store_true",
        help="Do not pass PLINK2 --double-id during VCF import.",
    )
    parser.add_argument(
        "--vcf-keep",
        help="Optional PLINK --keep file applied during VCF-to-PGEN conversion.",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite existing PGEN output files.")

    args = parser.parse_args()
    if args.chromosome < 1 or args.chromosome > 22:
        parser.error("--chromosome must be in 1..22")
    if args.new_id_max_allele_len <= 0:
        parser.error("--new-id-max-allele-len must be positive")
    return args


def resolve_existing_or_default(explicit_path: str | None, candidates: list[Path], label: str) -> Path:
    if explicit_path:
        path = Path(explicit_path).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(f"{label} does not exist: {path}")
        return path

    for path in candidates:
        if path.exists():
            return path.resolve()

    expected = ", ".join(str(path) for path in candidates)
    raise FileNotFoundError(f"Could not find {label}. Expected one of: {expected}")


def read_panel(path: Path) -> dict[str, dict[str, str]]:
    panel = {}
    with path.open() as fh:
        next(fh)
        for line in fh:
            if not line.strip():
                continue
            sample, pop, super_pop, gender, *_ = line.rstrip("\n").split("\t")
            panel[sample] = {"pop": pop, "super_pop": super_pop, "gender": gender}
    return panel


def read_sample_ids(path: Path | None, panel: dict[str, dict[str, str]]) -> list[str]:
    if path is None:
        return list(panel)

    sample_ids = []
    with path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            toks = line.split()
            sample_ids.append(toks[-1])
    return sample_ids


def hash_uniform(sample_id: str, seed: str) -> float:
    digest = hashlib.sha256(f"{sample_id}|{seed}".encode()).digest()
    return int.from_bytes(digest[:8], "big") / float(2**64)


def standardize(values: list[float]) -> list[float]:
    mean = sum(values) / len(values)
    var = sum((value - mean) ** 2 for value in values) / (len(values) - 1)
    sd = math.sqrt(var)
    return [(value - mean) / sd for value in values]


def cov_header(cov_mode: str) -> list[str]:
    if cov_mode == "superpop":
        return ["superpop_AFR", "superpop_AMR", "superpop_EAS", "superpop_EUR"]
    return ["sex_male", "pop_GWD", "pop_YRI", "pop_ESN"]


def cov_values(meta: dict[str, str], cov_mode: str) -> list[float]:
    if cov_mode == "superpop":
        return [1.0 if meta["super_pop"] == super_pop else 0.0 for super_pop in ("AFR", "AMR", "EAS", "EUR")]
    return [
        1.0 if meta["gender"] == "male" else 0.0,
        1.0 if meta["pop"] == "GWD" else 0.0,
        1.0 if meta["pop"] == "YRI" else 0.0,
        1.0 if meta["pop"] == "ESN" else 0.0,
    ]


def build_synthetic_rows(
    *,
    panel: dict[str, dict[str, str]],
    sample_ids: list[str],
    seed: str,
    cov_mode: str,
) -> tuple[list[str], list[tuple[str, float, list[float]]]]:
    missing = [sample_id for sample_id in sample_ids if sample_id not in panel]
    if missing:
        raise ValueError(f"{len(missing)} samples are missing from the panel; first missing sample: {missing[0]}")
    if len(sample_ids) < 2:
        raise ValueError("Need at least two samples to standardize the synthetic phenotype")

    pheno = standardize([hash_uniform(sample_id, seed) for sample_id in sample_ids])
    headers = cov_header(cov_mode)
    rows = []
    for sample_id, y in zip(sample_ids, pheno):
        rows.append((sample_id, y, cov_values(panel[sample_id], cov_mode)))
    return headers, rows


def write_split_tables(
    pheno_out: Path,
    cov_out: Path,
    headers: list[str],
    rows: list[tuple[str, float, list[float]]],
) -> None:
    pheno_out.parent.mkdir(parents=True, exist_ok=True)
    cov_out.parent.mkdir(parents=True, exist_ok=True)
    with pheno_out.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["person_id", "pheno"])
        for sample_id, y, _ in rows:
            writer.writerow([sample_id, f"{y:.12f}"])
    with cov_out.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["person_id", *headers])
        for sample_id, _, cov in rows:
            writer.writerow([sample_id, *[f"{value:.1f}" for value in cov]])
    print(f"Wrote split phenotype table for {len(rows)} samples to {pheno_out}")
    print(f"Wrote split covariate table for {len(rows)} samples to {cov_out}")


def main() -> int:
    args = parse_args()
    source_dir = Path(args.source_dir).expanduser().resolve()
    raw_dir = source_dir / "raw"
    pgen_prefix = (
        Path(args.pgen_prefix_out).expanduser().resolve()
        if args.pgen_prefix_out
        else source_dir / "pgen" / f"1000g.chr{args.chromosome}"
    )
    pheno_out = (
        Path(args.pheno_out).expanduser().resolve()
        if args.pheno_out
        else source_dir / "pheno" / "phenotype_data.csv"
    )
    cov_out = (
        Path(args.cov_out).expanduser().resolve()
        if args.cov_out
        else source_dir / "pheno" / "covariate_data.csv"
    )

    vcf_path = resolve_existing_or_default(
        args.vcf,
        [
            raw_dir / default_1000g_vcf_name(args.chromosome),
            source_dir / default_1000g_vcf_name(args.chromosome),
        ],
        "1000G VCF",
    )
    panel_file = resolve_existing_or_default(
        args.panel_file,
        [
            raw_dir / default_1000g_panel_name(),
            source_dir / default_1000g_panel_name(),
        ],
        "1000G panel file",
    )
    sample_file = Path(args.sample_file).expanduser().resolve() if args.sample_file else None
    keep_path = Path(args.vcf_keep).expanduser().resolve() if args.vcf_keep else None
    if keep_path and not keep_path.exists():
        raise FileNotFoundError(f"VCF keep file does not exist: {keep_path}")

    print(f"Source directory: {source_dir}")
    print(f"Input VCF:        {vcf_path}")
    print(f"Panel file:       {panel_file}")
    print(f"PGEN prefix out:  {pgen_prefix}")
    print(f"Phenotype out:    {pheno_out}")
    print(f"Covariate out:    {cov_out}")

    run_vcf_to_pgen(
        vcf_path=vcf_path,
        out_prefix=pgen_prefix,
        chromosome=args.chromosome,
        plink2=args.plink2,
        new_id_max_allele_len=args.new_id_max_allele_len,
        vcf_no_double_id=args.vcf_no_double_id,
        keep_path=keep_path,
        force=args.force,
    )

    panel = read_panel(panel_file)
    sample_ids = read_sample_ids(sample_file, panel)
    headers, rows = build_synthetic_rows(
        panel=panel,
        sample_ids=sample_ids,
        seed=args.seed,
        cov_mode=args.cov_mode,
    )
    write_split_tables(pheno_out, cov_out, headers, rows)

    print("\n1000G source inputs are ready.")
    print(f"PGEN prefix: {pgen_prefix}")
    print(f"Phenotype:   {pheno_out}")
    print(f"Covariates:  {cov_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
