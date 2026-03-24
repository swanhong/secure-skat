#!/usr/bin/env python3
import argparse
import json
import math
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np


@dataclass
class PartyBlockExport:
    party_name: str
    raw_path: Path
    psam_path: Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Simulate PRS-CS style phenotypes from a dataset root containing party1/party2 PGEN files."
    )
    parser.add_argument("--dataset-root", required=True, help="Dataset root containing party1/, party2/, ...")
    parser.add_argument("--num-causal-snps", required=True, type=int, help="Expected number of causal SNPs (K)")
    parser.add_argument("--h2", required=True, type=float, help="Target SNP heritability in [0, 1)")
    parser.add_argument("--seed", required=True, type=int, help="Random seed")
    parser.add_argument("--plink2", default="plink2", help="Path to plink2 executable")
    parser.add_argument(
        "--keep-cache",
        action="store_true",
        help="Keep exported .raw files under the simulation run directory",
    )
    return parser.parse_args()


def discover_party_dirs(dataset_root: Path):
    party_dirs = sorted(
        [p for p in dataset_root.iterdir() if p.is_dir() and p.name.startswith("party")],
        key=lambda p: p.name,
    )
    if not party_dirs:
        raise FileNotFoundError(f"No party directories found under {dataset_root}")
    return party_dirs


def detect_pfile_modifier(pfile_prefix: Path):
    if pfile_prefix.with_suffix(".pvar.zst").exists():
        return "vzs"
    return None


def read_block_sizes(party_dir: Path):
    chrom_sizes = []
    with (party_dir / "chrom_sizes.txt").open() as fh:
        for line in fh:
            line = line.strip()
            if line:
                chrom_sizes.append(int(line))
    if not chrom_sizes:
        raise ValueError(f"No block sizes found in {party_dir / 'chrom_sizes.txt'}")
    return chrom_sizes


def export_block_raw(plink2: str, party_dir: Path, block_idx: int, cache_dir: Path):
    geno_prefix = party_dir / "geno" / f"chr{block_idx}"
    if not geno_prefix.with_suffix(".pgen").exists():
        raise FileNotFoundError(f"Missing {geno_prefix}.pgen")
    sample_keep = party_dir / "sample_keep.txt"
    out_prefix = cache_dir / f"{party_dir.name}_chr{block_idx:02d}"
    raw_path = out_prefix.with_suffix(".raw")
    psam_path = geno_prefix.with_suffix(".psam")
    if raw_path.exists():
        return PartyBlockExport(party_name=party_dir.name, raw_path=raw_path, psam_path=psam_path)

    cmd = [plink2, "--pfile", str(geno_prefix)]
    modifier = detect_pfile_modifier(geno_prefix)
    if modifier:
        cmd.append(modifier)
    cmd.extend(["--keep", str(sample_keep), "--export", "A", "--out", str(out_prefix)])
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if not raw_path.exists():
        raise FileNotFoundError(f"Expected export output not found: {raw_path}")
    return PartyBlockExport(party_name=party_dir.name, raw_path=raw_path, psam_path=psam_path)


def read_psam_iids(psam_path: Path):
    iids = []
    with psam_path.open() as fh:
        header = next(fh)
        delim = "\t" if "\t" in header else None
        for line in fh:
            line = line.strip()
            if not line:
                continue
            toks = line.split(delim) if delim else line.split()
            if len(toks) < 2:
                raise ValueError(f"Unexpected PSAM row: {line}")
            iids.append(toks[1])
    return iids


def load_raw_matrix(raw_path: Path):
    with raw_path.open() as fh:
        header = fh.readline().strip().split("\t")
    if len(header) <= 6:
        raise ValueError(f"No genotype columns found in {raw_path}")

    variant_ids = np.array(header[6:], dtype=object)
    sample_ids = np.genfromtxt(
        raw_path,
        delimiter="\t",
        dtype=str,
        skip_header=1,
        usecols=(1,),
        encoding="utf-8",
    )
    geno = np.genfromtxt(
        raw_path,
        delimiter="\t",
        dtype=np.float64,
        skip_header=1,
        usecols=tuple(range(6, len(header))),
        missing_values="NA",
        filling_values=np.nan,
    )
    if geno.ndim == 1:
        geno = geno.reshape(1, -1)
    if sample_ids.ndim == 0:
        sample_ids = np.array([sample_ids.item()], dtype=object)
    return sample_ids.astype(object), variant_ids, geno


def first_pass_block(export_infos, expected_num_variants):
    block_variant_ids = None
    sample_ids_by_party = {}
    block_sum = None
    block_sumsq = None
    block_count = None

    for export in export_infos:
        raw_sample_ids, variant_ids, geno = load_raw_matrix(export.raw_path)
        psam_iids = np.array(read_psam_iids(export.psam_path), dtype=object)
        if raw_sample_ids.shape[0] != psam_iids.shape[0] or not np.array_equal(raw_sample_ids, psam_iids):
            raise ValueError(f"Sample order mismatch between {export.raw_path} and {export.psam_path}")

        if variant_ids.shape[0] != expected_num_variants:
            raise ValueError(
                f"Variant count mismatch in {export.raw_path}: expected {expected_num_variants}, got {variant_ids.shape[0]}"
            )
        if block_variant_ids is None:
            block_variant_ids = variant_ids
            block_sum = np.zeros(expected_num_variants, dtype=np.float64)
            block_sumsq = np.zeros(expected_num_variants, dtype=np.float64)
            block_count = np.zeros(expected_num_variants, dtype=np.int64)
        elif not np.array_equal(block_variant_ids, variant_ids):
            raise ValueError(f"Variant ID mismatch between party exports for {export.raw_path.name}")

        observed = ~np.isnan(geno)
        geno_zeroed = np.where(observed, geno, 0.0)
        block_sum += geno_zeroed.sum(axis=0)
        block_sumsq += np.square(geno_zeroed).sum(axis=0)
        block_count += observed.sum(axis=0)
        sample_ids_by_party[export.party_name] = raw_sample_ids.tolist()

    means = block_sum / np.maximum(block_count, 1)
    variances = block_sumsq / np.maximum(block_count, 1) - np.square(means)
    variances = np.where(variances < 0.0, 0.0, variances)
    stds = np.sqrt(variances)

    valid_mask = (block_count > 0) & (stds > 0.0)
    return block_variant_ids, sample_ids_by_party, means, stds, valid_mask


def second_pass_block(export_infos, means, stds, beta_block):
    y_parts = {}
    valid_stds = np.where(stds > 0.0, stds, 1.0)
    for export in export_infos:
        _raw_sample_ids, _variant_ids, geno = load_raw_matrix(export.raw_path)
        geno = np.where(np.isnan(geno), means, geno)
        z_block = (geno - means) / valid_stds
        if not np.isfinite(z_block).all():
            raise ValueError(f"Non-finite standardized genotype values in {export.raw_path}")
        y_parts[export.party_name] = np.einsum("ij,j->i", z_block, beta_block, dtype=np.float64)
    return y_parts


def main():
    args = parse_args()
    if not (0.0 <= args.h2 < 1.0):
        raise ValueError("--h2 must be in [0, 1)")
    if args.num_causal_snps <= 0:
        raise ValueError("--num-causal-snps must be positive")

    dataset_root = Path(args.dataset_root).resolve()
    party_dirs = discover_party_dirs(dataset_root)
    block_sizes = read_block_sizes(party_dirs[0])
    for party_dir in party_dirs[1:]:
        other_sizes = read_block_sizes(party_dir)
        if other_sizes != block_sizes:
            raise ValueError(f"Block sizes mismatch between parties: {party_dirs[0]} vs {party_dir}")

    total_variants = int(sum(block_sizes))
    if args.num_causal_snps > total_variants:
        raise ValueError("--num-causal-snps exceeds total variant count")

    timestamp = datetime.now().strftime("%y%m%d_%H%M%S")
    run_dir = dataset_root / "simulation" / f"prscs_{timestamp}_seed{args.seed}_k{args.num_causal_snps}_h2_{str(args.h2).replace('.', 'p')}"
    run_dir.mkdir(parents=True, exist_ok=False)
    cache_dir = run_dir / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    pi = args.num_causal_snps / total_variants
    beta_sd = math.sqrt(args.h2 / args.num_causal_snps)
    sigma_e = math.sqrt(1.0 - args.h2)

    sample_ids_global = []
    party_sample_ids = {}
    party_offsets = {}
    all_y_genetic = []
    all_variant_ids = []
    beta_blocks = []
    block_summaries = []

    running_sample_offset = 0
    for block_idx, expected_num_variants in enumerate(block_sizes, start=1):
        export_infos = [export_block_raw(args.plink2, party_dir, block_idx, cache_dir) for party_dir in party_dirs]
        variant_ids, block_sample_ids, means, stds, valid_mask = first_pass_block(export_infos, expected_num_variants)

        if not sample_ids_global:
            for party_dir in party_dirs:
                party_name = party_dir.name
                ids = block_sample_ids[party_name]
                party_sample_ids[party_name] = ids
                party_offsets[party_name] = (running_sample_offset, running_sample_offset + len(ids))
                running_sample_offset += len(ids)
                sample_ids_global.extend(ids)
                all_y_genetic.append(np.zeros(len(ids), dtype=np.float64))
        else:
            for party_dir in party_dirs:
                party_name = party_dir.name
                if block_sample_ids[party_name] != party_sample_ids[party_name]:
                    raise ValueError(f"Sample order changed across blocks for {party_name}")

        causal_draw = rng.random(expected_num_variants) <= pi
        causal_draw &= valid_mask
        beta_block = np.zeros(expected_num_variants, dtype=np.float64)
        if causal_draw.any():
            beta_block[causal_draw] = rng.normal(0.0, beta_sd, size=int(causal_draw.sum()))

        y_parts = second_pass_block(export_infos, means, stds, beta_block)
        for party_dir in party_dirs:
            party_name = party_dir.name
            offset_start, _offset_end = party_offsets[party_name]
            party_idx = list(party_offsets).index(party_name)
            all_y_genetic[party_idx] += y_parts[party_name]

        all_variant_ids.extend(variant_ids.tolist())
        beta_blocks.append(beta_block)
        block_summaries.append(
            {
                "block_index": block_idx,
                "num_variants": expected_num_variants,
                "num_valid_variants": int(valid_mask.sum()),
                "num_causal_variants": int(causal_draw.sum()),
                "mean_abs_genotype_mean": float(np.mean(np.abs(means[valid_mask]))) if valid_mask.any() else 0.0,
                "max_abs_genotype_mean": float(np.max(np.abs(means[valid_mask]))) if valid_mask.any() else 0.0,
            }
        )

    y_genetic = np.concatenate(all_y_genetic)
    epsilon = rng.normal(0.0, sigma_e, size=y_genetic.shape[0])
    y = y_genetic + epsilon
    if not np.isfinite(y).all():
        raise ValueError("Simulated phenotype contains non-finite values")

    beta_all = np.concatenate(beta_blocks)
    causal_mask_all = beta_all != 0.0
    realized_num_causal = int(causal_mask_all.sum())

    metadata = {
        "timestamp_local": timestamp,
        "dataset_root": str(dataset_root),
        "command": " ".join(map(str, sys.argv)),
        "seed": int(args.seed),
        "num_causal_snps_requested": int(args.num_causal_snps),
        "num_causal_snps_realized": realized_num_causal,
        "h2": float(args.h2),
        "pi": float(pi),
        "beta_sd": float(beta_sd),
        "sigma_e": float(sigma_e),
        "num_samples": int(y.shape[0]),
        "num_variants": total_variants,
        "y_mean": float(np.mean(y)),
        "y_var": float(np.var(y)),
        "genetic_component_mean": float(np.mean(y_genetic)),
        "genetic_component_var": float(np.var(y_genetic)),
        "noise_mean": float(np.mean(epsilon)),
        "noise_var": float(np.var(epsilon)),
        "blocks": block_summaries,
        "party_sample_counts": {party: len(ids) for party, ids in party_sample_ids.items()},
    }

    (run_dir / "simulation_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    with (run_dir / "causal_variant_ids.txt").open("w") as fh:
        for variant_id in np.array(all_variant_ids, dtype=object)[causal_mask_all]:
            fh.write(f"{variant_id}\n")

    with (run_dir / "sim_beta.tsv").open("w") as fh:
        fh.write("variant_id\tbeta\tis_causal\n")
        for variant_id, beta in zip(all_variant_ids, beta_all):
            fh.write(f"{variant_id}\t{beta:.18e}\t{int(beta != 0.0)}\n")

    y_offset = 0
    for party_dir in party_dirs:
        party_name = party_dir.name
        ids = party_sample_ids[party_name]
        n_party = len(ids)
        party_y = y[y_offset : y_offset + n_party]
        y_offset += n_party

        pheno_path = party_dir / "pheno.txt"
        backup_path = run_dir / f"{party_name}_pheno_before.txt"
        shutil.copy2(pheno_path, backup_path)
        with pheno_path.open("w") as fh:
            for value in party_y:
                fh.write(f"{value:.12f}\n")
        with (run_dir / f"{party_name}_pheno_after.txt").open("w") as fh:
            for value in party_y:
                fh.write(f"{value:.12f}\n")
        with (run_dir / f"{party_name}_sample_ids.txt").open("w") as fh:
            for sample_id in ids:
                fh.write(f"{sample_id}\n")

    if not args.keep_cache:
        shutil.rmtree(cache_dir)

    print(f"Simulation run directory: {run_dir}")
    print(f"Requested causal SNPs: {args.num_causal_snps}")
    print(f"Realized causal SNPs: {realized_num_causal}")
    print(f"Phenotype mean: {np.mean(y):.6f}")
    print(f"Phenotype variance: {np.var(y):.6f}")
    print(f"Genetic component variance: {np.var(y_genetic):.6f}")
    print(f"Noise variance: {np.var(epsilon):.6f}")


if __name__ == "__main__":
    main()
