from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

from cloud import cloud_cp_recursive
from utils import ensure_clean_dir, is_gcs_uri, resolve_path


def replace_toml_assignment(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"(?m)^({re.escape(key)}\s*=\s*).*$")
    updated, count = pattern.subn(lambda match: f"{match.group(1)}{value}", text, count=1)
    if count != 1:
        raise ValueError(f"Could not patch TOML key '{key}'")
    return updated


def patch_toml_file(path: Path, updates: dict[str, str]) -> None:
    text = path.read_text()
    for key, value in updates.items():
        text = replace_toml_assignment(text, key, value)
    path.write_text(text)


def emit_config_dir(
    args: argparse.Namespace,
    out_dataset: Path,
    config_out_dir: Path,
    n_party1: int,
    n_party2: int,
    num_covs: int,
    windows: list[dict[str, object]],
) -> None:
    template_dir = resolve_path(args.config_template_dir)
    if not template_dir.is_dir():
        raise FileNotFoundError(f"Missing config template directory: {template_dir}")
    ensure_clean_dir(config_out_dir, args.force)
    shutil.rmtree(config_out_dir)
    shutil.copytree(template_dir, config_out_dir)

    total_variants = sum(int(row["variant_count"]) for row in windows)
    num_blocks = len(windows)
    shared_keys_path = resolve_path(args.shared_keys_path)

    patch_toml_file(
        config_out_dir / "configGlobal.toml",
        {
            "num_inds": f"[0, {n_party1}, {n_party2}]",
            "num_snps": str(total_variants),
            "num_covs": str(num_covs),
            "geno_num_blocks": str(num_blocks),
            "use_precomputed_geno_count": "false",
            "skip_qc": "true",
            "skip_pca": "true",
            "imiss_ub": "1.0",
            "het_lb": "0.0",
            "het_ub": "1.0",
            "gmiss": "1.0",
            "maf_lb": "0.0",
            "hwe_ub": "1.0e18",
            "snp_dist_thres": "0",
        },
    )

    for party_id in (1, 2):
        party_dir = out_dataset / f"party{party_id}"
        patch_toml_file(
            config_out_dir / f"configLocal.Party{party_id}.toml",
            {
                "shared_keys_path": f'"{shared_keys_path}"',
                "geno_binary_file_prefix": f'"{party_dir / "geno" / "chr%d"}"',
                "geno_num_blocks": str(num_blocks),
                "geno_block_size_file": f'"{party_dir / "chrom_sizes.txt"}"',
                "pheno_file": f'"{party_dir / "pheno.txt"}"',
                "covar_file": f'"{party_dir / "cov.txt"}"',
                "snp_position_file": f'"{party_dir / "snp_pos.txt"}"',
                "sample_keep_file": f'"{party_dir / "sample_keep.txt"}"',
                "snp_ids_file": f'"{party_dir / "snp_ids.txt"}"',
                "geno_count_file": '""',
            },
        )

    patch_toml_file(
        config_out_dir / "configLocal.Party0.toml",
        {
            "shared_keys_path": f'"{shared_keys_path}"',
        },
    )
    print(f"Config directory written to {config_out_dir}")


def write_manifest(args: argparse.Namespace, out_dataset: Path, config_out_dir: Path, source_prefix: Path) -> None:
    if args.pheno_file and args.cov_file:
        phenotype_input_mode = "split_table"
    elif args.pheno_file:
        phenotype_input_mode = "merged_table"
    else:
        phenotype_input_mode = "aligned_vector_matrix"
    cov_selector = args.cov_cols or args.cov_col_indices or ("all_except_first_id_column" if args.cov_file else "")
    manifest = {
        "chromosome": str(args.chromosome),
        "input_type": "pgen" if args.pgen_prefix else "vcf",
        "pgen_prefix": args.pgen_prefix or "",
        "vcf": args.vcf or "",
        "vcf_keep": args.vcf_keep or "",
        "vcf_double_id": str(not args.vcf_no_double_id),
        "source_prefix": str(source_prefix),
        "phenotype_input_mode": phenotype_input_mode,
        "pheno_file": str(resolve_path(args.pheno_file)) if args.pheno_file else "",
        "cov_file": str(resolve_path(args.cov_file)) if args.cov_file else "",
        "id_col": args.id_col or "",
        "pheno_col": args.pheno_col or "",
        "pheno_col_index": str(args.pheno_col_index or ""),
        "cov_cols": args.cov_cols or "",
        "cov_col_indices": args.cov_col_indices or "",
        "cov_selector": cov_selector,
        "pheno_sep": args.pheno_sep or "",
        "cov_sep": args.cov_sep or "",
        "pheno_vector_file": str(resolve_path(args.pheno_vector_file)) if args.pheno_vector_file else "",
        "cov_matrix_file": str(resolve_path(args.cov_matrix_file)) if args.cov_matrix_file else "",
        "normalize_covariates": args.normalize_covariates,
        "normalize_phenotype": args.normalize_phenotype,
        "n_samples": str(args.n_samples),
        "party1_frac": str(args.party1_frac),
        "seed": str(args.seed),
        "maf_threshold": str(args.maf_threshold),
        "window_bp": str(args.window_bp),
        "min_rare_per_window": str(args.min_rare_per_window),
        "max_windows": str(args.max_windows),
        "block_mode": args.block_mode,
        "new_id_max_allele_len": str(args.new_id_max_allele_len),
        "out_dataset": str(out_dataset),
        "config_out_dir": str(config_out_dir),
    }
    lines = [f"{key}={value}" for key, value in manifest.items()]
    (out_dataset / "pgen_window_build_manifest.txt").write_text("\n".join(lines) + "\n")


def backup_outputs(args: argparse.Namespace, out_dataset: Path, config_out_dir: Path) -> None:
    if not args.backup_gcs_uri:
        return
    if not is_gcs_uri(args.backup_gcs_uri):
        raise ValueError("--backup-gcs-uri must start with gs://")
    dst = args.backup_gcs_uri.rstrip("/") + "/"
    cloud_cp_recursive(out_dataset, dst, args)
    cloud_cp_recursive(config_out_dir, dst, args)
    print(f"Backed up dataset/config to {dst}")
