from __future__ import annotations

import argparse
import os


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a secure-skat dataset from a cohort PGEN shard or VCF, "
            "phenotype/covariate inputs, and MAF-window variant sets."
        )
    )
    parser.add_argument("--chromosome", type=int, default=22, help="Chromosome label to write in snp_pos.txt.")
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--pgen-prefix",
        help=(
            "Input PGEN prefix, without extension. Can be local or gs://. "
            "For example: /data/cohort.chr22 or gs://bucket/path/cohort.chr22."
        ),
    )
    input_group.add_argument(
        "--vcf",
        help="Input VCF/BCF path. Can be local or gs://. Converted to PGEN with plink2 before windowing.",
    )
    parser.add_argument(
        "--raw-dir",
        default=".local/skat/raw",
        help="Local cache directory for gs:// input files (default: .local/skat/raw).",
    )
    parser.add_argument(
        "--work-dir",
        default=".local/skat/work",
        help="Local working directory for temporary source/filter files (default: .local/skat/work).",
    )
    parser.add_argument("--pheno-file", help="Phenotype/covariate table path for table mode.")
    parser.add_argument("--id-col", help="Column in --pheno-file matching the PGEN IID.")
    parser.add_argument("--pheno-col", help="Phenotype column to write as pheno.txt.")
    parser.add_argument(
        "--cov-cols",
        help="Comma-separated covariate columns to write as cov.txt. Do not include an intercept.",
    )
    parser.add_argument("--pheno-sep", default=None, help="Optional delimiter for --pheno-file.")
    parser.add_argument(
        "--pheno-vector-file",
        help="Phenotype vector already aligned to PSAM row order. Use with --cov-matrix-file.",
    )
    parser.add_argument(
        "--cov-matrix-file",
        help="Covariate matrix already aligned to PSAM row order. Use with --pheno-vector-file.",
    )
    parser.add_argument(
        "--out-dataset",
        required=True,
        help="Output secure-skat dataset root, e.g. dataset/pgen_chr22_windows.",
    )
    parser.add_argument(
        "--config-template-dir",
        default="config",
        help="Config template directory to copy and patch (default: config).",
    )
    parser.add_argument("--config-out-dir", required=True, help="Output config directory for the generated dataset.")
    parser.add_argument(
        "--shared-keys-path",
        default="example_data/keys",
        help="Shared key directory for local three-process runs (default: example_data/keys).",
    )
    parser.add_argument("--n-samples", type=int, default=2000, help="Smoke-test sample count; use 0 for all.")
    parser.add_argument("--party1-frac", type=float, default=0.5, help="Fraction assigned to party1 (default: 0.5).")
    parser.add_argument("--seed", type=int, default=1, help="Random seed for sample subset and party split.")
    parser.add_argument("--maf-threshold", type=float, default=0.01, help="Rare variant MAF upper bound.")
    parser.add_argument("--window-bp", type=int, default=50000, help="Fixed genomic bin/window width in bp.")
    parser.add_argument(
        "--min-rare-per-window",
        type=int,
        default=20,
        help="Minimum rare variants required for a window/block.",
    )
    parser.add_argument("--max-windows", type=int, default=8, help="Maximum number of selected windows; 0 keeps all.")
    parser.add_argument(
        "--block-mode",
        choices=("rare-only", "all-in-window"),
        default="rare-only",
        help=(
            "Which variants to materialize per selected window. 'rare-only' keeps only variants with "
            "MAF <= threshold; 'all-in-window' keeps every variant inside bins passing the rare-count threshold."
        ),
    )
    parser.add_argument("--plink2", default="plink2", help="Path to plink2 executable.")
    parser.add_argument(
        "--new-id-max-allele-len",
        type=int,
        default=1000,
        help=(
            "PLINK2 --new-id-max-allele-len value used with --set-all-var-ids. "
            "Needed for raw VCFs containing long indels (default: 1000)."
        ),
    )
    parser.add_argument(
        "--billing-project",
        default=os.environ.get("GOOGLE_PROJECT", ""),
        help="Billing/quota project for requester-pays GCS buckets (default: $GOOGLE_PROJECT).",
    )
    parser.add_argument(
        "--cloud-cli",
        choices=("gcloud", "gsutil"),
        default="gcloud",
        help="CLI used for gs:// copies (default: gcloud).",
    )
    parser.add_argument(
        "--backup-gcs-uri",
        default=None,
        help="Optional gs:// destination. The output dataset and config directory are copied there at the end.",
    )
    parser.add_argument(
        "--vcf-keep",
        default=None,
        help="Optional PLINK --keep file applied during VCF-to-PGEN conversion.",
    )
    parser.add_argument(
        "--vcf-no-double-id",
        action="store_true",
        help=(
            "Do not pass PLINK2 --double-id during VCF import. By default VCF sample IDs become both FID and IID, "
            "matching the 1000G keep-file convention used by this repo."
        ),
    )
    parser.add_argument("--force", action="store_true", help="Overwrite existing output dataset/config/work files.")

    args = parser.parse_args()
    if args.chromosome < 1 or args.chromosome > 22:
        parser.error("--chromosome must be in 1..22")
    if args.n_samples < 0:
        parser.error("--n-samples must be non-negative")
    if not (0.0 < args.party1_frac < 1.0):
        parser.error("--party1-frac must be in (0, 1)")
    if not (0.0 < args.maf_threshold < 0.5):
        parser.error("--maf-threshold must be in (0, 0.5)")
    if args.window_bp <= 0:
        parser.error("--window-bp must be positive")
    if args.min_rare_per_window <= 0:
        parser.error("--min-rare-per-window must be positive")
    if args.max_windows < 0:
        parser.error("--max-windows must be non-negative")
    if args.new_id_max_allele_len <= 0:
        parser.error("--new-id-max-allele-len must be positive")

    table_mode = args.pheno_file is not None
    aligned_mode = args.pheno_vector_file is not None or args.cov_matrix_file is not None
    if table_mode == aligned_mode:
        parser.error("choose exactly one phenotype mode: --pheno-file or --pheno-vector-file/--cov-matrix-file")
    if table_mode:
        missing = [
            flag
            for flag, value in (
                ("--id-col", args.id_col),
                ("--pheno-col", args.pheno_col),
                ("--cov-cols", args.cov_cols),
            )
            if not value
        ]
        if missing:
            parser.error(f"table mode requires {', '.join(missing)}")
    if aligned_mode and (not args.pheno_vector_file or not args.cov_matrix_file):
        parser.error("aligned mode requires both --pheno-vector-file and --cov-matrix-file")
    return args
