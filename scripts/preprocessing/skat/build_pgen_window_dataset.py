#!/usr/bin/env python3

from __future__ import annotations

import sys

from argument import parse_args
from config import backup_outputs, emit_config_dir, write_manifest
from plink import materialize_party_blocks, materialize_raw_pgen, run_plink_source
from samples import build_sample_files
from utils import ensure_clean_dir, repo_root, resolve_path
from windowing import read_source_variants, select_windows


def main() -> int:
    args = parse_args()
    root = repo_root()

    raw_dir = resolve_path(args.raw_dir)
    work_dir = resolve_path(args.work_dir)
    out_dataset = resolve_path(args.out_dataset)
    config_out_dir = resolve_path(args.config_out_dir)

    ensure_clean_dir(out_dataset, args.force)
    work_dir.mkdir(parents=True, exist_ok=True)

    raw_prefix = materialize_raw_pgen(args, raw_dir)
    n_party1, n_party2, num_covs = build_sample_files(args, raw_prefix, work_dir, out_dataset)
    source_prefix = run_plink_source(args, raw_prefix, work_dir)
    variants = read_source_variants(source_prefix)
    windows = select_windows(args, variants, work_dir, out_dataset)
    materialize_party_blocks(args, source_prefix, out_dataset, windows)
    emit_config_dir(args, out_dataset, config_out_dir, n_party1, n_party2, num_covs, windows)
    write_manifest(args, out_dataset, config_out_dir, source_prefix)
    backup_outputs(args, out_dataset, config_out_dir)

    print("")
    print("PGEN secure-skat dataset is ready.")
    print(f"Dataset: {out_dataset}")
    print(f"Config:  {config_out_dir}")
    print("")
    print("Run:")
    print(f"  cd {root}")
    print(f"  bash run_example.sh --mode skat --dataset {out_dataset} --config-dir {config_out_dir}")
    print("")
    print("Optional SKAT-O:")
    print(f"  bash run_example.sh --mode skato --skato-rho 0.5 --dataset {out_dataset} --config-dir {config_out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
