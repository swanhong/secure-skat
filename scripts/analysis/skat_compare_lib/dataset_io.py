"""Dataset-facing readers and PLINK raw export helpers."""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from .models import BlockGenotypeInput, CompareContext, DatasetInputs, ExportedBlock


def read_chrom_sizes(dataset_root: Path) -> list[int]:
    chrom_path = dataset_root / "party1" / "chrom_sizes.txt"
    if not chrom_path.exists():
        raise RuntimeError(f"Missing chrom_sizes.txt under {dataset_root}")
    values = [int(line.strip()) for line in chrom_path.read_text().splitlines() if line.strip()]
    if not values:
        raise RuntimeError(f"No blocks found in {chrom_path}")
    return values


def read_positions(party_dir: Path, total_dataset_variants: int) -> np.ndarray:
    pos_table = pd.read_csv(
        party_dir / "snp_pos.txt",
        sep="\t",
        header=None,
        dtype=int,
        engine="python",
    )
    if pos_table.shape[1] < 1:
        raise RuntimeError(f"No columns found in {party_dir / 'snp_pos.txt'}")
    pos_col = 1 if pos_table.shape[1] >= 2 else 0
    positions = pos_table.iloc[:, pos_col].to_numpy(dtype=int)
    if positions.size != total_dataset_variants:
        raise RuntimeError(
            "Position vector length mismatch: "
            f"expected {total_dataset_variants} variants but found {positions.size} "
            f"in {party_dir / 'snp_pos.txt'}"
        )
    return positions


def read_qc_filter(run_root: Path | None, total_dataset_variants: int) -> np.ndarray | None:
    if run_root is None:
        return None
    qc_path = run_root / "cache" / "party1" / "gkeep.txt"
    if not qc_path.exists():
        return None
    values = [line.strip() for line in qc_path.read_text().splitlines() if line.strip()]
    if len(values) != total_dataset_variants:
        raise RuntimeError(
            "QC filter length mismatch: "
            f"expected {total_dataset_variants} variants but found {len(values)} in {qc_path}"
        )
    return np.array([value == "1" for value in values], dtype=bool)


def read_pheno(party_dir: Path) -> np.ndarray:
    return np.atleast_1d(np.loadtxt(party_dir / "pheno.txt", dtype=float)).astype(float)


def read_cov(party_dir: Path) -> np.ndarray:
    return np.loadtxt(party_dir / "cov.txt", dtype=float, delimiter="\t", ndmin=2)


def export_block_matrix(party_dir: Path, block_index: int, cache_dir: Path, plink2: str) -> ExportedBlock:
    out_prefix = cache_dir / f"{party_dir.name}_block{block_index:02d}"
    raw_path = out_prefix.with_suffix(".raw")
    if not raw_path.exists():
        pfile_prefix = party_dir / "geno" / f"chr{block_index}"
        if not pfile_prefix.with_suffix(".pgen").exists():
            raise RuntimeError(f"Missing block genotype file: {pfile_prefix}.pgen")
        sample_keep = party_dir / "sample_keep.txt"
        cmd = [plink2, "--pfile", str(pfile_prefix)]
        if pfile_prefix.with_suffix(".pvar.zst").exists():
            cmd.append("vzs")
        cmd.extend(["--keep", str(sample_keep), "--export", "A", "--out", str(out_prefix)])
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(
                f"plink2 export failed for {party_dir.name} block {block_index}\n"
                f"{proc.stderr.strip()}"
            )

    raw_df = pd.read_csv(raw_path, sep=r"\s+", engine="python")
    if raw_df.shape[1] <= 6:
        raise RuntimeError(f"No genotype columns found in {raw_path}")

    geno = raw_df.iloc[:, 6:].to_numpy(dtype=float, copy=True)
    variant_ids = list(raw_df.columns[6:])
    return ExportedBlock(raw_path=raw_path, geno=geno, variant_ids=variant_ids)


def load_dataset_inputs(ctx: CompareContext) -> DatasetInputs:
    block_inputs: list[BlockGenotypeInput] = []

    for block_index in range(1, ctx.dataset.n_blocks + 1):
        party_exports = [
            export_block_matrix(party_dir, block_index, ctx.cache_dir, ctx.plink2)
            for party_dir in ctx.dataset.party_dirs
        ]
        if party_exports[0].variant_ids != party_exports[1].variant_ids:
            raise RuntimeError(f"Variant order mismatch at block {block_index}")

        local_alt_genotypes_by_party = [np.nan_to_num(export.geno.copy(), nan=0.0) for export in party_exports]
        variant_ids = list(party_exports[0].variant_ids)
        raw_paths_by_party = [export.raw_path for export in party_exports]

        block_start, block_end = ctx.dataset.block_offsets[block_index - 1]
        block_positions = np.asarray(ctx.all_positions[block_start:block_end], dtype=int)

        if ctx.secure_qc_filter is not None:
            block_filter = ctx.secure_qc_filter[block_start:block_end]
            variant_ids = [variant_id for variant_id, keep in zip(variant_ids, block_filter) if keep]
            block_positions = block_positions[block_filter]
            local_alt_genotypes_by_party = [geno[:, block_filter] for geno in local_alt_genotypes_by_party]

        block_inputs.append(
            BlockGenotypeInput(
                block_index=block_index,
                raw_paths_by_party=raw_paths_by_party,
                local_alt_genotypes_by_party=local_alt_genotypes_by_party,
                variant_ids=variant_ids,
                positions=block_positions,
            )
        )

    return DatasetInputs(block_inputs=block_inputs)
