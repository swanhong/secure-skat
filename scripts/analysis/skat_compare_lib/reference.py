"""R reference bridge for invoking the public SKAT package."""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from .common import require_executable
from .models import CompareContext, ManualResults, ReferenceResult


R_REFERENCE_HELPER = Path(__file__).resolve().parent.parent / "r_skat_reference.R"


def write_reference_inputs(ctx: CompareContext, manual: ManualResults) -> tuple[Path, Path, Path, Path]:
    selected_blocks = [manual.blocks[block - 1] for block in ctx.analysis_blocks]

    manifest_path = ctx.cache_dir / "reference_manifest.tsv"
    pheno_path = ctx.cache_dir / "reference_pheno.tsv"
    cov_path = ctx.cache_dir / "reference_cov.tsv"
    weights_path = ctx.cache_dir / "reference_weights.tsv"

    manifest_rows = []
    for block in selected_blocks:
        keep_path = ctx.cache_dir / f"reference_block{block.block_index:02d}_variant_ids.txt"
        keep_path.write_text("\n".join(block.variant_ids) + "\n")
        manifest_rows.append(
            {
                "block": block.block_index,
                "raw_path_party1": str(block.raw_paths_by_party[0]),
                "raw_path_party2": str(block.raw_paths_by_party[1]),
                "variant_ids_path": str(keep_path),
                "n_variants": block.n_variants,
            }
        )

    pd.DataFrame(manifest_rows).to_csv(manifest_path, sep="\t", index=False)
    pd.DataFrame({"y": ctx.model.y}).to_csv(pheno_path, sep="\t", index=False)
    cov_df = pd.DataFrame(ctx.model.x, columns=[f"X{i}" for i in range(1, ctx.model.x.shape[1] + 1)])
    cov_df.to_csv(cov_path, sep="\t", index=False)
    weights = np.concatenate([block.weight_vec for block in selected_blocks]).astype(float)
    pd.DataFrame({"weight": weights}).to_csv(weights_path, sep="\t", index=False)

    return manifest_path, pheno_path, cov_path, weights_path


def run_reference(ctx: CompareContext, manual: ManualResults) -> ReferenceResult:
    if ctx.skip_reference:
        return ReferenceResult(
            skat_q=float("nan"),
            burden_q=float("nan"),
            n_markers=0,
            summary_path=None,
            skipped_reason="skipped by --skip-reference",
        )

    rscript = require_executable("Rscript")
    manifest_path, pheno_path, cov_path, weights_path = write_reference_inputs(ctx, manual)
    out_path = ctx.cache_dir / "reference_summary.tsv"
    cmd = [
        rscript,
        str(R_REFERENCE_HELPER),
        "--manifest",
        str(manifest_path),
        "--pheno",
        str(pheno_path),
        "--cov",
        str(cov_path),
        "--weights",
        str(weights_path),
        "--out",
        str(out_path),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        stderr = proc.stderr.strip() or proc.stdout.strip()
        raise RuntimeError(f"R SKAT reference failed: {stderr}")
    summary_df = pd.read_csv(out_path, sep="\t")
    if summary_df.shape[0] != 1:
        raise RuntimeError(f"Expected exactly one row in {out_path}")
    row = summary_df.iloc[0]
    return ReferenceResult(
        skat_q=float(row["skat_q"]),
        burden_q=float(row["burden_q"]),
        n_markers=int(row["n_markers"]),
        summary_path=out_path,
    )
