"""R helper bridge for block-wise package-reference SKAT and burden statistics."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pandas as pd


R_REFERENCE_HELPER = Path(__file__).resolve().parent.parent / "r_skat_reference.R"


def require_executable(arg_name: str) -> str:
    path = shutil.which(arg_name)
    if not path:
        raise RuntimeError(f"{arg_name} not found in PATH")
    return path


def write_reference_inputs(arg_ctx: dict, arg_manual_result: dict) -> tuple[Path, Path, Path]:
    selected_blocks = sorted(arg_manual_result["blocks"], key=lambda block: int(block["block_index"]))

    manifest_path = arg_ctx["scratch_dir"] / "reference_manifest.tsv"
    pheno_path = arg_ctx["scratch_dir"] / "reference_pheno.tsv"
    cov_path = arg_ctx["scratch_dir"] / "reference_cov.tsv"

    manifest_rows = []
    for block in selected_blocks:
        keep_path = arg_ctx["scratch_dir"] / f"reference_block{block['block_index']:02d}_variant_ids.txt"
        keep_path.write_text("\n".join(block["variant_ids"]) + "\n")
        manifest_rows.append(
            {
                "block": block["block_index"],
                "raw_path_party1": str(block["raw_paths_by_party"][0]),
                "raw_path_party2": str(block["raw_paths_by_party"][1]),
                "variant_ids_path": str(keep_path),
                "n_variants": block["n_variants"],
            }
        )

    pd.DataFrame(manifest_rows).to_csv(manifest_path, sep="\t", index=False)
    pd.DataFrame({"y": arg_ctx["model"]["y"]}).to_csv(pheno_path, sep="\t", index=False)
    X = arg_ctx["model"]["X"]
    cov_df = pd.DataFrame(X, columns=[f"X{i}" for i in range(1, X.shape[1] + 1)])
    cov_df.to_csv(cov_path, sep="\t", index=False)

    return manifest_path, pheno_path, cov_path


def run_reference(arg_ctx: dict, arg_manual_result: dict) -> dict:
    if arg_ctx["skip_reference"]:
        return {
            "blocks": {},
            "n_blocks": 0,
            "summary_path": None,
            "skipped_reason": "skipped by --skip-reference",
        }

    rscript = require_executable("Rscript")
    manifest_path, pheno_path, cov_path = write_reference_inputs(arg_ctx, arg_manual_result)
    out_path = arg_ctx["output_dir"] / "reference_block_summary.tsv"
    cmd = [
        rscript,
        str(R_REFERENCE_HELPER),
        "--manifest",
        str(manifest_path),
        "--pheno",
        str(pheno_path),
        "--cov",
        str(cov_path),
        "--out",
        str(out_path),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        stderr = proc.stderr.strip() or proc.stdout.strip()
        raise RuntimeError(f"R SKAT reference failed: {stderr}")

    summary_df = pd.read_csv(out_path, sep="\t")
    required_cols = {"block", "skat_q", "burden_q", "n_markers"}
    missing_cols = required_cols.difference(summary_df.columns)
    if missing_cols:
        raise RuntimeError(f"Missing expected reference columns in {out_path}: {sorted(missing_cols)}")

    block_rows = {}
    for _, row in summary_df.iterrows():
        block_index = int(row["block"])
        block_rows[block_index] = {
            "skat_q": float(row["skat_q"]),
            "burden_q": float(row["burden_q"]),
            "n_markers": int(row["n_markers"]),
        }

    return {
        "blocks": block_rows,
        "n_blocks": len(block_rows),
        "summary_path": out_path,
        "skipped_reason": None,
    }
