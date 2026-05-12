"""Simplified plain-vs-secure SKAT comparison pipeline."""

from __future__ import annotations

import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from . import compute, plotting, reference


def require_executable(arg_name: str) -> str:
    path = shutil.which(arg_name)
    if not path:
        raise RuntimeError(f"{arg_name} not found in PATH")
    return path


def safe_rel_diff(arg_x: float, arg_y: float) -> float:
    if not (math.isfinite(arg_x) and math.isfinite(arg_y)):
        return float("nan")
    return abs(arg_x - arg_y) / max(abs(arg_y), 1e-12)


def safe_corr(arg_x: pd.Series | np.ndarray, arg_y: pd.Series | np.ndarray) -> float:
    x = np.asarray(arg_x, dtype=float)
    y = np.asarray(arg_y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    if int(np.sum(keep)) < 2:
        return float("nan")
    return float(np.corrcoef(x[keep], y[keep])[0, 1])


def safe_r2(arg_x: pd.Series | np.ndarray, arg_y: pd.Series | np.ndarray) -> float:
    corr = safe_corr(arg_x, arg_y)
    if not math.isfinite(corr):
        return float("nan")
    return float(corr * corr)


def count_finite_pairs(arg_x: pd.Series | np.ndarray, arg_y: pd.Series | np.ndarray) -> int:
    x = np.asarray(arg_x, dtype=float)
    y = np.asarray(arg_y, dtype=float)
    return int(np.sum(np.isfinite(x) & np.isfinite(y)))


def format_float(arg_value: float) -> str:
    return "NA" if not math.isfinite(arg_value) else f"{arg_value:.10e}"

def read_kv_file(arg_path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    if not arg_path.exists():
        return out
    for line in arg_path.read_text().splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        out[key.strip()] = value.strip()
    return out


def count_nonblank_lines(arg_path: Path) -> int:
    return sum(1 for line in arg_path.read_text().splitlines() if line.strip())


def parse_block_spec(arg_spec: str | None, arg_n_blocks: int) -> list[int]:
    if arg_spec is None:
        return list(range(1, arg_n_blocks + 1))

    out: set[int] = set()
    for token in [part.strip() for part in arg_spec.split(",") if part.strip()]:
        if token == "all":
            out.update(range(1, arg_n_blocks + 1))
            continue
        if token == "last":
            out.add(arg_n_blocks)
            continue
        if re.fullmatch(r"\d+-\d+", token):
            start_text, end_text = token.split("-", 1)
            start_idx = int(start_text)
            end_idx = int(end_text)
            if start_idx > end_idx:
                raise RuntimeError(f"Invalid block range: {token}")
            out.update(range(start_idx, end_idx + 1))
            continue
        try:
            out.add(int(token))
        except ValueError as exc:
            raise RuntimeError(f"Invalid block specification token: {token}") from exc

    keep = sorted(block_index for block_index in out if 1 <= block_index <= arg_n_blocks)
    if not keep:
        raise RuntimeError("No valid blocks remain after applying --blocks")
    return keep


def compute_block_offsets(arg_chrom_sizes: list[int]) -> list[tuple[int, int]]:
    offsets = []
    start_idx = 0
    for block_size in arg_chrom_sizes:
        end_idx = start_idx + int(block_size)
        offsets.append((start_idx, end_idx))
        start_idx = end_idx
    return offsets


def read_chrom_sizes(arg_dataset_root: Path) -> list[int]:
    chrom_path = arg_dataset_root / "party1" / "chrom_sizes.txt"
    if not chrom_path.exists():
        raise RuntimeError(f"Missing chrom_sizes.txt under {arg_dataset_root}")
    chrom_sizes = [int(line.strip()) for line in chrom_path.read_text().splitlines() if line.strip()]
    if not chrom_sizes:
        raise RuntimeError(f"No blocks found in {chrom_path}")
    return chrom_sizes


def read_positions(arg_party_dir: Path, arg_total_variants: int) -> np.ndarray:
    pos_df = pd.read_csv(arg_party_dir / "snp_pos.txt", sep="\t", header=None, dtype=int, engine="python")
    if pos_df.shape[1] < 1:
        raise RuntimeError(f"No columns found in {arg_party_dir / 'snp_pos.txt'}")
    pos_col = 1 if pos_df.shape[1] >= 2 else 0
    positions = pos_df.iloc[:, pos_col].to_numpy(dtype=int)
    if positions.size != arg_total_variants:
        raise RuntimeError(
            "Position vector length mismatch: "
            f"expected {arg_total_variants} variants but found {positions.size}"
        )
    return positions


def read_qc_filter(arg_run_root: Path | None, arg_total_variants: int) -> np.ndarray | None:
    if arg_run_root is None:
        return None
    qc_path = arg_run_root / "cache" / "party1" / "gkeep.txt"
    if not qc_path.exists():
        return None
    values = [line.strip() for line in qc_path.read_text().splitlines() if line.strip()]
    if len(values) != arg_total_variants:
        raise RuntimeError(
            "QC filter length mismatch: "
            f"expected {arg_total_variants} variants but found {len(values)} in {qc_path}"
        )
    return np.asarray([value == "1" for value in values], dtype=bool)


def read_pheno(arg_party_dir: Path) -> np.ndarray:
    return np.atleast_1d(np.loadtxt(arg_party_dir / "pheno.txt", dtype=float)).astype(float)


def read_cov(arg_party_dir: Path) -> np.ndarray:
    return np.loadtxt(arg_party_dir / "cov.txt", dtype=float, delimiter="\t", ndmin=2)


def resolve_run_root(arg_repo_root: Path, arg_run_id: str) -> Path:
    run_roots = [path for path in (arg_repo_root / "out").glob(f"output_*_{arg_run_id}") if path.is_dir()]
    if not run_roots:
        raise RuntimeError(f"No secure output directory found for run id: {arg_run_id}")
    run_roots.sort(key=lambda path: path.stat().st_mtime, reverse=True)
    return run_roots[0].resolve()


def resolve_run_root_arg(arg_repo_root: Path, arg_run_id: str | None, arg_run_root: str | None) -> Path:
    if arg_run_root:
        run_root = Path(arg_run_root)
        if not run_root.is_absolute():
            run_root = arg_repo_root / run_root
        run_root = run_root.resolve()
        if not run_root.is_dir():
            raise RuntimeError(f"Run root does not exist: {run_root}")
        return run_root
    if not arg_run_id:
        raise RuntimeError("Provide --run-id or --run-root")
    return resolve_run_root(arg_repo_root, arg_run_id)


def resolve_dataset_path(arg_repo_root: Path, arg_dataset_value: str) -> Path:
    dataset_root = Path(arg_dataset_value)
    if not dataset_root.is_absolute():
        dataset_root = arg_repo_root / dataset_root
    dataset_root = dataset_root.resolve()
    if not (dataset_root / "party1" / "chrom_sizes.txt").exists():
        raise RuntimeError(f"Missing dataset directory or chrom_sizes.txt: {dataset_root}")
    return dataset_root


def candidate_dataset_roots(arg_repo_root: Path) -> list[Path]:
    roots = []
    example_root = arg_repo_root / "example_data"
    if (example_root / "party1" / "chrom_sizes.txt").exists():
        roots.append(example_root)

    datasets_root = arg_repo_root / ".local" / "datasets"
    if datasets_root.exists():
        for child in sorted(datasets_root.iterdir()):
            if (child / "party1" / "chrom_sizes.txt").exists():
                roots.append(child)
    return roots


def infer_run_shape(arg_repo_root: Path, arg_run_root: Path, arg_run_metadata: dict[str, str]) -> tuple[int | None, int | None]:
    block_count = None
    total_variants = None

    if arg_run_metadata.get("dataset"):
        dataset_root = resolve_dataset_path(arg_repo_root, arg_run_metadata["dataset"])
        chrom_sizes = read_chrom_sizes(dataset_root)
        block_count = len(chrom_sizes)
        total_variants = sum(chrom_sizes)

    if block_count is None:
        qblock_files = sorted((arg_run_root / "party1").glob("qBlock_block*.txt"))
        if qblock_files:
            block_count = len(qblock_files)
        else:
            cache_files = sorted((arg_run_root / "cache" / "party1").glob("assoc_cache_dos_sum.skat.*.txt"))
            if cache_files:
                block_count = len(cache_files)

    if total_variants is None:
        qc_path = arg_run_root / "cache" / "party1" / "gkeep.txt"
        if qc_path.exists():
            total_variants = count_nonblank_lines(qc_path)

    return block_count, total_variants


def resolve_dataset(arg_repo_root: Path, arg_run_root: Path, arg_run_metadata: dict[str, str], arg_dataset_arg: str | None) -> dict:
    if arg_dataset_arg:
        dataset_root = resolve_dataset_path(arg_repo_root, arg_dataset_arg)
    elif arg_run_metadata.get("dataset"):
        dataset_root = resolve_dataset_path(arg_repo_root, arg_run_metadata["dataset"])
    else:
        block_count, total_variants = infer_run_shape(arg_repo_root, arg_run_root, arg_run_metadata)
        if block_count is None or total_variants is None:
            raise RuntimeError("Unable to infer dataset; provide --dataset explicitly")
        matches = []
        for candidate_root in candidate_dataset_roots(arg_repo_root):
            chrom_sizes = read_chrom_sizes(candidate_root)
            if len(chrom_sizes) == block_count and sum(chrom_sizes) == total_variants:
                matches.append(candidate_root.resolve())
        if len(matches) != 1:
            raise RuntimeError(
                "Dataset inference failed: expected exactly one local dataset match for "
                f"{block_count} blocks and {total_variants} variants, found {len(matches)}"
            )
        dataset_root = matches[0]

    chrom_sizes = read_chrom_sizes(dataset_root)
    total_variants = sum(chrom_sizes)
    n_blocks = len(chrom_sizes)
    run_block_count, run_total_variants = infer_run_shape(arg_repo_root, arg_run_root, arg_run_metadata)
    if run_block_count is not None and run_block_count != n_blocks:
        raise RuntimeError(
            f"Dataset/run block mismatch: dataset has {n_blocks} blocks, run implies {run_block_count} blocks"
        )
    if run_total_variants is not None and run_total_variants != total_variants:
        raise RuntimeError(
            f"Dataset/run variant-count mismatch: dataset has {total_variants} variants, run implies {run_total_variants} variants"
        )

    return {
        "dataset_root": dataset_root.resolve(),
        "party_dirs": [(dataset_root / "party1").resolve(), (dataset_root / "party2").resolve()],
        "chrom_sizes": chrom_sizes,
        "total_variants": total_variants,
        "n_blocks": n_blocks,
        "block_offsets": compute_block_offsets(chrom_sizes),
    }


def export_block_matrix(arg_party_dir: Path, arg_block_index: int, arg_cache_dir: Path, arg_plink2: str) -> dict:
    out_prefix = arg_cache_dir / f"{arg_party_dir.name}_block{arg_block_index:02d}"
    raw_path = out_prefix.with_suffix(".raw")
    if not raw_path.exists():
        pfile_prefix = arg_party_dir / "geno" / f"chr{arg_block_index}"
        if not pfile_prefix.with_suffix(".pgen").exists():
            raise RuntimeError(f"Missing block genotype file: {pfile_prefix}.pgen")
        cmd = [arg_plink2, "--pfile", str(pfile_prefix)]
        if pfile_prefix.with_suffix(".pvar.zst").exists():
            cmd.append("vzs")
        cmd.extend(["--keep", str(arg_party_dir / "sample_keep.txt"), "--export", "A", "--out", str(out_prefix)])
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(
                f"plink2 export failed for {arg_party_dir.name} block {arg_block_index}\n{proc.stderr.strip()}"
            )

    raw_df = pd.read_csv(raw_path, sep=r"\s+", engine="python")
    if raw_df.shape[1] <= 6:
        raise RuntimeError(f"No genotype columns found in {raw_path}")
    return {
        "raw_path": raw_path,
        "geno": raw_df.iloc[:, 6:].to_numpy(dtype=float, copy=True),
        "variant_ids": list(raw_df.columns[6:]),
    }


def raw_alt_to_secure_oriented(arg_raw_alt: np.ndarray) -> np.ndarray:
    raw_alt = np.asarray(arg_raw_alt, dtype=float)
    return np.where(np.isfinite(raw_alt), 2.0 - raw_alt, 0.0)


def load_plain_blocks(arg_ctx: dict) -> list[dict]:
    block_inputs = []
    for block_index in arg_ctx["analysis_blocks"]:
        party_exports = [
            export_block_matrix(party_dir, block_index, arg_ctx["scratch_dir"], arg_ctx["plink2"])
            for party_dir in arg_ctx["party_dirs"]
        ]
        if party_exports[0]["variant_ids"] != party_exports[1]["variant_ids"]:
            raise RuntimeError(f"Variant order mismatch at block {block_index}")

        variant_ids = list(party_exports[0]["variant_ids"])
        raw_paths_by_party = [party_export["raw_path"] for party_export in party_exports]
        local_secure_genotypes_by_party = [
            raw_alt_to_secure_oriented(party_export["geno"]) for party_export in party_exports
        ]

        start_idx, end_idx = arg_ctx["block_offsets"][block_index - 1]
        positions = np.asarray(arg_ctx["all_positions"][start_idx:end_idx], dtype=int)
        if arg_ctx["secure_qc_filter"] is not None:
            keep = arg_ctx["secure_qc_filter"][start_idx:end_idx]
            variant_ids = [variant_id for variant_id, keep_flag in zip(variant_ids, keep) if keep_flag]
            positions = positions[keep]
            local_secure_genotypes_by_party = [G_part[:, keep] for G_part in local_secure_genotypes_by_party]

        block_inputs.append(
            {
                "block_index": block_index,
                "raw_paths_by_party": raw_paths_by_party,
                "variant_ids": variant_ids,
                "positions": positions,
                "local_secure_genotypes_by_party": local_secure_genotypes_by_party,
            }
        )

    return block_inputs


def read_secure_scalar(arg_path: Path) -> float | None:
    if not arg_path.exists():
        return None
    value = np.loadtxt(arg_path, dtype=float, max_rows=1)
    return float(np.asarray(value).reshape(-1)[0])


def load_secure_compare(arg_ctx: dict) -> dict:
    party1_dir = arg_ctx["run_root"] / "party1"
    rare_variant_scale = arg_ctx["model"]["rare_variant_scale"]

    block_raw = {}
    selected_q_skat_raw_total = 0.0
    selected_q_burden_raw_total = 0.0
    selected_available = True
    n_blocks_available = 0
    missing_blocks = []

    for block_index in arg_ctx["analysis_blocks"]:
        q_skat_block_raw = read_secure_scalar(party1_dir / f"qBlock_block{block_index - 1}.txt")
        q_burden_block_raw = read_secure_scalar(party1_dir / f"qBurdenBlock_block{block_index - 1}.txt")
        block_raw[block_index] = {
            "q_skat_block_raw": q_skat_block_raw,
            "q_burden_block_raw": q_burden_block_raw,
        }
        if q_skat_block_raw is None or q_burden_block_raw is None:
            selected_available = False
            missing_blocks.append(block_index)
            continue
        n_blocks_available += 1
        selected_q_skat_raw_total += q_skat_block_raw
        selected_q_burden_raw_total += q_burden_block_raw

    selected_skat_q = selected_q_skat_raw_total * rare_variant_scale if selected_available else float("nan")
    selected_burden_q = (selected_q_burden_raw_total**2) * rare_variant_scale if selected_available else float("nan")
    run_skat_q = read_secure_scalar(party1_dir / "skat_out.txt")
    run_burden_q = read_secure_scalar(party1_dir / "burden_out.txt")
    used_run_scalars = not selected_available

    return {
        "blocks": block_raw,
        "selected_available": selected_available,
        "selected_q_skat_raw_total": selected_q_skat_raw_total,
        "selected_q_burden_raw_total": selected_q_burden_raw_total,
        "selected_skat_q": selected_skat_q,
        "selected_burden_q": selected_burden_q,
        "n_blocks_available": n_blocks_available,
        "n_blocks_requested": len(arg_ctx["analysis_blocks"]),
        "missing_blocks": missing_blocks,
        "run_skat_q": float("nan") if run_skat_q is None else run_skat_q,
        "run_burden_q": float("nan") if run_burden_q is None else run_burden_q,
        "summary_skat_q": (float("nan") if run_skat_q is None else run_skat_q) if used_run_scalars else selected_skat_q,
        "summary_burden_q": (float("nan") if run_burden_q is None else run_burden_q) if used_run_scalars else selected_burden_q,
        "used_run_scalars": used_run_scalars,
    }


def build_block_compare_frame(arg_ctx: dict, arg_manual_result: dict, arg_secure_result: dict) -> pd.DataFrame:
    rows = []
    rare_variant_scale = arg_ctx["model"]["rare_variant_scale"]
    reference_blocks = arg_ctx.get("reference_blocks", {})
    for block in arg_manual_result["blocks"]:
        block_index = int(block["block_index"])
        secure_block = arg_secure_result["blocks"].get(block_index, {})
        reference_block = reference_blocks.get(block_index, {})
        secure_q_skat_raw = secure_block.get("q_skat_block_raw")
        secure_q_burden_raw = secure_block.get("q_burden_block_raw")
        reference_skat_q = reference_block.get("skat_q", float("nan"))
        reference_burden_q = reference_block.get("burden_q", float("nan"))

        plain_skat_q = block["q_skat_block_raw"] * rare_variant_scale
        plain_burden_q = (block["q_burden_block_raw"]**2) * rare_variant_scale
        secure_skat_q = float("nan") if secure_q_skat_raw is None else secure_q_skat_raw * rare_variant_scale
        secure_burden_q = float("nan") if secure_q_burden_raw is None else (secure_q_burden_raw**2) * rare_variant_scale

        rows.append(
            {
                "block": block["block_index"],
                "n_variants": block["n_variants"],
                "block_start_bp": int(block["positions"][0]) if block["n_variants"] else np.nan,
                "block_end_bp": int(block["positions"][-1]) if block["n_variants"] else np.nan,
                "start_variant_id": block["variant_ids"][0] if block["n_variants"] else "",
                "end_variant_id": block["variant_ids"][-1] if block["n_variants"] else "",
                "plain_skat_q": plain_skat_q,
                "reference_skat_q": reference_skat_q,
                "secure_skat_q": secure_skat_q,
                "skat_abs_diff": abs(plain_skat_q - secure_skat_q) if math.isfinite(secure_skat_q) else np.nan,
                "skat_rel_diff": safe_rel_diff(plain_skat_q, secure_skat_q),
                "plain_burden_q": plain_burden_q,
                "reference_burden_q": reference_burden_q,
                "secure_burden_q": secure_burden_q,
                "burden_abs_diff": abs(plain_burden_q - secure_burden_q) if math.isfinite(secure_burden_q) else np.nan,
                "burden_rel_diff": safe_rel_diff(plain_burden_q, secure_burden_q),
                "plain_burden_sum": block["q_burden_block_raw"],
                "secure_burden_sum": float("nan") if secure_q_burden_raw is None else secure_q_burden_raw,
                "plain_vs_reference_skat_abs_diff": abs(plain_skat_q - reference_skat_q) if math.isfinite(reference_skat_q) else np.nan,
                "plain_vs_reference_skat_rel_diff": safe_rel_diff(plain_skat_q, reference_skat_q),
                "secure_vs_reference_skat_abs_diff": abs(secure_skat_q - reference_skat_q) if math.isfinite(secure_skat_q) and math.isfinite(reference_skat_q) else np.nan,
                "secure_vs_reference_skat_rel_diff": safe_rel_diff(secure_skat_q, reference_skat_q),
                "plain_vs_reference_burden_abs_diff": abs(plain_burden_q - reference_burden_q) if math.isfinite(reference_burden_q) else np.nan,
                "plain_vs_reference_burden_rel_diff": safe_rel_diff(plain_burden_q, reference_burden_q),
                "secure_vs_reference_burden_abs_diff": abs(secure_burden_q - reference_burden_q) if math.isfinite(secure_burden_q) and math.isfinite(reference_burden_q) else np.nan,
                "secure_vs_reference_burden_rel_diff": safe_rel_diff(secure_burden_q, reference_burden_q),
            }
        )
    return pd.DataFrame(rows)


def summarize_pairwise_errors(arg_block_df: pd.DataFrame, arg_metric: str, arg_left_prefix: str, arg_right_prefix: str) -> dict | None:
    left_col = f"{arg_left_prefix}_{arg_metric}_q"
    right_col = f"{arg_right_prefix}_{arg_metric}_q"
    if left_col not in arg_block_df.columns or right_col not in arg_block_df.columns:
        return None

    left = arg_block_df[left_col].to_numpy(dtype=float)
    right = arg_block_df[right_col].to_numpy(dtype=float)
    keep = np.isfinite(left) & np.isfinite(right)
    if int(np.sum(keep)) == 0:
        return None

    blocks = arg_block_df.loc[keep, "block"].to_numpy(dtype=int)
    left = left[keep]
    right = right[keep]
    abs_diff = np.abs(left - right)
    rel_diff = np.asarray([safe_rel_diff(left_val, right_val) for left_val, right_val in zip(left, right)], dtype=float)

    max_abs_idx = int(np.nanargmax(abs_diff))
    max_rel_idx = int(np.nanargmax(rel_diff))

    return {
        "metric": arg_metric,
        "comparison": f"{arg_left_prefix}_vs_{arg_right_prefix}",
        "n_blocks_compared": int(blocks.size),
        "max_abs_diff": float(abs_diff[max_abs_idx]),
        "max_abs_diff_block": int(blocks[max_abs_idx]),
        "max_rel_diff": float(rel_diff[max_rel_idx]),
        "max_rel_diff_block": int(blocks[max_rel_idx]),
        "mean_abs_diff": float(np.mean(abs_diff)),
        "mean_rel_diff": float(np.mean(rel_diff)),
    }


def build_summary_frame(arg_block_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for metric in ("skat", "burden"):
        for left_prefix, right_prefix in (
            ("plain", "reference"),
            ("secure", "reference"),
            ("plain", "secure"),
        ):
            row = summarize_pairwise_errors(arg_block_df, metric, left_prefix, right_prefix)
            if row is not None:
                rows.append(row)

    return pd.DataFrame(rows)


def write_summary_csv(arg_ctx: dict, arg_block_df: pd.DataFrame) -> tuple[Path, pd.DataFrame]:
    summary_path = arg_ctx["output_dir"] / "summary.csv"
    summary_df = build_summary_frame(arg_block_df)
    summary_df.to_csv(summary_path, index=False)
    return summary_path, summary_df


def print_preflight(arg_ctx: dict) -> None:
    print("\n--- Preflight ---")
    print(f"Command: {arg_ctx['command']}")
    print(f"Resolved dataset: {arg_ctx['dataset_root']}")
    print(f"Run root: {arg_ctx['run_root']}")
    print(f"Blocks: {arg_ctx['n_blocks']}")
    print(f"Variants: {arg_ctx['total_variants']}")
    print(
        "Analysis blocks: "
        f"{arg_ctx['analysis_blocks'][0]}..{arg_ctx['analysis_blocks'][-1]} "
        f"({len(arg_ctx['analysis_blocks'])} selected)"
    )
    if arg_ctx["command"] == "compare":
        print(f"Reference step: {'skipped' if arg_ctx['skip_reference'] else 'enabled'}")
        print(f"Plain mode: {arg_ctx['plain_mode']}")
        if arg_ctx["plain_mode"] == compute.PLAIN_MODE_LOCAL_WEIGHT_BURDEN:
            print(f"Local weight mode: {arg_ctx['local_weight_mode']}")
    else:
        print("Reference step: enabled")
    print(f"Output directory: {arg_ctx['output_dir']}")
    sys.stdout.flush()


def print_block_comparison_summary(arg_block_df: pd.DataFrame, arg_block_csv_path: Path, arg_png_paths: list[Path]) -> None:
    skat_pair_count = count_finite_pairs(arg_block_df["plain_skat_q"], arg_block_df["secure_skat_q"])
    burden_pair_count = count_finite_pairs(arg_block_df["plain_burden_q"], arg_block_df["secure_burden_q"])
    skat_corr = safe_corr(arg_block_df["plain_skat_q"], arg_block_df["secure_skat_q"])
    burden_corr = safe_corr(arg_block_df["plain_burden_q"], arg_block_df["secure_burden_q"])
    skat_r2 = safe_r2(arg_block_df["plain_skat_q"], arg_block_df["secure_skat_q"])
    burden_r2 = safe_r2(arg_block_df["plain_burden_q"], arg_block_df["secure_burden_q"])

    print("\n--- Block Comparison ---")
    print(f"Blocks retained: {len(arg_block_df)}")
    print(f"Block summary CSV: {arg_block_csv_path}")
    if len(arg_png_paths) >= 1:
        print(f"Block-level SKAT scatter plot: {arg_png_paths[0]}")
    else:
        print(f"Block-level SKAT scatter plot: not generated ({skat_pair_count} finite plain/secure pairs)")
    if len(arg_png_paths) >= 2:
        print(f"Block-level Burden scatter plot: {arg_png_paths[1]}")
    else:
        print(f"Block-level Burden scatter plot: not generated ({burden_pair_count} finite plain/secure pairs)")
    print(f"Block SKAT corr: {skat_corr:.10f}")
    print(f"Block SKAT r2: {skat_r2:.10f}")
    print(f"Block Burden corr: {burden_corr:.10f}")
    print(f"Block Burden r2: {burden_r2:.10f}")


def print_aggregate_comparison_summary(arg_ctx: dict, arg_manual_result: dict, arg_secure_result: dict) -> None:
    print("\n--- Aggregate Plain vs Secure ---")
    print(
        "Secure aggregate source: "
        f"{'run-level scalars' if arg_secure_result['used_run_scalars'] else 'per-block raw sums'}"
    )
    print(
        "Secure block raw values available: "
        f"{arg_secure_result['n_blocks_available']} / {arg_secure_result['n_blocks_requested']}"
    )
    if arg_secure_result["n_blocks_available"] == 0:
        print(
            "No block-level secure outputs were found under "
            f"{arg_ctx['run_root'] / 'party1'} "
            "(expected qBlock_block*.txt and qBurdenBlock_block*.txt)."
        )
    elif arg_secure_result["missing_blocks"]:
        preview = ",".join(str(block) for block in arg_secure_result["missing_blocks"][:8])
        suffix = "..." if len(arg_secure_result["missing_blocks"]) > 8 else ""
        print(f"Missing secure block raw outputs for blocks: {preview}{suffix}")

    comparisons = [
        ("SKAT", float(arg_manual_result["analysis_skat_q"]), float(arg_secure_result["summary_skat_q"])),
        ("Burden", float(arg_manual_result["analysis_burden_q"]), float(arg_secure_result["summary_burden_q"])),
    ]
    for label, plain_value, secure_value in comparisons:
        abs_diff = abs(plain_value - secure_value) if math.isfinite(plain_value) and math.isfinite(secure_value) else float("nan")
        rel_diff = safe_rel_diff(plain_value, secure_value)
        print(
            f"{label}: plain={format_float(plain_value)} "
            f"secure={format_float(secure_value)} "
            f"abs_diff={format_float(abs_diff)} "
            f"rel_diff={format_float(rel_diff)}"
        )


def print_compare_summary(
    arg_ctx: dict,
    arg_manual_result: dict,
    arg_secure_result: dict,
    arg_reference_result: dict,
    arg_summary_df: pd.DataFrame,
    arg_summary_path: Path,
    arg_png_paths: list[Path],
) -> None:
    print("\n--- Final Summary ---")
    if arg_reference_result["skipped_reason"] is None:
        print(f"Reference mode: block-wise SKAT package run ({arg_reference_result['n_blocks']} blocks)")
        print(f"Reference block summary TSV: {arg_reference_result['summary_path']}")
    else:
        print(f"Reference result: {arg_reference_result['skipped_reason']}")

    print_aggregate_comparison_summary(arg_ctx, arg_manual_result, arg_secure_result)

    if arg_summary_df.empty:
        print("No finite block-wise pairwise metrics were available to summarize.")

    for _, row in arg_summary_df.iterrows():
        print(
            f"{row['metric']} {row['comparison']}: "
            f"max_abs={format_float(float(row['max_abs_diff']))} "
            f"(block {int(row['max_abs_diff_block'])}), "
            f"max_rel={format_float(float(row['max_rel_diff']))} "
            f"(block {int(row['max_rel_diff_block'])}), "
            f"mean_abs={format_float(float(row['mean_abs_diff']))}, "
            f"mean_rel={format_float(float(row['mean_rel_diff']))}"
        )
    print(f"Summary CSV: {arg_summary_path}")
    for png_path in arg_png_paths:
        print(f"PNG: {png_path}")


def build_context(arg_ns: object) -> dict:
    repo_root = Path(arg_ns.repo_root).resolve()
    run_root = resolve_run_root_arg(repo_root, getattr(arg_ns, "run_id", None), getattr(arg_ns, "run_root", None))
    run_metadata = read_kv_file(run_root / "run_metadata.txt")
    dataset_info = resolve_dataset(repo_root, run_root, run_metadata, arg_ns.dataset)
    analysis_blocks = parse_block_spec(arg_ns.blocks, dataset_info["n_blocks"])
    plink2 = require_executable("plink2")
    all_positions = read_positions(dataset_info["party_dirs"][0], dataset_info["total_variants"])
    secure_qc_filter = read_qc_filter(run_root, dataset_info["total_variants"])

    y_parts = [read_pheno(party_dir) for party_dir in dataset_info["party_dirs"]]
    X_parts = []
    cov_n_cols = None
    for party_dir in dataset_info["party_dirs"]:
        X_part = np.asarray(read_cov(party_dir), dtype=float)
        if X_part.ndim == 1:
            X_part = X_part.reshape(1, -1)
        if cov_n_cols is None:
            cov_n_cols = int(X_part.shape[1])
        elif X_part.shape[1] != cov_n_cols:
            raise RuntimeError(
                "Covariate column-count mismatch across parties: "
                f"expected {cov_n_cols} columns but found {X_part.shape[1]} in {party_dir / 'cov.txt'}"
            )
        X_parts.append(X_part)
    model = compute.fit_null_model(np.vstack(X_parts).astype(float), np.concatenate(y_parts).astype(float))

    output_dir = run_root / "analysis"
    scratch_dir = output_dir / ".tmp"
    output_dir.mkdir(parents=True, exist_ok=True)
    scratch_dir.mkdir(parents=True, exist_ok=True)

    return {
        "command": arg_ns.command,
        "repo_root": repo_root,
        "run_id": getattr(arg_ns, "run_id", None) or run_metadata.get("run_id", run_root.name),
        "run_root": run_root,
        "run_metadata": run_metadata,
        "dataset_root": dataset_info["dataset_root"],
        "party_dirs": dataset_info["party_dirs"],
        "chrom_sizes": dataset_info["chrom_sizes"],
        "total_variants": dataset_info["total_variants"],
        "n_blocks": dataset_info["n_blocks"],
        "block_offsets": dataset_info["block_offsets"],
        "analysis_blocks": analysis_blocks,
        "plain_mode": getattr(arg_ns, "plain_mode", compute.PLAIN_MODE_STANDARD),
        "local_weight_mode": getattr(arg_ns, "local_weight_mode", compute.LOCAL_WEIGHT_MODE_DIRECT_TOTAL),
        "skip_reference": bool(getattr(arg_ns, "skip_reference", False)),
        "output_dir": output_dir,
        "scratch_dir": scratch_dir,
        "plink2": plink2,
        "all_positions": all_positions,
        "secure_qc_filter": secure_qc_filter,
        "model": model,
    }


def run_compare(arg_ctx: dict) -> int:
    # 1) Load the selected blocks plus the shared null-model context.
    print_preflight(arg_ctx)
    block_inputs = load_plain_blocks(arg_ctx)

    # 2) Compute the plain/manual statistics, then load the secure and R-package references.
    manual_result = compute.compute_manual_results(arg_ctx, block_inputs)
    secure_result = load_secure_compare(arg_ctx)
    reference_result = reference.run_reference(arg_ctx, manual_result)
    arg_ctx["reference_blocks"] = reference_result["blocks"]

    # 3) Assemble block-level comparison tables and write the main CSV artifact.
    block_df = build_block_compare_frame(arg_ctx, manual_result, secure_result)
    block_csv_path = arg_ctx["output_dir"] / "block_compare.csv"
    block_df.to_csv(block_csv_path, index=False)

    # 4) Render the block-level scatter plots used for quick visual comparison.
    png_paths = []
    skat_png_path = arg_ctx["output_dir"] / "block_compare_skat_scatter.png"
    burden_png_path = arg_ctx["output_dir"] / "block_compare_burden_scatter.png"
    if plotting.draw_scatter_png(
        block_df["plain_skat_q"].to_numpy(dtype=float),
        block_df["secure_skat_q"].to_numpy(dtype=float),
        skat_png_path,
        title="Block-Level SKAT Comparison",
        subtitle=(
            f"n = {count_finite_pairs(block_df['plain_skat_q'], block_df['secure_skat_q'])}, "
            f"corr = {safe_corr(block_df['plain_skat_q'], block_df['secure_skat_q']):.6f}, "
            f"r2 = {safe_r2(block_df['plain_skat_q'], block_df['secure_skat_q']):.6f}"
        ),
    ):
        png_paths.append(skat_png_path)
    if plotting.draw_scatter_png(
        block_df["plain_burden_q"].to_numpy(dtype=float),
        block_df["secure_burden_q"].to_numpy(dtype=float),
        burden_png_path,
        title="Block-Level Burden Comparison",
        subtitle=(
            f"n = {count_finite_pairs(block_df['plain_burden_q'], block_df['secure_burden_q'])}, "
            f"corr = {safe_corr(block_df['plain_burden_q'], block_df['secure_burden_q']):.6f}, "
            f"r2 = {safe_r2(block_df['plain_burden_q'], block_df['secure_burden_q']):.6f}"
        ),
    ):
        png_paths.append(burden_png_path)

    # 5) Write the compact summary table and print the final console report.
    summary_path, summary_df = write_summary_csv(arg_ctx, block_df)
    print_block_comparison_summary(block_df, block_csv_path, png_paths)
    print_compare_summary(arg_ctx, manual_result, secure_result, reference_result, summary_df, summary_path, png_paths)
    return 0


def run_reference_only(arg_ctx: dict) -> int:
    # Load the selected blocks, build the plain inputs, and hand them to the R package helper.
    print_preflight(arg_ctx)
    block_inputs = load_plain_blocks(arg_ctx)
    manual_result = compute.compute_manual_results(arg_ctx, block_inputs)
    reference_result = reference.run_reference(arg_ctx, manual_result)
    if reference_result["skipped_reason"] is None:
        print("\n--- Reference Summary ---")
        print(f"Reference blocks computed: {reference_result['n_blocks']}")
        print(f"Reference block summary TSV: {reference_result['summary_path']}")
        return 0
    raise RuntimeError(reference_result["skipped_reason"])
