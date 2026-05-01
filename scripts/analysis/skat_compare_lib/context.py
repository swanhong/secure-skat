"""Context building, dataset resolution, and model-input loading."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Sequence

import numpy as np

from .common import count_nonblank_lines, read_kv_file, require_executable, sanitize_path_tag
from .dataset_io import read_chrom_sizes, read_cov, read_pheno, read_positions, read_qc_filter
from .models import CompareContext, DatasetInfo, ModelInputs, RunHints


def parse_block_spec(spec: str | None, n_blocks: int, *, default_all: bool) -> list[int]:
    if spec is None:
        values = list(range(1, n_blocks + 1)) if default_all else [1, n_blocks]
        return sorted(set(values))

    raw_parts = [part.strip() for part in spec.split(",") if part.strip()]
    if not raw_parts:
        raise RuntimeError("Empty block specification")

    out: set[int] = set()
    for part in raw_parts:
        if part == "all":
            out.update(range(1, n_blocks + 1))
            continue
        if part == "last":
            out.add(n_blocks)
            continue
        if re.fullmatch(r"\d+-\d+", part):
            start_s, end_s = part.split("-", 1)
            start_v = int(start_s)
            end_v = int(end_s)
            if start_v > end_v:
                raise RuntimeError(f"Invalid block range: {part}")
            out.update(range(start_v, end_v + 1))
            continue
        try:
            out.add(int(part))
        except ValueError as exc:
            raise RuntimeError(f"Invalid block specification token: {part}") from exc

    filtered = sorted(block for block in out if 1 <= block <= n_blocks)
    if not filtered:
        raise RuntimeError("No valid blocks remain after applying the block specification")
    return filtered


def compute_block_offsets(chrom_sizes: Sequence[int]) -> list[tuple[int, int]]:
    offsets: list[tuple[int, int]] = []
    start = 0
    for size in chrom_sizes:
        end = start + int(size)
        offsets.append((start, end))
        start = end
    return offsets


def candidate_dataset_roots(repo_root: Path) -> list[Path]:
    candidates: list[Path] = []
    example_data = repo_root / "example_data"
    if (example_data / "party1" / "chrom_sizes.txt").exists():
        candidates.append(example_data)

    datasets_dir = repo_root / ".local" / "datasets"
    if datasets_dir.exists():
        for child in sorted(datasets_dir.iterdir()):
            if (child / "party1" / "chrom_sizes.txt").exists():
                candidates.append(child)
    return candidates


def resolve_run_root(repo_root: Path, run_id: str | None) -> Path | None:
    if not run_id:
        return None
    candidates = [
        path
        for path in (repo_root / "out").glob(f"output_*_{run_id}")
        if path.is_dir()
    ]
    if not candidates:
        raise RuntimeError(f"No secure output directory found for run id: {run_id}")
    candidates.sort(key=lambda path: path.stat().st_mtime, reverse=True)
    return candidates[0].resolve()


def resolve_dataset_path(repo_root: Path, dataset_value: str) -> Path:
    dataset_path = Path(dataset_value)
    if not dataset_path.is_absolute():
        dataset_path = repo_root / dataset_path
    dataset_path = dataset_path.resolve()
    if not (dataset_path / "party1" / "chrom_sizes.txt").exists():
        raise RuntimeError(f"Missing dataset directory or chrom_sizes.txt: {dataset_path}")
    return dataset_path


def infer_run_hints(repo_root: Path, run_root: Path | None, run_metadata: dict[str, str]) -> RunHints:
    if run_root is None:
        return RunHints(block_count=None, total_variants=None)

    block_count: int | None = None
    total_variants: int | None = None

    metadata_dataset = run_metadata.get("dataset")
    if metadata_dataset:
        dataset_root = resolve_dataset_path(repo_root, metadata_dataset)
        chrom_sizes = read_chrom_sizes(dataset_root)
        block_count = len(chrom_sizes)
        total_variants = sum(chrom_sizes)

    if block_count is None:
        qblock_files = sorted((run_root / "party1").glob("qBlock_block*.txt"))
        if qblock_files:
            block_count = len(qblock_files)
        else:
            cache_files = sorted((run_root / "cache" / "party1").glob("assoc_cache_dos_sum.skat.*.txt"))
            if cache_files:
                block_count = len(cache_files)

    if total_variants is None:
        qc_path = run_root / "cache" / "party1" / "gkeep.txt"
        if qc_path.exists():
            total_variants = count_nonblank_lines(qc_path)

    return RunHints(block_count=block_count, total_variants=total_variants)


def validate_dataset_against_run(dataset: DatasetInfo, run_hints: RunHints, run_root: Path | None) -> None:
    if run_root is None:
        return
    if run_hints.block_count is not None and dataset.n_blocks != run_hints.block_count:
        raise RuntimeError(
            f"Dataset/run block mismatch: dataset has {dataset.n_blocks} blocks, "
            f"but {run_root.name} implies {run_hints.block_count} blocks"
        )
    if run_hints.total_variants is not None and dataset.total_variants != run_hints.total_variants:
        raise RuntimeError(
            f"Dataset/run variant-count mismatch: dataset has {dataset.total_variants} variants, "
            f"but {run_root.name} implies {run_hints.total_variants} variants"
        )


def resolve_dataset(
    repo_root: Path,
    dataset_arg: str | None,
    run_root: Path | None,
    run_metadata: dict[str, str],
    run_hints: RunHints,
) -> DatasetInfo:
    dataset_root: Path | None = None

    if dataset_arg:
        dataset_root = resolve_dataset_path(repo_root, dataset_arg)
    elif run_metadata.get("dataset"):
        dataset_root = resolve_dataset_path(repo_root, run_metadata["dataset"])
    else:
        if run_root is None:
            raise RuntimeError("Manual/reference mode without --dataset requires an explicit dataset path")
        if run_hints.block_count is None or run_hints.total_variants is None:
            raise RuntimeError(
                "Unable to infer a unique dataset from the secure run; "
                "provide --dataset explicitly"
            )
        matches = []
        for candidate in candidate_dataset_roots(repo_root):
            chrom_sizes = read_chrom_sizes(candidate)
            if len(chrom_sizes) == run_hints.block_count and sum(chrom_sizes) == run_hints.total_variants:
                matches.append(candidate.resolve())
        if len(matches) != 1:
            raise RuntimeError(
                "Dataset inference failed: expected exactly one local dataset match for "
                f"{run_hints.block_count} blocks and {run_hints.total_variants} variants, "
                f"found {len(matches)}"
            )
        dataset_root = matches[0]

    chrom_sizes = read_chrom_sizes(dataset_root)
    dataset = DatasetInfo(
        root=dataset_root.resolve(),
        party_dirs=[(dataset_root / "party1").resolve(), (dataset_root / "party2").resolve()],
        chrom_sizes=chrom_sizes,
        total_variants=sum(chrom_sizes),
        n_blocks=len(chrom_sizes),
        block_offsets=compute_block_offsets(chrom_sizes),
    )
    validate_dataset_against_run(dataset, run_hints, run_root)
    return dataset


def load_model_inputs(party_dirs: Sequence[Path], cov_num_cols: int = 4) -> ModelInputs:
    party_pheno = [read_pheno(path) for path in party_dirs]
    party_cov = [read_cov(path) for path in party_dirs]

    party_sample_counts = [int(len(pheno)) for pheno in party_pheno]
    y = np.concatenate(party_pheno).astype(float)

    cov_rows: list[np.ndarray] = []
    for cov in party_cov:
        cov = np.asarray(cov, dtype=float)
        if cov.ndim == 1:
            cov = cov.reshape(1, -1)
        if cov.shape[1] < cov_num_cols:
            raise RuntimeError("Covariate file has fewer than 4 columns")
        cov_rows.append(cov[:, :cov_num_cols])
    x = np.vstack(cov_rows).astype(float)

    design = np.column_stack([np.ones(x.shape[0], dtype=float), x])
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    y_resid = y - design @ beta
    if np.isnan(y_resid).any():
        raise RuntimeError("Residual calculation produced NA values")

    q_matrix, _ = np.linalg.qr(design, mode="reduced")
    null_rss = float(np.sum(y_resid**2))
    dof = int(design.shape[0] - design.shape[1])
    if dof <= 0:
        raise RuntimeError("Null-model degrees of freedom must be positive")
    null_model_s2 = float(null_rss / dof)
    rare_variant_scale = float(dof / (2.0 * null_rss))

    return ModelInputs(
        y=y,
        x=x,
        design=design,
        y_resid=y_resid,
        q_matrix=q_matrix,
        null_rss=null_rss,
        dof=dof,
        rare_variant_scale=rare_variant_scale,
        null_model_s2=null_model_s2,
        n_total=int(y_resid.size),
        party_sample_counts=party_sample_counts,
    )


def build_context(args: argparse.Namespace) -> CompareContext:
    repo_root = Path(args.repo_root).resolve()
    run_root = resolve_run_root(repo_root, args.run_id)
    run_metadata = read_kv_file(run_root / "run_metadata.txt") if run_root else {}
    run_hints = infer_run_hints(repo_root, run_root, run_metadata)
    dataset = resolve_dataset(repo_root, args.dataset, run_root, run_metadata, run_hints)

    analysis_blocks = parse_block_spec(args.blocks, dataset.n_blocks, default_all=True)
    detail_blocks = parse_block_spec(args.detail_blocks, dataset.n_blocks, default_all=False)

    if args.window_bp is not None and args.window_bp <= 0:
        raise RuntimeError("--window-bp must be a positive integer")
    if args.step_bp is None and args.window_bp is not None:
        args.step_bp = args.window_bp
    if args.step_bp is not None and args.step_bp <= 0:
        raise RuntimeError("--step-bp must be a positive integer")
    if args.min_window_variants <= 0:
        raise RuntimeError("--min-window-variants must be a positive integer")
    if args.window_limit is not None and args.window_limit <= 0:
        raise RuntimeError("--window-limit must be a positive integer")

    needs_plink = args.command in {"compare", "manual", "reference"}
    plink2 = require_executable("plink2") if needs_plink else ""
    all_positions = read_positions(dataset.party_dirs[0], dataset.total_variants)
    secure_qc_filter = read_qc_filter(run_root, dataset.total_variants)
    model = load_model_inputs(dataset.party_dirs)

    run_tag = args.run_id or (run_root.name if run_root else "manual")
    cache_dir = (
        repo_root
        / ".local"
        / "tmp"
        / "skat_compare"
        / sanitize_path_tag(dataset.root)
        / sanitize_path_tag(run_tag)
    )
    cache_dir.mkdir(parents=True, exist_ok=True)
    csv_dir = cache_dir / "variant_debug_csv"
    if args.debug:
        csv_dir.mkdir(parents=True, exist_ok=True)

    return CompareContext(
        command=args.command,
        repo_root=repo_root,
        run_id=args.run_id,
        run_root=run_root,
        run_metadata=run_metadata,
        dataset=dataset,
        analysis_blocks=analysis_blocks,
        detail_blocks=detail_blocks,
        debug_mode=args.debug,
        window_bp=args.window_bp,
        step_bp=args.step_bp,
        min_window_variants=args.min_window_variants,
        window_limit=args.window_limit,
        window_output_tag=args.window_output_tag,
        skip_reference=args.skip_reference,
        cache_dir=cache_dir,
        csv_dir=csv_dir,
        plink2=plink2,
        all_positions=all_positions,
        secure_qc_filter=secure_qc_filter,
        model=model,
        run_hints=run_hints,
    )
