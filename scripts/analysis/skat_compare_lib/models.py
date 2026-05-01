"""Dataclasses shared across the SKAT comparison modules."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class RunHints:
    block_count: int | None
    total_variants: int | None


@dataclass
class DatasetInfo:
    root: Path
    party_dirs: list[Path]
    chrom_sizes: list[int]
    total_variants: int
    n_blocks: int
    block_offsets: list[tuple[int, int]]


@dataclass
class ModelInputs:
    y: np.ndarray
    x: np.ndarray
    design: np.ndarray
    y_resid: np.ndarray
    q_matrix: np.ndarray
    null_rss: float
    dof: int
    rare_variant_scale: float
    null_model_s2: float
    n_total: int
    party_sample_counts: list[int]


@dataclass
class CompareContext:
    command: str
    repo_root: Path
    run_id: str | None
    run_root: Path | None
    run_metadata: dict[str, str]
    dataset: DatasetInfo
    analysis_blocks: list[int]
    detail_blocks: list[int]
    debug_mode: bool
    window_bp: int | None
    step_bp: int | None
    min_window_variants: int
    window_limit: int | None
    window_output_tag: str | None
    skip_reference: bool
    cache_dir: Path
    csv_dir: Path
    plink2: str
    all_positions: np.ndarray
    secure_qc_filter: np.ndarray | None
    model: ModelInputs
    run_hints: RunHints


@dataclass
class ExportedBlock:
    raw_path: Path
    geno: np.ndarray
    variant_ids: list[str]


@dataclass
class BlockGenotypeInput:
    block_index: int
    raw_paths_by_party: list[Path]
    local_alt_genotypes_by_party: list[np.ndarray]
    variant_ids: list[str]
    positions: np.ndarray


@dataclass
class DatasetInputs:
    block_inputs: list[BlockGenotypeInput]


@dataclass
class SecureBlockArtifacts:
    local_dosage_sum: list[np.ndarray | None]
    secure_global_dosage_sum: np.ndarray | None
    p_vec: np.ndarray | None
    p_vec_imag: np.ndarray | None
    p_bar_vec: np.ndarray | None
    p_bar_vec_imag: np.ndarray | None
    score_vec: np.ndarray | None
    weight_vec: np.ndarray | None
    weight_vec_imag: np.ndarray | None
    score_sq_vec: np.ndarray | None
    weight_sq_vec: np.ndarray | None
    weight_sq_vec_imag: np.ndarray | None
    w2s2_vec: np.ndarray | None
    w2s2_vec_imag: np.ndarray | None
    wS_vec: np.ndarray | None
    q_skat_block_raw: float | None
    q_burden_block_raw: float | None


@dataclass
class ManualBlockResult:
    block_index: int
    raw_paths_by_party: list[Path]
    n_variants: int
    variant_ids: list[str]
    positions: np.ndarray
    local_alt_dosage_sum_vecs: list[np.ndarray]
    local_secure_dosage_sum_vecs: list[np.ndarray]
    alt_dosage_sum_vec: np.ndarray
    secure_dosage_sum_vec: np.ndarray
    p_vec: np.ndarray
    p_bar_vec: np.ndarray
    beta_base_vec: np.ndarray
    weight_vec: np.ndarray
    score_vec: np.ndarray
    score_sq_vec: np.ndarray
    weight_sq_vec: np.ndarray
    w2s2_vec: np.ndarray
    wS_vec: np.ndarray
    q_skat_block_raw: float
    q_burden_block_raw: float
    secure_artifacts: SecureBlockArtifacts | None = None


@dataclass
class ManualResults:
    blocks: list[ManualBlockResult]
    variant_total: int
    analysis_variant_total: int
    analysis_q_skat_raw_total: float
    analysis_skat_q: float
    analysis_q_burden_raw_total: float
    analysis_burden_q: float
    all_q_skat_raw_total: float
    all_skat_q: float
    all_q_burden_raw_total: float
    all_burden_q: float


@dataclass
class SecureSummary:
    secure_run_skat_q: float
    secure_run_burden_q: float
    selected_q_skat_raw_total: float
    selected_q_burden_raw_total: float
    selected_skat_q: float
    selected_burden_q: float
    selected_available: bool
    secure_skat_q_for_summary: float
    secure_burden_q_for_summary: float
    used_run_scalars: bool


@dataclass
class ReferenceResult:
    skat_q: float
    burden_q: float
    n_markers: int
    summary_path: Path | None
    skipped_reason: str | None = None

    @property
    def available(self) -> bool:
        return math.isfinite(self.skat_q) and math.isfinite(self.burden_q)


@dataclass
class BlockComparisonRow:
    block: int
    n_variants: int
    block_start_bp: int | float
    block_end_bp: int | float
    start_variant_id: str
    end_variant_id: str
    plain_skat_q: float
    secure_skat_q: float
    skat_abs_diff: float
    skat_rel_diff: float
    plain_burden_q: float
    secure_burden_q: float
    burden_abs_diff: float
    burden_rel_diff: float
    plain_burden_sum: float
    secure_burden_sum: float


@dataclass
class WindowComparisonRow:
    window_index: int
    block: int
    block_window_index: int
    window_start_bp: int
    window_end_bp: int
    n_variants: int
    start_variant_id: str
    end_variant_id: str
    plain_skat_q: float
    secure_skat_q: float
    skat_abs_diff: float
    skat_rel_diff: float
    plain_burden_q: float
    secure_burden_q: float
    burden_abs_diff: float
    burden_rel_diff: float
    plain_burden_sum: float
    secure_burden_sum: float
