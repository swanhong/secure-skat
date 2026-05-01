"""Readers for secure-output intermediates and run summaries."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np

from .common import trim_or_none
from .models import CompareContext, SecureBlockArtifacts, SecureSummary


def secure_party_dir(run_root: Path, party_idx: int) -> Path:
    return run_root / f"party{party_idx}"


def secure_cache_party_dir(run_root: Path, party_idx: int) -> Path:
    return run_root / "cache" / f"party{party_idx}"


def read_secure_matrix(path: Path) -> np.ndarray | None:
    if not path.exists():
        return None
    rows: list[list[float]] = []
    for line in path.read_text().splitlines():
        vals = [token.strip() for token in line.split(",")]
        row = [float(token) for token in vals if token]
        if row:
            rows.append(row)
    if not rows:
        return None
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise RuntimeError(f"Inconsistent secure matrix row widths in {path}")
    return np.asarray(rows, dtype=float)


def read_secure_vector(path: Path) -> np.ndarray | None:
    matrix = read_secure_matrix(path)
    if matrix is None:
        return None
    return matrix.reshape(-1)


def read_secure_vector_any(paths: Iterable[Path]) -> np.ndarray | None:
    for path in paths:
        vec = read_secure_vector(path)
        if vec is not None:
            return vec
    return None


def read_comma_numeric_vector(path: Path) -> np.ndarray | None:
    if not path.exists():
        return None
    values: list[float] = []
    for line in path.read_text().splitlines():
        for token in line.split(","):
            token = token.strip()
            if token:
                values.append(float(token))
    if not values:
        return None
    return np.asarray(values, dtype=float)


def load_secure_block_artifacts(
    ctx: CompareContext,
    block_index: int,
    n_variants: int,
) -> SecureBlockArtifacts | None:
    if ctx.run_root is None:
        return None

    secure_block_idx = block_index - 1
    run_root = ctx.run_root

    secure_local_dosage = [
        trim_or_none(
            read_secure_vector(
                secure_cache_party_dir(run_root, party_idx) / f"assoc_cache_dos_sum.skat.{secure_block_idx}.txt"
            ),
            n_variants,
        )
        for party_idx in (1, 2)
    ]

    secure_global_dosage_sum = None
    if all(vec is not None for vec in secure_local_dosage):
        secure_global_dosage_sum = secure_local_dosage[0] + secure_local_dosage[1]

    p_vec = trim_or_none(
        read_secure_vector_any(
            [
                secure_party_dir(run_root, 1) / f"p_block{secure_block_idx}.txt",
                secure_party_dir(run_root, 1) / f"p_enc_block{secure_block_idx}.txt",
                secure_party_dir(run_root, 1) / f"p_block{secure_block_idx}_real.txt",
                secure_party_dir(run_root, 1) / f"p_enc_block{secure_block_idx}_real.txt",
            ]
        ),
        n_variants,
    )
    p_vec_imag = trim_or_none(
        read_secure_vector_any(
            [
                secure_party_dir(run_root, 1) / f"p_block{secure_block_idx}_imag.txt",
                secure_party_dir(run_root, 1) / f"p_enc_block{secure_block_idx}_imag.txt",
            ]
        ),
        n_variants,
    )
    p_bar_vec = trim_or_none(
        read_secure_vector_any(
            [
                secure_party_dir(run_root, 1) / f"p_bar_block{secure_block_idx}.txt",
                secure_party_dir(run_root, 1) / f"p_bar_block{secure_block_idx}_real.txt",
            ]
        ),
        n_variants,
    )
    p_bar_vec_imag = trim_or_none(
        read_secure_vector_any(
            [secure_party_dir(run_root, 1) / f"p_bar_block{secure_block_idx}_imag.txt"]
        ),
        n_variants,
    )
    score_vec = trim_or_none(
        read_secure_vector(secure_party_dir(run_root, 1) / f"S_vec_block{secure_block_idx}.txt"),
        n_variants,
    )
    weight_vec = trim_or_none(
        read_secure_vector_any(
            [
                secure_party_dir(run_root, 1) / f"w_enc_block{secure_block_idx}.txt",
                secure_party_dir(run_root, 1) / f"w_enc_block{secure_block_idx}_real.txt",
            ]
        ),
        n_variants,
    )
    weight_vec_imag = trim_or_none(
        read_secure_vector(secure_party_dir(run_root, 1) / f"w_enc_block{secure_block_idx}_imag.txt"),
        n_variants,
    )
    score_sq_vec = trim_or_none(
        read_secure_vector(secure_party_dir(run_root, 1) / f"S2_block{secure_block_idx}.txt"),
        n_variants,
    )
    weight_sq_vec = trim_or_none(
        read_secure_vector_any(
            [
                secure_party_dir(run_root, 1) / f"w2_block{secure_block_idx}.txt",
                secure_party_dir(run_root, 1) / f"w2_block{secure_block_idx}_real.txt",
            ]
        ),
        n_variants,
    )
    weight_sq_vec_imag = trim_or_none(
        read_secure_vector(secure_party_dir(run_root, 1) / f"w2_block{secure_block_idx}_imag.txt"),
        n_variants,
    )
    w2s2_vec = trim_or_none(
        read_secure_vector_any(
            [
                secure_party_dir(run_root, 1) / f"w2S2_block{secure_block_idx}.txt",
                secure_party_dir(run_root, 1) / f"w2S2_block{secure_block_idx}_real.txt",
            ]
        ),
        n_variants,
    )
    w2s2_vec_imag = trim_or_none(
        read_secure_vector(secure_party_dir(run_root, 1) / f"w2S2_block{secure_block_idx}_imag.txt"),
        n_variants,
    )
    wS_vec = trim_or_none(
        read_secure_vector(secure_party_dir(run_root, 1) / f"wS_block{secure_block_idx}.txt"),
        n_variants,
    )
    q_block_vec = trim_or_none(
        read_secure_vector(secure_party_dir(run_root, 1) / f"qBlock_block{secure_block_idx}.txt"),
        1,
    )
    burden_vec = trim_or_none(
        read_secure_vector(secure_party_dir(run_root, 1) / f"qBurdenBlock_block{secure_block_idx}.txt"),
        1,
    )

    return SecureBlockArtifacts(
        local_dosage_sum=secure_local_dosage,
        secure_global_dosage_sum=secure_global_dosage_sum,
        p_vec=p_vec,
        p_vec_imag=p_vec_imag,
        p_bar_vec=p_bar_vec,
        p_bar_vec_imag=p_bar_vec_imag,
        score_vec=score_vec,
        weight_vec=weight_vec,
        weight_vec_imag=weight_vec_imag,
        score_sq_vec=score_sq_vec,
        weight_sq_vec=weight_sq_vec,
        weight_sq_vec_imag=weight_sq_vec_imag,
        w2s2_vec=w2s2_vec,
        w2s2_vec_imag=w2s2_vec_imag,
        wS_vec=wS_vec,
        q_skat_block_raw=None if q_block_vec is None else float(q_block_vec[0]),
        q_burden_block_raw=None if burden_vec is None else float(burden_vec[0]),
    )


def load_secure_summary(ctx: CompareContext) -> SecureSummary | None:
    if ctx.run_root is None:
        return None

    party1_dir = secure_party_dir(ctx.run_root, 1)
    skat_out_path = party1_dir / "skat_out.txt"
    burden_out_path = party1_dir / "burden_out.txt"

    secure_run_skat_q = (
        float(np.loadtxt(skat_out_path, dtype=float, max_rows=1))
        if skat_out_path.exists()
        else float("nan")
    )
    secure_run_burden_q = (
        float(np.loadtxt(burden_out_path, dtype=float, max_rows=1))
        if burden_out_path.exists()
        else float("nan")
    )

    selected_q_skat_raw_total = 0.0
    selected_q_burden_raw_total = 0.0
    selected_available = True
    for block_index in ctx.analysis_blocks:
        secure_block_idx = block_index - 1
        q_vec = trim_or_none(read_secure_vector(party1_dir / f"qBlock_block{secure_block_idx}.txt"), 1)
        burden_vec = trim_or_none(read_secure_vector(party1_dir / f"qBurdenBlock_block{secure_block_idx}.txt"), 1)
        if q_vec is None or burden_vec is None:
            selected_available = False
            break
        selected_q_skat_raw_total += float(q_vec[0])
        selected_q_burden_raw_total += float(burden_vec[0])

    selected_skat_q = (
        selected_q_skat_raw_total * ctx.model.rare_variant_scale if selected_available else float("nan")
    )
    selected_burden_q = (
        (selected_q_burden_raw_total**2) * ctx.model.rare_variant_scale
        if selected_available
        else float("nan")
    )

    used_run_scalars = not selected_available
    secure_skat_q_for_summary = secure_run_skat_q if used_run_scalars else selected_skat_q
    secure_burden_q_for_summary = secure_run_burden_q if used_run_scalars else selected_burden_q

    return SecureSummary(
        secure_run_skat_q=secure_run_skat_q,
        secure_run_burden_q=secure_run_burden_q,
        selected_q_skat_raw_total=selected_q_skat_raw_total,
        selected_q_burden_raw_total=selected_q_burden_raw_total,
        selected_skat_q=selected_skat_q,
        selected_burden_q=selected_burden_q,
        selected_available=selected_available,
        secure_skat_q_for_summary=secure_skat_q_for_summary,
        secure_burden_q_for_summary=secure_burden_q_for_summary,
        used_run_scalars=used_run_scalars,
    )
