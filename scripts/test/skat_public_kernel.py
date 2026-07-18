#!/usr/bin/env python3
"""Shared public-kernel loader for the SKAT profiling scripts.

Only public-list genotype blocks and covariates are read.  No phenotype,
private-variant block, manifest, or variant metadata is accessed here.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import numpy as np


RIDGE_REL_DEFAULT = 1e-6


@dataclass(frozen=True)
class PublicBlockStats:
    oriented_dosage: np.ndarray
    gtg: np.ndarray
    gtx: np.ndarray


@dataclass(frozen=True)
class PublicKernelSource:
    root: Path
    block_sizes: np.ndarray
    x_a: np.ndarray
    x_b: np.ndarray
    x_a_sum: np.ndarray
    x_b_sum: np.ndarray
    xtx_ridged: np.ndarray

    @classmethod
    def load(cls, root: Path, ridge_rel: float) -> "PublicKernelSource":
        block_sizes = np.loadtxt(root / "block_sizes.txt", dtype=np.int64, ndmin=1)
        if block_sizes.ndim != 1 or block_sizes.size == 0 or np.any(block_sizes <= 0):
            raise ValueError("invalid block_sizes.txt")

        x_a = load_covariates(root / "A" / "cov.txt")
        x_b = load_covariates(root / "B" / "cov.txt")
        if x_a.shape[1] != x_b.shape[1]:
            raise ValueError("covariate dimensions do not match")
        return cls(
            root,
            block_sizes,
            x_a,
            x_b,
            np.sum(x_a, axis=0),
            np.sum(x_b, axis=0),
            ridged_xtx(x_a, x_b, ridge_rel),
        )

    def selected_genes(self, spec: str | None) -> list[int]:
        return parse_gene_indices(spec, int(self.block_sizes.size))

    def kernels(
        self, genes: list[int], chunk_rows: int
    ) -> Iterator[tuple[int, int, np.ndarray]]:
        for gene_index in genes:
            size = int(self.block_sizes[gene_index])
            stats_a = public_block_stats(
                self.root / "A" / f"geno.{gene_index}.bin",
                self.x_a.shape[0],
                size,
                self.x_a,
                chunk_rows,
                self.x_a_sum,
            )
            stats_b = public_block_stats(
                self.root / "B" / f"geno.{gene_index}.bin",
                self.x_b.shape[0],
                size,
                self.x_b,
                chunk_rows,
                self.x_b_sum,
            )
            kernel = build_public_kernel(
                stats_a,
                stats_b,
                self.x_a.shape[0] + self.x_b.shape[0],
                self.xtx_ridged,
            )
            del stats_a, stats_b
            yield gene_index, size, kernel


def add_public_kernel_arguments(
    parser: argparse.ArgumentParser, *, default_genes: str
) -> None:
    parser.add_argument(
        "--fed-out",
        type=Path,
        default=Path(os.environ.get("FED_OUT", "~/fed_prep_out")).expanduser(),
        help="fed_prep.py output directory (default: FED_OUT or ~/fed_prep_out)",
    )
    parser.add_argument(
        "--genes",
        default=default_genes,
        help=(
            "zero-based public block indices, e.g. 0,2-5 or all "
            f"(default: {default_genes})"
        ),
    )
    parser.add_argument("--ridge-rel", type=float, default=RIDGE_REL_DEFAULT)
    parser.add_argument(
        "--chunk-rows",
        type=int,
        default=512,
        help="rows per streaming genotype chunk (default: 512)",
    )


def validate_public_kernel_arguments(args: argparse.Namespace) -> None:
    if not math.isfinite(args.ridge_rel) or args.ridge_rel < 0.0:
        raise ValueError("ridge-rel must be finite and nonnegative")
    if args.chunk_rows < 1:
        raise ValueError("chunk-rows must be at least one")


def emit(record: dict) -> None:
    """Emit one deterministic JSON record and nothing else on stdout."""
    print(json.dumps(record, sort_keys=True, separators=(",", ":"), allow_nan=False))


def finite(value: float) -> float | None:
    value = float(value)
    return value if math.isfinite(value) else None


def parse_gene_indices(spec: str | None, gene_count: int) -> list[int]:
    if spec is None or spec.strip().lower() == "all":
        return list(range(gene_count))

    selected: set[int] = set()
    for piece in spec.split(","):
        piece = piece.strip()
        if not piece:
            continue
        if "-" in piece:
            left, right = piece.split("-", 1)
            start, stop = int(left), int(right)
            if stop < start:
                raise ValueError("gene ranges must be increasing")
            selected.update(range(start, stop + 1))
        else:
            selected.add(int(piece))

    result = sorted(selected)
    if not result or result[0] < 0 or result[-1] >= gene_count:
        raise ValueError("gene index is outside block_sizes.txt")
    return result


def load_covariates(path: Path) -> np.ndarray:
    cov = np.loadtxt(path, dtype=np.float64, ndmin=2)
    if cov.ndim != 2 or cov.shape[0] == 0 or not np.all(np.isfinite(cov)):
        raise ValueError("invalid covariate matrix")
    return np.column_stack((np.ones(cov.shape[0], dtype=np.float64), cov))


def public_block_stats(
    path: Path,
    rows: int,
    columns: int,
    design: np.ndarray,
    chunk_rows: int,
    design_sum: np.ndarray | None = None,
) -> PublicBlockStats:
    """Compute oriented GtG/GtX in one file pass using sufficient statistics."""
    expected = rows * columns
    if path.stat().st_size != expected:
        raise ValueError("public genotype block has an unexpected size")
    raw = np.memmap(path, dtype=np.int8, mode="r", shape=(rows, columns))

    dosage = np.zeros(columns, dtype=np.float64)
    gtg = np.zeros((columns, columns), dtype=np.float64)
    gtx = np.zeros((columns, design.shape[1]), dtype=np.float64)
    for start in range(0, rows, chunk_rows):
        stop = min(start + chunk_rows, rows)
        chunk = np.asarray(raw[start:stop], dtype=np.float64)
        np.maximum(chunk, 0.0, out=chunk)
        dosage += np.sum(chunk, axis=0)
        gtg += chunk.T @ chunk
        gtx += chunk.T @ design[start:stop]
    del raw

    # Match orientGenotypeLocal: G' = G*S + 1*b^T, with S=-1 and b=2
    # exactly on columns whose local dosage is strictly greater than n.
    flip = dosage > rows
    sign = np.where(flip, -1.0, 1.0)
    offset = np.where(flip, 2.0, 0.0)
    signed_dosage = sign * dosage
    oriented_dosage = signed_dosage + rows * offset

    gtg *= sign[:, None]
    gtg *= sign[None, :]
    gtx *= sign[:, None]
    flipped = np.flatnonzero(flip)
    if flipped.size:
        gtg[:, flipped] += 2.0 * signed_dosage[:, None]
        gtg[flipped, :] += 2.0 * signed_dosage[None, :]
        gtg[np.ix_(flipped, flipped)] += 4.0 * rows
        if design_sum is None:
            design_sum = np.sum(design, axis=0)
        gtx[flipped, :] += 2.0 * design_sum[None, :]
    return PublicBlockStats(oriented_dosage=oriented_dosage, gtg=gtg, gtx=gtx)


def ridged_xtx(x_a: np.ndarray, x_b: np.ndarray, ridge_rel: float) -> np.ndarray:
    """Match the sum of the per-data-party ridges in computeBetaHatEnc."""
    xtx = x_a.T @ x_a + x_b.T @ x_b
    covariate_count = xtx.shape[0]
    if covariate_count > 1:
        epsilon = ridge_rel * float(np.trace(xtx[1:, 1:])) / covariate_count
        indices = np.arange(1, covariate_count)
        xtx[indices, indices] += epsilon
    return xtx


def build_public_kernel(
    stats_a: PublicBlockStats,
    stats_b: PublicBlockStats,
    n_total: int,
    xtx_ridged: np.ndarray,
) -> np.ndarray:
    """Construct the normalized public-public K used by the secure path."""
    dosage = stats_a.oriented_dosage + stats_b.oriented_dosage
    allele_frequency = dosage / (2.0 * n_total)
    weight = 25.0 * np.power(1.0 - allele_frequency, 24)

    # The per-gene sufficient statistics are consumed here, so aggregate the
    # large Gram matrix in place instead of allocating another m-by-m copy.
    gtx = stats_a.gtx.copy()
    gtx += stats_b.gtx
    projected = np.linalg.solve(xtx_ridged, gtx.T)
    residual_gram = stats_a.gtg
    residual_gram += stats_b.gtg
    residual_gram -= gtx @ projected
    residual_gram /= n_total
    residual_gram *= weight[:, None]
    residual_gram *= weight[None, :]
    # Current K has a leading 1/2.  Symmetrize away solve/BLAS roundoff.
    np.add(residual_gram, residual_gram.T, out=residual_gram)
    residual_gram *= 0.25
    return residual_gram
