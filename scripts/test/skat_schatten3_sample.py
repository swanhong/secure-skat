#!/usr/bin/env python3
"""Validate a few-entry Schatten-3 estimator on prepared federated AoU data.

The script intentionally reads only public-variant genotype blocks and
covariates produced by scripts/preprocessing/fed_prep.py.  It never reads
phenotypes, private-variant blocks, variant identifiers, positions, or the
private fields in manifest.json.  Normal stdout is JSONL containing only
gene-level aggregate trace components and estimator-error summaries.

Run this inside the AoU Controlled Tier.  The aggregate output is still a
research result and must go through the applicable AoU export review.

Example:
    FED_OUT=/path/to/fed_prep_out python3 scripts/test/skat_schatten3_sample.py \
        --genes 0-19 --constants 1,2,5 --seeds 200
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

try:
    from scipy import sparse
except ImportError as exc:  # pragma: no cover - exercised on the server
    raise SystemExit("scipy is required") from exc


RIDGE_REL_DEFAULT = 1e-6
TINY = np.finfo(np.float64).tiny


@dataclass(frozen=True)
class ExactMoments:
    s1: float
    s2: float
    s3: float
    d3: float
    w3: float
    t3: float


@dataclass(frozen=True)
class PublicBlockStats:
    oriented_dosage: np.ndarray
    gtg: np.ndarray
    gtx: np.ndarray


def emit(record: dict) -> None:
    """Emit one deterministic JSON record and nothing else on stdout."""
    print(json.dumps(record, sort_keys=True, separators=(",", ":"), allow_nan=False))


def finite(value: float) -> float | None:
    value = float(value)
    return value if math.isfinite(value) else None


def parse_positive_floats(spec: str) -> list[float]:
    values = sorted(
        set(float(piece.strip()) for piece in spec.split(",") if piece.strip())
    )
    if not values or any(value <= 0.0 for value in values):
        raise ValueError("constants must be positive")
    return values


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
    if cov.ndim != 2 or cov.shape[0] == 0:
        raise ValueError("invalid covariate matrix")
    return np.column_stack((np.ones(cov.shape[0], dtype=np.float64), cov))


def public_block_stats(
    path: Path,
    rows: int,
    columns: int,
    design: np.ndarray,
    chunk_rows: int,
) -> PublicBlockStats:
    """Stream one public block and retain only aggregate sufficient statistics."""
    expected = rows * columns
    if path.stat().st_size != expected:
        raise ValueError("public genotype block has an unexpected size")
    raw = np.memmap(path, dtype=np.int8, mode="r", shape=(rows, columns))

    dosage = np.zeros(columns, dtype=np.int64)
    for start in range(0, rows, chunk_rows):
        stop = min(start + chunk_rows, rows)
        chunk = np.maximum(raw[start:stop], 0)
        dosage += np.sum(chunk, axis=0, dtype=np.int64)

    # Match orientGenotypeLocal: flip when local dosage is strictly > n.
    flip = dosage > rows
    oriented_dosage = np.where(flip, 2 * rows - dosage, dosage).astype(np.float64)
    gtg = np.zeros((columns, columns), dtype=np.float64)
    gtx = np.zeros((columns, design.shape[1]), dtype=np.float64)

    for start in range(0, rows, chunk_rows):
        stop = min(start + chunk_rows, rows)
        chunk = np.maximum(raw[start:stop], 0).astype(np.float64)
        chunk[:, flip] = 2.0 - chunk[:, flip]
        gtg += chunk.T @ chunk
        gtx += chunk.T @ design[start:stop]

    del raw
    return PublicBlockStats(oriented_dosage=oriented_dosage, gtg=gtg, gtx=gtx)


def ridged_xtx(x_a: np.ndarray, x_b: np.ndarray, ridge_rel: float) -> np.ndarray:
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

    gtg = stats_a.gtg + stats_b.gtg
    gtx = stats_a.gtx + stats_b.gtx
    projected = np.linalg.solve(xtx_ridged, gtx.T)
    residual_gram = (gtg - gtx @ projected) / n_total
    kernel = 0.5 * (weight[:, None] * residual_gram) * weight[None, :]
    return 0.5 * (kernel + kernel.T)


def exact_moments(kernel: np.ndarray) -> ExactMoments:
    diagonal = np.diag(kernel)
    upper_i, upper_j = np.triu_indices(kernel.shape[0], k=1)
    upper = kernel[upper_i, upper_j]

    d3 = float(np.sum(diagonal**3))
    w3 = float(
        3.0 * np.sum((diagonal[upper_i] + diagonal[upper_j]) * upper**2)
    )
    s1 = float(np.sum(diagonal))
    s2 = float(np.sum(kernel * kernel))

    squared = kernel @ kernel
    s3 = float(np.sum(squared * kernel.T))
    t3 = s3 - d3 - w3
    return ExactMoments(s1=s1, s2=s2, s3=s3, d3=d3, w3=w3, t3=t3)


def q(values: Sequence[float], probability: float) -> float:
    return float(np.quantile(np.asarray(values, dtype=np.float64), probability))


def summarize_errors(estimates: Sequence[float], reference: float, scale: float) -> dict:
    errors = np.asarray(estimates, dtype=np.float64) - reference
    normalized = errors / max(abs(scale), TINY)
    absolute_normalized = np.abs(normalized)
    return {
        "bias_over_scale": finite(float(np.mean(normalized))),
        "error_over_scale_abs_max": finite(float(np.max(absolute_normalized))),
        "error_over_scale_abs_p50": finite(q(absolute_normalized, 0.50)),
        "error_over_scale_abs_p95": finite(q(absolute_normalized, 0.95)),
        "estimate_mean": finite(float(np.mean(estimates))),
        "estimate_std": (
            finite(float(np.std(estimates, ddof=1))) if len(estimates) > 1 else 0.0
        ),
    }


def sampled_triangle_trace(
    size: int,
    upper_i: np.ndarray,
    upper_j: np.ndarray,
    upper_values: np.ndarray,
    selected: np.ndarray,
) -> float:
    i = upper_i[selected]
    j = upper_j[selected]
    values = upper_values[selected]
    rows = np.concatenate((i, j))
    columns = np.concatenate((j, i))
    data = np.concatenate((values, values))
    sampled = sparse.csr_matrix((data, (rows, columns)), shape=(size, size))
    sampled_squared = sampled @ sampled
    return float(sampled_squared.multiply(sampled).sum())


def estimate_grid(
    kernel: np.ndarray,
    exact: ExactMoments,
    constants: Sequence[float],
    seeds: int,
    base_seed: int,
) -> Iterable[dict]:
    size = kernel.shape[0]
    diagonal = np.diag(kernel)
    upper_i, upper_j = np.triu_indices(size, k=1)
    upper_values = kernel[upper_i, upper_j]
    upper_squares = upper_values**2
    upper_wedges = (diagonal[upper_i] + diagonal[upper_j]) * upper_squares
    edge_count = upper_values.size
    d2 = float(np.sum(diagonal**2))

    probabilities = [
        min(1.0, constant * size ** (-2.0 / 3.0)) for constant in constants
    ]
    accumulators = {
        constant: {"edges": [], "s2": [], "w3": [], "t3": [], "s3": []}
        for constant in constants
    }

    seed_sequence = np.random.SeedSequence(base_seed)
    for child_seed in seed_sequence.spawn(seeds):
        rng = np.random.default_rng(child_seed)
        uniforms = rng.random(edge_count)
        for constant, probability in zip(constants, probabilities):
            accumulator = accumulators[constant]
            if probability >= 1.0:
                selected = np.ones(edge_count, dtype=bool)
                s2_hat = exact.s2
                w3_hat = exact.w3
                t3_hat = exact.t3
            else:
                selected = uniforms < probability
                s2_hat = d2 + (2.0 / probability) * float(np.sum(upper_squares[selected]))
                w3_hat = (3.0 / probability) * float(np.sum(upper_wedges[selected]))
                triangle_trace = sampled_triangle_trace(
                    size, upper_i, upper_j, upper_values, selected
                )
                t3_hat = triangle_trace / probability**3

            accumulator["edges"].append(int(np.count_nonzero(selected)))
            accumulator["s2"].append(s2_hat)
            accumulator["w3"].append(w3_hat)
            accumulator["t3"].append(t3_hat)
            accumulator["s3"].append(exact.d3 + w3_hat + t3_hat)

    triangle_count = math.comb(size, 3) if size >= 3 else 0
    for constant, probability in zip(constants, probabilities):
        accumulator = accumulators[constant]
        yield {
            "constant": constant,
            "expected_sampled_edges": finite(edge_count * probability),
            "expected_sampled_triangles": finite(triangle_count * probability**3),
            "negative_s3_rate": finite(
                float(np.mean(np.asarray(accumulator["s3"], dtype=np.float64) <= 0.0))
            ),
            "p": probability,
            "sampled_edges_mean": finite(float(np.mean(accumulator["edges"]))),
            "s2": summarize_errors(accumulator["s2"], exact.s2, exact.s2),
            "s3": summarize_errors(accumulator["s3"], exact.s3, exact.s3),
            "t3": summarize_errors(accumulator["t3"], exact.t3, exact.s3),
            "w3": summarize_errors(accumulator["w3"], exact.w3, exact.s3),
        }


def ratio(value: float, denominator: float) -> float | None:
    if denominator == 0.0:
        return None
    return finite(value / denominator)


def run_self_test() -> None:
    rng = np.random.default_rng(7)

    raw_geno = rng.integers(-1, 3, size=(17, 9), dtype=np.int8)
    design = np.column_stack((np.ones(17), rng.normal(size=(17, 3))))
    with tempfile.TemporaryDirectory() as temporary_directory:
        block_path = Path(temporary_directory) / "geno.bin"
        raw_geno.tofile(block_path)
        stats = public_block_stats(block_path, 17, 9, design, chunk_rows=5)

    oriented = np.maximum(raw_geno, 0).astype(np.float64)
    flip = oriented.sum(axis=0) > oriented.shape[0]
    oriented[:, flip] = 2.0 - oriented[:, flip]
    if not np.allclose(stats.oriented_dosage, oriented.sum(axis=0)):
        raise AssertionError("streamed dosage/orientation failed")
    if not np.allclose(stats.gtg, oriented.T @ oriented):
        raise AssertionError("streamed GtG failed")
    if not np.allclose(stats.gtx, oriented.T @ design):
        raise AssertionError("streamed GtX failed")

    factor = rng.normal(size=(13, 7))
    kernel = factor @ factor.T / factor.shape[1]
    exact = exact_moments(kernel)

    direct_s3 = float(np.trace(kernel @ kernel @ kernel))
    if not np.isclose(exact.s3, direct_s3, rtol=1e-12, atol=1e-12):
        raise AssertionError("exact S3 computation failed")
    if not np.isclose(
        exact.s3, exact.d3 + exact.w3 + exact.t3, rtol=1e-12, atol=1e-12
    ):
        raise AssertionError("D3/W3/T3 decomposition failed")

    # Exhaust all edge subsets of a 3x3 kernel to verify unbiased scaling.
    small = kernel[:3, :3]
    small_exact = exact_moments(small)
    diagonal = np.diag(small)
    upper_i, upper_j = np.triu_indices(3, k=1)
    upper = small[upper_i, upper_j]
    probability = 0.37
    expected = np.zeros(3, dtype=np.float64)
    for mask_bits in range(1 << upper.size):
        selected = np.array(
            [bool(mask_bits & (1 << edge)) for edge in range(upper.size)]
        )
        selected_count = int(np.count_nonzero(selected))
        mass = probability**selected_count * (1.0 - probability) ** (
            upper.size - selected_count
        )
        s2_hat = float(np.sum(diagonal**2)) + (2.0 / probability) * float(
            np.sum(upper[selected] ** 2)
        )
        w3_hat = (3.0 / probability) * float(
            np.sum(
                (diagonal[upper_i[selected]] + diagonal[upper_j[selected]])
                * upper[selected] ** 2
            )
        )
        t3_hat = (
            sampled_triangle_trace(3, upper_i, upper_j, upper, selected)
            / probability**3
        )
        expected += mass * np.array((s2_hat, w3_hat, t3_hat))
    if not np.allclose(
        expected,
        (small_exact.s2, small_exact.w3, small_exact.t3),
        rtol=1e-12,
        atol=1e-12,
    ):
        raise AssertionError("few-entry estimator is not unbiased")

    full_constant = kernel.shape[0] ** (2.0 / 3.0)
    result = next(estimate_grid(kernel, exact, [full_constant], 1, 11))
    for component in ("s2", "s3", "w3", "t3"):
        if result[component]["error_over_scale_abs_max"] > 1e-12:
            raise AssertionError("p=1 estimator failed")
    emit({"record": "self_test", "status": "ok"})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Test a few-entry Schatten-3 estimator on public SKAT blocks."
    )
    parser.add_argument(
        "--fed-out",
        type=Path,
        default=Path(os.environ.get("FED_OUT", "~/fed_prep_out")).expanduser(),
        help="fed_prep.py output directory (default: FED_OUT or ~/fed_prep_out)",
    )
    parser.add_argument(
        "--genes",
        default="0",
        help="zero-based public block indices, e.g. 0,2-5 or all (default: 0)",
    )
    parser.add_argument(
        "--constants",
        default="1,2,5",
        help="C values in p=min(1,C*m^(-2/3)) (default: 1,2,5)",
    )
    parser.add_argument(
        "--seeds", type=int, default=50, help="Monte Carlo repeats per gene"
    )
    parser.add_argument("--seed", type=int, default=20260717, help="base random seed")
    parser.add_argument("--ridge-rel", type=float, default=RIDGE_REL_DEFAULT)
    parser.add_argument(
        "--chunk-rows",
        type=int,
        default=512,
        help="rows per streaming genotype chunk (default: 512)",
    )
    parser.add_argument("--self-test", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.self_test:
        run_self_test()
        return
    if args.seeds < 1:
        raise ValueError("seeds must be at least one")
    if args.ridge_rel < 0.0:
        raise ValueError("ridge-rel must be nonnegative")
    if args.chunk_rows < 1:
        raise ValueError("chunk-rows must be at least one")

    constants = parse_positive_floats(args.constants)
    root = args.fed_out
    block_sizes = np.loadtxt(root / "block_sizes.txt", dtype=np.int64, ndmin=1)
    if block_sizes.ndim != 1 or block_sizes.size == 0 or np.any(block_sizes <= 0):
        raise ValueError("invalid block_sizes.txt")
    genes = parse_gene_indices(args.genes, int(block_sizes.size))

    x_a = load_covariates(root / "A" / "cov.txt")
    x_b = load_covariates(root / "B" / "cov.txt")
    if x_a.shape[1] != x_b.shape[1]:
        raise ValueError("covariate dimensions do not match")
    xtx_ridged = ridged_xtx(x_a, x_b, args.ridge_rel)

    emit(
        {
            "constants": constants,
            "kernel": "normalized_public_public_K",
            "private_blocks_read": False,
            "record": "config",
            "schema": "skat_few_entry_schatten3_v1",
            "seeds": args.seeds,
        }
    )

    for gene_index in genes:
        size = int(block_sizes[gene_index])
        stats_a = public_block_stats(
            root / "A" / f"geno.{gene_index}.bin",
            x_a.shape[0],
            size,
            x_a,
            args.chunk_rows,
        )
        stats_b = public_block_stats(
            root / "B" / f"geno.{gene_index}.bin",
            x_b.shape[0],
            size,
            x_b,
            args.chunk_rows,
        )
        kernel = build_public_kernel(
            stats_a,
            stats_b,
            x_a.shape[0] + x_b.shape[0],
            xtx_ridged,
        )
        del stats_a, stats_b
        exact = exact_moments(kernel)

        emit(
            {
                "d3": finite(exact.d3),
                "d3_over_s3": ratio(exact.d3, exact.s3),
                "gene_index": gene_index,
                "m_public": size,
                "record": "exact",
                "s1": finite(exact.s1),
                "s2": finite(exact.s2),
                "s3": finite(exact.s3),
                "t3": finite(exact.t3),
                "t3_over_s3": ratio(exact.t3, exact.s3),
                "w3": finite(exact.w3),
                "w3_over_s3": ratio(exact.w3, exact.s3),
            }
        )

        for summary in estimate_grid(
            kernel, exact, constants, args.seeds, args.seed + gene_index
        ):
            emit(
                {
                    "gene_index": gene_index,
                    "m_public": size,
                    "record": "estimate",
                    **summary,
                }
            )


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, np.linalg.LinAlgError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc
