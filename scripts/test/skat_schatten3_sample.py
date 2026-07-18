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
import itertools
import math
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

from skat_public_kernel import (
    PublicKernelSource,
    add_public_kernel_arguments,
    build_public_kernel,
    emit,
    finite,
    public_block_stats,
    ridged_xtx,
    validate_public_kernel_arguments,
)

SCHEMA = "skat_few_entry_schatten3_v3"
TINY = np.finfo(np.float64).tiny


@dataclass(frozen=True)
class ExactMoments:
    s1: float
    s2: float
    s3: float
    d3: float
    w3: float
    t3: float


def output(record: dict) -> None:
    emit({"schema": SCHEMA, **record})


def parse_positive_floats(spec: str) -> list[float]:
    values = sorted(
        set(float(piece.strip()) for piece in spec.split(",") if piece.strip())
    )
    if not values or any(not math.isfinite(value) or value <= 0.0 for value in values):
        raise ValueError("constants must be finite and positive")
    return values


def exact_moments(kernel: np.ndarray) -> ExactMoments:
    diagonal = np.diag(kernel)
    diagonal_squared = diagonal * diagonal
    row_squared_norms = np.einsum("ij,ij->i", kernel, kernel, optimize=True)

    d3 = float(np.dot(diagonal_squared, diagonal))
    w3 = 3.0 * float(np.dot(diagonal, row_squared_norms - diagonal_squared))
    s1 = float(np.sum(diagonal))
    s2 = float(np.sum(row_squared_norms))

    squared = kernel @ kernel
    s3 = float(np.einsum("ij,ji->", squared, kernel, optimize=True))
    t3 = s3 - d3 - w3
    return ExactMoments(s1=s1, s2=s2, s3=s3, d3=d3, w3=w3, t3=t3)


def q(values: Sequence[float], probability: float) -> float:
    return float(np.quantile(np.asarray(values, dtype=np.float64), probability))


def summarize_errors(
    estimates: Sequence[float], reference: float, scale: float
) -> dict:
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
                selected_count = edge_count
                s2_hat = exact.s2
                w3_hat = exact.w3
                t3_hat = exact.t3
            else:
                selected = uniforms < probability
                selected_count = int(np.count_nonzero(selected))
                s2_hat = d2 + (2.0 / probability) * float(
                    np.sum(upper_squares, where=selected, initial=0.0)
                )
                w3_hat = (3.0 / probability) * float(
                    np.sum(upper_wedges, where=selected, initial=0.0)
                )
                triangle_trace = sampled_triangle_trace(
                    size, upper_i, upper_j, upper_values, selected
                )
                t3_hat = triangle_trace / probability**3

            accumulator["edges"].append(selected_count)
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


@dataclass
class ContributionAccumulator:
    """Streaming concentration summary retaining only the largest contributions."""

    total: float = 0.0
    squared_total: float = 0.0
    count: int = 0
    largest: np.ndarray | None = None

    def add(self, values: np.ndarray) -> None:
        if values.size == 0:
            return
        self.total += float(np.sum(values))
        self.squared_total += float(np.dot(values, values))
        self.count += int(values.size)
        candidates = (
            values
            if self.largest is None
            else np.concatenate((self.largest, values))
        )
        if candidates.size > 1000:
            candidates = np.partition(candidates, candidates.size - 1000)[-1000:]
        self.largest = candidates

    def summary(self) -> dict:
        largest = (
            np.sort(self.largest)[::-1]
            if self.largest is not None
            else np.empty(0, dtype=np.float64)
        )
        result = {
            "effective_count": (
                finite(self.total * self.total / self.squared_total)
                if self.squared_total > 0.0
                else None
            )
        }
        for top in (1, 10, 100, 1000):
            used = min(top, self.count)
            result[f"top_{top}_fraction"] = (
                finite(float(np.sum(largest[:used])) / self.total)
                if self.total != 0.0
                else None
            )
        return result


def exact_concentration(kernel: np.ndarray, exact: ExactMoments) -> dict:
    diagonal = np.diag(kernel)
    offdiag_s2 = ContributionAccumulator()
    w3 = ContributionAccumulator()
    for row in range(kernel.shape[0] - 1):
        upper_squares = np.square(kernel[row, row + 1 :])
        offdiag_s2.add(2.0 * upper_squares)
        w3.add(3.0 * (diagonal[row] + diagonal[row + 1 :]) * upper_squares)

    d2 = float(np.dot(diagonal, diagonal))
    return {
        "d2": finite(d2),
        "d2_over_s2": ratio(d2, exact.s2),
        "offdiag_s2": finite(exact.s2 - d2),
        "offdiag_s2_concentration": offdiag_s2.summary(),
        "offdiag_s2_over_s2": ratio(exact.s2 - d2, exact.s2),
        "w3_concentration": w3.summary(),
    }


def draw_distinct_pairs(
    rng: np.random.Generator,
    count: int,
    first_probability: np.ndarray,
    second_probability: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    first_parts: list[np.ndarray] = []
    second_parts: list[np.ndarray] = []
    accepted = 0
    size = first_probability.size
    acceptance = 1.0 - float(np.dot(first_probability, second_probability))
    if acceptance <= 0.0:
        raise ValueError("pair proposal cannot draw distinct indices")
    while accepted < count:
        remaining = count - accepted
        batch = min(1_000_000, max(1024, int(math.ceil(1.1 * remaining / acceptance))))
        first = rng.choice(size, size=batch, p=first_probability)
        second = rng.choice(size, size=batch, p=second_probability)
        keep = first != second
        first_parts.append(first[keep])
        second_parts.append(second[keep])
        accepted += int(np.count_nonzero(keep))
    first = np.concatenate(first_parts)[:count]
    second = np.concatenate(second_parts)[:count]
    return np.minimum(first, second), np.maximum(first, second)


def draw_distinct_triples(
    rng: np.random.Generator,
    count: int,
    probability: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    parts: list[np.ndarray] = []
    accepted = 0
    size = probability.size
    p2 = float(np.dot(probability, probability))
    p3 = float(np.sum(probability**3))
    acceptance = 1.0 - 3.0 * p2 + 2.0 * p3
    if acceptance <= 0.0:
        raise ValueError("triangle proposal cannot draw distinct indices")
    while accepted < count:
        remaining = count - accepted
        batch = min(1_000_000, max(1024, int(math.ceil(1.1 * remaining / acceptance))))
        triples = np.column_stack(
            (
                rng.choice(size, size=batch, p=probability),
                rng.choice(size, size=batch, p=probability),
                rng.choice(size, size=batch, p=probability),
            )
        )
        keep = (
            (triples[:, 0] != triples[:, 1])
            & (triples[:, 0] != triples[:, 2])
            & (triples[:, 1] != triples[:, 2])
        )
        parts.append(np.sort(triples[keep], axis=1))
        accepted += int(np.count_nonzero(keep))
    triples = np.concatenate(parts, axis=0)[:count]
    return triples[:, 0], triples[:, 1], triples[:, 2]


@dataclass(frozen=True)
class ImportancePlan:
    constant: float
    probability: float
    requested_budget: int
    budget: int
    pair_samples: int
    triangle_samples: int

    @property
    def evaluations(self) -> int:
        return self.pair_samples + 3 * self.triangle_samples


def importance_plans(
    size: int, constants: Sequence[float], has_triangles: bool
) -> list[ImportancePlan]:
    edge_count = math.comb(size, 2)
    plans: list[ImportancePlan] = []
    for constant in constants:
        probability = min(1.0, constant * size ** (-2.0 / 3.0))
        requested = max(1, int(math.ceil(edge_count * probability)))
        budget = max(4, requested) if has_triangles else requested
        if has_triangles:
            pair_samples = max(1, budget // 2)
            triangle_samples = (budget - pair_samples) // 3
            if triangle_samples == 0:
                pair_samples, triangle_samples = budget - 3, 1
        else:
            pair_samples, triangle_samples = budget, 0
        plans.append(
            ImportancePlan(
                constant,
                probability,
                requested,
                budget,
                pair_samples,
                triangle_samples,
            )
        )
    return plans


def diagonal_normalizers(diagonal: np.ndarray) -> tuple[float, float, float]:
    """Stable sums for unordered diagonal pairs, wedges, and triples."""
    prefix_1 = 0.0
    prefix_2 = 0.0
    pair_sum = 0.0
    wedge_sum = 0.0
    triple_sum = 0.0
    for value in diagonal:
        value = float(value)
        triple_sum += value * pair_sum
        wedge_sum += value * value * prefix_1 + value * prefix_2
        pair_sum += value * prefix_1
        prefix_1 += value
        prefix_2 += value * value
    return pair_sum, wedge_sum, triple_sum


def importance_grid(
    kernel: np.ndarray,
    exact: ExactMoments,
    constants: Sequence[float],
    seeds: int,
    base_seed: int,
) -> Iterable[dict]:
    """Diagonal-aware importance estimator with nested entry-evaluation budgets."""
    size = kernel.shape[0]
    diagonal = np.maximum(np.diag(kernel), 0.0)
    d1 = float(np.sum(diagonal))
    diagonal_squared = diagonal * diagonal
    d2 = float(np.sum(diagonal_squared))
    positive = int(np.count_nonzero(diagonal))
    if positive < 2 or d1 == 0.0 or d2 == 0.0:
        raise ValueError("importance sampling requires at least two positive diagonals")

    probability_d = diagonal / d1
    probability_d2 = diagonal_squared / d2
    z2, zw, zt = diagonal_normalizers(diagonal)
    has_triangles = positive >= 3 and zt > 0.0
    plans = importance_plans(size, constants, has_triangles)
    max_pair_samples = max(plan.pair_samples for plan in plans)
    max_triangle_samples = max(plan.triangle_samples for plan in plans)
    accumulators = {
        plan.constant: {"s2": [], "w3": [], "t3": [], "s3": []}
        for plan in plans
    }

    seed_sequence = np.random.SeedSequence(base_seed + 1_000_000_007)
    for child_seed in seed_sequence.spawn(seeds):
        rng = np.random.default_rng(child_seed)
        use_q2 = rng.random(max_pair_samples) < 0.5
        pair_i = np.empty(max_pair_samples, dtype=np.int64)
        pair_j = np.empty(max_pair_samples, dtype=np.int64)
        for selected_q2, first_probability in (
            (True, probability_d),
            (False, probability_d2),
        ):
            positions = np.flatnonzero(use_q2 == selected_q2)
            if positions.size:
                sampled_i, sampled_j = draw_distinct_pairs(
                    rng, int(positions.size), first_probability, probability_d
                )
                pair_i[positions], pair_j[positions] = sampled_i, sampled_j

        di = diagonal[pair_i]
        dj = diagonal[pair_j]
        q2 = di * dj / z2
        qw = (di + dj) * di * dj / zw
        qmix = 0.5 * (q2 + qw)
        kij_squared = np.square(kernel[pair_i, pair_j])
        s2_prefix = np.cumsum(2.0 * kij_squared / qmix)
        w3_prefix = np.cumsum(3.0 * (di + dj) * kij_squared / qmix)

        if max_triangle_samples:
            tri_i, tri_j, tri_k = draw_distinct_triples(
                rng, max_triangle_samples, probability_d
            )
            qtri = diagonal[tri_i] * diagonal[tri_j] * diagonal[tri_k] / zt
            triangle = (
                kernel[tri_i, tri_j]
                * kernel[tri_j, tri_k]
                * kernel[tri_k, tri_i]
            )
            t3_prefix = np.cumsum(6.0 * triangle / qtri)
        else:
            t3_prefix = np.empty(0, dtype=np.float64)

        for plan in plans:
            estimates = accumulators[plan.constant]
            s2_hat = d2 + float(s2_prefix[plan.pair_samples - 1] / plan.pair_samples)
            w3_hat = float(w3_prefix[plan.pair_samples - 1] / plan.pair_samples)
            t3_hat = (
                float(t3_prefix[plan.triangle_samples - 1] / plan.triangle_samples)
                if plan.triangle_samples
                else 0.0
            )
            estimates["s2"].append(s2_hat)
            estimates["w3"].append(w3_hat)
            estimates["t3"].append(t3_hat)
            estimates["s3"].append(exact.d3 + w3_hat + t3_hat)

    for plan in plans:
        estimates = accumulators[plan.constant]
        yield {
            "constant": plan.constant,
            "diagonal_entries": size,
            "entry_evaluation_budget": plan.budget,
            "entry_evaluation_budget_requested": plan.requested_budget,
            "entry_evaluations": plan.evaluations,
            "entry_evaluations_with_diagonal": size + plan.evaluations,
            "negative_s3_rate": finite(
                float(np.mean(np.asarray(estimates["s3"]) <= 0.0))
            ),
            "pair_samples": plan.pair_samples,
            "p": plan.probability,
            "s2": summarize_errors(estimates["s2"], exact.s2, exact.s2),
            "s3": summarize_errors(estimates["s3"], exact.s3, exact.s3),
            "t3": summarize_errors(estimates["t3"], exact.t3, exact.s3),
            "triangle_samples": plan.triangle_samples,
            "w3": summarize_errors(estimates["w3"], exact.w3, exact.s3),
        }


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

    split = 8
    with tempfile.TemporaryDirectory() as temporary_directory:
        root = Path(temporary_directory)
        raw_geno[:split].tofile(root / "a.bin")
        raw_geno[split:].tofile(root / "b.bin")
        stats_a = public_block_stats(
            root / "a.bin", split, raw_geno.shape[1], design[:split], chunk_rows=3
        )
        stats_b = public_block_stats(
            root / "b.bin",
            raw_geno.shape[0] - split,
            raw_geno.shape[1],
            design[split:],
            chunk_rows=4,
        )
    xtx_ridged = ridged_xtx(design[:split], design[split:], 1e-6)
    built_kernel = build_public_kernel(stats_a, stats_b, raw_geno.shape[0], xtx_ridged)
    stacked_oriented = np.vstack(
        (
            np.maximum(raw_geno[:split], 0).astype(np.float64),
            np.maximum(raw_geno[split:], 0).astype(np.float64),
        )
    )
    for cohort in (stacked_oriented[:split], stacked_oriented[split:]):
        cohort_flip = cohort.sum(axis=0) > cohort.shape[0]
        cohort[:, cohort_flip] = 2.0 - cohort[:, cohort_flip]
    residual_gram = (
        stacked_oriented.T @ stacked_oriented
        - stacked_oriented.T
        @ design
        @ np.linalg.solve(xtx_ridged, design.T @ stacked_oriented)
    ) / raw_geno.shape[0]
    frequency = stacked_oriented.sum(axis=0) / (2.0 * raw_geno.shape[0])
    weight = 25.0 * np.power(1.0 - frequency, 24)
    direct_kernel = 0.5 * weight[:, None] * residual_gram * weight[None, :]
    direct_kernel = 0.5 * (direct_kernel + direct_kernel.T)
    if not np.allclose(built_kernel, direct_kernel, rtol=1e-12, atol=1e-12):
        raise AssertionError("public kernel does not match direct GtPG construction")

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

    importance_diagonal = np.diag(kernel)
    importance_i, importance_j = np.triu_indices(kernel.shape[0], k=1)
    importance_d1 = float(np.sum(importance_diagonal))
    importance_d2 = float(np.sum(importance_diagonal**2))
    importance_d3 = float(np.sum(importance_diagonal**3))
    importance_z2 = 0.5 * (importance_d1**2 - importance_d2)
    importance_zw = importance_d1 * importance_d2 - importance_d3
    stable_z2, stable_zw, stable_zt = diagonal_normalizers(importance_diagonal)
    closed_zt = (
        importance_d1**3
        - 3.0 * importance_d1 * importance_d2
        + 2.0 * importance_d3
    ) / 6.0
    if not np.allclose(
        (stable_z2, stable_zw, stable_zt),
        (importance_z2, importance_zw, closed_zt),
        rtol=1e-12,
        atol=1e-12,
    ):
        raise AssertionError("stable importance normalizers failed")
    importance_q2 = (
        importance_diagonal[importance_i]
        * importance_diagonal[importance_j]
        / importance_z2
    )
    importance_qw = (
        (importance_diagonal[importance_i] + importance_diagonal[importance_j])
        * importance_diagonal[importance_i]
        * importance_diagonal[importance_j]
        / importance_zw
    )
    importance_qmix = 0.5 * (importance_q2 + importance_qw)
    if not np.isclose(np.sum(importance_qmix), 1.0):
        raise AssertionError("pair importance probabilities do not sum to one")
    importance_upper_squared = kernel[importance_i, importance_j] ** 2
    expected_offdiag_s2 = float(
        np.sum(importance_qmix * 2.0 * importance_upper_squared / importance_qmix)
    )
    expected_w3 = float(
        np.sum(
            importance_qmix
            * 3.0
            * (importance_diagonal[importance_i] + importance_diagonal[importance_j])
            * importance_upper_squared
            / importance_qmix
        )
    )
    if not np.isclose(expected_offdiag_s2, exact.s2 - importance_d2):
        raise AssertionError("S2 importance estimator is not unbiased")
    if not np.isclose(expected_w3, exact.w3):
        raise AssertionError("W3 importance estimator is not unbiased")

    triples = np.asarray(
        list(itertools.combinations(range(kernel.shape[0]), 3)), dtype=np.int64
    )
    importance_zt = (
        importance_d1**3
        - 3.0 * importance_d1 * importance_d2
        + 2.0 * importance_d3
    ) / 6.0
    importance_qt = (
        importance_diagonal[triples[:, 0]]
        * importance_diagonal[triples[:, 1]]
        * importance_diagonal[triples[:, 2]]
        / importance_zt
    )
    triangle_values = (
        6.0
        * kernel[triples[:, 0], triples[:, 1]]
        * kernel[triples[:, 1], triples[:, 2]]
        * kernel[triples[:, 2], triples[:, 0]]
    )
    expected_t3 = float(np.sum(importance_qt * triangle_values / importance_qt))
    if not np.isclose(np.sum(importance_qt), 1.0):
        raise AssertionError("triangle importance probabilities do not sum to one")
    if not np.isclose(expected_t3, exact.t3, rtol=1e-10, atol=1e-12):
        raise AssertionError("T3 importance estimator is not unbiased")

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
    importance_result = next(importance_grid(kernel, exact, [1.0], 3, 13))
    for component in ("s2", "s3", "w3", "t3"):
        if importance_result[component]["estimate_mean"] is None:
            raise AssertionError("importance sampling returned a non-finite estimate")
    if (
        importance_result["entry_evaluations"]
        > importance_result["entry_evaluation_budget"]
    ):
        raise AssertionError("importance estimator exceeded its entry budget")
    tiny_plan = importance_plans(3, [0.01], has_triangles=True)[0]
    if tiny_plan.evaluations > tiny_plan.budget or tiny_plan.triangle_samples < 1:
        raise AssertionError("minimum unbiased importance budget failed")
    output({"record": "self_test", "status": "ok"})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Test a few-entry Schatten-3 estimator on public SKAT blocks."
    )
    add_public_kernel_arguments(parser, default_genes="0")
    parser.add_argument(
        "--constants",
        default="1,2,5",
        help="C values in p=min(1,C*m^(-2/3)) (default: 1,2,5)",
    )
    parser.add_argument(
        "--seeds", type=int, default=50, help="Monte Carlo repeats per gene"
    )
    parser.add_argument("--seed", type=int, default=20260717, help="base random seed")
    parser.add_argument(
        "--estimator",
        choices=("importance", "uniform", "both"),
        default="importance",
        help="estimator to run (default: importance)",
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
    validate_public_kernel_arguments(args)

    constants = parse_positive_floats(args.constants)
    source = PublicKernelSource.load(args.fed_out, args.ridge_rel)
    genes = source.selected_genes(args.genes)

    output(
        {
            "constants": constants,
            "kernel": "normalized_public_public_K",
            "private_blocks_read": False,
            "record": "config",
            "seeds": args.seeds,
            "estimator": args.estimator,
        }
    )

    for gene_index, size, kernel in source.kernels(genes, args.chunk_rows):
        exact = exact_moments(kernel)
        concentration = exact_concentration(kernel, exact)

        output(
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
                "t3_subtraction_condition": finite(
                    (abs(exact.s3) + abs(exact.d3) + abs(exact.w3))
                    / max(abs(exact.t3), TINY)
                ),
                "w3": finite(exact.w3),
                "w3_over_s3": ratio(exact.w3, exact.s3),
                **concentration,
            }
        )

        if args.estimator in ("uniform", "both"):
            for summary in estimate_grid(
                kernel, exact, constants, args.seeds, args.seed + gene_index
            ):
                output(
                    {
                        "gene_index": gene_index,
                        "m_public": size,
                        "record": "uniform_estimate",
                        **summary,
                    }
                )
        if args.estimator in ("importance", "both"):
            for summary in importance_grid(
                kernel, exact, constants, args.seeds, args.seed + gene_index
            ):
                output(
                    {
                        "gene_index": gene_index,
                        "m_public": size,
                        "record": "importance_estimate",
                        **summary,
                    }
                )


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, np.linalg.LinAlgError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc
