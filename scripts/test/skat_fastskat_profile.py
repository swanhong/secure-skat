#!/usr/bin/env python3
"""Profile whether FastSKAT-style spectral deflation can replace public tau_3.

This script reads only the public-list genotype blocks and covariates produced by
scripts/preprocessing/fed_prep.py.  It never reads phenotypes, private genotype
blocks, variant identifiers, positions, or private manifest fields.

Stdout is deterministic JSONL containing only gene-level aggregate spectral
concentration and approximation-error summaries.  It never prints eigenvalues,
eigenvectors, genotype entries, covariates, or variant metadata.  The output is
still an AoU research result and must go through the applicable export review.

The kernel being profiled is K_pp, not the full public/private block kernel.  The
result therefore answers the implementation question: can the public Hutchinson
tau_3 pass be replaced by a FastSKAT-style top-k plus Satterthwaite-tail model?

Example:
    FED_OUT=~/runs/out260713034114 python3 scripts/test/skat_fastskat_profile.py \
        --genes all --ks 1,2,5,10,20,50,100,200
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path
from typing import Sequence

import numpy as np

# Running this file directly places scripts/test on sys.path.
from skat_schatten3_sample import (  # noqa: E402
    RIDGE_REL_DEFAULT,
    build_public_kernel,
    emit,
    finite,
    load_covariates,
    parse_gene_indices,
    public_block_stats,
    ridged_xtx,
)


SCHEMA = "skat_fastskat_profile_v1"
SQRT_2 = math.sqrt(2.0)


def parse_positive_ints(spec: str) -> list[int]:
    values = sorted(set(int(piece.strip()) for piece in spec.split(",") if piece.strip()))
    if not values or values[0] < 1:
        raise ValueError("ks must contain positive integers")
    return values


def parse_probabilities(spec: str, name: str, include_one: bool) -> list[float]:
    values = sorted(
        set(float(piece.strip()) for piece in spec.split(",") if piece.strip())
    )
    upper_ok = (lambda value: value <= 1.0) if include_one else (lambda value: value < 1.0)
    if not values or any(value <= 0.0 or not upper_ok(value) for value in values):
        suffix = "(0,1]" if include_one else "(0,1)"
        raise ValueError(f"{name} must contain probabilities in {suffix}")
    return values


def safe_ratio(numerator: float, denominator: float) -> float | None:
    if denominator <= 0.0:
        return None
    return finite(numerator / denominator)


def generalized_saddle_sf(
    statistic: float, weights: np.ndarray, dfs: np.ndarray
) -> float:
    """Lugannani-Rice upper tail for sum_j weights[j] * chi2(dfs[j])."""
    weights = np.asarray(weights, dtype=np.float64)
    dfs = np.asarray(dfs, dtype=np.float64)
    keep = (weights > 0.0) & (dfs > 0.0)
    weights = weights[keep]
    dfs = dfs[keep]
    if weights.size == 0:
        raise ValueError("saddlepoint mixture is empty")

    mean = float(np.dot(dfs, weights))
    if statistic <= mean:
        raise ValueError("saddlepoint evaluator is restricted to the upper tail")

    upper = 0.5 / float(np.max(weights))
    low, high = 0.0, upper * (1.0 - 1e-14)

    def first_derivative(value: float) -> float:
        denominator = 1.0 - 2.0 * weights * value
        return float(np.sum(dfs * weights / denominator))

    if first_derivative(high) <= statistic:
        raise ValueError("failed to bracket saddlepoint root")

    for _ in range(160):
        middle = 0.5 * (low + high)
        if first_derivative(middle) < statistic:
            low = middle
        else:
            high = middle

    saddle = 0.5 * (low + high)
    denominator = 1.0 - 2.0 * weights * saddle
    cumulant = -0.5 * float(np.sum(dfs * np.log(denominator)))
    second_derivative = float(
        np.sum(2.0 * dfs * weights * weights / (denominator * denominator))
    )
    w_squared = max(2.0 * (saddle * statistic - cumulant), 0.0)
    w = math.sqrt(w_squared)
    u = saddle * math.sqrt(second_derivative)
    if w <= 1e-12 or u <= 0.0:
        raise ValueError("unstable saddlepoint evaluation")

    normal_sf = 0.5 * math.erfc(w / SQRT_2)
    normal_pdf = math.exp(-0.5 * w * w) / math.sqrt(2.0 * math.pi)
    probability = normal_sf + normal_pdf * (1.0 / u - 1.0 / w)
    return min(1.0, max(np.finfo(np.float64).tiny, float(probability)))


def saddle_quantile(weights: np.ndarray, target: float) -> float:
    dfs = np.ones(weights.size, dtype=np.float64)
    mean = float(np.sum(weights))
    variance = 2.0 * float(np.sum(weights * weights))
    low = mean * (1.0 + 1e-12)
    high = mean + 8.0 * math.sqrt(max(variance, np.finfo(np.float64).tiny))
    while generalized_saddle_sf(high, weights, dfs) > target:
        high = mean + 2.0 * (high - mean)

    for _ in range(100):
        middle = 0.5 * (low + high)
        if generalized_saddle_sf(middle, weights, dfs) > target:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def wh_sf(statistic: float, s1: float, s2: float, s3: float) -> float:
    if s2 <= 0.0 or s3 <= 0.0:
        return 1.0
    u = (statistic - s1) * s3 / (s2 * s2)
    h = 2.0 * s3 * s3 / (9.0 * s2 * s2 * s2)
    argument = 1.0 + u
    if argument <= 0.0 or h <= 0.0:
        return 1.0
    z_value = (math.cbrt(argument) - 1.0 + h) / math.sqrt(h)
    return 0.5 * math.erfc(z_value / SQRT_2)


def minimum_rank_for_mass(
    descending_eigenvalues: np.ndarray, power: int, threshold: float
) -> int:
    contributions = descending_eigenvalues**power
    total = float(np.sum(contributions))
    if total <= 0.0:
        return 0
    index = int(np.searchsorted(np.cumsum(contributions), threshold * total, side="left"))
    return min(index + 1, descending_eigenvalues.size)


def spectrum_from_kernel(kernel: np.ndarray) -> tuple[np.ndarray, dict]:
    eigenvalues = np.linalg.eigvalsh(kernel)
    maximum = float(np.max(eigenvalues))
    if maximum <= 0.0:
        raise ValueError("kernel has no positive eigenvalue")

    numerical_floor = max(1, kernel.shape[0]) * np.finfo(np.float64).eps * maximum
    material_negative = eigenvalues < -100.0 * numerical_floor
    if np.any(material_negative):
        raise ValueError("kernel has a material negative eigenvalue")

    clipped_count = int(np.count_nonzero(eigenvalues < 0.0))
    eigenvalues = np.maximum(eigenvalues, 0.0)[::-1]
    rank = {
        "negative_eigenvalues_clipped": clipped_count,
        "rank_rel_1e_8": int(np.count_nonzero(eigenvalues > maximum * 1e-8)),
        "rank_rel_1e_10": int(np.count_nonzero(eigenvalues > maximum * 1e-10)),
        "rank_rel_1e_12": int(np.count_nonzero(eigenvalues > maximum * 1e-12)),
    }
    return eigenvalues, rank


def tail_model(eigenvalues: np.ndarray, count: int) -> dict:
    count = min(count, eigenvalues.size)
    powers = {power: eigenvalues**power for power in (1, 2, 3)}
    totals = {power: float(np.sum(values)) for power, values in powers.items()}
    heads = {
        power: float(np.sum(values[:count])) for power, values in powers.items()
    }
    tails = {power: totals[power] - heads[power] for power in (1, 2, 3)}

    if count == eigenvalues.size or tails[1] <= 0.0 or tails[2] <= 0.0:
        tail3_approx = tails[3]
        scale = 0.0
        degrees_of_freedom = 0.0
    else:
        tail3_approx = tails[2] * tails[2] / tails[1]
        scale = tails[2] / tails[1]
        degrees_of_freedom = tails[1] * tails[1] / tails[2]

    hybrid_s3 = heads[3] + tail3_approx
    return {
        "count": count,
        "degrees_of_freedom": degrees_of_freedom,
        "head_s1_fraction": safe_ratio(heads[1], totals[1]),
        "head_s2_fraction": safe_ratio(heads[2], totals[2]),
        "head_s3_fraction": safe_ratio(heads[3], totals[3]),
        "hybrid_s3_error_over_s3": safe_ratio(hybrid_s3 - totals[3], totals[3]),
        "scale": scale,
        "tail_s1": tails[1],
        "tail_s2": tails[2],
        "tail_s3": tails[3],
        "tail3_approx": tail3_approx,
        "tail3_error_over_tail3": safe_ratio(tail3_approx - tails[3], tails[3]),
        "totals": totals,
        "heads": heads,
    }


def pvalue_summary(
    eigenvalues: np.ndarray,
    model: dict,
    statistic: float,
    target: float,
) -> dict:
    count = int(model["count"])
    fast_weights = eigenvalues[:count].copy()
    fast_dfs = np.ones(count, dtype=np.float64)
    if model["degrees_of_freedom"] > 0.0:
        fast_weights = np.append(fast_weights, model["scale"])
        fast_dfs = np.append(fast_dfs, model["degrees_of_freedom"])

    fast_p = generalized_saddle_sf(statistic, fast_weights, fast_dfs)
    totals = model["totals"]
    hybrid_s3 = model["heads"][3] + model["tail3_approx"]
    wh_exact = wh_sf(statistic, totals[1], totals[2], totals[3])
    wh_hybrid = wh_sf(statistic, totals[1], totals[2], hybrid_s3)

    def log10_error(probability: float) -> float:
        return math.log10(probability) - math.log10(target)

    return {
        "fastskat_log10_p_error": finite(log10_error(fast_p)),
        "fastskat_p_over_reference": finite(fast_p / target),
        "hybrid_wh_log10_p_error": finite(log10_error(wh_hybrid)),
        "hybrid_wh_p_over_reference": finite(wh_hybrid / target),
        "reference": "full_spectrum_lugannani_rice",
        "target_p": target,
        "wh_exact_moments_log10_p_error": finite(log10_error(wh_exact)),
        "wh_exact_moments_p_over_reference": finite(wh_exact / target),
    }


def profile_gene(
    eigenvalues: np.ndarray,
    ks: Sequence[int],
    mass_thresholds: Sequence[float],
    p_targets: Sequence[float],
) -> tuple[dict, list[dict]]:
    s1 = float(np.sum(eigenvalues))
    s2 = float(np.sum(eigenvalues**2))
    s3 = float(np.sum(eigenvalues**3))
    summary = {
        "effective_rank_s1_s2": finite(s1 * s1 / s2),
        "effective_rank_s2_s3": finite(s2**3 / (s3 * s3)),
        "mass_ranks": {
            str(threshold): {
                "s1": minimum_rank_for_mass(eigenvalues, 1, threshold),
                "s2": minimum_rank_for_mass(eigenvalues, 2, threshold),
                "s3": minimum_rank_for_mass(eigenvalues, 3, threshold),
            }
            for threshold in mass_thresholds
        },
    }

    quantiles = {
        target: saddle_quantile(eigenvalues, target) for target in p_targets
    }
    records: list[dict] = []
    for requested_count in ks:
        model = tail_model(eigenvalues, requested_count)
        record = {
            key: finite(value) if isinstance(value, float) else value
            for key, value in model.items()
            if key not in ("totals", "heads", "tail_s1", "tail_s2", "tail_s3", "tail3_approx")
        }
        record["pvalue_checks"] = [
            pvalue_summary(eigenvalues, model, quantiles[target], target)
            for target in p_targets
        ]
        records.append(record)
    return summary, records


def run_self_test() -> None:
    from scipy.stats import chi2

    equal = np.full(12, 0.7, dtype=np.float64)
    statistic = float(0.7 * chi2.isf(1e-3, df=12))
    got = generalized_saddle_sf(
        statistic, np.asarray([0.7]), np.asarray([12.0])
    )
    if abs(math.log10(got) - math.log10(1e-3)) > 0.05:
        raise AssertionError("generalized saddlepoint failed equal-weight check")

    model = tail_model(equal, 3)
    if abs(model["tail3_error_over_tail3"]) > 1e-12:
        raise AssertionError("Satterthwaite tail3 must be exact for equal eigenvalues")

    unequal = np.asarray([5.0, 2.0, 1.0, 0.5, 0.1], dtype=np.float64)
    model = tail_model(unequal, 1)
    if model["tail3_approx"] > model["tail_s3"] * (1.0 + 1e-12):
        raise AssertionError("tail3 approximation violated Cauchy-Schwarz lower bound")
    if minimum_rank_for_mass(unequal, 3, 0.9) > minimum_rank_for_mass(unequal, 1, 0.9):
        raise AssertionError("higher moments should not need more head eigenvalues here")
    emit({"record": "self_test", "schema": SCHEMA, "status": "ok"})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Profile FastSKAT top-k and Satterthwaite-tail accuracy on public SKAT kernels."
    )
    parser.add_argument(
        "--fed-out",
        type=Path,
        default=Path(os.environ.get("FED_OUT", "~/fed_prep_out")).expanduser(),
        help="fed_prep.py output directory (default: FED_OUT or ~/fed_prep_out)",
    )
    parser.add_argument(
        "--genes",
        default="all",
        help="zero-based block indices, e.g. 0,2-5 or all (default: all)",
    )
    parser.add_argument(
        "--ks", default="1,2,5,10,20,50,100,200", help="top-k values"
    )
    parser.add_argument(
        "--mass-thresholds",
        default="0.9,0.95,0.99,0.999",
        help="head mass thresholds (default: 0.9,0.95,0.99,0.999)",
    )
    parser.add_argument(
        "--p-targets",
        default="0.001,0.0001,0.0000025",
        help="upper-tail probabilities for LR comparisons",
    )
    parser.add_argument("--ridge-rel", type=float, default=RIDGE_REL_DEFAULT)
    parser.add_argument("--chunk-rows", type=int, default=512)
    parser.add_argument("--self-test", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.self_test:
        run_self_test()
        return
    if args.ridge_rel < 0.0:
        raise ValueError("ridge-rel must be nonnegative")
    if args.chunk_rows < 1:
        raise ValueError("chunk-rows must be at least one")

    ks = parse_positive_ints(args.ks)
    mass_thresholds = parse_probabilities(
        args.mass_thresholds, "mass-thresholds", include_one=True
    )
    p_targets = parse_probabilities(args.p_targets, "p-targets", include_one=False)

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
            "genes_requested": len(genes),
            "kernel": "normalized_public_public_K",
            "ks": ks,
            "mass_thresholds": mass_thresholds,
            "private_blocks_read": False,
            "pvalue_reference": "full_spectrum_lugannani_rice",
            "p_targets": p_targets,
            "phenotypes_read": False,
            "record": "config",
            "schema": SCHEMA,
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

        eigenvalues, rank = spectrum_from_kernel(kernel)
        del kernel
        summary, topk_records = profile_gene(
            eigenvalues, ks, mass_thresholds, p_targets
        )
        emit(
            {
                "gene_index": gene_index,
                "m_public": size,
                "record": "spectrum",
                "schema": SCHEMA,
                **rank,
                **summary,
            }
        )
        for record in topk_records:
            emit(
                {
                    "gene_index": gene_index,
                    "m_public": size,
                    "record": "topk",
                    "schema": SCHEMA,
                    **record,
                }
            )


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, np.linalg.LinAlgError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc
