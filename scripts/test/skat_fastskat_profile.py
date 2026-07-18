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
import math
import sys
from dataclasses import dataclass
from typing import Sequence

import numpy as np

# Running this file directly places scripts/test on sys.path.
from skat_public_kernel import (  # noqa: E402
    PublicKernelSource,
    add_public_kernel_arguments,
    emit,
    finite,
    validate_public_kernel_arguments,
)


SCHEMA = "skat_fastskat_profile_v2"
SQRT_2 = math.sqrt(2.0)
SECURE_WH_ARGUMENT_BOUNDS = (0.1, 9.0)
SECURE_WH_S2_BOUNDS = (1e-8, 1e12)


def output(record: dict) -> None:
    emit({"schema": SCHEMA, **record})


def parse_positive_ints(spec: str) -> list[int]:
    values = sorted(
        set(int(piece.strip()) for piece in spec.split(",") if piece.strip())
    )
    if not values or values[0] < 1:
        raise ValueError("ks must contain positive integers")
    return values


def parse_probabilities(spec: str, name: str, include_one: bool) -> list[float]:
    values = sorted(
        set(float(piece.strip()) for piece in spec.split(",") if piece.strip())
    )
    upper_ok = (
        (lambda value: value <= 1.0)
        if include_one
        else (lambda value: value < 1.0)
    )
    if not values or any(
        not math.isfinite(value) or value <= 0.0 or not upper_ok(value)
        for value in values
    ):
        suffix = "(0,1]" if include_one else "(0,1)"
        raise ValueError(f"{name} must contain probabilities in {suffix}")
    return values


def safe_ratio(numerator: float, denominator: float) -> float | None:
    if denominator <= 0.0:
        return None
    return finite(numerator / denominator)


@dataclass(frozen=True)
class SaddleMixture:
    """Cached positive mixture for Lugannani-Rice upper-tail diagnostics."""

    weights: np.ndarray
    dfs: np.ndarray
    upper: float

    @classmethod
    def build(cls, weights: np.ndarray, dfs: np.ndarray) -> "SaddleMixture":
        weights = np.asarray(weights, dtype=np.float64)
        dfs = np.asarray(dfs, dtype=np.float64)
        keep = (weights > 0.0) & (dfs > 0.0)
        weights, dfs = weights[keep], dfs[keep]
        if weights.size == 0:
            raise ValueError("saddlepoint mixture is empty")
        upper = 0.5 / float(np.max(weights)) * (1.0 - 1e-14)
        return cls(weights, dfs, upper)

    def first_derivative(self, saddle: float) -> float:
        denominator = 1.0 - 2.0 * self.weights * saddle
        return float(np.sum(self.dfs * self.weights / denominator))

    def sf_at_saddle(
        self, saddle: float, statistic: float | None = None
    ) -> tuple[float, float]:
        denominator = 1.0 - 2.0 * self.weights * saddle
        if statistic is None:
            statistic = float(np.sum(self.dfs * self.weights / denominator))
        cumulant = -0.5 * float(np.sum(self.dfs * np.log(denominator)))
        second_derivative = float(
            np.sum(
                2.0
                * self.dfs
                * self.weights
                * self.weights
                / (denominator * denominator)
            )
        )
        w = math.sqrt(max(2.0 * (saddle * statistic - cumulant), 0.0))
        u = saddle * math.sqrt(second_derivative)
        if w <= 1e-12 or u <= 0.0:
            return statistic, 0.5
        normal_sf = 0.5 * math.erfc(w / SQRT_2)
        normal_pdf = math.exp(-0.5 * w * w) / math.sqrt(2.0 * math.pi)
        probability = normal_sf + normal_pdf * (1.0 / u - 1.0 / w)
        probability = min(
            1.0, max(np.finfo(np.float64).tiny, float(probability))
        )
        return statistic, probability

    def sf(self, statistic: float) -> float:
        mean = float(np.dot(self.dfs, self.weights))
        if statistic <= mean:
            raise ValueError("saddlepoint evaluator is restricted to the upper tail")
        if self.first_derivative(self.upper) <= statistic:
            raise ValueError("failed to bracket saddlepoint root")

        low, high = 0.0, self.upper
        for _ in range(64):
            middle = 0.5 * (low + high)
            if middle == low or middle == high:
                break
            if self.first_derivative(middle) < statistic:
                low = middle
            else:
                high = middle
        return self.sf_at_saddle(0.5 * (low + high), statistic)[1]

    def quantile(self, target: float) -> float:
        """Invert LR directly in saddle space, avoiding a nested root search."""
        if target <= 0.0 or target > 0.1:
            raise ValueError("LR upper-tail targets must lie in (0,0.1]")
        low, high = 0.0, self.upper
        for _ in range(64):
            middle = 0.5 * (low + high)
            if middle == low or middle == high:
                break
            _, probability = self.sf_at_saddle(middle)
            if probability > target:
                low = middle
            else:
                high = middle
        return self.sf_at_saddle(0.5 * (low + high))[0]


def generalized_saddle_sf(
    statistic: float, weights: np.ndarray, dfs: np.ndarray
) -> float:
    return SaddleMixture.build(weights, dfs).sf(statistic)


def saddle_quantile(weights: np.ndarray, target: float) -> float:
    return SaddleMixture.build(weights, np.ones(weights.size)).quantile(target)


def wh_parameters(
    statistic: float, s1: float, s2: float, s3: float
) -> tuple[float, float]:
    if s2 <= 0.0 or s3 <= 0.0:
        return 0.0, 0.0
    u = (statistic - s1) * s3 / (s2 * s2)
    h = 2.0 * s3 * s3 / (9.0 * s2 * s2 * s2)
    return 1.0 + u, h


def wh_sf(
    statistic: float,
    s1: float,
    s2: float,
    s3: float,
    argument_bounds: tuple[float, float] | None = None,
) -> float:
    argument, h = wh_parameters(statistic, s1, s2, s3)
    if argument_bounds is not None:
        argument = min(argument_bounds[1], max(argument_bounds[0], argument))
    if argument <= 0.0 or h <= 0.0:
        return 1.0
    z_value = (math.cbrt(argument) - 1.0 + h) / math.sqrt(h)
    return 0.5 * math.erfc(z_value / SQRT_2)


def secure_clamped_wh(
    statistic: float, s1: float, s2: float, s3: float
) -> tuple[float, float, bool, bool]:
    secure_s2 = min(SECURE_WH_S2_BOUNDS[1], max(SECURE_WH_S2_BOUNDS[0], s2))
    argument, _ = wh_parameters(statistic, s1, secure_s2, s3)
    s2_clamp_active = secure_s2 != s2
    argument_clamp_active = not (
        SECURE_WH_ARGUMENT_BOUNDS[0]
        <= argument
        <= SECURE_WH_ARGUMENT_BOUNDS[1]
    )
    probability = wh_sf(
        statistic, s1, secure_s2, s3, SECURE_WH_ARGUMENT_BOUNDS
    )
    return probability, argument, s2_clamp_active, argument_clamp_active


@dataclass(frozen=True)
class SpectrumCache:
    eigenvalues: np.ndarray
    prefix: np.ndarray
    suffix: np.ndarray
    totals: np.ndarray

    @classmethod
    def build(cls, eigenvalues: np.ndarray) -> "SpectrumCache":
        eigenvalues = np.asarray(eigenvalues, dtype=np.float64)
        powers = np.vstack((eigenvalues, eigenvalues**2, eigenvalues**3))
        prefix = np.cumsum(powers, axis=1)
        suffix = np.cumsum(powers[:, ::-1], axis=1)[:, ::-1]
        return cls(eigenvalues, prefix, suffix, suffix[:, 0])

    def mass_rank(self, power: int, threshold: float) -> int:
        total = float(self.totals[power - 1])
        if total <= 0.0:
            return 0
        index = int(
            np.searchsorted(
                self.prefix[power - 1], threshold * total, side="left"
            )
        )
        return min(index + 1, self.eigenvalues.size)

    def head_tail(self, count: int) -> tuple[dict[int, float], dict[int, float]]:
        count = min(count, self.eigenvalues.size)
        heads = {
            power: (float(self.prefix[power - 1, count - 1]) if count else 0.0)
            for power in (1, 2, 3)
        }
        tails = {
            power: (
                float(self.suffix[power - 1, count])
                if count < self.eigenvalues.size
                else 0.0
            )
            for power in (1, 2, 3)
        }
        return heads, tails


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


def tail_model(spectrum: SpectrumCache, count: int) -> dict:
    count = min(count, spectrum.eigenvalues.size)
    totals = {power: float(spectrum.totals[power - 1]) for power in (1, 2, 3)}
    heads, tails = spectrum.head_tail(count)

    if count == spectrum.eigenvalues.size or tails[1] <= 0.0 or tails[2] <= 0.0:
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


def log_error_fields(prefix: str, probability: float, target: float) -> dict:
    return {
        f"{prefix}_log10_p_error": finite(
            math.log10(probability) - math.log10(target)
        ),
        f"{prefix}_p_over_reference": finite(probability / target),
    }


def fastskat_mixture(spectrum: SpectrumCache, model: dict) -> SaddleMixture:
    count = int(model["count"])
    fast_weights = spectrum.eigenvalues[:count].copy()
    fast_dfs = np.ones(count, dtype=np.float64)
    if model["degrees_of_freedom"] > 0.0:
        fast_weights = np.append(fast_weights, model["scale"])
        fast_dfs = np.append(fast_dfs, model["degrees_of_freedom"])
    return SaddleMixture.build(fast_weights, fast_dfs)


def pvalue_summary(
    mixture: SaddleMixture, model: dict, statistic: float, target: float
) -> dict:
    fast_p = mixture.sf(statistic)
    totals = model["totals"]
    hybrid_s3 = model["heads"][3] + model["tail3_approx"]
    wh_ideal = wh_sf(statistic, totals[1], totals[2], hybrid_s3)
    wh_secure, secure_argument, s2_clamped, argument_clamped = secure_clamped_wh(
        statistic, totals[1], totals[2], hybrid_s3
    )

    return {
        **log_error_fields("fastskat", fast_p, target),
        **log_error_fields("hybrid_wh_ideal", wh_ideal, target),
        **log_error_fields("hybrid_wh_secure_clamped", wh_secure, target),
        "hybrid_wh_secure_argument": finite(secure_argument),
        "hybrid_wh_argument_clamp_active": argument_clamped,
        "hybrid_wh_s2_clamp_active": s2_clamped,
        "target_p": target,
    }


def reference_summary(
    spectrum: SpectrumCache, statistic: float, target: float
) -> dict:
    s1, s2, s3 = (float(value) for value in spectrum.totals)
    ideal = wh_sf(statistic, s1, s2, s3)
    secure, secure_argument, s2_clamped, argument_clamped = secure_clamped_wh(
        statistic, s1, s2, s3
    )
    return {
        **log_error_fields("wh_ideal_exact_moments", ideal, target),
        **log_error_fields("wh_secure_clamped_exact_moments", secure, target),
        "target_p": target,
        "wh_secure_argument": finite(secure_argument),
        "wh_argument_clamp_active": argument_clamped,
        "wh_s2_clamp_active": s2_clamped,
    }


def profile_gene(
    eigenvalues: np.ndarray,
    ks: Sequence[int],
    mass_thresholds: Sequence[float],
    p_targets: Sequence[float],
) -> tuple[dict, list[dict]]:
    spectrum = SpectrumCache.build(eigenvalues)
    s1, s2, s3 = (float(value) for value in spectrum.totals)
    full_mixture = SaddleMixture.build(eigenvalues, np.ones(eigenvalues.size))
    quantiles = {target: full_mixture.quantile(target) for target in p_targets}
    summary = {
        "effective_rank_s1_s2": finite(s1 * s1 / s2),
        "effective_rank_s2_s3": finite(s2**3 / (s3 * s3)),
        "mass_ranks": {
            str(threshold): {
                "s1": spectrum.mass_rank(1, threshold),
                "s2": spectrum.mass_rank(2, threshold),
                "s3": spectrum.mass_rank(3, threshold),
            }
            for threshold in mass_thresholds
        },
        "reference_checks": [
            reference_summary(spectrum, quantiles[target], target)
            for target in p_targets
        ],
    }

    records: list[dict] = []
    requested_by_count: dict[int, list[int]] = {}
    for requested_count in ks:
        count = min(requested_count, eigenvalues.size)
        requested_by_count.setdefault(count, []).append(requested_count)
    for count, requested_counts in requested_by_count.items():
        model = tail_model(spectrum, count)
        mixture = fastskat_mixture(spectrum, model)
        record = {
            key: finite(value) if isinstance(value, float) else value
            for key, value in model.items()
            if key
            not in (
                "totals",
                "heads",
                "tail_s1",
                "tail_s2",
                "tail_s3",
                "tail3_approx",
            )
        }
        record["requested_counts"] = requested_counts
        record["pvalue_checks"] = [
            pvalue_summary(mixture, model, quantiles[target], target)
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
    mixture = SaddleMixture.build(equal, np.ones(equal.size))
    quantile = mixture.quantile(1e-4)
    if abs(mixture.sf(quantile) / 1e-4 - 1.0) > 1e-10:
        raise AssertionError("direct saddle-space quantile inversion failed")
    edge_quantile = mixture.quantile(0.1)
    if abs(mixture.sf(edge_quantile) / 0.1 - 1.0) > 1e-10:
        raise AssertionError("LR quantile tail-boundary check failed")

    equal_spectrum = SpectrumCache.build(equal)
    model = tail_model(equal_spectrum, 3)
    if abs(model["tail3_error_over_tail3"]) > 1e-12:
        raise AssertionError("Satterthwaite tail3 must be exact for equal eigenvalues")

    unequal = np.asarray([5.0, 2.0, 1.0, 0.5, 0.1], dtype=np.float64)
    unequal_spectrum = SpectrumCache.build(unequal)
    model = tail_model(unequal_spectrum, 1)
    if model["tail3_approx"] > model["tail_s3"] * (1.0 + 1e-12):
        raise AssertionError("tail3 approximation violated Cauchy-Schwarz lower bound")
    if unequal_spectrum.mass_rank(3, 0.9) > unequal_spectrum.mass_rank(1, 0.9):
        raise AssertionError(
            "higher moments should not need more head eigenvalues here"
        )

    full_model = tail_model(unequal_spectrum, unequal.size)
    unequal_mixture = SaddleMixture.build(unequal, np.ones(unequal.size))
    full_quantile = unequal_mixture.quantile(1e-4)
    full_check = pvalue_summary(
        fastskat_mixture(unequal_spectrum, full_model),
        full_model,
        full_quantile,
        1e-4,
    )
    if abs(full_check["fastskat_p_over_reference"] - 1.0) > 1e-10:
        raise AssertionError("full-k FastSKAT mixture must match the LR reference")
    deep_target = 2.5e-6
    deep_quantile = unequal_mixture.quantile(deep_target)
    if not reference_summary(
        unequal_spectrum, deep_quantile, deep_target
    )["wh_argument_clamp_active"]:
        raise AssertionError("secure WH clamp diagnostic failed")

    known_factor = np.arange(1.0, 13.0).reshape(4, 3)
    known_kernel = known_factor @ known_factor.T
    known_eigenvalues, rank = spectrum_from_kernel(known_kernel)
    if rank["rank_rel_1e_10"] != 2 or np.any(known_eigenvalues < 0.0):
        raise AssertionError("known PSD spectrum/rank check failed")
    output({"record": "self_test", "status": "ok"})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Profile FastSKAT top-k and Satterthwaite-tail accuracy on public "
            "SKAT kernels."
        )
    )
    add_public_kernel_arguments(parser, default_genes="all")
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
    parser.add_argument("--self-test", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.self_test:
        run_self_test()
        return
    validate_public_kernel_arguments(args)

    ks = parse_positive_ints(args.ks)
    mass_thresholds = parse_probabilities(
        args.mass_thresholds, "mass-thresholds", include_one=True
    )
    p_targets = parse_probabilities(args.p_targets, "p-targets", include_one=False)
    if p_targets[-1] > 0.1:
        raise ValueError("p-targets must be upper-tail probabilities at most 0.1")

    source = PublicKernelSource.load(args.fed_out, args.ridge_rel)
    genes = source.selected_genes(args.genes)

    output(
        {
            "genes_requested": len(genes),
            "kernel": "normalized_public_public_K",
            "ks": ks,
            "mass_thresholds": mass_thresholds,
            "private_blocks_read": False,
            "pvalue_reference": "full_spectrum_lugannani_rice",
            "pvalue_reference_is_exact": False,
            "p_targets": p_targets,
            "phenotypes_read": False,
            "record": "config",
            "secure_wh_argument_bounds": list(SECURE_WH_ARGUMENT_BOUNDS),
            "secure_wh_emulation": "s2_and_argument_clamps_only",
            "secure_wh_s2_bounds": list(SECURE_WH_S2_BOUNDS),
        }
    )

    for gene_index, size, kernel in source.kernels(genes, args.chunk_rows):
        eigenvalues, rank = spectrum_from_kernel(kernel)
        del kernel
        summary, topk_records = profile_gene(
            eigenvalues, ks, mass_thresholds, p_targets
        )
        output(
            {
                "gene_index": gene_index,
                "m_public": size,
                "record": "spectrum",
                **rank,
                **summary,
            }
        )
        for record in topk_records:
            output(
                {
                    "gene_index": gene_index,
                    "m_public": size,
                    "record": "topk",
                    **record,
                }
            )


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, np.linalg.LinAlgError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc
