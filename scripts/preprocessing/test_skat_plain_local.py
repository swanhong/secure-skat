#!/usr/bin/env python3
"""Numerical regression tests for the plaintext SKAT mixture p-value references."""

import math
import unittest

import numpy as np

from skat_plain_local import _chi2_cdf, skat_davies_p, skat_davies_p_batch


class TestSKATDavies(unittest.TestCase):
    def test_rank_one_is_exact_and_scale_invariant(self):
        want = math.erfc(math.sqrt(5.0 / 2.0))
        for scale in (1e-6, 1.0, 1e3, 1e6):
            with self.subTest(scale=scale):
                self.assertAlmostEqual(skat_davies_p([scale], 5.0 * scale), want, delta=1e-14)

    def test_equal_rank_two_is_exact_and_scale_invariant(self):
        # lambda*(chi-square_1 + chi-square_1) = lambda*chi-square_2.
        want = math.exp(-5.0)
        for scale in (1e-6, 1.0, 1e3, 1e6):
            with self.subTest(scale=scale):
                self.assertAlmostEqual(
                    skat_davies_p([scale, scale], 10.0 * scale), want, delta=1e-14)

    def test_tbr1_rank_one_regression(self):
        # The retired Python integrator returned about 0.5 after scaling this chr2-like case.
        q_over_lambda = 11.0537975
        want = math.erfc(math.sqrt(q_over_lambda / 2.0))
        self.assertAlmostEqual(
            skat_davies_p([1e6], q_over_lambda * 1e6), want, delta=1e-14)

    def test_bounded_ruben_low_tail_regressions(self):
        # Independent R::SKAT qfc values; these were failure/wrong-answer cases for direct Imhof.
        cases = (
            ([1.0, 0.7], 1e-6, 0.9999994026326544),
            ([1.0, 0.001], 1e-6, 0.9999841905898519),
            ([1.0, 0.2, 0.05], 2.1e-5, 0.9999997432422292),
        )
        for weights, q, want in cases:
            with self.subTest(weights=weights, q=q):
                self.assertAlmostEqual(skat_davies_p(weights, q), want, delta=1e-8)

    def test_high_rank_low_tail_never_returns_spurious_half(self):
        weights = np.random.default_rng(7).uniform(0.01, 1.0, 50)
        self.assertGreater(skat_davies_p(weights, 2e-5), 1.0 - 1e-12)

    def test_mixed_spectrum_batch_matches_high_accuracy_qfc(self):
        # Independent R::SKAT qfc/Davies references (lim=1e7).
        weights = np.array([1.0, 0.2, 0.05])
        q_and_want = (
            (3.0, 0.09960883351601679),
            (10.0, 0.001820259189572404),
            (30.0, 4.979889034473217e-8),
            (50.0, 1.7683632336229493e-12),
        )
        cases, expected = [], []
        for scale in (1e-6, 1.0, 1e6):
            for q, want in q_and_want:
                cases.append((weights * scale, q * scale))
                expected.append(want)
        got = skat_davies_p_batch(cases)
        for value, want in zip(got, expected):
            self.assertAlmostEqual(value, want, delta=max(2e-10, want * 1e-6))
        for offset in range(0, len(got), len(q_and_want)):
            self.assertTrue(all(got[offset + i] > got[offset + i + 1]
                                for i in range(len(q_and_want) - 1)))

    def test_unresolved_extreme_tail_fails_instead_of_clipping(self):
        with self.assertRaisesRegex(RuntimeError, "did not resolve"):
            skat_davies_p([1.0], 2000.0)
        with self.assertWarnsRegex(RuntimeWarning, "2/3 cases unresolved"):
            values = skat_davies_p_batch(
                [([1.0, 0.2, 0.05], 3.0),
                 ([1.0, 0.2, 0.05], 60.0),
                 ([1.0], 2000.0)],
                on_unresolved="nan")
        self.assertAlmostEqual(values[0], 0.09960883351601679, delta=2e-10)
        self.assertTrue(math.isnan(values[1]))
        self.assertTrue(math.isnan(values[2]))

    def test_incomplete_gamma_nonconvergence_is_explicit(self):
        with self.assertRaisesRegex(ArithmeticError, "did not converge"):
            _chi2_cdf(1e6, 1e6)

    def test_many_small_positive_modes_are_not_pruned(self):
        weights = [1.0] + [0.9e-10] * 1000
        self.assertAlmostEqual(
            skat_davies_p(weights, 1.0), 0.31731052964950024, delta=2e-10)

    def test_tiny_negative_roundoff_is_dropped(self):
        want = math.erfc(math.sqrt(5.0 / 2.0))
        self.assertAlmostEqual(skat_davies_p([1.0, -1e-13], 5.0), want, delta=1e-14)

    def test_materially_negative_eigenvalue_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "negative eigenvalue"):
            skat_davies_p([1.0, -1e-4], 1.0)

    def test_normalized_statistic_overflow_is_explicit(self):
        with self.assertRaisesRegex(OverflowError, "not finite"):
            skat_davies_p([1e-300, 5e-301], 1e308)


if __name__ == "__main__":
    unittest.main()
