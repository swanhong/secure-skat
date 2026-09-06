import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import Barrier
from unittest.mock import Mock

import numpy as np

from .model import GenePlan, GeneRef, GeneVariants, PhenoCovRows, VariantRef
from .pipeline import assign_roles, extract_genotypes, prepare_blocks


class AssignRolesTest(unittest.TestCase):
    def test_shared_rate(self):
        variants = tuple(
            VariantRef(f"1:{position}:A:G", position, "PASS")
            for position in range(1, 6)
        )
        groups = (
            GeneVariants(GeneRef("ENSG1", "G1", "1", 0), variants),
            GeneVariants(GeneRef("ENSG2", "G2", "1", 1), variants[1:4]),
        )

        default = assign_roles(groups, role_seed=42)
        self.assertEqual(
            tuple(role for _, role in default[0].variant_roles),
            ("private", "shared", "shared", "shared", "public_only"),
        )
        self.assertEqual(tuple(variant for variant, _ in default[0].variant_roles), variants)

        expected = tuple(
            tuple((variant, "shared") for variant in group.variants)
            for group in groups
        )
        for seed in (1, 999):
            plans = assign_roles(groups, role_seed=seed, shared_rate=1.0)
            self.assertEqual(tuple(plan.variant_roles for plan in plans), expected)

    def test_prepare_blocks_skips_completed_output(self):
        with TemporaryDirectory() as directory:
            out_dir = Path(directory)
            (out_dir / "pos.txt").touch()
            extractor = Mock(side_effect=AssertionError("extractor called"))
            rows = PhenoCovRows((), (), ())

            result = prepare_blocks(Path("unused"), (), rows, rows, 42, out_dir, extractor=extractor)

            self.assertEqual(result, out_dir)
            extractor.assert_not_called()

    def test_extract_genotypes_runs_cohorts_in_parallel(self):
        variant = VariantRef("1:1:A:G", 1, "PASS")
        plan = GenePlan(GeneRef("ENSG1", "G1", "1", 0), ((variant, "shared"),))
        rows_a = PhenoCovRows(("A",), (), ())
        rows_b = PhenoCovRows(("B",), (), ())
        barrier = Barrier(2)

        def extractor(_, sample_ids, variant_keys):
            barrier.wait(timeout=5)
            value = 1 if sample_ids == ("A",) else 2
            return np.full((1, 1), value, dtype=np.int8), variant_keys

        result = extract_genotypes(Path("unused"), rows_a, rows_b, (plan,), extractor)

        self.assertEqual(result[1][0, 0], 1)
        self.assertEqual(result[3][0, 0], 2)


if __name__ == "__main__":
    unittest.main()
