import unittest

from .model import GeneRef, GeneVariants, VariantRef
from .pipeline import assign_roles


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


if __name__ == "__main__":
    unittest.main()
