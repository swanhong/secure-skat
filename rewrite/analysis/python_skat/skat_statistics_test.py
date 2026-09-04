import unittest

from skat_statistics import noncentral_chi_square_sf


class NoncentralChiSquareTest(unittest.TestCase):
    def test_survival_probability_does_not_exceed_one(self) -> None:
        probability = noncentral_chi_square_sf(1e-15, 1.0, 50.0)

        self.assertEqual(probability, 1.0)


if __name__ == "__main__":
    unittest.main()
