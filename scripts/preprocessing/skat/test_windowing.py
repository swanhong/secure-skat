#!/usr/bin/env python3

from __future__ import annotations

import math
import tempfile
import unittest
from pathlib import Path

from windowing import read_source_variants


class WindowingTests(unittest.TestCase):
    def test_read_source_variants_keeps_nonfinite_afreq_as_nan(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            prefix = Path(td) / "source"
            prefix.with_suffix(".pvar").write_text(
                "\n".join(
                    [
                        "#CHROM POS ID REF ALT",
                        "22 100 v1 A C",
                        "22 200 v2 G T",
                    ]
                )
                + "\n"
            )
            prefix.with_suffix(".afreq").write_text(
                "\n".join(
                    [
                        "#CHROM ID REF ALT ALT_FREQS OBS_CT",
                        "22 v1 A C 0.01 20",
                        "22 v2 G T NA 0",
                    ]
                )
                + "\n"
            )

            variants = read_source_variants(prefix)

            self.assertEqual(len(variants), 2)
            self.assertAlmostEqual(float(variants[0]["maf"]), 0.01)
            self.assertTrue(math.isnan(float(variants[1]["maf"])))


if __name__ == "__main__":
    unittest.main()
