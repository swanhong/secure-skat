import csv
import os
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import fed_prep


class MultiPhenotypePrepTest(unittest.TestCase):
    def test_lipid_phenotype_order(self):
        self.assertEqual(fed_prep.LIPID_PHENO_COLS, (
            "LDLC_final_mgdl_6sd_masked",
            "HDLC_mgdl_6sd_masked",
            "TotChol_corrected_mvp_explicit_duration_mgdl_6sd_masked",
            "ln_Trig_6sd_masked",
            "nonHDL_corrected_mvp_explicit_duration_mgdl_6sd_masked",
        ))

    def test_load_phenos_keeps_only_complete_finite_rows(self):
        columns = ("p0", "p1", "p2")
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "pheno.csv")
            with open(path, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(("person_id", *columns))
                writer.writerow(("a", "1", "2", "3"))
                writer.writerow(("b", "4", "", "6"))
                writer.writerow(("c", "7", "nan", "9"))
                writer.writerow(("d", "10", "11", "12"))

            phenos = fed_prep.load_phenos(path, columns)
            self.assertEqual(phenos, {"a": [1.0, 2.0, 3.0], "d": [10.0, 11.0, 12.0]})

            output = os.path.join(tmp, "pheno.txt")
            fed_prep.write_pheno(["d", "a"], phenos, output)
            np.testing.assert_allclose(np.loadtxt(output), [[10, 11, 12], [1, 2, 3]])


if __name__ == "__main__":
    unittest.main()
