import unittest

import pandas as pd

from summarize_metrics import make_accuracy_table


class AccuracyTableTest(unittest.TestCase):
    def test_reports_all_and_non_failed_r_squared(self) -> None:
        rows = []
        for chromosome, secure_values, reference_values, converged in (
            (21, (0.1, 0.1), (0.1, 0.01), ("1", "0")),
            (22, (0.1, 0.01), (0.1, 0.01), ("1", "1")),
        ):
            for gene_index in range(2):
                rows.append(
                    {
                        "ancestry": "EUR",
                        "phenotype_index": "0",
                        "phenotype_name": "phenotype1",
                        "chromosome": str(chromosome),
                        "secure_burden_p": secure_values[gene_index],
                        "r_burden_p": reference_values[gene_index],
                        "secure_skat_wh_p": secure_values[gene_index],
                        "r_skat_liu_p": reference_values[gene_index],
                        "r_skat_davies_p": reference_values[gene_index],
                        "r_skat_davies_converged": converged[gene_index],
                    }
                )

        table = make_accuracy_table(pd.DataFrame(rows))
        davies = table.loc[
            table["Comparison"] == "SKAT WH vs Davies+fallback"
        ].iloc[0]

        self.assertEqual(davies["#gene (total(failed))"], "4(1)")
        self.assertEqual(davies["R^2 (all, R::SKAT)"], "0.000000")
        self.assertEqual(davies["R^2 (non-failed)"], "1.000000")
        self.assertEqual(davies["Worst R^2 (all)"], "-1.000000")
        self.assertEqual(davies["Worst phenotype"], "0:phenotype1")
        self.assertEqual(davies["Worst Chr"], "21")
        self.assertEqual(davies["Worst #gene (total(failed))"], "2(1)")


if __name__ == "__main__":
    unittest.main()
