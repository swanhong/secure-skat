import tempfile
import unittest
from pathlib import Path

import compare_results
import plot_results


class PlotResultsTest(unittest.TestCase):
    def test_r_squared_uses_negative_log10_p(self) -> None:
        count, score = compare_results.r_squared(
            [
                {"secure": "0.01", "reference": "0.001"},
                {"secure": "0.1", "reference": "0.1"},
            ],
            "secure",
            "reference",
        )

        self.assertEqual(count, 2)
        self.assertAlmostEqual(score, 0.5)

    def test_scatter_values_use_negative_log10_p(self) -> None:
        secure, reference = plot_results.negative_log10_p_value_pairs(
            [
                {"secure": "0.01", "reference": "0.001"},
                {"secure": "0", "reference": "1"},
                {"secure": "-0.1", "reference": "0.1"},
            ],
            "secure",
            "reference",
        )

        self.assertEqual(len(secure), 2)
        self.assertAlmostEqual(secure[0], 2)
        self.assertAlmostEqual(reference[0], 3)
        self.assertAlmostEqual(secure[1], 300)
        self.assertAlmostEqual(reference[1], 0)

    def test_writes_scatter_plots(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            temp_dir = Path(directory)
            run_dir = temp_dir / "run"
            comparison_root = run_dir / "comparison"
            comparison_dir = comparison_root / "EUR"
            comparison_dir.mkdir(parents=True)
            chromosome_dir = run_dir / "prepared" / "EUR" / "chr21"
            chromosome_dir.mkdir(parents=True)

            (chromosome_dir / "genes.txt").write_text(
                "GENE1\nGENE2\n",
                encoding="utf-8",
            )
            (chromosome_dir / "block_sizes.txt").write_text(
                "1\n1\n",
                encoding="utf-8",
            )
            (chromosome_dir / "pos.txt").write_text(
                "21 100\n21 200\n",
                encoding="utf-8",
            )

            config_path = temp_dir / "run.conf"
            config_path.write_text(
                f'run_dir = "{run_dir.as_posix()}"\n'
                'ancestries = ["EUR"]\n',
                encoding="utf-8",
            )

            (comparison_dir / "all_comparison.csv").write_text(
                "chromosome,gene_index,gene_id,"
                "phenotype_index,phenotype_name,"
                "secure_burden_p,r_burden_p,"
                "secure_skat_wh_p,r_skat_liu_p\n"
                "21,0,GENE1,0,phenotype1,0.1,0.1,0.2,0.21\n"
                "21,1,GENE2,0,phenotype1,0.8,0.79,0.9,0.88\n",
                encoding="utf-8",
            )

            plot_results.plot_results(config_path)

            plots_dir = comparison_dir / "plots"
            self.assertTrue(
                (plots_dir / "scatter_burden_pheno0_chr21.png").exists()
            )
            self.assertTrue(
                (plots_dir / "scatter_skat_liu_pheno0_chr21.png").exists()
            )
            self.assertTrue((plots_dir / "_SUCCESS").exists())
            self.assertTrue(
                (comparison_root / "_PLOTS_SUCCESS").exists()
            )


if __name__ == "__main__":
    unittest.main()
