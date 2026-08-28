import tempfile
import unittest
from pathlib import Path

import plot_results


class PlotResultsTest(unittest.TestCase):
    def test_writes_scatter_plots(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            temp_dir = Path(directory)
            run_dir = temp_dir / "run"
            comparison_dir = run_dir / "comparison"
            comparison_dir.mkdir(parents=True)

            config_path = temp_dir / "run.conf"
            config_path.write_text(
                f'run_dir = "{run_dir.as_posix()}"\n',
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


if __name__ == "__main__":
    unittest.main()