import csv
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import run_reference


R_HEADER = (
    "gene_index\tgene_id\tphenotype_index\tburden_p\t"
    "skat_davies_p\tskat_davies_converged\tskat_liu_p\n"
)


class RunReferenceTest(unittest.TestCase):
    def test_writes_chromosome_and_combined_results(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            temp_dir = Path(directory)
            run_dir = temp_dir / "run"
            config_path = temp_dir / "run.conf"
            config_path.write_text(
                f'run_dir = "{run_dir.as_posix()}"\n'
                "chromosomes = [21, 22]\n"
                'ancestries = ["EUR"]\n'
                'phenotype_columns = ["phenotype1", "phenotype2"]\n',
                encoding="utf-8",
            )

            r_outputs = [
                R_HEADER
                + "0\tGENE21\t0\t0.1\t0.2\t1\t0.3\n"
                + "0\tGENE21\t1\t0.4\t0.5\t1\t0.6\n",
                R_HEADER
                + "0\tGENE22\t0\t1\t1\tNA\t1\n"
                + "0\tGENE22\t1\t0.7\t0.8\t0\t0.9\n",
            ]
            completed = [
                subprocess.CompletedProcess([], 0, stdout=output)
                for output in r_outputs
            ]

            with patch.object(
                run_reference.subprocess,
                "run",
                side_effect=completed,
            ) as run_r:
                run_reference.run_reference(config_path)

            reference_root = run_dir / "reference"
            reference_dir = reference_root / "EUR"
            self.assertTrue((reference_root / "_SUCCESS").exists())
            self.assertEqual(run_r.call_count, 2)
            self.assertEqual(
                run_r.call_args_list[0].args[0][-1],
                str(run_dir / "prepared" / "EUR" / "chr21"),
            )

            with (reference_dir / "all_r_results.csv").open(
                newline="",
                encoding="utf-8",
            ) as result_file:
                rows = list(csv.DictReader(result_file))

            self.assertEqual(len(rows), 4)
            self.assertEqual(
                [row["chromosome"] for row in rows],
                ["21", "21", "22", "22"],
            )
            self.assertEqual(
                [row["phenotype_name"] for row in rows],
                ["phenotype1", "phenotype2"] * 2,
            )
            self.assertEqual(rows[2]["r_skat_davies_converged"], "NA")


if __name__ == "__main__":
    unittest.main()
