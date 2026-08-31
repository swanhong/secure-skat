import csv
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


KEY_COLUMNS = ("gene_index", "gene_id", "phenotype_index")
TOLERANCES = {
    "burden_p": 1e-10,
    "skat_davies_p": 1e-10,
    "skat_liu_p": 1e-8,
}


def write_fixture(directory: Path) -> None:
    sample_counts = {"A": 4, "B": 4}
    covariates = np.array([[0], [1], [0], [1], [1], [0], [1], [0]])
    phenotypes = np.array(
        [
            [0.2, 1], [-0.4, 0.5], [1.1, -0.3], [0.7, 1.2],
            [-1, 0.1], [0.3, -0.7], [1.4, -0.4], [-0.2, 0.8],
        ]
    )
    public = (
        np.array(
            [[9, 2], [0, 2], [1, 1], [2, 2], [0, 0], [1, 1], [2, 2], [1, 1]]
        ),
        np.empty((8, 0), dtype=np.int8),
        np.zeros((8, 1), dtype=np.int8),
    )
    private_b = (
        np.array([[0], [1], [2], [0]]),
        np.array([[0], [1], [2], [0]]),
        np.empty((4, 0), dtype=np.int8),
    )

    (directory / "genes.txt").write_text("MULTI\nRANK1\nEMPTY\n")
    (directory / "block_sizes.txt").write_text("2\n0\n1\n")
    offset = 0
    for cohort, count in sample_counts.items():
        cohort_dir = directory / cohort
        (cohort_dir / "geno").mkdir(parents=True)
        (cohort_dir / "private").mkdir()
        cohort_slice = slice(offset, offset + count)
        np.savetxt(cohort_dir / "cov.txt", covariates[cohort_slice])
        np.savetxt(cohort_dir / "pheno.txt", phenotypes[cohort_slice])
        for gene_index, genotype in enumerate(public):
            genotype[cohort_slice].astype(np.int8).tofile(
                cohort_dir / "geno" / f"block.{gene_index}.bin"
            )
            if cohort == "B":
                private_b[gene_index].tofile(
                    cohort_dir / "private" / f"block.{gene_index}.bin"
                )
        offset += count


def run_reference(command: list[str]) -> list[dict[str, str]]:
    environment = os.environ | {"OPENBLAS_NUM_THREADS": "1"}
    completed = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
        env=environment,
    )
    return list(csv.DictReader(completed.stdout.splitlines(), delimiter="\t"))


class ReferenceParityTest(unittest.TestCase):
    def test_python_matches_r_skat(self) -> None:
        scripts = Path(__file__).resolve().parent
        with tempfile.TemporaryDirectory() as directory:
            fixture = Path(directory)
            write_fixture(fixture)
            r_rows = run_reference(
                ["Rscript", str(scripts / "r_skat/compute_reference.R"), str(fixture)]
            )
            python_rows = run_reference(
                [
                    sys.executable,
                    str(scripts / "python_skat/compute_reference.py"),
                    str(fixture),
                ]
            )

        self.assertEqual(
            [tuple(row[column] for column in KEY_COLUMNS) for row in r_rows],
            [
                ("0", "MULTI", "0"), ("0", "MULTI", "1"),
                ("1", "RANK1", "0"), ("1", "RANK1", "1"),
                ("2", "EMPTY", "0"), ("2", "EMPTY", "1"),
            ],
        )
        self.assertEqual(
            [tuple(row[column] for column in KEY_COLUMNS) for row in r_rows],
            [tuple(row[column] for column in KEY_COLUMNS) for row in python_rows],
        )
        expected_convergence = ["1", "1", "1", "1", "NA", "NA"]
        self.assertEqual(
            [row["skat_davies_converged"] for row in r_rows],
            expected_convergence,
        )
        self.assertEqual(
            [row["skat_davies_converged"] for row in python_rows],
            expected_convergence,
        )
        for column, tolerance in TOLERANCES.items():
            error = max(
                abs(float(r_row[column]) - float(python_row[column]))
                for r_row, python_row in zip(r_rows, python_rows)
            )
            self.assertLess(error, tolerance, column)


if __name__ == "__main__":
    unittest.main()
