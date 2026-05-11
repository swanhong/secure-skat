#!/usr/bin/env python3

from __future__ import annotations

import argparse
import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from argument import parse_args
from samples import build_sample_files


def quiet_build_sample_files(*args: object) -> tuple[int, int, int]:
    with contextlib.redirect_stdout(io.StringIO()):
        return build_sample_files(*args)


def make_args(tmp: Path, **overrides: object) -> argparse.Namespace:
    values = {
        "pheno_file": str(tmp / "phenotype.csv"),
        "cov_file": str(tmp / "covariate.csv"),
        "id_col": None,
        "pheno_col": "pheno",
        "pheno_col_index": None,
        "cov_cols": None,
        "cov_col_indices": None,
        "pheno_sep": None,
        "cov_sep": None,
        "n_samples": 0,
        "party1_frac": 0.5,
        "seed": 7,
    }
    values.update(overrides)
    return argparse.Namespace(**values)


def write_psam(prefix: Path) -> None:
    prefix.with_suffix(".psam").write_text(
        "\n".join(
            [
                "#FID IID",
                "S1 S1",
                "S2 S2",
                "S3 S3",
                "S4 S4",
            ]
        )
        + "\n"
    )


class SplitInputTests(unittest.TestCase):
    def test_split_tables_join_by_first_column_and_use_all_covariates(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            raw_prefix = tmp / "raw"
            write_psam(raw_prefix)
            (tmp / "phenotype.csv").write_text(
                "\n".join(
                    [
                        "person_id,pheno,other_trait",
                        "S1,1.5,10",
                        "S2,,20",
                        "S3,3.5,30",
                        "SX,9.0,90",
                    ]
                )
                + "\n"
            )
            (tmp / "covariate.csv").write_text(
                "\n".join(
                    [
                        "person_id,age,pc1,pc2",
                        "S1,51,0.1,0.2",
                        "S2,52,0.3,0.4",
                        "S3,53,0.5,0.6",
                        "S4,54,0.7,0.8",
                    ]
                )
                + "\n"
            )

            out_dataset = tmp / "dataset"
            n_party1, n_party2, num_covs = quiet_build_sample_files(
                make_args(tmp),
                raw_prefix,
                tmp / "work",
                out_dataset,
            )

            self.assertEqual(n_party1 + n_party2, 2)
            self.assertEqual(num_covs, 3)
            pheno_lines = (out_dataset / "party1" / "pheno.txt").read_text().splitlines()
            pheno_lines += (out_dataset / "party2" / "pheno.txt").read_text().splitlines()
            self.assertEqual(sorted(float(value) for value in pheno_lines), [1.5, 3.5])
            for party in ("party1", "party2"):
                cov_path = out_dataset / party / "cov.txt"
                for line in cov_path.read_text().splitlines():
                    self.assertEqual(len(line.split("\t")), 3)

    def test_split_tables_support_1_based_indices(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            raw_prefix = tmp / "raw"
            write_psam(raw_prefix)
            (tmp / "phenotype.csv").write_text("person_id,pheno,other_trait\nS1,1.5,10\nS3,3.5,30\n")
            (tmp / "covariate.csv").write_text("person_id,age,pc1,pc2\nS1,51,0.1,0.2\nS3,53,0.5,0.6\n")

            out_dataset = tmp / "dataset"
            _, _, num_covs = quiet_build_sample_files(
                make_args(tmp, pheno_col=None, pheno_col_index=3, cov_col_indices="2,4"),
                raw_prefix,
                tmp / "work",
                out_dataset,
            )

            self.assertEqual(num_covs, 2)
            pheno_lines = (out_dataset / "party1" / "pheno.txt").read_text().splitlines()
            pheno_lines += (out_dataset / "party2" / "pheno.txt").read_text().splitlines()
            self.assertEqual(sorted(float(value) for value in pheno_lines), [10.0, 30.0])
            all_cov_lines = (out_dataset / "party1" / "cov.txt").read_text().splitlines()
            all_cov_lines += (out_dataset / "party2" / "cov.txt").read_text().splitlines()
            self.assertEqual(sorted(all_cov_lines), ["51\t0.2", "53\t0.6"])

    def test_split_tables_fail_on_nonnumeric_covariate(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            raw_prefix = tmp / "raw"
            write_psam(raw_prefix)
            (tmp / "phenotype.csv").write_text("person_id,pheno\nS1,1.5\nS3,3.5\n")
            (tmp / "covariate.csv").write_text("person_id,age\nS1,51\nS3,not_numeric\n")

            with self.assertRaisesRegex(ValueError, "Non-numeric covariate value"):
                quiet_build_sample_files(make_args(tmp), raw_prefix, tmp / "work", tmp / "dataset")

    def test_cli_rejects_ambiguous_phenotype_selector(self) -> None:
        argv = [
            "build_pgen_window_dataset.py",
            "--pgen-prefix",
            "dummy",
            "--pheno-file",
            "phenotype.csv",
            "--cov-file",
            "covariate.csv",
            "--pheno-col",
            "pheno",
            "--pheno-col-index",
            "2",
            "--out-dataset",
            "dataset",
            "--config-out-dir",
            "config",
        ]
        with patch.object(sys, "argv", argv), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            parse_args()

    def test_cli_keeps_legacy_merged_table_mode(self) -> None:
        argv = [
            "build_pgen_window_dataset.py",
            "--pgen-prefix",
            "dummy",
            "--pheno-file",
            "pheno_cov.tsv",
            "--id-col",
            "sample_id",
            "--pheno-col",
            "pheno",
            "--cov-cols",
            "age,pc1",
            "--out-dataset",
            "dataset",
            "--config-out-dir",
            "config",
        ]
        with patch.object(sys, "argv", argv):
            args = parse_args()
        self.assertIsNone(args.cov_file)
        self.assertEqual(args.id_col, "sample_id")
        self.assertEqual(args.cov_cols, "age,pc1")

    def test_cli_rejects_zero_phenotype_index(self) -> None:
        argv = [
            "build_pgen_window_dataset.py",
            "--pgen-prefix",
            "dummy",
            "--pheno-file",
            "phenotype.csv",
            "--cov-file",
            "covariate.csv",
            "--pheno-col-index",
            "0",
            "--out-dataset",
            "dataset",
            "--config-out-dir",
            "config",
        ]
        with patch.object(sys, "argv", argv), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            parse_args()


if __name__ == "__main__":
    unittest.main()
