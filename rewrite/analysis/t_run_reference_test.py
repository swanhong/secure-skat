import csv
import io
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import run_reference


REFERENCE_HEADER = (
    "gene_index\tgene_id\tphenotype_index\tburden_p\t"
    "skat_davies_p\tskat_davies_converged\tskat_liu_p\n"
)


class RunReferenceTest(unittest.TestCase):
    def test_runs_r_and_python_engines(self) -> None:
        for engine, launcher, script_name, blas_threads in (
            ("r", "Rscript", "compute_reference.R", "1"),
            ("python", sys.executable, "compute_reference.py", "1"),
        ):
            with self.subTest(engine=engine), tempfile.TemporaryDirectory() as directory:
                run_dir = Path(directory) / "run"
                config_path = Path(directory) / "run.conf"
                config_path.write_text(
                    f'run_dir = "{run_dir.as_posix()}"\n'
                    "chromosomes = [21, 22]\n"
                    'ancestries = ["EUR"]\n'
                    'phenotype_columns = ["phenotype1", "phenotype2"]\n',
                    encoding="utf-8",
                )

                def run_backend(
                    command: list[str],
                    **kwargs: object,
                ) -> subprocess.CompletedProcess[str]:
                    self.assertEqual(command[0], launcher)
                    self.assertEqual(Path(command[1]).name, script_name)
                    self.assertEqual(
                        kwargs["env"]["OPENBLAS_NUM_THREADS"],
                        blas_threads,
                    )
                    chromosome = Path(command[-1]).name.removeprefix("chr")
                    output = (
                        REFERENCE_HEADER
                        + f"0\tGENE{chromosome}\t0\t0.1\t0.2\t1\t0.3\n"
                        + f"0\tGENE{chromosome}\t1\t0.4\t0.5\tNA\t0.6\n"
                    )
                    return subprocess.CompletedProcess(command, 0, stdout=output)

                with (
                    patch("os.cpu_count", return_value=8),
                    patch.dict(run_reference.os.environ, {}, clear=True),
                    patch.object(
                        run_reference,
                        "monotonic",
                        side_effect=(100.0, 130.0, 160.0),
                    ),
                    patch.object(
                        run_reference.subprocess,
                        "run",
                        side_effect=run_backend,
                    ),
                    patch("sys.stdout", new_callable=io.StringIO) as output,
                ):
                    run_reference.run_reference(config_path, engine)

                progress = output.getvalue()
                self.assertIn(
                    "Reference progress: 1 done / 1 remaining "
                    "(0.5 min elapsed)",
                    progress,
                )
                self.assertIn(
                    "Reference progress: 2 done / 0 remaining "
                    "(1.0 min elapsed)",
                    progress,
                )

                reference_root = run_dir / "reference"
                with (
                    reference_root / "EUR" / "all_r_results.csv"
                ).open(newline="", encoding="utf-8") as result_file:
                    rows = list(csv.DictReader(result_file))

                self.assertTrue((reference_root / "_SUCCESS").exists())
                self.assertEqual(
                    [row["chromosome"] for row in rows],
                    ["21", "21", "22", "22"],
                )
                self.assertEqual(
                    [row["phenotype_name"] for row in rows],
                    ["phenotype1", "phenotype2"] * 2,
                )
                self.assertEqual(
                    [row["r_skat_davies_converged"] for row in rows],
                    ["1", "NA", "1", "NA"],
                )

    def test_worker_defaults(self) -> None:
        with (
            patch("os.cpu_count", return_value=32),
            patch.dict(run_reference.os.environ, {}, clear=True),
        ):
            self.assertEqual(run_reference.worker_settings("r", 6), (6, 1))
            self.assertEqual(
                run_reference.worker_settings("python", 66),
                (16, 1),
            )

    def test_worker_overrides(self) -> None:
        with (
            patch("os.cpu_count", return_value=32),
            patch.dict(
                run_reference.os.environ,
                {
                    "REFERENCE_WORKERS": "8",
                    "REFERENCE_BLAS_THREADS": "2",
                },
                clear=True,
            ),
        ):
            self.assertEqual(
                run_reference.worker_settings("r", 66),
                (8, 2),
            )


if __name__ == "__main__":
    unittest.main()
