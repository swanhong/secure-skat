import csv
import subprocess
import tempfile
import unittest
from pathlib import Path


class SummarizeMetricsTest(unittest.TestCase):
    def test_prints_chromosome_and_ancestry_totals(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            run_dir = Path(temporary_directory)
            metrics_dir = run_dir / "metrics" / "EUR"
            metrics_dir.mkdir(parents=True)
            metrics_path = metrics_dir / "metrics_party1.csv"

            with metrics_path.open("w", newline="") as output:
                writer = csv.DictWriter(
                    output,
                    fieldnames=[
                        "process",
                        "ancestry",
                        "stage",
                        "parent_stage",
                        "measurement_kind",
                        "chromosome",
                        "count",
                        "duration_seconds",
                    ],
                )
                writer.writeheader()
                for chromosome, total in ((21, 100), (22, 200)):
                    for stage, duration in (
                        ("compute_weights", 2),
                        ("packed_statistics", total - 10),
                        ("first_gtg_action", 20),
                        ("second_gtg_action", 10),
                        ("private_trace_correction", 30),
                        ("finalize", 5),
                        ("release", 1),
                        ("chromosome_total", total),
                    ):
                        writer.writerow(
                            {
                                "process": "party1",
                                "ancestry": "EUR",
                                "stage": stage,
                                "parent_stage": "",
                                "measurement_kind": "",
                                "chromosome": chromosome,
                                "count": 1,
                                "duration_seconds": duration,
                            }
                        )

            script = Path(__file__).resolve().parents[2] / "summarize_metrics.sh"
            completed = subprocess.run(
                ["bash", str(script), str(run_dir)],
                check=True,
                capture_output=True,
                text=True,
            )

            self.assertIn("Metrics timing summary (party1)", completed.stdout)
            self.assertIn("Ancestry", completed.stdout)
            self.assertIn("GtG1", completed.stdout)
            self.assertIn("EUR", completed.stdout)
            self.assertIn("1m40.000s", completed.stdout)
            self.assertIn("3m20.000s", completed.stdout)
            self.assertIn("5m00.000s", completed.stdout)


if __name__ == "__main__":
    unittest.main()
