import csv
import subprocess
import tempfile
import unittest
from pathlib import Path


class SummarizeMetricsTest(unittest.TestCase):
    def test_prints_configured_timing_communication_and_r_squared(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            run_dir = Path(temporary_directory)
            config_path = run_dir / "run.aou.conf"
            config_path.write_text(
                f'run_dir = "{run_dir.as_posix()}"\n'
                'ancestries = ["EUR"]\n',
                encoding="utf-8",
            )

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
                        "sent_bytes",
                        "received_bytes",
                        "sent_message_count",
                        "received_message_count",
                    ],
                )
                writer.writeheader()
                for stage, sent, received, sent_messages, received_messages in (
                    ("sample_count_exchange", 100, 200, 1, 2),
                    ("collective_setup", 200, 300, 2, 3),
                    ("null_model", 300, 400, 3, 4),
                ):
                    writer.writerow(
                        {
                            "process": "party1",
                            "ancestry": "EUR",
                            "stage": stage,
                            "chromosome": "",
                            "count": 1,
                            "duration_seconds": 1,
                            "sent_bytes": sent,
                            "received_bytes": received,
                            "sent_message_count": sent_messages,
                            "received_message_count": received_messages,
                        }
                    )

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
                        is_total = stage == "chromosome_total"
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
                                "sent_bytes": chromosome * 1024
                                if is_total
                                else 999_999,
                                "received_bytes": chromosome * 2048
                                if is_total
                                else 999_999,
                                "sent_message_count": chromosome
                                if is_total
                                else 999,
                                "received_message_count": chromosome * 2
                                if is_total
                                else 999,
                            }
                        )

            comparison_dir = run_dir / "comparison" / "EUR"
            comparison_dir.mkdir(parents=True)
            comparison_path = comparison_dir / "all_comparison.csv"
            with comparison_path.open("w", newline="") as output:
                fieldnames = [
                    "chromosome",
                    "gene_index",
                    "gene_id",
                    "phenotype_index",
                    "phenotype_name",
                    "secure_burden_p",
                    "r_burden_p",
                    "secure_skat_wh_p",
                    "r_skat_liu_p",
                    "r_skat_davies_p",
                    "r_skat_davies_converged",
                ]
                writer = csv.DictWriter(output, fieldnames=fieldnames)
                writer.writeheader()
                for chromosome, phenotype_index, phenotype_name, secure_values in (
                    (21, 0, "trait-a", (0.1, 0.01)),
                    (22, 1, "trait-b", (0.1, 0.001)),
                ):
                    for gene_index, (secure_value, reference_value) in enumerate(
                        zip(secure_values, (0.1, 0.01))
                    ):
                        writer.writerow(
                            {
                                "chromosome": chromosome,
                                "gene_index": gene_index,
                                "gene_id": f"gene-{gene_index}",
                                "phenotype_index": phenotype_index,
                                "phenotype_name": phenotype_name,
                                "secure_burden_p": secure_value,
                                "r_burden_p": reference_value,
                                "secure_skat_wh_p": secure_value,
                                "r_skat_liu_p": reference_value,
                                "r_skat_davies_p": reference_value,
                                "r_skat_davies_converged": 1,
                            }
                        )

            script = Path(__file__).resolve().parents[2] / "summarize_metrics.sh"
            completed = subprocess.run(
                ["bash", str(script), str(config_path)],
                check=True,
                capture_output=True,
                text=True,
            )

            self.assertIn(f"Config: {config_path}", completed.stdout)
            self.assertIn(f"Run directory: {run_dir}", completed.stdout)
            self.assertIn("Timing summary (party1)", completed.stdout)
            self.assertIn("Ancestry", completed.stdout)
            self.assertIn("GtG1", completed.stdout)
            self.assertIn("EUR", completed.stdout)
            self.assertIn("1m40.000s", completed.stdout)
            self.assertIn("3m20.000s", completed.stdout)
            self.assertIn("5m00.000s", completed.stdout)
            self.assertIn("Communication summary (party1)", completed.stdout)
            self.assertIn("Setup", completed.stdout)
            self.assertIn("21.00 KiB", completed.stdout)
            self.assertIn("42.00 KiB", completed.stdout)
            self.assertIn("63.00 KiB", completed.stdout)
            self.assertNotIn("976.56 KiB", completed.stdout)
            self.assertIn("R^2 summary on -log10(p)", completed.stdout)
            self.assertIn("SKAT WH vs Davies", completed.stdout)
            self.assertIn("trait-b", completed.stdout)
            self.assertIn("0.000000", completed.stdout)
            self.assertIn("-1.000000", completed.stdout)


if __name__ == "__main__":
    unittest.main()
