import io
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

from rewrite.workbench.timing import (
    TimingRecorder,
    aggregate_timing_events,
    empty_timing_row,
    print_timing_summary,
    read_timing_rows,
    write_timing_rows,
)


class TimingTest(unittest.TestCase):
    def test_recorder_uses_clock_and_records_failure(self) -> None:
        times = iter((10.0, 11.25, 20.0, 20.5))
        recorder = TimingRecorder(
            component="preprocessing",
            chromosome="chr21",
            clock=lambda: next(times),
        )

        with recorder.measure(scope="chromosome", phase="success"):
            pass
        with self.assertRaisesRegex(ValueError, "failed phase"):
            with recorder.measure(scope="chromosome", phase="failure"):
                raise ValueError("failed phase")

        self.assertEqual(recorder.rows[0]["elapsed_seconds"], "1.250000000")
        self.assertEqual(recorder.rows[0]["status"], "success")
        self.assertEqual(recorder.rows[1]["elapsed_seconds"], "0.500000000")
        self.assertEqual(recorder.rows[1]["status"], "failure")

    def test_aggregate_and_print_timing_summary(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            chromosome_root = root / "chromosomes"
            final = root / "final"
            for chromosome in (21, 22):
                rows = [
                    empty_timing_row(
                        component="workflow",
                        scope="chromosome",
                        chromosome=f"chr{chromosome}",
                        phase="preprocessing",
                        event_index=0,
                        elapsed_seconds="2.000000000",
                    ),
                    empty_timing_row(
                        component="workflow",
                        scope="chromosome",
                        chromosome=f"chr{chromosome}",
                        phase="chromosome_total",
                        event_index=1,
                        elapsed_seconds="10.000000000",
                    ),
                    empty_timing_row(
                        component="secure",
                        scope="gene",
                        chromosome=f"chr{chromosome}",
                        party="1",
                        phase="batch_total",
                        event_index=2,
                        gene_index="0",
                        gene_id=f"gene-{chromosome}",
                        measurement_kind="amortized",
                        elapsed_seconds="3.000000000",
                    ),
                    empty_timing_row(
                        component="secure",
                        scope="gene",
                        chromosome=f"chr{chromosome}",
                        party="2",
                        phase="batch_total",
                        event_index=3,
                        gene_index="0",
                        gene_id=f"gene-{chromosome}",
                        measurement_kind="amortized",
                        elapsed_seconds="4.000000000",
                    ),
                    empty_timing_row(
                        component="secure",
                        scope="batch",
                        chromosome=f"chr{chromosome}",
                        party="1",
                        phase="prepare_gene_batch",
                        event_index=4,
                        trace_mode="exact",
                        elapsed_seconds="2.000000000",
                    ),
                    empty_timing_row(
                        component="secure",
                        scope="batch",
                        chromosome=f"chr{chromosome}",
                        party="1",
                        phase="kernel_statistics",
                        event_index=5,
                        trace_mode="exact",
                        elapsed_seconds="5.000000000",
                    ),
                    empty_timing_row(
                        component="r_reference",
                        scope="gene",
                        chromosome=f"chr{chromosome}",
                        phase="r_gene_total",
                        event_index=6,
                        gene_index="0",
                        gene_id=f"gene-{chromosome}",
                        elapsed_seconds="1.500000000",
                    ),
                ]
                write_timing_rows(
                    chromosome_root / f"chr{chromosome}" / "timing_events.csv",
                    rows,
                )

            coordinator = root / "coordinator_timing.csv"
            write_timing_rows(coordinator, [empty_timing_row(
                component="coordinator",
                scope="run",
                phase="run_total",
                event_index=0,
                elapsed_seconds="12.000000000",
            )])

            summary = aggregate_timing_events(
                chromosome_root,
                final,
                chromosomes=(21, 22),
                coordinator_timing=coordinator,
            )
            output = io.StringIO()
            with redirect_stdout(output):
                print_timing_summary(summary)

            self.assertEqual(len(read_timing_rows(final / "all_timing_events.csv")), 15)
            self.assertTrue((final / "timing_summary.csv").exists())
            rollups = [
                row
                for row in read_timing_rows(final / "timing_summary.csv")
                if row["phase"] == "phenotype_independent_gene_work"
            ]
            self.assertEqual([row["elapsed_seconds"] for row in rollups], [
                "7.000000000",
                "7.000000000",
            ])
            rendered = output.getvalue()
            self.assertIn("chr21", rendered)
            self.assertIn("gene-21 4.000000 1.500000", rendered)
            self.assertIn("RUN WALL TIME: 12.000 seconds", rendered)


if __name__ == "__main__":
    unittest.main()
