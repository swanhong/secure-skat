import csv
import tempfile
import unittest
from pathlib import Path

from rewrite.analysis.plots import write_all_plots
from rewrite.analysis.results import (
    RESULT_FIELDS,
    aggregate_results,
    join_results,
    read_rows,
    write_rows,
)


def write_table(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class AnalysisTest(unittest.TestCase):
    def test_join_handles_nonconvergence_and_zero_reference(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            secure = root / "secure.csv"
            reference = root / "reference.csv"
            output = root / "joined.csv"
            write_table(secure, [self.secure_row("gene1", "0.2", "0.3")])
            write_table(reference, [{
                "chromosome": 22,
                "gene_index": 0,
                "gene_id": "gene1",
                "phenotype_index": 0,
                "r_burden_p": 0,
                "r_skat_davies_p": "NA",
                "r_skat_davies_converged": 0,
            }])

            rows = join_results(secure, reference, output)

            self.assertEqual(rows[0]["burden_abs_error"], "0.20000000000000001")
            self.assertEqual(rows[0]["burden_rel_error"], "NA")
            self.assertEqual(rows[0]["skat_abs_error"], "NA")
            self.assertEqual(rows[0]["skat_rel_error"], "NA")

    def test_two_chromosome_aggregate_and_plots(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            chromosomes = root / "chromosomes"
            final = root / "final"

            for chromosome, converged in ((21, "true"), (22, "false")):
                row = self.joined_row(chromosome, converged)
                chromosome_dir = chromosomes / f"chr{chromosome}"
                write_rows(chromosome_dir / "gene_results.csv", [row], RESULT_FIELDS)
                (chromosome_dir / "_SUCCESS").touch()

            rows = aggregate_results(chromosomes, final, chromosomes=(21, 22))
            write_all_plots(final / "all_gene_results.csv", final / "plots")

            self.assertEqual([row["global_gene_index"] for row in rows], ["0", "1"])
            self.assertEqual(len(read_rows(final / "error_summary.csv")), 6)
            self.assertEqual(len(list((final / "plots").glob("scatter_*.png"))), 4)
            self.assertEqual(len(list((final / "plots").glob("manhattan_*.png"))), 4)

    @staticmethod
    def secure_row(gene_id: str, burden: str, skat: str) -> dict[str, object]:
        return {
            "chromosome": 22,
            "gene_index": 0,
            "gene_id": gene_id,
            "gene_symbol": "GENE1",
            "gene_order": 0,
            "phenotype_index": 0,
            "phenotype_name": "trait",
            "secure_burden_p": burden,
            "secure_skat_wh_p": skat,
        }

    @staticmethod
    def joined_row(chromosome: int, converged: str) -> dict[str, object]:
        return {
            "chromosome": chromosome,
            "gene_index": 0,
            "gene_id": f"gene{chromosome}",
            "gene_symbol": f"GENE{chromosome}",
            "gene_order": 0,
            "phenotype_index": 0,
            "phenotype_name": "trait",
            "secure_burden_p": 0.2,
            "r_burden_p": 0.21,
            "burden_abs_error": 0.01,
            "burden_rel_error": 0.01 / 0.21,
            "secure_skat_wh_p": 0.3,
            "r_skat_davies_p": 0.31,
            "r_skat_davies_converged": converged,
            "skat_abs_error": 0.01 if converged == "true" else "NA",
            "skat_rel_error": 0.01 / 0.31 if converged == "true" else "NA",
        }


if __name__ == "__main__":
    unittest.main()
