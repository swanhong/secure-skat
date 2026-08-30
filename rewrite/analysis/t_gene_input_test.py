import csv
import tempfile
import unittest
from pathlib import Path

from rewrite.preprocessing.input import read_gene_panel
from rewrite.testdata.aou.prepare_aou import normalize_annotation


class GeneInputTest(unittest.TestCase):
    def test_normalize_annotation_omits_empty_gene_ids(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            source = root / "source.tsv"
            pvar = root / "genotype.pvar"
            annotation = root / "annotation.tsv"
            gene_panel = root / "gene_panel.tsv"

            source.write_text(
                "variant_key\tgene_id\tgene_symbol\tannotation\tgnomad_af\n"
                "21:100:A:G\t\t\tpLoF\t0.001\n"
                "21:200:C:T\tENSG1\tGENE1\tpLoF\t0.002\n"
            )
            pvar.write_text(
                "#CHROM\tPOS\tID\tREF\tALT\n"
                "21\t100\t21:100:A:G\tA\tG\n"
                "21\t200\t21:200:C:T\tC\tT\n"
            )
            annotation.write_text("stale annotation\n")
            gene_panel.write_text(
                "gene_id\tgene_symbol\tchromosome\torder_index\n"
                "\t\t21\t0\n"
            )

            normalize_annotation(
                source=source,
                pvar_path=pvar,
                annotation_output=annotation,
                gene_panel_output=gene_panel,
                chromosome=21,
                maf_column="gnomad_af",
            )

            with annotation.open(newline="") as input_file:
                annotation_rows = list(
                    csv.DictReader(input_file, delimiter="\t")
                )
            with gene_panel.open(newline="") as input_file:
                gene_rows = list(csv.DictReader(input_file, delimiter="\t"))

            self.assertEqual(
                [row["gene_id"] for row in annotation_rows],
                ["ENSG1"],
            )
            self.assertEqual(
                [row["gene_id"] for row in gene_rows],
                ["ENSG1"],
            )
            self.assertEqual(gene_rows[0]["order_index"], "0")

    def test_read_gene_panel_rejects_empty_gene_id(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            panel = Path(temporary_directory) / "gene_panel.tsv"
            panel.write_text(
                "gene_id\tgene_symbol\tchromosome\torder_index\n"
                "\t\t21\t0\n"
            )

            with self.assertRaisesRegex(ValueError, "gene_id is required"):
                read_gene_panel(panel)


if __name__ == "__main__":
    unittest.main()
