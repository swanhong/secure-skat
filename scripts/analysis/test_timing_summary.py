import io
import json
import os
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))
import timing_summary


class TimingSummaryPhenotypeTest(unittest.TestCase):
    def test_parses_and_labels_phenotype_timing(self):
        row = timing_summary.parse_timing_line(
            "[timing] scope=secure mode=skat_p party=2 stage=phenotype_score_ql "
            "parent=packed_first_pass kind=phenotype status=done milliseconds=12.5 "
            "count=1 phenotype=1",
            "secure",
        )
        with tempfile.TemporaryDirectory() as directory:
            manifest = Path(directory) / "manifest.json"
            manifest.write_text(json.dumps({"phenotypes": ["LDL", "HDL"]}))
            timing_summary.add_phenotype_names([row], manifest)

        self.assertEqual(row["phenotype"], "1")
        self.assertEqual(row["phenotype_name"], "HDL")
        self.assertEqual(row["kind"], "phenotype")

    def test_prints_named_party_two_summary(self):
        rows = []
        for scope, milliseconds in (("phenotype_null_rss", "12.5"),
                                    ("phenotype_score_ql", "37.5")):
            row = {field: "" for field in timing_summary.FIELDS}
            row.update({"scope": scope, "party": "2", "kind": "phenotype",
                        "phenotype": "0", "phenotype_name": "LDL",
                        "milliseconds": milliseconds})
            rows.append(row)
        output = io.StringIO()
        with redirect_stdout(output):
            timing_summary.print_phenotype_summary(rows)
        self.assertIn("LDL", output.getvalue())
        self.assertIn("0.050s", output.getvalue())


if __name__ == "__main__":
    unittest.main()
