import csv
import json
import os
import subprocess
import tarfile
import tempfile
import unittest
from pathlib import Path

from rewrite.workbench.plots import write_all_plots
from rewrite.preprocessing.cli import (
    normalize_annotation_variant_ids,
    normalize_array_sample_table,
    validate_public_widths,
)
from rewrite.workbench.results import (
    RESULT_FIELDS,
    aggregate_results,
    join_results,
    read_rows,
    write_rows,
)


def write_table(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=list(rows[0]),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


class WorkbenchTest(unittest.TestCase):
    def test_coordinator_runs_chromosomes_in_parallel_before_aggregation(self) -> None:
        root = Path(__file__).parents[2]
        with tempfile.TemporaryDirectory() as directory:
            temporary = Path(directory)
            code_root = temporary / "code"
            r_root = temporary / "r"
            bin_dir = temporary / "bin"
            annotations = temporary / "annotations"
            results = temporary / "results"
            for path in (
                code_root / "rewrite/workbench",
                r_root / "bin",
                bin_dir,
                annotations,
            ):
                path.mkdir(parents=True)

            self.write_executable(code_root / "secure-skat", "#!/bin/bash\n")
            self.write_executable(code_root / "plink2", "#!/bin/bash\n")
            self.write_executable(
                code_root / "rewrite/workbench/run_chromosome_task.sh",
                """#!/usr/bin/env bash
printf '%s\t%s\n' "$CHROMOSOME" "$PORT_BASE" >> "$COORDINATOR_CALLS"
case "$CHROMOSOME" in
  chr20)
    while ! grep -q '^chr21' "$COORDINATOR_CALLS"; do sleep 0.01; done
    sleep 0.2
    touch "$COORDINATOR_STATE/chr20.done"
    ;;
  chr21)
    sleep 0.4
    ;;
  chr22)
    [ -f "$COORDINATOR_STATE/chr20.done" ] || exit 9
    ;;
esac
mkdir -p "$CHROMOSOME_RESULTS"
touch "$CHROMOSOME_RESULTS/_SUCCESS"
""",
            )
            self.write_executable(r_root / "bin/conda-unpack", "#!/bin/bash\n")
            self.write_executable(r_root / "bin/Rscript", "#!/bin/bash\n")
            code_archive = temporary / "code.tar.gz"
            r_archive = temporary / "r.tar.gz"
            self.write_archive(code_root, code_archive)
            self.write_archive(r_root, r_archive)

            for chromosome in (20, 21, 22):
                (annotations / f"chr{chromosome}_annotation.tsv").write_text(
                    "variant_key\tgene_id\tgene_symbol\tannotation\n"
                )
            phenotype = temporary / "phenotype.csv"
            covariate = temporary / "covariate.tsv"
            phenotype.touch()
            covariate.touch()

            self.write_executable(
                bin_dir / "gcloud",
                """#!/usr/bin/env bash
destination=${!#}
mkdir -p "$destination"
for argument in "$@"; do
  case "$argument" in
    gs://*.pgen|gs://*.pvar|gs://*.psam)
      touch "$destination/$(basename "$argument")"
      ;;
  esac
done
""",
            )
            self.write_executable(
                bin_dir / "python3",
                """#!/usr/bin/env bash
printf '%s\n' "$@" > "$AGGREGATE_ARGUMENTS"
while [ "$#" -gt 0 ]; do
  if [ "$1" = "--output-dir" ]; then
    mkdir -p "$2"
    exit 0
  fi
  shift
done
exit 1
""",
            )

            coordinator_calls = temporary / "coordinator-calls.txt"
            coordinator_state = temporary / "coordinator-state"
            coordinator_state.mkdir()
            aggregate_arguments = temporary / "aggregate-arguments.txt"
            environment = os.environ.copy()
            environment.update({
                "PATH": f"{bin_dir}:{environment['PATH']}",
                "CODE_BUNDLE": str(code_archive),
                "R_ENV_ARCHIVE": str(r_archive),
                "PHENOTYPE": str(phenotype),
                "COVARIATE": str(covariate),
                "ANNOTATIONS": str(annotations),
                "COORDINATOR_RESULTS": str(results),
                "WORK": str(temporary / "work"),
                "LOG": "-",
                "project": "test-project",
                "genotype_prefix": "gs://test/exome",
                "chromosomes": "20,21,22",
                "max_parallel_chromosomes": "2",
                "COORDINATOR_CALLS": str(coordinator_calls),
                "COORDINATOR_STATE": str(coordinator_state),
                "AGGREGATE_ARGUMENTS": str(aggregate_arguments),
            })

            result = subprocess.run(
                [
                    "bash",
                    str(root / "rewrite/workbench/run_coordinator_task.sh"),
                ],
                cwd=root,
                env=environment,
                capture_output=True,
                text=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                sorted(coordinator_calls.read_text().splitlines()),
                ["chr20\t18000", "chr21\t18010", "chr22\t18000"],
            )
            arguments = aggregate_arguments.read_text().splitlines()
            chromosome_start = arguments.index("--chromosomes") + 1
            self.assertEqual(arguments[chromosome_start:], ["20", "21", "22"])
            self.assertTrue((results / "_SUCCESS").exists())

    def test_submit_uses_home_plink2_and_returns_after_one_job(self) -> None:
        root = Path(__file__).parents[2]
        with tempfile.TemporaryDirectory() as directory:
            temporary = Path(directory)
            bin_dir = temporary / "bin"
            home_dir = temporary / "home"
            annotation_dir = temporary / "annotations"
            bin_dir.mkdir()
            home_dir.mkdir()
            annotation_dir.mkdir()
            for chromosome in (21, 22):
                (annotation_dir / f"chr{chromosome}_annotation.tsv").write_text(
                    "variant_key\tgene_id\tgene_symbol\tannotation\n"
                )

            plink2 = home_dir / "plink2"
            plink2.write_text("#!/usr/bin/env bash\nexit 0\n")
            plink2.chmod(0o755)
            dsub_arguments = temporary / "dsub-arguments.txt"

            self.write_executable(
                bin_dir / "go",
                """#!/usr/bin/env bash
while [ "$#" -gt 0 ]; do
  if [ "$1" = "-o" ]; then
    output=$2
    mkdir -p "$(dirname "$output")"
    printf '#!/usr/bin/env bash\\nexit 0\\n' > "$output"
    chmod +x "$output"
    exit 0
  fi
  shift
done
exit 1
""",
            )
            self.write_executable(
                bin_dir / "gcloud",
                "#!/usr/bin/env bash\nexit 0\n",
            )
            self.write_executable(
                bin_dir / "dsub",
                """#!/usr/bin/env bash
printf 'CALL\\n' >> "$DSUB_ARGUMENTS"
printf '%s\\n' "$@" >> "$DSUB_ARGUMENTS"
echo fake-coordinator-job
""",
            )

            config = temporary / "run.conf"
            config.write_text(
                f"annotation_dir={annotation_dir}\n"
                "data_bits=70\n"
                "chromosomes=21,22\n"
                "project=test-project\n"
                "service_account=test@example.com\n"
                "output_root=gs://test/secure-skat\n"
                "r_env_archive=gs://test/secure-skat/runtime/r.tar.gz\n"
            )
            environment = os.environ.copy()
            environment["PATH"] = f"{bin_dir}:/usr/bin:/bin"
            environment["HOME"] = str(home_dir)
            environment["DSUB_ARGUMENTS"] = str(dsub_arguments)

            result = subprocess.run(
                ["bash", str(root / "submit_main.sh"), str(config)],
                cwd=root,
                env=environment,
                capture_output=True,
                text=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            arguments = dsub_arguments.read_text().splitlines()
            self.assertEqual(arguments.count("CALL"), 1)
            self.assertNotIn("--wait", arguments)
            self.assertNotIn("--tasks", arguments)
            self.assertIn("max_parallel_chromosomes=2", arguments)
            self.assertIn(
                str(root / "rewrite/workbench/run_coordinator_task.sh"),
                arguments,
            )

    @staticmethod
    def write_executable(path: Path, contents: str) -> None:
        path.write_text(contents)
        path.chmod(0o755)

    @staticmethod
    def write_archive(source: Path, output: Path) -> None:
        with tarfile.open(output, "w:gz") as archive:
            for path in source.rglob("*"):
                archive.add(path, arcname=path.relative_to(source))

    def test_maps_annotation_coordinates_to_pvar_ids(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            annotation = root / "annotation.tsv"
            pvar = root / "genotype.pvar"
            output = root / "normalized.tsv"
            annotation.write_text(
                "variant_key\tgene_id\tgene_symbol\tannotation\n"
                "100:A:G\tg1\tGENE1\tpLoF\n"
                "200:C:T\tg2\tGENE2\tpLoF\n"
            )
            pvar.write_text(
                "##fileformat=PVARv1\n"
                "#CHROM\tPOS\tID\tREF\tALT\tFILTER\n"
                "chr21\t100\tchr21:100:opaque1\tA\tG\tPASS\n"
            )

            positions = normalize_annotation_variant_ids(
                annotation,
                pvar,
                output,
            )

            with output.open(newline="") as file:
                rows = list(csv.DictReader(file, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["variant_key"], "chr21:100:opaque1")
            self.assertEqual(positions, {"chr21:100:opaque1": 100})

    def test_expands_aou_ancestry_pcs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "ancestry.tsv"
            output = root / "covariates.tsv"
            columns = [f"PC{index}" for index in range(1, 17)]
            source.write_text(
                "research_id\tpca_features\tancestry_pred\n"
                "s1\t[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]\teur\n"
                "s2\t[21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37]\tafr\n"
            )

            normalize_array_sample_table(
                source,
                output,
                delimiter="\t",
                id_column="research_id",
                array_column="pca_features",
                value_columns=columns,
            )

            with output.open(newline="") as file:
                rows = list(csv.DictReader(file, delimiter="\t"))
            self.assertEqual(list(rows[0]), ["IID", *columns])
            self.assertEqual(rows[0]["IID"], "s1")
            self.assertEqual(rows[0]["PC1"], "1.0")
            self.assertEqual(rows[0]["PC16"], "16.0")
            self.assertEqual(rows[1]["PC16"], "36.0")

    def test_rejects_too_few_ancestry_pcs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "ancestry.tsv"
            source.write_text(
                "research_id\tpca_features\n"
                "s1\t[0.1,0.2]\n"
            )

            with self.assertRaisesRegex(
                ValueError,
                "pca_features for s1 has 2 values; expected at least 3",
            ):
                normalize_array_sample_table(
                    source,
                    root / "covariates.tsv",
                    delimiter="\t",
                    id_column="research_id",
                    array_column="pca_features",
                    value_columns=["PC1", "PC2", "PC3"],
                )

    def test_public_width_boundary(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            (output / "genes.txt").write_text("allowed\ntoo_wide\n")
            (output / "block_sizes.txt").write_text("4096\n4097\n")
            (output / "metadata.json").write_text(json.dumps({"chromosome": "22"}))

            errors = validate_public_widths(output)

            self.assertEqual(len(errors), 1)
            self.assertEqual(errors[0]["gene_id"], "too_wide")
            self.assertEqual(errors[0]["public_variant_count"], 4097)
            self.assertTrue((output / "validation_errors.csv").exists())

    def test_join_nonconverged_and_zero_reference(self) -> None:
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
