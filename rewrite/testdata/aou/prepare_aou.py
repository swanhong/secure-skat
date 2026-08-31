from __future__ import annotations

import argparse
import csv
import os
import shlex
import shutil
import subprocess
from pathlib import Path


DEFAULT_GENOTYPE = (
    "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/"
    "snpindel/exome/pgen/exome.chr{chromosome}"
)
DEFAULT_PHENOTYPE = (
    "gs://gwas-data-wgs-wb-jaunty-blueberry-8679/pheno/"
    "v9_final_lipid_med_corrected_short_read_tot.csv"
)
DEFAULT_ANCESTRY = (
    "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/"
    "snpindel/aux/ancestry/ancestry_preds.tsv"
)
REQUIRED_ANNOTATION_COLUMNS = (
    "variant_key",
    "gene_id",
    "gene_symbol",
)
MISSING_FREQUENCIES = {"", ".", "NA", "NaN", "nan"}


def run(command: list[str]) -> None:
    print("+", shlex.join(command), flush=True)
    subprocess.run(command, check=True)


def billing_project() -> str:
    project = os.environ.get("GOOGLE_PROJECT") or os.environ.get(
        "GOOGLE_CLOUD_PROJECT",
        "",
    )
    if not project:
        raise ValueError(
            "--billing-project is required when GOOGLE_PROJECT and "
            "GOOGLE_CLOUD_PROJECT are unset"
        )
    return project


def copy_gcs_files(
    sources: tuple[str, ...],
    destination: Path,
    project: str,
) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    missing = tuple(
        source
        for source in sources
        if not (destination / source.rsplit("/", 1)[-1]).exists()
    )
    if not missing:
        for source in sources:
            print("exists:", destination / source.rsplit("/", 1)[-1])
        return

    run([
        "gsutil",
        "-u",
        project,
        "-m",
        "cp",
        *missing,
        f"{destination}/",
    ])


def copy_input(source: str, destination: Path, project: str) -> Path:
    if destination.exists():
        print("exists:", destination)
        return destination

    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.startswith("gs://"):
        run([
            "gsutil",
            "-u",
            project,
            "cp",
            source,
            str(destination),
        ])
    else:
        local_source = Path(os.path.expandvars(source)).expanduser()
        shutil.copy2(local_source, destination)
        print(f"copied: {local_source} -> {destination}")
    return destination


def localize_genotype(
    source_template: str,
    chromosome: int,
    output_dir: Path,
    project: str,
) -> Path:
    source_prefix = source_template.format(chromosome=chromosome)
    if not source_prefix.startswith("gs://"):
        raise ValueError("AoU genotype input must be a gs:// prefix")

    genotype_dir = output_dir / "genotype"
    copy_gcs_files(
        tuple(f"{source_prefix}.{extension}" for extension in ("pgen", "pvar", "psam")),
        genotype_dir,
        project,
    )
    return genotype_dir / source_prefix.rsplit("/", 1)[-1]


def key_pvar(pgen_prefix: Path, plink2_bin: str) -> None:
    pvar_path = Path(f"{pgen_prefix}.pvar")
    with pvar_path.open() as file:
        for line in file:
            if line.startswith("#"):
                continue
            chromosome, position, variant_id, ref, alt, *_ = line.rstrip().split("\t")
            expected_id = f"{chromosome}:{position}:{ref}:{alt}"
            if variant_id == expected_id:
                print("exists: keyed", pvar_path)
                return
            break

    keyed_prefix = pgen_prefix.with_name(f"{pgen_prefix.name}.keyed")
    run([
        plink2_bin,
        "--pfile",
        str(pgen_prefix),
        "--set-all-var-ids",
        "@:#:$r:$a",
        "--new-id-max-allele-len",
        "1000",
        "--make-just-pvar",
        "--out",
        str(keyed_prefix),
    ])
    keyed_pvar = Path(f"{keyed_prefix}.pvar")
    keyed_pvar.replace(pvar_path)
    Path(f"{keyed_prefix}.log").unlink(missing_ok=True)


def variant_alleles(variant_key: str) -> tuple[int, str, str]:
    fields = variant_key.split(":")
    if len(fields) < 3:
        raise ValueError(f"invalid annotation variant key: {variant_key}")
    return int(fields[-3]), fields[-2], fields[-1]


def annotation_keys(
    annotation_path: Path,
    maf_column: str,
) -> set[tuple[int, str, str]]:
    with annotation_path.open(newline="") as file:
        reader = csv.DictReader(file, delimiter="\t")
        columns = set(reader.fieldnames or ())
        required = {*REQUIRED_ANNOTATION_COLUMNS, maf_column}
        missing = required - columns
        if missing:
            raise ValueError(
                f"{annotation_path} is missing columns: {sorted(missing)}"
            )
        return {
            variant_alleles(row["variant_key"])
            for row in reader
        }


def pvar_ids(
    pvar_path: Path,
    wanted: set[tuple[int, str, str]],
) -> dict[tuple[int, str, str], str]:
    with pvar_path.open() as file:
        for line in file:
            if line.startswith("#CHROM"):
                columns = line.lstrip("#").rstrip().split("\t")
                break
        else:
            raise ValueError(f"PVAR header not found: {pvar_path}")

        position_column = columns.index("POS")
        id_column = columns.index("ID")
        ref_column = columns.index("REF")
        alt_column = columns.index("ALT")
        ids = {}
        for line in file:
            fields = line.rstrip().split("\t")
            key = (
                int(fields[position_column]),
                fields[ref_column],
                fields[alt_column],
            )
            if key in wanted:
                ids[key] = fields[id_column]
        return ids


def minor_allele_frequency(value: str) -> float:
    if value.strip() in MISSING_FREQUENCIES:
        return 0.0
    frequency = float(value)
    if frequency < 0 or frequency > 1:
        raise ValueError(f"allele frequency must be between 0 and 1: {value}")
    return min(frequency, 1.0 - frequency)


def normalize_annotation(
    source: Path,
    pvar_path: Path,
    annotation_output: Path,
    gene_panel_output: Path,
    chromosome: int,
    maf_column: str,
) -> None:
    if annotation_output.exists() and gene_panel_output.exists():
        with gene_panel_output.open(newline="") as file:
            reader = csv.DictReader(file, delimiter="\t")
            has_empty_gene_id = any(
                not row["gene_id"].strip() for row in reader
            )
        if not has_empty_gene_id:
            print("exists:", annotation_output)
            print("exists:", gene_panel_output)
            return
        print("regenerating gene panel with empty gene_id:", gene_panel_output)

    wanted = annotation_keys(source, maf_column)
    ids = pvar_ids(pvar_path, wanted)
    annotation_output.parent.mkdir(parents=True, exist_ok=True)
    gene_panel_output.parent.mkdir(parents=True, exist_ok=True)

    genes: dict[str, tuple[str, int]] = {}
    written = 0
    with (
        source.open(newline="") as source_file,
        annotation_output.open("w", newline="") as output_file,
    ):
        reader = csv.DictReader(source_file, delimiter="\t")
        source_columns = reader.fieldnames or []
        extra_columns = [
            column
            for column in source_columns
            if column not in REQUIRED_ANNOTATION_COLUMNS and column != "MAF"
        ]
        output_columns = [*REQUIRED_ANNOTATION_COLUMNS, *extra_columns, "MAF"]
        writer = csv.DictWriter(
            output_file,
            fieldnames=output_columns,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()

        for row in reader:
            alleles = variant_alleles(row["variant_key"])
            variant_id = ids.get(alleles)
            if variant_id is None:
                continue

            gene_id = row["gene_id"].strip()
            if not gene_id:
                continue

            row["variant_key"] = variant_id
            row["gene_id"] = gene_id
            row["MAF"] = format(
                minor_allele_frequency(row[maf_column]),
                ".17g",
            )
            writer.writerow(row)
            written += 1

            gene_symbol = row["gene_symbol"]
            position = alleles[0]
            previous = genes.get(gene_id)
            if previous is None or position < previous[1]:
                genes[gene_id] = (gene_symbol, position)

    if written == 0:
        annotation_output.unlink()
        raise ValueError(
            f"no annotation variants from {source} matched {pvar_path}"
        )

    ordered_genes = sorted(
        genes.items(),
        key=lambda item: (item[1][1], item[0]),
    )
    with gene_panel_output.open("w", newline="") as file:
        writer = csv.writer(file, delimiter="\t", lineterminator="\n")
        writer.writerow(("gene_id", "gene_symbol", "chromosome", "order_index"))
        for order_index, (gene_id, (gene_symbol, _)) in enumerate(ordered_genes):
            writer.writerow((gene_id, gene_symbol, chromosome, order_index))

    print(f"created: {annotation_output}, ({written} annotations)")
    print(f"created: {gene_panel_output}, ({len(ordered_genes)} genes)")


def prepare_aou(args: argparse.Namespace) -> None:
    output_dir = args.output_dir
    project = args.billing_project or billing_project()
    plink2_bin = os.path.expandvars(args.plink2_bin)
    pgen_prefixes = {}

    for chromosome in args.chromosome:
        prefix = localize_genotype(
            source_template=args.genotype,
            chromosome=chromosome,
            output_dir=output_dir,
            project=project,
        )
        key_pvar(prefix, plink2_bin)
        pgen_prefixes[chromosome] = prefix

    copy_input(
        args.phenotype,
        output_dir / "phenotype.csv",
        project,
    )
    copy_input(
        args.ancestry,
        output_dir / "ancestry_preds.tsv",
        project,
    )

    for chromosome in args.chromosome:
        annotation_source = Path(
            os.path.expandvars(
                args.annotation.format(chromosome=chromosome)
            )
        ).expanduser()
        normalize_annotation(
            source=annotation_source,
            pvar_path=Path(f"{pgen_prefixes[chromosome]}.pvar"),
            annotation_output=output_dir / "annotation" / f"chr{chromosome}.tsv",
            gene_panel_output=output_dir / "gene_panel" / f"chr{chromosome}.tsv",
            chromosome=chromosome,
            maf_column=args.maf_column,
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--chromosome",
        nargs="+",
        type=int,
        required=True,
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "generated",
    )
    parser.add_argument("--genotype", default=DEFAULT_GENOTYPE)
    parser.add_argument("--phenotype", default=DEFAULT_PHENOTYPE)
    parser.add_argument("--ancestry", default=DEFAULT_ANCESTRY)
    parser.add_argument(
        "--annotation",
        default="~/fed_prep_out/chr{chromosome}_annotation.tsv",
    )
    parser.add_argument("--maf-column", default="gnomad_af")
    parser.add_argument("--billing-project")
    parser.add_argument(
        "--plink2-bin",
        default=os.environ.get("PLINK2", "~/plink2"),
    )
    args = parser.parse_args()

    if len(set(args.chromosome)) != len(args.chromosome):
        parser.error("chromosomes must not contain duplicates")
    if any(chromosome < 1 or chromosome > 22 for chromosome in args.chromosome):
        parser.error("chromosomes must be between 1 and 22")

    args.output_dir = args.output_dir.expanduser()
    args.plink2_bin = str(Path(args.plink2_bin).expanduser())
    prepare_aou(args)


if __name__ == "__main__":
    main()
