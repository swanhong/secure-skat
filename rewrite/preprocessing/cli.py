from __future__ import annotations

import argparse
import csv
import json
import sys
import tempfile
from collections import defaultdict
from collections.abc import Mapping
from pathlib import Path

from rewrite.preprocessing.input import load_inputs
from rewrite.preprocessing.model import PrepOptions
from rewrite.preprocessing.pipeline import prepare_blocks
from rewrite.workbench.timing import TimingRecorder


MAX_PUBLIC_VARIANTS = 4096


def normalize_chromosome(value: str) -> str:
    value = value.strip()
    return value[3:] if value.lower().startswith("chr") else value


def variant_position(variant_key: str) -> int:
    fields = variant_key.split(":")
    if len(fields) < 2:
        raise ValueError(f"variant key has no position: {variant_key}")
    return int(fields[1])


def build_gene_panel(
    annotation_path: Path,
    chromosome: str,
    output_path: Path,
    variant_positions: Mapping[str, int] | None = None,
) -> None:
    target_chromosome = normalize_chromosome(chromosome)
    genes = {}

    with annotation_path.open(newline="") as file:
        reader = csv.DictReader(file, delimiter="\t")
        required = {"variant_key", "gene_id", "gene_symbol"}
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(
                "annotation columns not found: " + ", ".join(sorted(missing))
            )

        for row in reader:
            key = row["variant_key"]
            gene_id = row["gene_id"]
            if variant_positions is None:
                key_chromosome = normalize_chromosome(key.split(":", 1)[0])
                if key_chromosome != target_chromosome:
                    continue
                position = variant_position(key)
            else:
                position = variant_positions.get(key)
                if position is None:
                    continue
            current = genes.get(gene_id)
            if current is None or position < current[1]:
                genes[gene_id] = (row["gene_symbol"], position)

    ordered = sorted(
        genes.items(),
        key=lambda item: (item[1][1], item[0]),
    )

    with output_path.open("w", newline="") as file:
        writer = csv.writer(file, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", "gene_symbol", "chromosome", "order_index"])
        for order_index, (gene_id, (gene_symbol, _)) in enumerate(ordered):
            writer.writerow([
                gene_id,
                gene_symbol,
                target_chromosome,
                order_index,
            ])


def read_gene_selection(value: str):
    if value == "all":
        return "all"

    path = Path(value)
    return {
        line.strip().split()[0]
        for line in path.read_text().splitlines()
        if line.strip()
    }


def parse_mask(values: list[str]) -> dict[str, str | set[str]]:
    grouped = defaultdict(set)
    for value in values:
        column, separator, accepted = value.partition("=")
        if not separator or not column or not accepted:
            raise ValueError(f"mask must be COLUMN=VALUE: {value}")
        grouped[column].add(accepted)

    return {
        column: next(iter(accepted)) if len(accepted) == 1 else accepted
        for column, accepted in grouped.items()
    }


def parse_sample_count(value: str):
    return "all" if value == "all" else int(value)


def normalize_sample_table(
    input_path: Path,
    output_path: Path,
    delimiter: str,
    id_column: str,
    value_columns: list[str],
) -> None:
    with input_path.open(newline="") as source:
        reader = csv.DictReader(source, delimiter=delimiter)
        required = {id_column, *value_columns}
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(
                f"{input_path} columns not found: " + ", ".join(sorted(missing))
            )

        with output_path.open("w", newline="") as output:
            writer = csv.DictWriter(
                output,
                fieldnames=["IID", *value_columns],
                delimiter=delimiter,
                lineterminator="\n",
            )
            writer.writeheader()
            for row in reader:
                writer.writerow({
                    "IID": row[id_column],
                    **{column: row[column] for column in value_columns},
                })


def normalize_array_sample_table(
    input_path: Path,
    output_path: Path,
    delimiter: str,
    id_column: str,
    array_column: str,
    value_columns: list[str],
) -> None:
    with input_path.open(newline="") as source:
        reader = csv.DictReader(source, delimiter=delimiter)
        required = {id_column, array_column}
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(
                f"{input_path} columns not found: " + ", ".join(sorted(missing))
            )

        with output_path.open("w", newline="") as output:
            writer = csv.DictWriter(
                output,
                fieldnames=["IID", *value_columns],
                delimiter=delimiter,
                lineterminator="\n",
            )
            writer.writeheader()
            for row in reader:
                values = [
                    float(value)
                    for value in row[array_column].strip().strip("[]").split(",")
                ]
                if len(values) < len(value_columns):
                    raise ValueError(
                        f"{array_column} for {row[id_column]} has {len(values)} values; "
                        f"expected at least {len(value_columns)}"
                    )
                writer.writerow({
                    "IID": row[id_column],
                    **dict(zip(value_columns, values)),
                })


def normalize_annotation_variant_ids(
    annotation_path: Path,
    pvar_path: Path,
    output_path: Path,
) -> dict[str, int]:
    with annotation_path.open(newline="") as file:
        reader = csv.DictReader(file, delimiter="\t")
        if not reader.fieldnames or "variant_key" not in reader.fieldnames:
            raise ValueError(f"{annotation_path} has no variant_key column")
        fieldnames = reader.fieldnames
        requested_keys = {row["variant_key"] for row in reader}

    normalized_ids = {}
    positions = {}
    with pvar_path.open(newline="") as file:
        reader = csv.reader(file, delimiter="\t")
        for header in reader:
            if header and header[0] == "#CHROM":
                break
        else:
            raise ValueError("PVAR header not found")

        position_column = header.index("POS")
        id_column = header.index("ID")
        reference_column = header.index("REF")
        alternate_column = header.index("ALT")
        for row in reader:
            if not row:
                continue
            position = int(row[position_column])
            variant_id = row[id_column]
            coordinate_key = (
                f"{position}:{row[reference_column]}:{row[alternate_column]}"
            )
            for source_key in (variant_id, coordinate_key):
                if source_key not in requested_keys:
                    continue
                previous = normalized_ids.get(source_key)
                if previous is not None and previous != variant_id:
                    raise ValueError(
                        f"annotation key maps to multiple PVAR IDs: {source_key}"
                    )
                normalized_ids[source_key] = variant_id
                positions[variant_id] = position

    if not normalized_ids:
        raise ValueError(f"no {annotation_path} variants matched {pvar_path}")

    with annotation_path.open(newline="") as source, output_path.open(
        "w", newline=""
    ) as output:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.DictWriter(
            output,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in reader:
            variant_id = normalized_ids.get(row["variant_key"])
            if variant_id is None:
                continue
            row["variant_key"] = variant_id
            writer.writerow(row)

    return positions


def count_rows(path: Path) -> int:
    with path.open() as file:
        return sum(1 for line in file if line.strip())


def write_workbench_metadata(
    output_dir: Path,
    inputs,
    options: PrepOptions,
) -> None:
    gene_by_id = {gene.gene_id: gene for gene in inputs.gene_panel}
    genes = [line.strip() for line in (output_dir / "genes.txt").read_text().splitlines()]

    with (output_dir / "gene_metadata.tsv").open("w", newline="") as file:
        writer = csv.writer(file, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "gene_index",
            "gene_id",
            "gene_symbol",
            "chromosome",
            "gene_order",
        ])
        for gene_index, gene_id in enumerate(genes):
            gene = gene_by_id[gene_id]
            writer.writerow([
                gene_index,
                gene.gene_id,
                gene.gene_symbol,
                gene.chromosome,
                gene.order_index,
            ])

    (output_dir / "phenotypes.txt").write_text(
        "".join(f"{column}\n" for column in options.phenotype_columns)
    )

    mask = {
        column: [value] if isinstance(value, str) else sorted(value)
        for column, value in options.mask.items()
    }
    metadata = {
        "chromosome": normalize_chromosome(options.chromosome),
        "sample_count_a": count_rows(output_dir / "A" / "cov.txt"),
        "sample_count_b": count_rows(output_dir / "B" / "cov.txt"),
        "phenotype_columns": list(options.phenotype_columns),
        "covariate_columns": list(options.covariate_columns),
        "mask": mask,
        "sample_seed": options.sample_seed,
        "role_seed": options.role_seed,
    }
    (output_dir / "metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n"
    )


def validate_public_widths(output_dir: Path) -> list[dict[str, str | int]]:
    genes = [line.strip() for line in (output_dir / "genes.txt").read_text().splitlines()]
    counts = [
        int(line)
        for line in (output_dir / "block_sizes.txt").read_text().splitlines()
        if line.strip()
    ]
    if len(genes) != len(counts):
        raise ValueError("genes.txt and block_sizes.txt have different lengths")

    chromosome = json.loads((output_dir / "metadata.json").read_text())["chromosome"]
    errors = []
    for gene_index, (gene_id, count) in enumerate(zip(genes, counts)):
        if count > MAX_PUBLIC_VARIANTS:
            errors.append({
                "chromosome": chromosome,
                "gene_index": gene_index,
                "gene_id": gene_id,
                "public_variant_count": count,
                "limit": MAX_PUBLIC_VARIANTS,
                "error": "public variant count exceeds secure limit",
            })

    if errors:
        path = output_dir / "validation_errors.csv"
        with path.open("w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=list(errors[0]))
            writer.writeheader()
            writer.writerows(errors)

    return errors


def prepare(args: argparse.Namespace) -> int:
    chromosome = normalize_chromosome(args.chromosome)
    mask = parse_mask(args.mask)
    timing = TimingRecorder("preprocessing", chromosome=f"chr{chromosome}")

    try:
        with timing.measure(
            scope="chromosome",
            phase="preprocessing_total",
        ):
            with tempfile.TemporaryDirectory() as directory:
                temporary = Path(directory)
                gene_panel_path = temporary / "gene_panel.tsv"
                annotation_path = temporary / "annotation.tsv"
                phenotype_path = temporary / "phenotype.csv"
                covariate_path = temporary / "covariate.tsv"

                with timing.measure(
                    scope="chromosome",
                    phase="normalize_annotation_variant_ids",
                    parent_phase="preprocessing_total",
                ):
                    variant_positions = normalize_annotation_variant_ids(
                        args.annotation,
                        Path(f"{args.pgen_prefix}.pvar"),
                        annotation_path,
                    )
                with timing.measure(
                    scope="chromosome",
                    phase="build_gene_panel",
                    parent_phase="preprocessing_total",
                ):
                    build_gene_panel(
                        annotation_path,
                        chromosome,
                        gene_panel_path,
                        variant_positions,
                    )
                with timing.measure(
                    scope="chromosome",
                    phase="normalize_phenotype",
                    parent_phase="preprocessing_total",
                ):
                    normalize_sample_table(
                        args.phenotype,
                        phenotype_path,
                        ",",
                        args.phenotype_id_column,
                        args.phenotype_columns,
                    )
                with timing.measure(
                    scope="chromosome",
                    phase="normalize_covariate",
                    parent_phase="preprocessing_total",
                ):
                    if args.covariate_array_column:
                        normalize_array_sample_table(
                            args.covariate,
                            covariate_path,
                            "\t",
                            args.covariate_id_column,
                            args.covariate_array_column,
                            args.covariate_columns,
                        )
                    else:
                        normalize_sample_table(
                            args.covariate,
                            covariate_path,
                            "\t",
                            args.covariate_id_column,
                            args.covariate_columns,
                        )
                with timing.measure(
                    scope="chromosome",
                    phase="load_inputs",
                    parent_phase="preprocessing_total",
                ):
                    inputs = load_inputs(
                        pgen_prefix=args.pgen_prefix,
                        gene_panel_path=gene_panel_path,
                        annotation_path=annotation_path,
                        phenotype_path=phenotype_path,
                        covariate_path=covariate_path,
                    )

                options = PrepOptions(
                    chromosome=chromosome,
                    gene_selection=read_gene_selection(args.genes),
                    mask=mask,
                    phenotype_columns=tuple(args.phenotype_columns),
                    covariate_columns=tuple(args.covariate_columns),
                    samples_per_cohort=parse_sample_count(args.samples_per_cohort),
                    sample_seed=args.sample_seed,
                    role_seed=args.role_seed,
                    out_dir=args.out,
                )
                with timing.measure(
                    scope="chromosome",
                    phase="prepare_blocks",
                    parent_phase="preprocessing_total",
                ):
                    output_dir = prepare_blocks(
                        inputs,
                        options,
                        timing=timing,
                    )
                with timing.measure(
                    scope="chromosome",
                    phase="write_workbench_metadata",
                    parent_phase="preprocessing_total",
                ):
                    write_workbench_metadata(output_dir, inputs, options)
                timing.set_defaults(
                    sample_count_a=count_rows(output_dir / "A" / "cov.txt"),
                    sample_count_b=count_rows(output_dir / "B" / "cov.txt"),
                )

            with timing.measure(
                scope="chromosome",
                phase="validate_public_widths",
                parent_phase="preprocessing_total",
            ):
                errors = validate_public_widths(output_dir)
            for error in errors:
                print(
                    f"{error['gene_id']}: public variants "
                    f"{error['public_variant_count']} exceed {error['limit']}",
                    file=sys.stderr,
                )
            return 1 if errors else 0
    finally:
        timing.write(args.timing_output)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m rewrite.preprocessing")
    subparsers = parser.add_subparsers(dest="command", required=True)

    command = subparsers.add_parser("prepare")
    command.add_argument("--pgen-prefix", type=Path, required=True)
    command.add_argument("--annotation", type=Path, required=True)
    command.add_argument("--phenotype", type=Path, required=True)
    command.add_argument("--covariate", type=Path, required=True)
    command.add_argument("--phenotype-id-column", default="IID")
    command.add_argument("--covariate-id-column", default="IID")
    command.add_argument("--covariate-array-column")
    command.add_argument("--phenotype-columns", nargs="+", required=True)
    command.add_argument("--covariate-columns", nargs="+", required=True)
    command.add_argument("--chromosome", required=True)
    command.add_argument("--genes", default="all")
    command.add_argument("--mask", action="append", default=[], metavar="COLUMN=VALUE")
    command.add_argument("--samples-per-cohort", default="all")
    command.add_argument("--sample-seed", type=int, default=42)
    command.add_argument("--role-seed", type=int, default=42)
    command.add_argument("--out", type=Path, required=True)
    command.add_argument("--timing-output", type=Path)
    command.set_defaults(handler=prepare)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    return args.handler(args)
