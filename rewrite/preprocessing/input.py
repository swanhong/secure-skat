import csv
from pathlib import Path

from .model import (
    AnnotationRef,
    GeneRef,
    VariantRef,
    PrepInputs,
)

MISSING_VALUES = {"", "NA", "NaN", "nan", ".", None}


def read_pvar(
    pvar_path: Path,
) -> tuple[VariantRef, ...]:
    with pvar_path.open(newline="") as file:
        reader = csv.reader(
            file,
            delimiter="\t",
        )

        for header in reader:
            if header and header[0] == "#CHROM":
                break
        else:
            raise ValueError("PVAR header not found")

        position_column = header.index("POS")
        id_column = header.index("ID")
        filter_column = (
            header.index("FILTER") if "FILTER" in header else None
        )

        return tuple(
            VariantRef(
                key=row[id_column],
                position=int(row[position_column]),
                filter_value=(
                    "PASS" if filter_column is None else row[filter_column]
                ),
            )
            for row in reader
            if row
        )

def read_psam(
    psam_path: Path,
) -> tuple[str, ...]:
    with psam_path.open(newline="") as file:
        reader = csv.reader(
            file,
            delimiter="\t",
        )
        header = next(reader)
        header[0] = header[0].lstrip("#")
        iid_column = header.index("IID")

        return tuple(
            row[iid_column]
            for row in reader
            if row
        )

def read_gene_panel(
    panel_path: Path,
) -> tuple[GeneRef, ...]:
    with panel_path.open(newline="") as file:
        reader = csv.DictReader(
            file,
            delimiter="\t",
        )

        return tuple(
            GeneRef(
                gene_id=row["gene_id"],
                gene_symbol=row["gene_symbol"],
                chromosome=row["chromosome"],
                order_index=int(row["order_index"]),
            )
            for row in reader
        )

def read_annotations(
    annotation_path: Path,
) -> tuple[
    tuple[AnnotationRef, ...],
    tuple[str, ...],
]:
    with annotation_path.open(newline="") as file:
        reader = csv.DictReader(
            file,
            delimiter="\t",
        )
        annotation_columns = tuple(
            reader.fieldnames[3:]
        )

        annotations = tuple(
            AnnotationRef(
                variant_key=row["variant_key"],
                gene_id=row["gene_id"],
                gene_symbol=row["gene_symbol"],
                values={
                    column: row[column]
                    for column in annotation_columns
                },
            )
            for row in reader
        )

        return annotations, annotation_columns

def read_sample_table(
    table_path: Path,
    delimiter: str,
) -> dict[str, dict[str, float]]:
    with table_path.open(newline="") as file:
        reader = csv.DictReader(
            file,
            delimiter=delimiter,
        )
        if not reader.fieldnames:
            raise ValueError(f"no header in {table_path}")

        value_columns = tuple(
            column
            for column in reader.fieldnames
            if column != "IID"
        )

        return {
            row["IID"]: {
                column: (
                    float("nan")
                    if row[column] in MISSING_VALUES
                    else float(row[column])
                )
                for column in value_columns
            }
            for row in reader
        }

def load_inputs(
    pgen_prefix: Path,
    gene_panel_path: Path,
    annotation_path: Path,
    phenotype_path: Path,
    covariate_path: Path,
) -> PrepInputs:
    pvar_path = Path(f"{pgen_prefix}.pvar")
    psam_path = Path(f"{pgen_prefix}.psam")

    annotations, annotation_columns = read_annotations(
        annotation_path
    )

    return PrepInputs(
        pgen_prefix=pgen_prefix,
        psam_ids=read_psam(psam_path),
        gene_panel=read_gene_panel(gene_panel_path),
        variants=read_pvar(pvar_path),
        annotations=annotations,
        annotation_columns=annotation_columns,
        phenotypes=read_sample_table(
            phenotype_path,
            delimiter=",",
        ),
        covariates=read_sample_table(
            covariate_path,
            delimiter="\t",
        ),
    )