from collections.abc import Collection, Mapping
from dataclasses import dataclass
from functools import partial
from pathlib import Path
from typing import Literal, TextIO
import json
import sys

from .input import (
    load_inputs,
    load_sample_inputs,
    read_gene_panel,
)
from .model import (
    GeneRef,
    PrepInputs,
    GeneVariants,
)
from .pipeline import (
    GenotypeExtractor,
    prepare_blocks,
    select_gene_variants,
    select_rows,
)
from .plink import plink_extract
from .output import write_selected_genes
from .selection import (
    select_file_genes,
    select_random_gene_groups,
)

@dataclass(frozen=True)
class GeneSelectionRequest:
    mode: Literal["random", "file", "all"]
    per_chromosome: int
    seed: int
    path: Path | None


@dataclass(frozen=True)
class PrepareRequest:
    run_dir: Path
    chromosomes: tuple[int, ...]

    genotype: str
    gene_panel: str
    annotation: str
    phenotype: Path
    covariate: Path
    ancestry: Path
    plink2_bin: str

    phenotype_id_column: str
    covariate_id_column: str
    covariate_column: str
    ancestry_id_column: str
    ancestry_column: str
    phenotype_columns: tuple[str, ...]
    ancestries: tuple[str, ...]
    num_cov: int
    mask: Mapping[str, str | Collection[str]]
    max_maf: float | None

    samples_per_cohort: int
    sample_seed: int
    role_seed: int

    gene_selection: GeneSelectionRequest


def chromosome_path(
    template: str,
    chromosome: int,
) -> Path:
    return Path(
        template.replace(
            "{chromosome}",
            str(chromosome),
        )
    )


def select_gene_groups(
    request: GeneSelectionRequest,
    inputs: PrepInputs,
    chromosome: int,
    mask: Mapping[str, str | Collection[str]],
    max_maf: float | None,
    file_genes: Collection[GeneRef],
) -> tuple[GeneVariants, ...]:
    chromosome_name = str(chromosome)

    if request.mode == "random":
        all_groups = select_gene_variants(
            gene_panel=inputs.gene_panel,
            variants=inputs.variants,
            annotations=inputs.annotations,
            annotation_columns=inputs.annotation_columns,
            chromosome=chromosome_name,
            gene_selection="all",
            mask=mask,
            max_maf=max_maf,
        )
        return select_random_gene_groups(
            gene_groups=all_groups,
            gene_count=request.per_chromosome,
            seed=request.seed,
        )

    if request.mode == "file":
        genes = select_file_genes(
            inputs=inputs,
            chromosome=chromosome_name,
            selected_genes=file_genes,
        )
        gene_selection = tuple(gene.gene_id for gene in genes)
    elif request.mode == "all":
        gene_selection = "all"
    else:
        raise ValueError(
            f"unsupported gene selection mode: {request.mode}"
        )
    return select_gene_variants(
        gene_panel=inputs.gene_panel,
        variants=inputs.variants,
        annotations=inputs.annotations,
        annotation_columns=inputs.annotation_columns,
        chromosome=chromosome_name,
        gene_selection=gene_selection,
        mask=mask,
    )


def prepare_chromosomes(
    request: PrepareRequest,
    extractor: GenotypeExtractor | None = None,
) -> tuple[GeneRef, ...]:
    if not request.chromosomes:
        raise ValueError("at least one chromosome is required")

    selected_genes_path = request.run_dir / "selected_genes.tsv"
    if request.gene_selection.mode != "all" and selected_genes_path.exists():
        raise ValueError(
            f"selected genes file already exists: {selected_genes_path}"
        )

    if extractor is None:
        extractor = partial(
            plink_extract,
            plink2_bin=request.plink2_bin,
        )
    sample_inputs = load_sample_inputs(
        phenotype_path=request.phenotype,
        covariate_path=request.covariate,
        ancestry_path=request.ancestry,
        phenotype_id_column=request.phenotype_id_column,
        phenotype_columns=request.phenotype_columns,
        covariate_id_column=request.covariate_id_column,
        covariate_column=request.covariate_column,
        ancestry_id_column=request.ancestry_id_column,
        ancestry_column=request.ancestry_column,
        num_cov=request.num_cov,
    )

    if request.samples_per_cohort == 0:
        samples_per_cohort: int | Literal["all"] = "all"
    else:
        samples_per_cohort = request.samples_per_cohort

    file_genes: tuple[GeneRef, ...] = ()
    if request.gene_selection.mode == "file":
        if request.gene_selection.path is None:
            raise ValueError(
                "gene selection path is required in file mode"
            )
        file_genes = read_gene_panel(
            request.gene_selection.path
        )

    selected_genes = []
    rows_by_ancestry = {}
    print(
        f"Running preprocessing for "
        f"{len(request.chromosomes)} chromosomes and "
        f"{len(request.ancestries)} ancestries...",
        flush=True,
    )
    for index, chromosome in enumerate(request.chromosomes):
        inputs = load_inputs(
            pgen_prefix=chromosome_path(
                request.genotype,
                chromosome,
            ),
            gene_panel_path=chromosome_path(
                request.gene_panel,
                chromosome,
            ),
            annotation_path=chromosome_path(
                request.annotation,
                chromosome,
            ),
        )

        if index == 0:
            reference_chromosome = chromosome
            reference_psam_ids = inputs.psam_ids

            for ancestry in request.ancestries:
                rows_by_ancestry[ancestry] = select_rows(
                    psam_ids=reference_psam_ids,
                    phenotypes=sample_inputs.phenotypes,
                    covariates=sample_inputs.covariates,
                    ancestries=sample_inputs.ancestries,
                    ancestry=ancestry,
                    phenotype_columns=request.phenotype_columns,
                    samples_per_cohort=samples_per_cohort,
                    sample_seed=request.sample_seed,
                )
        elif inputs.psam_ids != reference_psam_ids:
            raise ValueError(
                f"ordered PSAM IDs differ for chromosome "
                f"{reference_chromosome} and chromosome {chromosome}"
            )

        chromosome_groups = select_gene_groups(
            request=request.gene_selection,
            inputs=inputs,
            chromosome=chromosome,
            mask=request.mask,
            max_maf=request.max_maf,
            file_genes=file_genes,
        )
        chromosome_genes = tuple(
            group.gene for group in chromosome_groups
        )
        selected_genes.extend(chromosome_genes)

        for ancestry in request.ancestries:
            rows_a, rows_b = rows_by_ancestry[ancestry]
            prepare_blocks(
                pgen_prefix=inputs.pgen_prefix,
                gene_variants=chromosome_groups,
                rows_a=rows_a,
                rows_b=rows_b,
                role_seed=request.role_seed,
                out_dir=(
                    request.run_dir
                    / "prepared"
                    / ancestry
                    / f"chr{chromosome}"
                ),
                extractor=extractor,
            )

    selected_genes = tuple(selected_genes)
    if request.gene_selection.mode != "all":
        write_selected_genes(
            path=selected_genes_path,
            genes=selected_genes,
        )
        print("Wrote selected genes to", selected_genes_path)
    return selected_genes 

def read_prepare_request(
    stream: TextIO,
) -> PrepareRequest:
    payload = json.load(stream)
    gene_selection = payload["gene_selection"]
    gene_selection_path = gene_selection["path"]

    return PrepareRequest(
        run_dir=Path(payload["run_dir"]),
        chromosomes=tuple(payload["chromosomes"]),
        genotype=payload["genotype"],
        gene_panel=payload["gene_panel"],
        annotation=payload["annotation"],
        phenotype=Path(payload["phenotype"]),
        covariate=Path(payload["covariate"]),
        ancestry=Path(payload["ancestry"]),
        plink2_bin=payload["plink2_bin"],
        phenotype_id_column=payload["phenotype_id_column"],
        covariate_id_column=payload["covariate_id_column"],
        covariate_column=payload["covariate_column"],
        ancestry_id_column=payload["ancestry_id_column"],
        ancestry_column=payload["ancestry_column"],
        phenotype_columns=tuple(payload["phenotype_columns"]),
        ancestries=tuple(
            ancestry.strip().upper()
            for ancestry in payload["ancestries"]
        ),
        num_cov=payload["num_cov"],
        mask=payload["mask"],
        max_maf=payload["max_maf"],
        samples_per_cohort=payload["samples_per_cohort"],
        sample_seed=payload["sample_seed"],
        role_seed=payload["role_seed"],
        gene_selection=GeneSelectionRequest(
            mode=gene_selection["mode"],
            per_chromosome=gene_selection["per_chromosome"],
            seed=gene_selection["seed"],
            path=(
                Path(gene_selection_path)
                if gene_selection_path
                else None
            ),
        ),
    )


def main() -> None:
    request = read_prepare_request(sys.stdin)
    prepare_chromosomes(request)


if __name__ == "__main__":
    main()
