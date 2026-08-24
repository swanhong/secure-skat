import math
import random
from collections.abc import Collection, Mapping, Sequence, Callable
from pathlib import Path
from typing import Literal
import numpy as np

from .model import (
    AnnotationRef,
    GeneBlock,
    GeneRef,
    GeneVariants,
    PhenoCovRows,
    VariantRef,
    GenePlan,
    PrepInputs,
    PrepOptions,
)

from .output import write_outputs
from .plink import plink_extract


def require_columns(
        table: Mapping[str, Mapping[str, float]],
        columns: Sequence[str],
        label: str,
) -> None:
    row = next(iter(table.values()), None)
    if row is None:
        raise ValueError(f"{label} table is empty")

    missing = [column for column in columns if column not in row]
    if missing:
        raise ValueError(
            f"{label} columns not found: {', '.join(missing)}"
        )

def extract_genotypes(
        pgen_prefix: Path,
        rows_a: PhenoCovRows,
        rows_b: PhenoCovRows,
        plans: Sequence[GenePlan],
        extractor: Callable[
            [Path, tuple[str, ...], tuple[str, ...]],
            tuple[np.ndarray, tuple[str, ...]],
        ],
) -> tuple[
    tuple[GenePlan, ...],
    np.ndarray,
    dict[str, int],
    np.ndarray,
    dict[str, int],
]:
    """Extract party-specific genotypes and filter unavailable variants."""

    # from plan, collect variants to extract, and their roles
    keys_a = []
    keys_b = []
    seen_a = set()
    seen_b = set()
    for plan in plans:
        for variant, role in plan.variant_roles:
            if role in {"shared", "public_only"}:
                if variant.key not in seen_a:
                    keys_a.append(variant.key)
                    seen_a.add(variant.key)
            if role in {"shared", "private"}:
                if variant.key not in seen_b:
                    keys_b.append(variant.key)
                    seen_b.add(variant.key)

    # extract genotypes for the variants in the plans
    # extractor is for PLINK2 pgen/pvar/psam files
    geno_a, emitted_a = extractor(
        pgen_prefix,
        rows_a.sample_ids,
        tuple(keys_a),
    )
    geno_b, emitted_b = extractor(
        pgen_prefix,
        rows_b.sample_ids,
        tuple(keys_b),
    )

    key_column_a = {
        key: column for column, key in enumerate(emitted_a)
    }
    key_column_b = {
        key: column for column, key in enumerate(emitted_b)
    }

    emitted_a_set = set(emitted_a)
    emitted_b_set = set(emitted_b)
    filtered_plans = []
    for plan in plans:
        retained = []
        for variant, role in plan.variant_roles:
            if role == "shared":
                keep = (variant.key in emitted_a_set and variant.key in emitted_b_set)
            elif role == "public_only":
                keep = (variant.key in emitted_a_set)
            else:
                keep = (variant.key in emitted_b_set)
            if keep:
                retained.append((variant, role))
        filtered_plans.append(
            GenePlan(
                gene=plan.gene,
                variant_roles=tuple(retained),
            )
        )
    return (
        tuple(filtered_plans),
        geno_a,
        key_column_a,
        geno_b,
        key_column_b,
    )

def build_blocks(
        plans: Sequence[GenePlan],
        geno_a: np.ndarray,
        key_column_a: Mapping[str, int],
        geno_b: np.ndarray,
        key_column_b: Mapping[str, int],
) -> tuple[GeneBlock, ...]:
    """Build gene-local public and private genotype blocks."""
    blocks = []
    for plan in plans:
        public = [
            (variant, role) for variant, role in plan.variant_roles
            if role in {"shared", "public_only"}
        ]
        private = [
            variant for variant, role in plan.variant_roles
            if role == "private"
        ]

        public_columns_a = [
            key_column_a[variant.key] for variant, _ in public
        ]
        private_columns_b = [
            key_column_b[variant.key] for variant in private
        ]

        public_a = geno_a[:, public_columns_a].copy()
        public_b = np.zeros((geno_b.shape[0], len(public)), dtype=np.int8)

        for block_column, (variant, role) in enumerate(public):
            if role == "shared":
                public_b[:, block_column] = (
                    geno_b[:, key_column_b[variant.key]]
                )
        private_b = geno_b[:, private_columns_b].copy()
        blocks.append(
            GeneBlock(
                gene=plan.gene,
                public_variants=tuple(variant for variant, _ in public),
                private_variants=tuple(private),
                public_a=public_a,
                public_b=public_b,
                private_b=private_b,
            )
        )
    return tuple(blocks)

def assign_roles(
        gene_variants: Sequence[GeneVariants],
        role_seed: int=42,
) -> tuple[GenePlan, ...]:
    """Assign roles to variants for each gene.

    Args:
        gene_variants: A sequence of GeneVariants objects.
        role_seed: Random seed for role assignment.
    Returns:
        A tuple of GenePlan objects, one for each gene in gene_variants.
    """
    rng = random.Random(role_seed)
    plans = []

    for gene_group in gene_variants:
        shuffled = list(gene_group.variants)
        rng.shuffle(shuffled)

        shared_count = round(0.6 * len(shuffled))
        public_only_count = round(0.2 * len(shuffled))
        public_only_end = shared_count + public_only_count

        role_by_variant = {}

        for index, variant in enumerate(shuffled):
            if index < shared_count:
                role_by_variant[variant] = "shared"
            elif index < public_only_end:
                role_by_variant[variant] = "public_only"
            else:
                role_by_variant[variant] = "private"

        plans.append(
            GenePlan(
                gene=gene_group.gene,
                variant_roles=tuple(
                    (variant, role_by_variant[variant])
                    for variant in gene_group.variants
                ),
            )
        )
    return tuple(plans)

def select_gene_variants(
        gene_panel: Sequence[GeneRef],
        variants: Sequence[VariantRef],
        annotations: Sequence[AnnotationRef],
        annot_columns: Collection[str],
        chr: str,
        gene_selection: Collection[str] | Literal["all"],
        mask: Mapping[str, str | Collection[str]],
) -> tuple[GeneVariants, ...]:
    """Select variants for genes in a gene panel.

    Args:
        gene_panel: Genes to consider.
        variants: Variants to consider.
        annotations: Annotations for the variants.
        annot_columns: Annotation columns to use for filtering.
        chr: Chromosome to filter on.
        gene_selection: Gene IDs to select. If "all", all genes are selected.
        mask: chosen annotation values to filter on.
            input example: {"LoF": "HC", "consequence": {"missense_variant", "stop_gained"}}

    Returns:
        A tuple of GeneVariants objects
    """
    for column in mask:
        if column not in set(annot_columns):
            raise ValueError(f"mask column not found: {column}")

    accepted_values_by_column = {
        column: (
            {accepted_values}
            if isinstance(accepted_values, str)
            else set(accepted_values)
        ) for column, accepted_values in mask.items()
    }

    # A bare str is a Collection[str], so `in` below would substring-match.
    if isinstance(gene_selection, str) and gene_selection != "all":
        raise TypeError(
            "gene_selection must be 'all' or a collection of gene ids"
        )

    selected_genes = tuple(
        gene for gene in gene_panel
        if gene.chromosome == chr and (
            gene_selection == "all" or gene.gene_id in gene_selection
        )
    )
    selected_gene_ids = {
        gene.gene_id for gene in selected_genes
    }

    annotations_by_variant = {}
    seen_occurrences = set()
    for annotation in annotations:
        if annotation.gene_id not in selected_gene_ids:
            continue

        occurrence = (annotation.gene_id, annotation.variant_key)
        if occurrence in seen_occurrences:
            raise ValueError(
                f"duplicate annotation row for gene {annotation.gene_id} "
                f"and variant {annotation.variant_key}"
            )
        seen_occurrences.add(occurrence)

        annotations_by_variant.setdefault(
            annotation.variant_key,
            [],
        ).append(annotation)

    variants_by_gene = {
        gene.gene_id:[] for gene in selected_genes
    }

    for variant in variants:
        if variant.filter_value not in {"PASS", "."}:
            continue

        for annotation in annotations_by_variant.get(variant.key, ()):
            if any(
                annotation.values[column] not in accepted_values
                for column, accepted_values in accepted_values_by_column.items()
            ):
                continue
            variants_by_gene[annotation.gene_id].append(variant)

    return tuple(
        GeneVariants(
            gene=gene,
            variants=tuple(variants_by_gene[gene.gene_id]),
        )
        for gene in selected_genes
    )

def select_rows(
    psam_ids: Sequence[str],
    phenotypes: Mapping[str, Mapping[str, float]],
    covariates: Mapping[str, Mapping[str, float]],
    phenotype_columns: Sequence[str],
    covariate_columns: Sequence[str],
    samples_per_cohort: int | Literal["all"],
    sample_seed: int=42,
) -> tuple[PhenoCovRows, PhenoCovRows]:
    # filter samples
    #   - have both have both phenotypes and covariates,
    #   - all values are finite (not nan, inf, -inf, ...)
    require_columns(phenotypes, phenotype_columns, "phenotype")
    require_columns(covariates, covariate_columns, "covariate")

    eligible = []
    for sample_id in psam_ids:
        if sample_id not in phenotypes or sample_id not in covariates:
            continue

        phenotype_values = tuple(
            phenotypes[sample_id][column]
            for column in phenotype_columns
        )
        covariate_values = tuple(
            covariates[sample_id][column]
            for column in covariate_columns
        )

        if not all(
            math.isfinite(value)
            for value in phenotype_values
        ):
            continue

        if not all(
            math.isfinite(value)
            for value in covariate_values
        ):
            continue

        eligible.append(
            (
                sample_id,
                covariate_values,
                phenotype_values,
            )
        )

    random.Random(sample_seed).shuffle(eligible)

    if samples_per_cohort == "all":
        if len(eligible) < 2:
            raise ValueError(
                f"need 2 eligible samples, found {len(eligible)}"
            )

        n_b = len(eligible) // 2
        n_a = len(eligible) - n_b
    else:
        if samples_per_cohort <= 0:
            raise ValueError(
                "samples_per_cohort must be positive or 'all'"
            )

        required = 2 * samples_per_cohort
        if len(eligible) < required:
            raise ValueError(
                f"need {required} eligible samples, "
                f"found {len(eligible)}"
            )

        n_a = n_b = samples_per_cohort

    selected_a = eligible[:n_a]
    selected_b = eligible[n_a:n_a + n_b]

    rows_a = PhenoCovRows(
        sample_ids=tuple(row[0] for row in selected_a),
        covariates=tuple(row[1] for row in selected_a),
        phenotypes=tuple(row[2] for row in selected_a),
    )
    rows_b = PhenoCovRows(
        sample_ids=tuple(row[0] for row in selected_b),
        covariates=tuple(row[1] for row in selected_b),
        phenotypes=tuple(row[2] for row in selected_b),
    )

    return rows_a, rows_b

def prepare_blocks(
        inputs: PrepInputs,
        options: PrepOptions,
        extractor: Callable[
            [Path, tuple[str, ...], tuple[str, ...],],
            tuple[np.ndarray, tuple[str, ...]],
        ] = plink_extract,
) -> Path:
    gene_variants = select_gene_variants(
        gene_panel=inputs.gene_panel,
        variants=inputs.variants,
        annotations=inputs.annotations,
        annot_columns=inputs.annotation_columns,
        chr=options.chromosome,
        gene_selection=options.gene_selection,
        mask=options.mask,
    )

    rows_a, rows_b = select_rows(
        psam_ids=inputs.psam_ids,
        phenotypes=inputs.phenotypes,
        covariates=inputs.covariates,
        phenotype_columns=options.phenotype_columns,
        covariate_columns=options.covariate_columns,
        samples_per_cohort=options.samples_per_cohort,
        sample_seed=options.sample_seed,
    )

    plans = assign_roles(
        gene_variants=gene_variants,
        role_seed=options.role_seed,
    )

    (
        extracted_plans,
        geno_a,
        key_column_a,
        geno_b,
        key_column_b,
    ) = extract_genotypes(
        pgen_prefix=inputs.pgen_prefix,
        rows_a=rows_a,
        rows_b=rows_b,
        plans=plans,
        extractor=extractor,
    )

    blocks = build_blocks(
        plans=extracted_plans,
        geno_a=geno_a,
        key_column_a=key_column_a,
        geno_b=geno_b,
        key_column_b=key_column_b,
    )

    write_outputs(
        blocks=blocks,
        rows_a=rows_a,
        rows_b=rows_b,
        out_dir=options.out_dir,
    )

    return options.out_dir
