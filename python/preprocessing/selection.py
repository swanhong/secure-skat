
import random
from collections.abc import Collection, Mapping, Sequence

from .model import GeneRef, GeneVariants, PrepInputs


def split_sample_ids(
    sample_ids: Sequence[str], ancestries: Mapping[str, str],
    num_parties: int, party: int,
) -> tuple[str, ...]:
    """Partition a common source once, independently of analysis/sample filters."""
    if num_parties != 2 or party not in (1, 2):
        raise ValueError("debug sample splitting requires 2 partitions and party 1 or 2")
    if len(set(sample_ids)) != len(sample_ids):
        raise ValueError("duplicate source sample IDs")
    # Alternating across ancestry boundaries also balances odd-sized groups.
    ordered = sorted(sample_ids, key=lambda sample: (ancestries.get(sample, ""), sample))
    selected = set(ordered[party - 1::num_parties])
    return tuple(sample for sample in sample_ids if sample in selected)


def select_random_gene_groups(
        gene_groups: Sequence[GeneVariants],
        gene_count: int,
) -> tuple[GeneVariants, ...]:
    if gene_count < 1:
        raise ValueError("random count must be positive")
    candidates = tuple(
        group for group in gene_groups if group.variants
    )
    if not candidates:
        raise ValueError("no genes with variants available for random selection")
    if gene_count > len(candidates):
        print(f"Requested {gene_count} genes; only {len(candidates)} have variants. Using all available genes.")
        return candidates

    selected_ids = {
        group.gene.gene_id
        for group in random.Random().sample(
            candidates,
            gene_count,
        )
    }

    return tuple(
        group for group in candidates if group.gene.gene_id in selected_ids
    )

def select_file_genes(
    inputs: PrepInputs,
    chromosome: str,
    selected_genes: Collection[GeneRef],
) -> tuple[GeneRef, ...]:
    requested = tuple(
        gene
        for gene in selected_genes
        if gene.chromosome == chromosome
    )

    requested_ids = set()
    for gene in requested:
        if gene.gene_id in requested_ids:
            raise ValueError(
                f"duplicate selected gene {gene.gene_id} "
                f"on chromosome {chromosome}"
            )
        requested_ids.add(gene.gene_id)

    available = tuple(
        gene
        for gene in inputs.gene_panel
        if (
            gene.chromosome == chromosome
            and gene.gene_id in requested_ids
        )
    )

    available_ids = {
        gene.gene_id
        for gene in available
    }
    missing = requested_ids - available_ids
    if missing:
        raise ValueError(
            f"selected genes not found on chromosome {chromosome}: "
            f"{', '.join(sorted(missing))}"
        )

    return available
