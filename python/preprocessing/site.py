"""Single-site preparation using an ordered, allele-specific public variant list."""
import csv

from .model import GenePlan, GeneRef, GeneVariants, VariantRef
from .pipeline import build_blocks, extract_genotypes
from .output import write_outputs

FIELDS = ("chromosome", "gene_id", "gene_symbol", "order_index", "variant_key", "position")


def validate_key(variant, chromosome):
    fields = variant.key.split(":")
    if (len(fields) != 4 or fields[0].removeprefix("chr") != chromosome
            or fields[1] != str(variant.position) or not fields[2] or not fields[3]
            or "," in fields[3]):
        raise ValueError(f"expected chromosome:position:REF:ALT key: {variant.key}")


def write_public_var_list(path, groups):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(FIELDS)
        for group in groups:
            gene = group.gene
            for variant in group.variants or (None,):
                if variant is not None:
                    validate_key(variant, gene.chromosome)
                writer.writerow((gene.chromosome, gene.gene_id, gene.gene_symbol,
                                 gene.order_index, variant.key if variant else "",
                                 variant.position if variant else ""))


def read_public_var_list(path):
    groups = []
    seen_genes = set()
    variants = []
    seen_variants = set()
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if tuple(reader.fieldnames or ()) != FIELDS:
            raise ValueError("invalid public variant list header")
        for row in reader:
            gene = GeneRef(row["gene_id"], row["gene_symbol"], row["chromosome"], int(row["order_index"]))
            if gene.chromosome not in {str(chromosome) for chromosome in range(1, 23)}:
                raise ValueError("public variant list requires autosomes 1-22")
            if not groups or groups[-1].gene != gene:
                if gene.gene_id in seen_genes:
                    raise ValueError(f"duplicate/noncontiguous gene: {gene.gene_id}")
                seen_genes.add(gene.gene_id)
                groups.append(GeneVariants(gene, ()))
                variants.append([])
                seen_variants = set()
            if row["variant_key"]:
                variant = VariantRef(row["variant_key"], int(row["position"]), "PASS")
                validate_key(variant, gene.chromosome)
                if variant.key in seen_variants:
                    raise ValueError(f"duplicate public variant in gene {gene.gene_id}")
                variants[-1].append(variant)
                seen_variants.add(variant.key)
    if not groups:
        raise ValueError("public variant list has no genes")
    return tuple(GeneVariants(group.gene, tuple(values)) for group, values in zip(groups, variants))


def prepare_site_blocks(*, pgen_prefix, gene_variants, rows_a, rows_b, out_dir,
                        extractor, party, public_groups, available_variants):
    groups = gene_variants if party == 1 else public_groups
    # The variant list is written last, after all genotype and phenotype files.
    if (out_dir / "public_variants.tsv").is_file():
        if (read_public_var_list(out_dir / "public_variants.tsv") != tuple(groups)
                or not (out_dir / ("A" if party == 1 else "B")).is_dir()):
            raise ValueError("existing site output differs; use a new run_dir or prepare --clear")
        return out_dir
    local = {group.gene.gene_id: group for group in gene_variants}
    available = {variant.key for variant in available_variants}
    public_keys = {variant.key for group in groups for variant in group.variants}
    plans = []
    for group in groups:
        if party == 1:
            roles = [(variant, "public_only") for variant in group.variants]
        else:
            roles = [
                (variant, "shared" if variant.key in available else "public_only")
                for variant in group.variants
            ]
            roles.extend(
                (variant, "private") for variant in local[group.gene.gene_id].variants
                if variant.key not in public_keys
            )
        plans.append(GenePlan(group.gene, tuple(roles)))
    extracted_plans, geno_a, columns_a, geno_b, columns_b = extract_genotypes(
        pgen_prefix, rows_a, rows_b, plans, extractor, filter_unavailable=False,
    )
    blocks = build_blocks(extracted_plans, geno_a, columns_a, geno_b, columns_b)
    write_outputs(out_dir, blocks, rows_a, rows_b, party=party)
    write_public_var_list(out_dir / "public_variants.tsv", groups)
    return out_dir
