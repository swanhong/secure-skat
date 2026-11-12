from collections.abc import Callable, Collection, Mapping
from dataclasses import asdict, dataclass, replace
from functools import partial
from pathlib import Path
from typing import Literal, TextIO
import math
import hashlib
import json
from tempfile import TemporaryDirectory
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
    eligible_rows,
    prepare_blocks,
    select_gene_variants,
    select_rows,
)
from .plink import pgen_allele_counts, pgen_extract
from .output import write_selected_genes
from .site import (
    prepare_site_blocks, read_local_cut_list, read_public_var_list,
    write_local_cut_list, write_public_var_list,
)
from .selection import (
    select_file_genes,
    select_random_gene_groups,
    split_sample_ids,
)

@dataclass(frozen=True)
class GeneSelectionRequest:
    mode: Literal["random", "file", "all"]
    per_chromosome: int
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
    covariate_columns: tuple[str, ...]
    is_cov_single_column: bool
    ancestry_id_column: str
    ancestry_column: str
    phenotype_columns: tuple[str, ...]
    ancestries: tuple[str, ...]
    local_max_maf: float | None

    samples_per_cohort: int
    shared_rate: float

    gene_selection: GeneSelectionRequest
    debug_split_num_parties: int = 0
    prepared_cache_dir: Path | None = None
    party: int = 0
    public_var_list: Path | None = None
    mode: Literal["both", "joint", "local"] = "joint"
    global_max_maf: tuple[float, ...] = ()


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


def split_input_config(config: Mapping) -> dict:
    """Normalize the input settings shared by TOML configs and PrepareRequest."""
    fields = ("chromosomes", "ancestries", "phenotype_columns",
              "samples_per_cohort", "shared_rate",
              "phenotype_id_column", "covariate_id_column", "covariate_column",
              "covariate_columns", "is_cov_single_column", "ancestry_id_column", "ancestry_column")
    result = {name: config.get(name, [] if name == "covariate_columns" else "") for name in fields}
    selection = config["gene_selection"]
    result["gene_selection"] = {"mode": selection["mode"], "per_chromosome": selection.get("per_chromosome", 0),
                                "path": str(selection.get("path") or "")}
    for name in ("genotype", "gene_panel", "annotation", "phenotype", "covariate", "ancestry"):
        result[name] = str(Path(config[name]).resolve())
    return json.loads(json.dumps(result))


def validate_split_inputs(phenotype: Path, config: dict) -> dict | None:
    marker = phenotype.resolve().parent / "split.json"
    if not marker.exists():
        return None
    manifest = json.loads(marker.read_text())
    if manifest["config"] != split_input_config(config):
        raise ValueError(f"split settings changed: {marker}; regenerate both site inputs in empty directories")
    files = dict(manifest["outputs"])
    selection = config["gene_selection"]
    if selection["mode"] == "file":
        path = str(Path(selection["path"]).resolve())
        if path in manifest["split"]["sources"]:
            files[path] = manifest["split"]["sources"][path]
    for name, signature in files.items():
        path = Path(name)
        if not path.is_file() or [path.stat().st_size, path.stat().st_mtime_ns] != signature:
            raise ValueError(f"split input missing or changed: {path}; regenerate both site inputs")
    return manifest


def debug_split_role(key: str, shared_rate: float) -> str:
    """Role from the key alone, so both debug-split site processes agree without communicating."""
    u = int.from_bytes(hashlib.sha256(key.encode()).digest()[:8], "big") / 2**64
    if u < shared_rate:
        return "shared"
    return "public_only" if u < shared_rate + (1 - shared_rate) / 2 else "private"


def local_cut_list_path(public_var_list: Path, chromosome: int) -> Path:
    """Party 1 writes the local-cut list next to its public variant list."""
    return chromosome_path(str(public_var_list), chromosome).with_name("local_cut_variants.tsv")


def local_maf_passing(counts: Mapping[str, int], cohort_size: int, local_max_maf: float | None) -> set[str] | None:
    """Keys whose cohort MAF is at most local_max_maf; None disables the cut."""
    if local_max_maf is None:
        return None
    return {key for key, count in counts.items()
            if min(count, 2 * cohort_size - count) / (2 * cohort_size) <= local_max_maf}


def select_gene_groups(
    request: GeneSelectionRequest,
    inputs: PrepInputs,
    chromosome: int,
    allowed_keys: Collection[str] | None,
    file_genes: Collection[GeneRef],
) -> tuple[GeneVariants, ...]:
    gene_selection = "all"
    if request.mode == "file":
        genes = select_file_genes(inputs, str(chromosome), file_genes)
        gene_selection = tuple(gene.gene_id for gene in genes)
    elif request.mode not in {"all", "random"}:
        raise ValueError(f"unsupported gene selection mode: {request.mode}")
    groups = select_gene_variants(
        gene_panel=inputs.gene_panel,
        variants=inputs.variants,
        annotations=inputs.annotations,
        chromosome=str(chromosome),
        gene_selection=gene_selection,
        allowed_keys=allowed_keys,
    )
    if request.mode == "random":
        return select_random_gene_groups(groups, request.per_chromosome)
    return groups


def prepared_cache_config(request: PrepareRequest, chromosome: int, ancestry: str, scope: str = "joint") -> dict:
    """Identify prepared inputs without reading the genotype file contents."""
    def source(path: Path) -> dict:
        path = path.resolve()
        stat = path.stat()
        return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}

    config = asdict(request)
    if config["debug_split_num_parties"] == 0:
        del config["debug_split_num_parties"]
    for name in ("run_dir", "prepared_cache_dir", "chromosomes", "ancestries", "plink2_bin"):
        del config[name]
    prefix = chromosome_path(request.genotype, chromosome)
    config["genotype"] = {
        extension: source(Path(f"{prefix}.{extension}"))
        for extension in ("pgen", "pvar", "psam")
    }
    for name in ("annotation", "gene_panel"):
        config[name] = source(chromosome_path(getattr(request, name), chromosome))
    for name in ("phenotype", "ancestry"):
        config[name] = source(getattr(request, name))
    config["covariate"] = source(Path(str(request.covariate).replace("{anc}", ancestry.lower())))
    config["cohort_covariate"] = source(Path(str(request.covariate).replace("{anc}", "all")))
    selection = config["gene_selection"]
    selection["path"] = source(selection["path"]) if selection["mode"] == "file" else None
    if request.party and scope == "joint":
        config["public_var_list"] = ({"sha256": hashlib.sha256(chromosome_path(str(request.public_var_list), chromosome).read_bytes()).hexdigest()}
                                     if request.public_var_list else None)
        config["local_cut_list"] = ({"sha256": hashlib.sha256(local_cut_list_path(request.public_var_list, chromosome).read_bytes()).hexdigest()}
                                    if request.public_var_list else None)
    else:
        del config["party"], config["public_var_list"]
    marker = request.phenotype.resolve().parent / "split.json"
    if marker.exists():
        config["split"] = json.loads(marker.read_text())["split"]
    config.pop("mode")
    if scope == "joint":
        config.pop("global_max_maf")
    else:
        config["local_max_maf"] = max(request.global_max_maf)
        config["source_party"] = request.party
    config.update(format_version=6, chromosome=chromosome, ancestry_group=ancestry, scope=scope)
    return config


def prepare_cached_blocks(
    request: PrepareRequest, chromosome: int, ancestry: str, build: Callable[..., Path],
    selection: dict, scope: str = "joint",
) -> Path:
    cache_root = request.prepared_cache_dir.resolve()
    output_root = request.run_dir / "prepared"
    if scope == "local":
        cache_root = cache_root / "local"
        output_root = output_root / "local"
    output = output_root / ancestry / f"chr{chromosome}"
    if output.exists() and not output.is_symlink():
        raise ValueError(f"{output} is not a link; use a new run_dir or prepare --clear")
    cache_config = prepared_cache_config(request, chromosome, ancestry, scope)
    cache_config["selection_sha256"] = hashlib.sha256(
        json.dumps(selection, sort_keys=True).encode()
    ).hexdigest()
    config = json.dumps(cache_config, sort_keys=True, indent=2) + "\n"
    parent = cache_root / ancestry
    parent.mkdir(parents=True, exist_ok=True)
    cached = parent / f"chr{chromosome}"
    if cached.exists():
        config_path = cached / "config.json"
        if not config_path.is_file() or config_path.read_text() != config:
            raise ValueError(f"prepared cache config missing or changed: {cached}; use a new run_dir")
        print(f"Using existing prepared cache: {cached}", flush=True)
        print(f"    ancestry: {ancestry}, chromosome: {chromosome}, local_max_maf: {cache_config['local_max_maf']}, shared_rate: {request.shared_rate}", flush=True)
    else:
        print("Creating prepared cache:", cached, flush=True)
        with TemporaryDirectory(prefix=f".chr{chromosome}.", dir=parent) as temporary:
            directory = Path(temporary)
            build(out_dir=directory)
            (directory / "config.json").write_text(config)
            directory.rename(cached)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.is_symlink():
        output.unlink()
    output.symlink_to(cached, target_is_directory=True)
    return output


def prepare_local_site(*, out_dir: Path, max_maf: float, **kwargs) -> Path:
    result = prepare_site_blocks(out_dir=out_dir, joint_stat=False, **kwargs)
    (out_dir / "max_maf.txt").write_text(f"{max_maf!r}\n")
    return result


def prepare_chromosomes(
    request: PrepareRequest,
    extractor: GenotypeExtractor | None = None,
) -> tuple[GeneRef, ...]:
    if request.mode not in {"both", "joint", "local"}:
        raise ValueError("mode must be both, joint, or local")
    joint = request.mode != "local"
    local = request.mode != "joint"
    if local and request.party not in (1, 2):
        raise ValueError("local preparation requires party 1 or 2")
    if local and (not request.global_max_maf or any(
            not math.isfinite(t) or not 0 < t <= 0.5 for t in request.global_max_maf)):
        raise ValueError("global_max_maf must contain thresholds in (0, 0.5]")
    if request.debug_split_num_parties not in (0, 2):
        raise ValueError("debug_split_num_parties must be 0 or 2")
    if request.debug_split_num_parties and request.party not in (1, 2):
        raise ValueError("debug_split_num_parties requires --party 1 or 2")
    split_manifest = validate_split_inputs(request.phenotype, asdict(request))
    if request.debug_split_num_parties and split_manifest is not None:
        raise ValueError("debug_split_num_parties requires common source inputs, not already split inputs")
    gene_selection = (GeneSelectionRequest("all", 0, None) if split_manifest is not None
                      else request.gene_selection)
    if not request.chromosomes:
        raise ValueError("at least one chromosome is required")

    if request.party not in (0, 1, 2) or (joint and (request.party == 2) != (request.public_var_list is not None)):
        raise ValueError("party must be 0/1/2; public-var-list is required only for Party 2")
    public_groups = ()
    if joint and request.party == 2:
        paths = tuple(dict.fromkeys(chromosome_path(str(request.public_var_list), chromosome)
                                    for chromosome in request.chromosomes))
        public_groups = tuple(group for path in paths for group in read_public_var_list(path))
        if len({group.gene.gene_id for group in public_groups}) != len(public_groups):
            raise ValueError("duplicate gene across public variant lists")
        # Party 1's locally cut variants are common there, so they must not enter as private variants here.
        public_cut_keys = set().union(*(
            read_local_cut_list(path) for path in dict.fromkeys(
                local_cut_list_path(request.public_var_list, chromosome) for chromosome in request.chromosomes)
        ))
    if public_groups and {int(group.gene.chromosome) for group in public_groups} != set(request.chromosomes):
        raise ValueError("public variant list chromosomes must match configured chromosomes")
    selected_genes_path = request.run_dir / "selected_genes.tsv"

    if extractor is None:
        extractor = pgen_extract
    sample_inputs = load_sample_inputs(
        phenotype_path=request.phenotype,
        covariate_path=request.covariate,
        ancestry_path=request.ancestry,
        phenotype_id_column=request.phenotype_id_column,
        phenotype_columns=request.phenotype_columns,
        covariate_id_column=request.covariate_id_column,
        covariate_column=request.covariate_column,
        covariate_columns=request.covariate_columns,
        is_cov_single_column=request.is_cov_single_column,
        # The site cohort for the local MAF cut and allele counts is the ALL analysis sample.
        ancestry_groups=tuple(dict.fromkeys((*request.ancestries, "ALL"))),
        ancestry_id_column=request.ancestry_id_column,
        ancestry_column=request.ancestry_column,
    )

    samples_per_cohort = request.samples_per_cohort or "all"

    file_genes: tuple[GeneRef, ...] = ()
    if gene_selection.mode == "file" and (local or request.party != 2):
        if gene_selection.path is None:
            raise ValueError(
                "gene selection path is required in file mode"
            )
        file_genes = read_gene_panel(
            gene_selection.path
        )

    selected_genes = []
    export_groups = []
    local_genes = []
    local_cut_keys = set()
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
            sample_ids = reference_psam_ids
            if request.debug_split_num_parties:
                sample_ids = split_sample_ids(
                    reference_psam_ids, sample_inputs.ancestries,
                    request.debug_split_num_parties, request.party,
                )
                print(f"Debug sample split: party {request.party}/{request.debug_split_num_parties}, "
                      f"{len(sample_ids)}/{len(reference_psam_ids)} source samples", flush=True)
            cohort_ids = tuple(row[0] for row in eligible_rows(
                sample_ids, sample_inputs.phenotypes, sample_inputs.covariates["ALL"],
                sample_inputs.ancestries, "ALL", request.phenotype_columns,
            ))
            if not cohort_ids:
                raise ValueError("no ALL analysis samples for the local MAF cohort")
            print(f"Local MAF cohort (ALL analysis samples): {len(cohort_ids)}", flush=True)

            for ancestry in request.ancestries:
                rows_by_ancestry[ancestry] = select_rows(
                    psam_ids=sample_ids,
                    phenotypes=sample_inputs.phenotypes,
                    covariates=sample_inputs.covariates[ancestry],
                    ancestries=sample_inputs.ancestries,
                    ancestry=ancestry,
                    phenotype_columns=request.phenotype_columns,
                    samples_per_cohort=samples_per_cohort,
                    party=request.party,
                )
        elif inputs.psam_ids != reference_psam_ids:
            raise ValueError(
                f"ordered PSAM IDs differ for chromosome "
                f"{reference_chromosome} and chromosome {chromosome}"
            )
        # A non-PASS row counts as absent everywhere, including Party 2's public-list
        # variants, the same as a row removed by QC.
        inputs = replace(inputs, variants=tuple(
            variant for variant in inputs.variants if variant.filter_value in {"PASS", "."}
        ))
        if request.debug_split_num_parties:
            # Both sites read one source; hide B-only variants from A and A-only variants from B.
            hidden = "private" if request.party == 1 else "public_only"
            inputs = replace(inputs, variants=tuple(
                variant for variant in inputs.variants
                if debug_split_role(variant.key, request.shared_rate) != hidden
            ))

        chromosome_public = tuple(group for group in public_groups if group.gene.chromosome == str(chromosome))
        candidate_panel = (tuple(group.gene for group in chromosome_public)
                           if joint and request.party == 2 else inputs.gene_panel)
        select_candidates = partial(
            select_gene_variants, gene_panel=candidate_panel, variants=inputs.variants,
            annotations=inputs.annotations,
            chromosome=str(chromosome), gene_selection="all",
        )
        candidates = select_candidates()
        # Party 2 also counts every public variant, including ones it would cut locally.
        local_candidates = (select_gene_variants(
            gene_panel=inputs.gene_panel, variants=inputs.variants, annotations=inputs.annotations,
            chromosome=str(chromosome), gene_selection="all",
        ) if local and joint and request.party == 2 else candidates)
        count_keys = {variant.key for group in (*candidates, *chromosome_public, *local_candidates)
                      for variant in group.variants}
        count_keys &= {variant.key for variant in inputs.variants}
        cohort_counts = pgen_allele_counts(inputs.pgen_prefix, cohort_ids, tuple(count_keys))
        if local:
            local_groups = select_gene_groups(
                gene_selection, inputs, chromosome,
                local_maf_passing(cohort_counts, len(cohort_ids), max(request.global_max_maf)),
                file_genes,
            )
            local_genes.extend(group.gene for group in local_groups)
            for ancestry, (rows_a, rows_b) in rows_by_ancestry.items():
                # Store a standalone site in the existing A layout, even for AoU/B.
                site_rows, empty_rows = (rows_a, rows_b) if request.party == 1 else (rows_b, rows_a)
                build = partial(
                    prepare_local_site, max_maf=max(request.global_max_maf), pgen_prefix=inputs.pgen_prefix,
                    gene_variants=local_groups, rows_a=site_rows, rows_b=empty_rows,
                    extractor=extractor, party=1, public_groups=(),
                    available_variants=inputs.variants,
                    cohort_counts=cohort_counts, cohort_size=len(cohort_ids),
                )
                if request.prepared_cache_dir is not None:
                    prepare_cached_blocks(request, chromosome, ancestry, build, selection={
                        "samples": site_rows.sample_ids,
                        "genes": [(g.gene.gene_id, [v.key for v in g.variants]) for g in local_groups],
                    }, scope="local")
                else:
                    build(out_dir=request.run_dir / "prepared" / "local" / ancestry / f"chr{chromosome}")
        if not joint:
            continue
        allowed_keys = local_maf_passing(cohort_counts, len(cohort_ids), request.local_max_maf)
        if request.party == 2:
            if allowed_keys is None:
                allowed_keys = {variant.key for group in candidates for variant in group.variants}
            allowed_keys -= public_cut_keys
            chromosome_groups = select_candidates(allowed_keys=allowed_keys)
        else:
            chromosome_groups = select_gene_groups(
                request=gene_selection,
                inputs=inputs,
                chromosome=chromosome,
                allowed_keys=allowed_keys,
                file_genes=file_genes,
            )
        if request.party == 1 and allowed_keys is not None:
            selected_ids = {group.gene.gene_id for group in chromosome_groups}
            local_cut_keys.update(
                variant.key for group in candidates if group.gene.gene_id in selected_ids
                for variant in group.variants
                if variant.key in cohort_counts and variant.key not in allowed_keys
            )
        selected_genes.extend(group.gene for group in chromosome_groups)
        if request.party == 1:
            export_groups.extend(chromosome_groups)

        for ancestry in request.ancestries:
            rows_a, rows_b = rows_by_ancestry[ancestry]
            common = dict(pgen_prefix=inputs.pgen_prefix, gene_variants=chromosome_groups,
                          rows_a=rows_a, rows_b=rows_b, extractor=extractor)
            if request.party:
                build = partial(
                    prepare_site_blocks, **common, party=request.party,
                    public_groups=chromosome_public, available_variants=inputs.variants,
                    cohort_counts=cohort_counts, cohort_size=len(cohort_ids),
                )
            else:
                build = partial(prepare_blocks, **common, shared_rate=request.shared_rate)
            if request.prepared_cache_dir is not None:
                prepare_cached_blocks(request, chromosome, ancestry, build, selection={
                    "samples_a": rows_a.sample_ids,
                    "samples_b": rows_b.sample_ids,
                    "genes": [(group.gene.gene_id, [variant.key for variant in group.variants])
                              for group in chromosome_groups],
                })
            else:
                build(out_dir=request.run_dir / "prepared" / ancestry / f"chr{chromosome}")

    if joint and request.party == 1:
        path = request.run_dir / "public_variants.tsv"
        write_public_var_list(path, export_groups)
        print("Wrote public variant list to", path)
        path = local_cut_list_path(path, request.chromosomes[0])
        write_local_cut_list(path, local_cut_keys)
        print(f"Wrote {len(local_cut_keys)} locally cut variants to", path)

    selected_genes = tuple(selected_genes)
    if joint and gene_selection.mode != "all":
        write_selected_genes(
            path=selected_genes_path,
            genes=selected_genes,
        )
        print("Wrote selected genes to", selected_genes_path)
    return selected_genes if joint else tuple(local_genes)

def read_prepare_request(
    stream: TextIO,
) -> PrepareRequest:
    payload = json.load(stream)
    gene_selection = payload["gene_selection"]
    gene_selection_path = gene_selection["path"]

    return PrepareRequest(
        run_dir=Path(payload["run_dir"]),
        mode=payload.get("mode", "joint"),
        global_max_maf=tuple(payload.get("global_max_maf", ())),
        party=payload.get("party", 0),
        public_var_list=Path(payload["public_var_list"]) if payload.get("public_var_list") else None,
        prepared_cache_dir=Path(payload["prepared_cache_dir"]) if payload.get("prepared_cache_dir") else None,
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
        covariate_columns=tuple(payload["covariate_columns"] or ()),
        is_cov_single_column=payload["is_cov_single_column"],
        ancestry_id_column=payload["ancestry_id_column"],
        ancestry_column=payload["ancestry_column"],
        phenotype_columns=tuple(payload["phenotype_columns"]),
        ancestries=tuple(
            ancestry.strip().upper()
            for ancestry in payload["ancestries"]
        ),
        local_max_maf=payload["local_max_maf"],
        samples_per_cohort=payload["samples_per_cohort"],
        debug_split_num_parties=payload.get("debug_split_num_parties", 0),
        shared_rate=payload["shared_rate"],
        gene_selection=GeneSelectionRequest(
            mode=gene_selection["mode"],
            per_chromosome=gene_selection["per_chromosome"],
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
