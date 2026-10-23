from collections.abc import Callable, Collection, Mapping
from dataclasses import asdict, dataclass
from functools import partial
from pathlib import Path
from typing import Literal, TextIO
import fcntl
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
    prepare_blocks,
    select_gene_variants,
    select_rows,
)
from .plink import pgen_extract
from .output import write_selected_genes
from .site import prepare_site_blocks, read_public_var_list, write_public_var_list
from .selection import (
    select_file_genes,
    select_random_gene_groups,
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
    num_cov: int
    mask: Mapping[str, str | Collection[str]]
    max_maf: float | None

    samples_per_cohort: int
    shared_rate: float

    gene_selection: GeneSelectionRequest
    prepared_cache_dir: Path | None = None
    party: int = 0
    public_var_list: Path | None = None


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
    fields = ("chromosomes", "ancestries", "phenotype_columns", "num_cov", "max_maf",
              "samples_per_cohort", "shared_rate",
              "phenotype_id_column", "covariate_id_column", "covariate_column",
              "covariate_columns", "is_cov_single_column", "ancestry_id_column", "ancestry_column")
    result = {name: config.get(name, [] if name == "covariate_columns" else "") for name in fields}
    result["max_maf"] = config.get("max_maf")
    result["mask"] = {key: sorted({value} if isinstance(value, str) else value)
                      for key, value in config["mask"].items()}
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
        max_maf=max_maf,
    )


def prepared_cache_config(request: PrepareRequest, chromosome: int, ancestry: str) -> dict:
    """Identify prepared inputs without reading the genotype file contents."""
    def source(path: Path) -> dict:
        path = path.resolve()
        stat = path.stat()
        return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}

    config = asdict(request)
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
    selection = config["gene_selection"]
    selection["path"] = source(selection["path"]) if selection["mode"] == "file" else None
    if request.party:
        config["public_var_list"] = ({"sha256": hashlib.sha256(chromosome_path(str(request.public_var_list), chromosome).read_bytes()).hexdigest()}
                                     if request.public_var_list else None)
    else:
        del config["party"], config["public_var_list"]
    config["mask"] = {
        column: sorted({values} if isinstance(values, str) else set(values))
        for column, values in request.mask.items()
    }
    marker = request.phenotype.resolve().parent / "split.json"
    if marker.exists():
        config["split"] = json.loads(marker.read_text())["split"]
    config.update(format_version=1, chromosome=chromosome, ancestry_group=ancestry)
    return config


def prepare_cached_blocks(
    request: PrepareRequest, chromosome: int, ancestry: str, build: Callable[..., Path],
    selection: dict,
) -> Path:
    cache_root = request.prepared_cache_dir.resolve()
    output = request.run_dir / "prepared" / ancestry / f"chr{chromosome}"
    if output.exists() and not output.is_symlink():
        raise ValueError(f"{output} is not a link; use a new run_dir or prepare --clear")
    cache_config = prepared_cache_config(request, chromosome, ancestry)
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
        print(f"    ancestry: {ancestry}, chromosome: {chromosome}, MAF: {request.max_maf}, shared_rate: {request.shared_rate}", flush=True)
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


def prepare_chromosomes(
    request: PrepareRequest,
    extractor: GenotypeExtractor | None = None,
) -> tuple[GeneRef, ...]:
    validate_split_inputs(request.phenotype, asdict(request))
    if not request.chromosomes:
        raise ValueError("at least one chromosome is required")

    if request.party not in (0, 1, 2) or (request.party == 2) != (request.public_var_list is not None):
        raise ValueError("party must be 0/1/2; public-var-list is required only for Party 2")
    public_groups = ()
    if request.party == 2:
        paths = tuple(dict.fromkeys(chromosome_path(str(request.public_var_list), chromosome)
                                    for chromosome in request.chromosomes))
        public_groups = tuple(group for path in paths for group in read_public_var_list(path))
        if len({group.gene.gene_id for group in public_groups}) != len(public_groups):
            raise ValueError("duplicate gene across public variant lists")
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
        ancestry_groups=request.ancestries,
        ancestry_id_column=request.ancestry_id_column,
        ancestry_column=request.ancestry_column,
        num_cov=request.num_cov,
    )

    if request.samples_per_cohort == 0:
        samples_per_cohort: int | Literal["all"] = "all"
    else:
        samples_per_cohort = request.samples_per_cohort

    file_genes: tuple[GeneRef, ...] = ()
    if request.gene_selection.mode == "file" and request.party != 2:
        if request.gene_selection.path is None:
            raise ValueError(
                "gene selection path is required in file mode"
            )
        file_genes = read_gene_panel(
            request.gene_selection.path
        )

    selected_genes = []
    export_groups = []
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
                    party=request.party,
                )
        elif inputs.psam_ids != reference_psam_ids:
            raise ValueError(
                f"ordered PSAM IDs differ for chromosome "
                f"{reference_chromosome} and chromosome {chromosome}"
            )

        chromosome_public = tuple(group for group in public_groups if group.gene.chromosome == str(chromosome))
        if request.party == 2:
            chromosome_groups = select_gene_variants(
                gene_panel=tuple(group.gene for group in chromosome_public),
                variants=inputs.variants, annotations=inputs.annotations,
                annotation_columns=inputs.annotation_columns, chromosome=str(chromosome),
                gene_selection="all", mask=request.mask, max_maf=request.max_maf,
            )
        else:
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
        if request.party == 1:
            export_groups.extend(chromosome_groups)

        for ancestry in request.ancestries:
            rows_a, rows_b = rows_by_ancestry[ancestry]
            build = partial(
                prepare_blocks,
                pgen_prefix=inputs.pgen_prefix,
                gene_variants=chromosome_groups,
                rows_a=rows_a,
                rows_b=rows_b,
                shared_rate=request.shared_rate,
                extractor=extractor,
            )
            if request.party:
                build = partial(
                    prepare_site_blocks, pgen_prefix=inputs.pgen_prefix,
                    gene_variants=chromosome_groups, rows_a=rows_a, rows_b=rows_b,
                    extractor=extractor, party=request.party,
                    public_groups=chromosome_public, available_variants=inputs.variants,
                )
            if request.prepared_cache_dir is not None:
                prepare_cached_blocks(request, chromosome, ancestry, build, selection={
                    "samples_a": rows_a.sample_ids,
                    "samples_b": rows_b.sample_ids,
                    "genes": [(group.gene.gene_id, [variant.key for variant in group.variants])
                              for group in chromosome_groups],
                })
            else:
                build(out_dir=request.run_dir / "prepared" / ancestry / f"chr{chromosome}")

    if request.party == 1:
        path = request.run_dir / "public_variants.tsv"
        write_public_var_list(path, export_groups)
        print("Wrote public variant list to", path)

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
        num_cov=payload["num_cov"],
        mask=payload["mask"],
        max_maf=payload["max_maf"],
        samples_per_cohort=payload["samples_per_cohort"],
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
