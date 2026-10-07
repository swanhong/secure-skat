from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import shutil
import subprocess
import tomllib
from pathlib import Path
from tempfile import TemporaryDirectory
from python.preprocessing.input import load_inputs, load_sample_inputs, read_gene_panel
from python.preprocessing.model import GeneVariants
from python.preprocessing.output import write_selected_genes
from python.preprocessing.pipeline import assign_roles, select_rows
from python.preprocessing.prepare import (
    GeneSelectionRequest, select_gene_groups, split_input_config, validate_split_inputs,
)



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
DEFAULT_COVARIATE_ROOT = "gs://gwas-data-wgs-wb-jaunty-blueberry-8679/v9_intermediate_results/pca_aou"
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


def write_table(path: Path, header, rows, delimiter="\t") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter=delimiter, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def site_path(config: dict, field: str, chromosome=0, ancestry="") -> Path:
    return Path(config[field].replace("{chromosome}", str(chromosome))
                .replace("{anc}", ancestry.lower()))


def split_aou_inputs(args: argparse.Namespace, pgen_prefixes: dict[int, Path]) -> None:
    # step1: Read MVP/A and AoU/B configs and normalize ancestry names.
    configs = []
    for directory in (args.config_mvp, args.config):
        config = {"shared_rate": 0.6, "is_cov_single_column": True,
                  "sample_seed": 0, "role_seed": 0, "samples_per_cohort": 0}
        for name in ("configGlobal.toml", "configPrepare.toml"):
            with (directory / name).open("rb") as stream:
                config.update(tomllib.load(stream))
        config["ancestries"] = [value.strip().upper() for value in config["ancestries"]]
        config["mask"] = dict(tuple(part.strip() for part in item.split("=", 1))
                              for item in config["masks"])
        configs.append(config)
    mvp, aou = configs

    # step2: Check that both sites use matching sample selection and variant filters.
    for field in ("chromosomes", "ancestries", "phenotype_columns", "num_cov",
                  "samples_per_cohort", "sample_seed", "masks", "max_maf"):
        if mvp.get(field) != aou.get(field):
            raise ValueError(f"A/B splitting requires matching {field} in both configs")

    # step3: Separate source and site output paths and prevent overwriting existing files.
    roots = [site_path(config, "phenotype").resolve().parent for config in configs]
    source = args.output_dir.resolve()
    for index, (config, root) in enumerate(zip(configs, roots)):
        for other in (source, roots[1 - index]):
            if root.is_relative_to(other) or other.is_relative_to(root):
                raise ValueError("source, MVP and AoU input directories must be separate")
        for field in ("genotype", "gene_panel", "annotation", "covariate", "ancestry"):
            if not site_path(config, field).resolve().is_relative_to(root):
                raise ValueError(f"split {field} must be below {root}")

    source_files = [source / "phenotype.csv", source / "ancestry_preds.tsv"]
    source_files += [source / f"covariate/{ancestry.lower()}_pca.eigenvec" for ancestry in mvp["ancestries"]]
    for chromosome, prefix in pgen_prefixes.items():
        source_files += [Path(f"{prefix}.{extension}") for extension in ("pgen", "pvar", "psam")]
        source_files += [source / f"{folder}/chr{chromosome}.tsv" for folder in ("annotation", "gene_panel")]
    if mvp["gene_selection"].get("path"):
        source_files.append(Path(mvp["gene_selection"]["path"]))
    split = {"configs": [split_input_config(config) for config in configs],
             "sources": {str(path.resolve()): [path.stat().st_size, path.stat().st_mtime_ns]
                         for path in source_files}}
    manifests = [validate_split_inputs(site_path(config, "phenotype"), config) for config in configs]
    if all(manifests):
        if any(manifest["split"] != split for manifest in manifests):
            raise ValueError("split source or settings changed; regenerate both site inputs in empty directories")
        print(f"exists: split inputs in {roots[0]} and {roots[1]}")
        return
    if any(root.exists() and any(root.iterdir()) for root in roots):
        raise ValueError("incomplete or unmarked split; remove both site output directories or use empty paths")
    for chromosome in pgen_prefixes:
        with (source / f"annotation/chr{chromosome}.tsv").open() as stream:
            columns = next(csv.reader(stream, delimiter="\t"))
        required = set(mvp["mask"]) | ({"MAF"} if mvp.get("max_maf") is not None else set())
        missing = required - set(columns)
        if missing:
            raise ValueError(f"mask columns missing on chromosome {chromosome}: {sorted(missing)}; "
                             "provide these annotation columns or configure masks matching the source")

    # step4: Read phenotype, covariate, and ancestry tables from the AoU source.
    samples = load_sample_inputs(
        phenotype_path=source / "phenotype.csv",
        covariate_path=source / "covariate/{anc}_pca.eigenvec",
        ancestry_path=source / "ancestry_preds.tsv",
        phenotype_id_column=aou["phenotype_id_column"],
        phenotype_columns=tuple(aou["phenotype_columns"]),
        covariate_id_column=aou["covariate_id_column"],
        covariate_column=aou.get("covariate_column", ""),
        covariate_columns=tuple(aou.get("covariate_columns", ())),
        is_cov_single_column=aou["is_cov_single_column"],
        ancestry_groups=tuple(aou["ancestries"]),
        ancestry_id_column=aou["ancestry_id_column"],
        ancestry_column=aou["ancestry_column"], num_cov=aou["num_cov"],
    )

    # step5: Prepare gene selection and annotation masks from the MVP config.
    selection = mvp["gene_selection"]
    request = GeneSelectionRequest(selection["mode"], selection.get("per_chromosome", 0),
                                   selection.get("seed", 0),
                                   Path(selection["path"]) if selection.get("path") else None)
    file_genes = read_gene_panel(request.path) if request.mode == "file" else ()
    mask = mvp["mask"]
    reference_ids = None
    site_ids = []

    # step6: Read source sample IDs, variants, genes, and annotations for each chromosome.
    for chromosome, prefix in pgen_prefixes.items():
        inputs = load_inputs(prefix, source / f"gene_panel/chr{chromosome}.tsv",
                             source / f"annotation/chr{chromosome}.tsv")
        if reference_ids is None:
            # step7: Split samples into disjoint A/B groups by ancestry once, using the first chromosome.
            reference_ids = inputs.psam_ids
            if len(set(reference_ids)) != len(reference_ids):
                raise ValueError("duplicate source sample IDs")
            rows_by_site = [{}, {}]
            for ancestry in mvp["ancestries"]:
                rows = select_rows(reference_ids, samples.phenotypes, samples.covariates,
                                   samples.ancestries, ancestry, tuple(mvp["phenotype_columns"]),
                                   mvp["samples_per_cohort"] or "all", mvp["sample_seed"])
                for index in (0, 1):
                    rows_by_site[index][ancestry] = rows[index]

            # step8: Write selected sample tables using each site's paths and column names.
            for config, rows in zip(configs, rows_by_site):
                members = [(ancestry, group, i) for ancestry, group in rows.items()
                           for i in range(len(group.sample_ids))]
                ids = {group.sample_ids[i] for _, group, i in members}
                site_ids.append(tuple(sample for sample in reference_ids if sample in ids))
                write_table(site_path(config, "phenotype"),
                            (config["phenotype_id_column"], *config["phenotype_columns"]),
                            ((group.sample_ids[i], *group.phenotypes[i]) for _, group, i in members), ",")
                write_table(site_path(config, "ancestry"),
                            (config["ancestry_id_column"], config["ancestry_column"]),
                            ((group.sample_ids[i], ancestry) for ancestry, group, i in members))
                covariate_rows = {}
                for ancestry, group, i in members:
                    path = site_path(config, "covariate", ancestry=ancestry)
                    values = ((json.dumps(group.covariates[i]),) if config["is_cov_single_column"]
                              else group.covariates[i])
                    covariate_rows.setdefault(path, []).append((group.sample_ids[i], *values))
                columns = ((config["covariate_column"],) if config["is_cov_single_column"]
                           else tuple(config["covariate_columns"]))
                for path, rows in covariate_rows.items():
                    write_table(path, (config["covariate_id_column"], *columns), rows)
        elif inputs.psam_ids != reference_ids:
            raise ValueError(f"ordered source PSAM IDs differ on chromosome {chromosome}")

        # step9: Apply gene selection, annotation masks, and the MAF threshold for this chromosome.
        groups = select_gene_groups(request, inputs, chromosome, mask, mvp.get("max_maf"), file_genes)

        # step10: Assign each variant one site role, even when it belongs to multiple genes.
        seen, unique_groups = set(), []
        for group in groups:
            unique = tuple(variant for variant in group.variants if variant.key not in seen)
            seen.update(variant.key for variant in unique)
            unique_groups.append(GeneVariants(group.gene, unique))
        roles = {variant.key: role for plan in assign_roles(unique_groups, mvp["role_seed"], mvp["shared_rate"])
                 for variant, role in plan.variant_roles}
        if not groups or (request.mode == "random" and any(
                all(roles[v.key] == "private" for v in group.variants) for group in groups)):
            raise ValueError(f"selected genes need MVP public variants on chromosome {chromosome}")
        pairs = {(group.gene.gene_id, variant.key) for group in groups for variant in group.variants}

        # step11: Select shared/public_only variants for A and shared/private variants for B.
        for index, config in enumerate(configs):
            allowed = {"shared", "public_only"} if index == 0 else {"shared", "private"}
            keys = {variant.key for variant in inputs.variants if roles.get(variant.key) in allowed}
            if not keys or not site_ids[index]:
                raise ValueError(f"no samples or variants for site {index + 1}, chromosome {chromosome}")
            output = site_path(config, "genotype", chromosome)
            output.parent.mkdir(parents=True, exist_ok=True)

            # step12: Pass sample and variant lists to PLINK and write site inputs in PGEN format.
            with TemporaryDirectory() as temporary:
                keep, extract = Path(temporary) / "keep.tsv", Path(temporary) / "variants.txt"
                with Path(f"{prefix}.psam").open() as stream:
                    reader = csv.DictReader(stream, delimiter="\t")
                    reader.fieldnames = [name.lstrip("#") for name in reader.fieldnames]
                    fields = [name for name in ("FID", "IID", "SID") if name in reader.fieldnames]
                    members = set(site_ids[index])
                    write_table(keep, ("#" + fields[0], *fields[1:]),
                                (tuple(row[name] for name in fields) for row in reader if row["IID"] in members))
                extract.write_text("".join(f"{key}\n" for key in keys))
                with Path(f"{prefix}.pvar").open() as stream:
                    chromosome_name = next(line for line in stream if not line.startswith("#")).split()[0]
                run([os.path.expandvars(args.plink2_bin), "--pfile", str(prefix),
                     "--keep", str(keep), "--extract", str(extract),
                     "--output-chr", "chrMT" if chromosome_name.startswith("chr") else "MT",
                     "--make-pgen", "--out", str(output)])
                Path(f"{output}.log").unlink(missing_ok=True)

            # step13: Write the selected gene panel and annotations retained at this site.
            panel = site_path(config, "gene_panel", chromosome)
            panel.parent.mkdir(parents=True, exist_ok=True)
            write_selected_genes(panel, tuple(group.gene for group in groups))
            write_table(site_path(config, "annotation", chromosome),
                        ("variant_key", "gene_id", "gene_symbol", *inputs.annotation_columns),
                        ((row.variant_key, row.gene_id, row.gene_symbol,
                          *(row.values[column] for column in inputs.annotation_columns))
                         for row in inputs.annotations if row.variant_key in keys
                         and (row.gene_id, row.variant_key) in pairs))

    # step14: Record completion after writing both sites and all chromosomes.
    for root, config in zip(roots, configs):
        outputs = {str(path.resolve()): [path.stat().st_size, path.stat().st_mtime_ns]
                   for path in root.rglob("*") if path.is_file()}
        (root / "split.json").write_text(json.dumps(
            {"config": split_input_config(config), "split": split, "outputs": outputs},
            sort_keys=True, indent=2) + "\n")


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
    copy_gcs_files(
        tuple(f"{DEFAULT_COVARIATE_ROOT}/{ancestry}_pca.eigenvec" for ancestry in ("ALL", "afr", "amr", "eur")),
        output_dir / "covariate",
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
    if args.config_mvp is not None:
        split_aou_inputs(args, pgen_prefixes)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--config-mvp", type=Path, help="split normalized inputs into the configured MVP/A and AoU/B paths")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[3] / "data" / "aou" / "generated",
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

    with (args.config / "configGlobal.toml").open("rb") as config_file:
        config = tomllib.load(config_file)
    args.chromosome = config["chromosomes"]

    if len(set(args.chromosome)) != len(args.chromosome):
        parser.error("chromosomes must not contain duplicates")
    if any(chromosome < 1 or chromosome > 22 for chromosome in args.chromosome):
        parser.error("chromosomes must be between 1 and 22")

    args.output_dir = args.output_dir.expanduser()
    args.plink2_bin = str(Path(args.plink2_bin).expanduser())
    prepare_aou(args)


if __name__ == "__main__":
    main()
