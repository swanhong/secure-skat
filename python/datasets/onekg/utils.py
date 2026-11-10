import csv
import gzip
import json
import math
import random
import subprocess
from pathlib import Path


def download_if_missing(url: str, destination: Path) -> Path:
    if destination.exists():
        print(f"exists: {destination}")
        return destination

    destination.parent.mkdir(parents=True, exist_ok=True)
    # download to partial, then rename to destination if successful
    partial = destination.with_name(f"{destination.name}.part")
    subprocess.run(["curl", "-fL", "--retry", "5", "--retry-all-errors", "-C", "-", "-o", str(partial), url,
        ],
        check=True,
    )
    partial.replace(destination)
    print(f"downloaded: {destination}")
    return destination

def create_pgen(
        vcf_path: Path,
        panel_path: Path,
        out_prefix: Path,
        keep_path: Path,
) -> Path:
    output_files = [
        Path(f"{out_prefix}.pgen"),
        Path(f"{out_prefix}.pvar"),
        Path(f"{out_prefix}.psam"),
    ]

    if all(path.exists() for path in output_files):
        print(f"exists: {out_prefix}")
        return out_prefix

    # read sample IDs from panel file
    with panel_path.open() as panel:
        columns = next(panel).split()
        sample_column = columns.index("sample")

        sample_ids = []
        for line in panel:
            sample_id = line.split()[sample_column]
            sample_ids.append(sample_id)

    # create output directories and write keep file
    out_prefix.parent.mkdir(parents=True, exist_ok=True)    
    keep_path.parent.mkdir(parents=True, exist_ok=True)
    keep_path.write_text(
        "".join(f"{sample_id}\n" for sample_id in sample_ids)
    )

    subprocess.run(
        [
            "plink2", "--vcf", str(vcf_path), 
            "--keep", str(keep_path), "--make-pgen", 
            "--out", str(out_prefix),
            "--var-filter",
            "--import-max-alleles", "2",
            "--set-all-var-ids", "@:#:$r:$a",
            "--new-id-max-allele-len", "1000",
        ],
        check=True,
    )

    print("created: " + ", ".join(str(path) for path in output_files))
    print("number of samples: " + str(len(sample_ids)))
    return out_prefix

def create_pca(
    pgen_prefixes: list[Path],
    work_dir: Path,
    num_components: int = 16,
) -> Path:
    work_dir.mkdir(parents=True, exist_ok=True)

    if len(pgen_prefixes) == 1:
        pca_input = pgen_prefixes[0]
    else:
        merge_list = work_dir / "pmerge_list.txt"
        merge_list.write_text(
            "".join(f"{prefix.resolve()}\n" for prefix in pgen_prefixes)
        )
        pca_input = work_dir / "merged"
        merged_files = [
            Path(f"{pca_input}.pgen"),
            Path(f"{pca_input}.pvar"),
            Path(f"{pca_input}.psam"),
        ]
        complete = all(path.is_file() for path in merged_files)
        if complete:
            with merged_files[0].open("rb") as file:
                header = file.read(12)
            # The standard PGEN header's variant/sample counts are finalized last.
            complete = (len(header) == 12 and header[:3] == b"\x6c\x1b\x10"
                        and int.from_bytes(header[3:7], "little") > 0
                        and int.from_bytes(header[7:11], "little") > 0)
        if not complete:
            print(f"Rebuilding incomplete PCA merge: {pca_input}", flush=True)
            partial = work_dir / "merged.partial"
            for path in merged_files + [Path(f"{partial}{path.suffix}") for path in merged_files] + [
                work_dir / name for name in ("pruned.prune.in", "pruned.prune.out", "pca.eigenvec", "pca.eigenval")
            ]:
                path.unlink(missing_ok=True)
            subprocess.run(
                [
                    "plink2",
                    "--pmerge-list", str(merge_list),
                    "--merge-info-mode", "erase",
                    "--out", str(partial),
                ],
                check=True,
            )
            for path in merged_files:
                Path(f"{partial}{path.suffix}").replace(path)
            Path(f"{partial}.log").replace(work_dir / "merged.log")

    prune_prefix = work_dir / "pruned"
    prune_in = Path(f"{prune_prefix}.prune.in")
    if not prune_in.exists():
        subprocess.run(
            [
                "plink2", "--pfile", str(pca_input),
                "--maf", "0.05",
                "--geno", "0.02",
                "--rm-dup", "force-first",
                "--indep-pairwise", "200", "50", "0.2",
                "--out", str(prune_prefix),
            ],
            check=True,
        )

    pca_prefix = work_dir / "pca"
    eigenvec = Path(f"{pca_prefix}.eigenvec")
    eigenval = Path(f"{pca_prefix}.eigenval")
    if eigenvec.exists() and eigenval.exists():
        print(f"exists: {eigenvec}, {eigenval}")
        return eigenvec

    subprocess.run(
        [
            "plink2", "--pfile", str(pca_input),
            "--extract", str(prune_in),
            "--rm-dup", "force-first",
            "--pca", str(num_components),
            "--out", str(pca_prefix),
        ],
        check=True,
    )
    print(f"created: {eigenvec}, {eigenval}")
    return eigenvec


def create_ancestry_table(
    panel_path: Path,
    eigenvec_path: Path,
    sample_ids: tuple[str, ...],
    out_path: Path,
    num_components: int = 16,
) -> Path:
    if out_path.exists():
        print(f"exists: {out_path}")
        return out_path
    ancestries = {}
    with panel_path.open() as panel:
        columns = next(panel).split()
        sample_column = columns.index("sample")
        super_pop_column = columns.index("super_pop")
        for line in panel:
            fields = line.split()
            ancestries[fields[sample_column]] = fields[super_pop_column]

    pca_values = {}
    with eigenvec_path.open() as eigenvec:
        columns = next(eigenvec).lstrip("#").split()
        iid_column = columns.index("IID")
        pc_columns = [
            columns.index(f"PC{index}")
            for index in range(1, num_components + 1)
        ]
        for line in eigenvec:
            fields = line.split()
            values = tuple(float(fields[index]) for index in pc_columns)
            if not all(math.isfinite(value) for value in values):
                raise ValueError(f"non-finite PCA value for {fields[iid_column]}")
            pca_values[fields[iid_column]] = values

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "research_id",
                "ancestry_pred",
                "probabilities",
                "pca_features",
                "ancestry_pred_other",
            ]
        )
        for sample_id in sample_ids:
            ancestry = ancestries[sample_id]
            features = json.dumps(
                pca_values[sample_id],
                separators=(",", ":"),
            )
            writer.writerow([sample_id, ancestry, "", features, ancestry])

    print(f"created: {out_path}, ({len(sample_ids)} samples)")
    return out_path


def create_inputs_from_gencode(
        gtf_path: Path,
        pvar_path: Path,
        gene_panel_path: Path,
        annotation_path: Path,
        chromosome: str = "22",
        seed: int = 42,
) -> tuple[Path, Path]:
    if gene_panel_path.exists() and annotation_path.exists():
        with annotation_path.open(newline="") as stream:
            columns = next(csv.reader(stream, delimiter="\t"))
        if columns != ["variant_key", "gene_id", "gene_symbol"]:
            raise ValueError(f"{annotation_path}: old synthetic annotation schema; regenerate gene panel and annotation in new paths")
        print(f"create_inputs_from_gencode - exists: {gene_panel_path}, {annotation_path}")
        return gene_panel_path, annotation_path
    chromosome = chromosome.removeprefix("chr")
    gtf_chromosome = f"chr{chromosome}"

    rng = random.Random(seed)
    genes = []
    with gzip.open(gtf_path, "rt") as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")
            if fields[0] != gtf_chromosome or fields[2] != "gene":
                continue

            attributes = {}
            # fields[8] contains attributes like 'gene_id "ENSG00000186092.5"; gene_name "OR4F5"; ...'
            for item in fields[8].rstrip(";").split("; "):
                key, value = item.strip().split(" ", 1)
                attributes[key] = value.strip('"')

            # filter for protein-coding genes, to make it simpler
            if attributes["gene_type"] != "protein_coding":
                continue

            genes.append(
                (
                    int(fields[3]),  # start position
                    int(fields[4]),  # end position
                    attributes["gene_id"],
                    attributes["gene_name"]
                )
            )
    gene_panel_path.parent.mkdir(parents=True, exist_ok=True)
    with gene_panel_path.open("w") as gene_panel:
        gene_panel.write(
            "gene_id\tgene_symbol\tchromosome\torder_index\n"
        )
        for order_index, (_, _, gene_id, gene_symbol) in enumerate(genes):
            gene_panel.write(
                f"{gene_id}\t{gene_symbol}\t{chromosome}\t{order_index}\n"
            )

    annotation_path.parent.mkdir(parents=True, exist_ok=True)
    gene_index = 0
    active_genes = []
    annotation_count = 0
    hc_count = 0
    annotation_block = []

    with (
        pvar_path.open() as pvar,
        annotation_path.open("w") as annotation,
    ):
        annotation.write(
            "variant_key\tgene_id\tgene_symbol\n"
        )

        for line in pvar:
            if line.startswith("##"):
                continue

            if line.startswith("#"):
                columns = line.lstrip("#").rstrip().split()
                chr_column = columns.index("CHROM")
                pos_column = columns.index("POS")
                id_column = columns.index("ID")
                continue

            fields = line.rstrip().split()
            variant_key = fields[id_column]

            if fields[chr_column].removeprefix("chr") != chromosome:
                continue

            pos = int(fields[pos_column])
            while gene_index < len(genes) and pos >= genes[gene_index][0]:
                active_genes.append(genes[gene_index])
                gene_index += 1

            active_genes = [
                gene for gene in active_genes if pos <= gene[1]
            ]

            for _, _, gene_id, gene_symbol in active_genes:
                # Preserve the legacy RNG draw so the selected HC rows stay identical.
                rng.choice(("missense_variant", "synonymous_variant"))
                annotation_block.append((variant_key, gene_id, gene_symbol))
                annotation_count += 1
                if len(annotation_block) == 100:
                    row = annotation_block[rng.randrange(100)]
                    annotation.write("\t".join(row) + "\n")
                    annotation_block.clear()
                    hc_count += 1

        # The legacy generator marked the final incomplete block entirely LC.

    print(f"created: {gene_panel_path}, ({len(genes)} genes)")
    print(
        f"created: {annotation_path}, "
        f"({annotation_count} candidate variant-gene rows, "
        f"{hc_count} synthetic HC retained, {annotation_count - hc_count} discarded)"
    )
    return gene_panel_path, annotation_path


def create_phenotype(
    psam_path: Path,
    out_path: Path,
    num_pheno: int = 1,
    seed: int = 42,
) -> Path:
    if out_path.exists():
        print(f"create_phenotype - exists: {out_path}")
        return out_path
    with psam_path.open() as psam:
        columns = next(psam).lstrip("#").split()
        iid_column = columns.index("IID")
        sample_ids = [
            line.split()[iid_column]
            for line in psam
        ]

    rng = random.Random(seed)
    phenotype_columns = [f"phenotype{i+1}" for i in range(num_pheno)]
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w") as output:
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(["IID", *phenotype_columns])

        for sample_id in sample_ids:
            phenotypes = [rng.gauss(0.0, 1.0) for _ in range(num_pheno)]
            writer.writerow([sample_id, *[f"{phenotype:.12f}" for phenotype in phenotypes]])

    print(f"created: {out_path}, ({len(sample_ids)} samples)")
    return out_path
