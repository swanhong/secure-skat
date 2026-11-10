from __future__ import annotations

import argparse
from pathlib import Path
import sys
import tomllib

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from python.analysis.config import normalize_chromosomes
from python.datasets.onekg import utils


VCF_TEMPLATE = (
    "1kGP_high_coverage_Illumina.chr{chromosome}."
    "filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
)
VCF_BASE_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/"
    "data_collections/1000G_2504_high_coverage/working/"
    "20220422_3202_phased_SNV_INDEL_SV/"
)
PANEL_NAME = "integrated_call_samples_v3.20130502.ALL.panel"
PANEL_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/"
    f"release/20130502/{PANEL_NAME}"
)
GTF_NAME = "gencode.v50.annotation.gtf.gz"
GTF_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/"
    f"Gencode_human/release_50/{GTF_NAME}"
)


def read_psam_ids(path: Path) -> tuple[str, ...]:
    with path.open() as file:
        columns = next(file).lstrip("#").split()
        iid_column = columns.index("IID")
        return tuple(line.split()[iid_column] for line in file if line.strip())


def prepare_1kgenome(root: Path, args) -> None:
    chromosomes = args.chromosome
    generated = root / "generated"
    raw = generated / "raw"
    panel = utils.download_if_missing(PANEL_URL, raw / PANEL_NAME)
    gtf = utils.download_if_missing(GTF_URL, raw / GTF_NAME)

    psam_ids = ()
    pgen_prefixes = []
    for index, chromosome in enumerate(chromosomes):
        vcf_name = VCF_TEMPLATE.format(chromosome=chromosome)
        vcf = utils.download_if_missing(VCF_BASE_URL + vcf_name, raw / vcf_name)
        prefix = utils.create_pgen(
            vcf_path=vcf,
            panel_path=panel,
            out_prefix=generated / "genotype" / f"chr{chromosome}",
            keep_path=generated / "work" / "phase3.keep",
        )
        pgen_prefixes.append(prefix)
        current_ids = read_psam_ids(prefix.with_suffix(".psam"))
        if index == 0:
            psam_ids = current_ids
        elif current_ids != psam_ids:
            raise ValueError(
                f"ordered PSAM IDs differ for chr{chromosomes[0]} and chr{chromosome}"
            )

        gene_panel = generated / "gene_panel" / f"chr{chromosome}.tsv"
        annotation = generated / "annotation" / f"chr{chromosome}.tsv"
        utils.create_inputs_from_gencode(
            gtf_path=gtf,
            pvar_path=prefix.with_suffix(".pvar"),
            gene_panel_path=gene_panel,
            annotation_path=annotation,
            chromosome=str(chromosome),
        )

    chromosome_set = "chr" + "_".join(str(value) for value in chromosomes)
    eigenvec = utils.create_pca(
        pgen_prefixes=pgen_prefixes,
        work_dir=generated / "work" / "pca" / chromosome_set,
    )
    utils.create_ancestry_table(
        panel_path=panel,
        eigenvec_path=eigenvec,
        sample_ids=psam_ids,
        out_path=generated / "ancestry_pred.tsv",
    )
    sample_path = pgen_prefixes[0].with_suffix(".psam")
    utils.create_phenotype(
        sample_path,
        generated / "phenotype.csv",
        num_pheno=args.num_pheno,
        seed=42,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument(
        "--num-pheno",
        type=int,
        default=1,
        help="number of phenotypes to generate",
    )
    args = parser.parse_args()

    with (args.config / "configGlobal.toml").open("rb") as config_file:
        config = tomllib.load(config_file)
    args.chromosome = normalize_chromosomes(config["chromosomes"])

    if args.num_pheno < 1:
        raise ValueError("num-pheno must be at least 1")
    prepare_1kgenome(Path(__file__).resolve().parents[3] / "data" / "1kg", args)


if __name__ == "__main__":
    main()
