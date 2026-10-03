from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class PlainInput:
    directory_a: Path
    directory_b: Path
    genes: list[str]
    public_variant_counts: list[int]
    sample_count_a: int
    sample_count_b: int
    covariates: np.ndarray
    phenotypes: np.ndarray


def read_text_matrix(path: Path) -> np.ndarray:
    return np.loadtxt(path, ndmin=2)


def read_plain_input(directory_a: Path, directory_b: Path | None = None) -> PlainInput:
    directory_a = directory_a.resolve(strict=True)
    directory_b = directory_a if directory_b is None else directory_b.resolve(strict=True)
    for name in ("genes.txt", "block_sizes.txt", "public_variants.tsv"):
        if directory_b != directory_a and (directory_a / name).read_bytes() != (directory_b / name).read_bytes():
            raise ValueError(f"A/B input mismatch: {name}")
    genes = (directory_a / "genes.txt").read_text().splitlines()
    public_variant_counts = [
        int(value)
        for value in (directory_a / "block_sizes.txt").read_text().split()
    ]
    if len(genes) != len(public_variant_counts):
        raise ValueError("genes.txt and block_sizes.txt have different lengths")

    covariates_a = read_text_matrix(directory_a / "A/cov.txt")
    covariates_b = read_text_matrix(directory_b / "B/cov.txt")
    phenotypes_a = read_text_matrix(directory_a / "A/pheno.txt")
    phenotypes_b = read_text_matrix(directory_b / "B/pheno.txt")
    return PlainInput(
        directory_a=directory_a,
        directory_b=directory_b,
        genes=genes,
        public_variant_counts=public_variant_counts,
        sample_count_a=covariates_a.shape[0],
        sample_count_b=covariates_b.shape[0],
        covariates=np.vstack((covariates_a, covariates_b)),
        phenotypes=np.vstack((phenotypes_a, phenotypes_b)),
    )


def read_dosages(path: Path, rows: int, columns: int) -> np.ndarray:
    values = np.fromfile(path, dtype=np.int8)
    return values.reshape(rows, columns).astype(np.float64)


def read_gene_genotype(input_data: PlainInput, gene_index: int) -> np.ndarray:
    block = f"block.{gene_index}.bin"
    public_count = input_data.public_variant_counts[gene_index]
    public_a = read_dosages(
        input_data.directory_a / "A/geno" / block,
        input_data.sample_count_a,
        public_count,
    )
    public_b = read_dosages(
        input_data.directory_b / "B/geno" / block,
        input_data.sample_count_b,
        public_count,
    )
    private_b = read_dosages(
        input_data.directory_b / "B/private" / block,
        input_data.sample_count_b,
        -1,
    )
    private_a = np.zeros(
        (input_data.sample_count_a, private_b.shape[1]),
        dtype=np.float64,
    )
    return np.vstack(
        (
            np.hstack((public_a, private_a)),
            np.hstack((public_b, private_b)),
        )
    )
