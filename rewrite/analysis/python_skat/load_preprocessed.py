from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class PlainInput:
    directory: Path
    genes: list[str]
    public_variant_counts: list[int]
    sample_count_a: int
    sample_count_b: int
    covariates: np.ndarray
    phenotypes: np.ndarray


def read_text_matrix(path: Path) -> np.ndarray:
    return np.loadtxt(path, ndmin=2)


def read_plain_input(directory: Path) -> PlainInput:
    directory = directory.resolve(strict=True)
    genes = (directory / "genes.txt").read_text().splitlines()
    public_variant_counts = [
        int(value)
        for value in (directory / "block_sizes.txt").read_text().split()
    ]
    if len(genes) != len(public_variant_counts):
        raise ValueError("genes.txt and block_sizes.txt have different lengths")

    covariates_a = read_text_matrix(directory / "A/cov.txt")
    covariates_b = read_text_matrix(directory / "B/cov.txt")
    phenotypes_a = read_text_matrix(directory / "A/pheno.txt")
    phenotypes_b = read_text_matrix(directory / "B/pheno.txt")
    return PlainInput(
        directory=directory,
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
        input_data.directory / "A/geno" / block,
        input_data.sample_count_a,
        public_count,
    )
    public_b = read_dosages(
        input_data.directory / "B/geno" / block,
        input_data.sample_count_b,
        public_count,
    )
    private_b = read_dosages(
        input_data.directory / "B/private" / block,
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
