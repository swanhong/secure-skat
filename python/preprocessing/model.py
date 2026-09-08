from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Literal
import numpy as np

VariantRole = Literal[
    "shared",
    "public_only",
    "private",
]


@dataclass(frozen=True)
class GeneRef:
    gene_id: str
    gene_symbol: str
    chromosome: str
    order_index: int


@dataclass(frozen=True)
class VariantRef:
    key: str
    position: int
    filter_value: str


@dataclass(frozen=True)
class AnnotationRef:
    variant_key: str
    gene_id: str
    gene_symbol: str
    values: Mapping[str, str]


@dataclass(frozen=True)
class GeneVariants:
    gene: GeneRef
    variants: tuple[VariantRef, ...]


@dataclass(frozen=True)
class GenePlan:
    gene: GeneRef
    variant_roles: tuple[
        tuple[VariantRef, VariantRole],
        ...,
    ]

@dataclass(frozen=True)
class GeneBlock:
    gene: GeneRef
    public_variants: tuple[VariantRef, ...]
    private_variants: tuple[VariantRef, ...]
    public_a: np.ndarray
    public_b: np.ndarray
    private_b: np.ndarray


@dataclass(frozen=True)
class PhenoCovRows:
    sample_ids: tuple[str, ...]
    covariates: tuple[tuple[float, ...], ...]
    phenotypes: tuple[tuple[float, ...], ...]

@dataclass(frozen=True)
class PrepInputs:
    pgen_prefix: Path
    psam_ids: tuple[str, ...]
    gene_panel: tuple[GeneRef, ...]
    variants: tuple[VariantRef, ...]
    annotations: tuple[AnnotationRef, ...]
    annotation_columns: tuple[str, ...]

@dataclass(frozen=True)
class SampleInputs:
    phenotypes: Mapping[str, Mapping[str, float]]
    covariates: Mapping[str, tuple[float, ...]]
    ancestries: Mapping[str, str]
