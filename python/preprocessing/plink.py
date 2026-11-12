from pathlib import Path

import numpy as np
import pgenlib

from .input import read_psam


def sample_subset(pgen_prefix: Path, sample_ids: tuple[str, ...]) -> tuple[int, np.ndarray]:
    psam_ids = read_psam(Path(f"{pgen_prefix}.psam"))
    sample_index = {sample_id: index for index, sample_id in enumerate(psam_ids)}
    if len(sample_index) != len(psam_ids) or len(set(sample_ids)) != len(sample_ids):
        raise ValueError("duplicate sample IDs")
    try:
        indices = np.array([sample_index[sample_id] for sample_id in sample_ids], dtype=np.uint32)
    except KeyError:
        raise ValueError("sample set mismatch") from None
    return len(psam_ids), indices


def biallelic_variants(pvar, wanted: set[str]) -> tuple[list[int], list[str]]:
    variant_indices, emitted_keys = [], []
    for index in range(pvar.get_variant_ct()):
        key = pvar.get_variant_id(index).decode()
        if key not in wanted or pvar.get_allele_ct(index) != 2:
            continue
        fields = key.split(":", 3)
        if (len(fields) != 4 or fields[0].removeprefix("chr") != pvar.get_variant_chrom(index).decode().removeprefix("chr")
                or fields[1] != str(pvar.get_variant_pos(index))):
            raise ValueError(f"{key}: variant coordinates do not match PVAR")
        if pvar.get_allele_code(index, 0).decode() != fields[2]:
            raise ValueError(f"{key}: reference allele is not the key's REF")
        if pvar.get_allele_code(index, 1).decode() != fields[-1]:
            raise ValueError(f"{key}: counted allele is not the key's ALT")
        variant_indices.append(index)
        emitted_keys.append(key)
    if len(set(emitted_keys)) != len(emitted_keys):
        raise ValueError("duplicate variant IDs")
    return variant_indices, emitted_keys


def read_dosage_batch(reader, variant_indices: list[int], sample_count: int) -> np.ndarray:
    """Read variants × samples (PSAM order) dosages; missing calls stay -9."""
    batch = np.array(variant_indices, dtype=np.uint32)
    dosages = np.empty((len(batch), sample_count), dtype=np.float32)
    reader.read_dosages_list(batch, dosages)
    if not np.isin(dosages, (-9, 0, 1, 2)).all():
        raise ValueError("fractional dosages are unsupported")
    return dosages


def pgen_extract(
    pgen_prefix: Path,
    sample_ids: tuple[str, ...],
    variant_keys: tuple[str, ...],
) -> tuple[np.ndarray, tuple[str, ...]]:
    """Read biallelic ALT dosages, retaining -9, then restore requested row order."""
    if not variant_keys:
        return np.zeros((len(sample_ids), 0), dtype=np.int8), ()

    raw_sample_count, indices = sample_subset(pgen_prefix, sample_ids)
    order = np.argsort(indices)
    restore = np.argsort(order)

    with pgenlib.PvarReader(str(pgen_prefix).encode() + b".pvar") as pvar:
        variant_indices, emitted_keys = biallelic_variants(pvar, set(variant_keys))
        matrix = np.empty((len(sample_ids), len(emitted_keys)), dtype=np.int8)
        if not emitted_keys or not sample_ids:
            return matrix, tuple(emitted_keys)

        with pgenlib.PgenReader(
            str(pgen_prefix).encode() + b".pgen", pvar=pvar,
            raw_sample_ct=raw_sample_count, sample_subset=indices[order],
        ) as reader:
            for start in range(0, len(emitted_keys), 256):
                end = start + 256
                dosages = read_dosage_batch(reader, variant_indices[start:end], len(sample_ids))
                matrix[:, start:end] = dosages[:, restore].T
    return matrix, tuple(emitted_keys)


def pgen_allele_counts(
    pgen_prefix: Path,
    sample_ids: tuple[str, ...],
    variant_keys: tuple[str, ...],
) -> dict[str, int]:
    """ALT allele count per biallelic variant over the samples; equals the column sum of extracted blocks."""
    if not variant_keys or not sample_ids:
        return {}

    raw_sample_count, indices = sample_subset(pgen_prefix, sample_ids)
    counts = {}
    with pgenlib.PvarReader(str(pgen_prefix).encode() + b".pvar") as pvar:
        variant_indices, emitted_keys = biallelic_variants(pvar, set(variant_keys))
        with pgenlib.PgenReader(
            str(pgen_prefix).encode() + b".pgen", pvar=pvar,
            raw_sample_ct=raw_sample_count, sample_subset=np.sort(indices),
        ) as reader:
            for start in range(0, len(emitted_keys), 256):
                end = start + 256
                dosages = read_dosage_batch(reader, variant_indices[start:end], len(sample_ids))
                counts.update(zip(emitted_keys[start:end], dosages.clip(min=0).sum(axis=1).astype(int).tolist()))
    return counts
