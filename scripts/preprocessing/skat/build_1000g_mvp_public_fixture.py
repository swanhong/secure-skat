#!/usr/bin/env python3
r"""Build a 1000G chr22 AoU/MVP public-variant two-party fixture.

The generated public dataset keeps the existing secure-skat PGEN layout:

  dataset/party1
  dataset/party2

where party1 is AoU and party2 is MVP.  MVP publishes V_M per block.  AoU's
public block uses the same V_M order, but variants in V_M \ V_A are masked to
all-reference before PGEN materialization.  Since the secure runtime later
converts ALT dosages by 2 - ALT, this reference fill becomes secure value 2.

The AoU-only hidden set H = V_A \ V_M is materialized separately under
dataset/party1_hidden for the future two-party runtime path.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int
    vid: str
    ref: str
    alt: str


@dataclass(frozen=True)
class SelectedBlock:
    block_id: int
    window_id: int
    overlap: list[Variant]
    mvp_only: list[Variant]
    aou_hidden: list[Variant]

    @property
    def public(self) -> list[Variant]:
        return sorted([*self.overlap, *self.mvp_only], key=lambda v: (v.pos, v.vid))

    @property
    def hidden(self) -> list[Variant]:
        return sorted(self.aou_hidden, key=lambda v: (v.pos, v.vid))


def repo_root_from_script() -> Path:
    return Path(__file__).resolve().parents[3]


def resolve_path(path_text: str | Path) -> Path:
    path = Path(path_text)
    if path.is_absolute():
        return path
    return (repo_root_from_script() / path).resolve()


def pfile_path(prefix: Path, ext: str) -> Path:
    return Path(str(prefix) + ext)


def run(cmd: list[str]) -> None:
    print("+ " + " ".join(map(str, cmd)), flush=True)
    proc = subprocess.run(cmd, text=True, capture_output=True)
    if proc.returncode != 0:
        if proc.stdout:
            print(proc.stdout, file=sys.stderr)
        if proc.stderr:
            print(proc.stderr, file=sys.stderr)
        raise RuntimeError(f"Command failed with exit code {proc.returncode}: {' '.join(map(str, cmd))}")


def ensure_clean_dir(path: Path, force: bool) -> None:
    if path.exists():
        if not force:
            raise FileExistsError(f"{path} already exists; pass --force to overwrite")
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def read_psam_samples(psam_path: Path) -> list[str]:
    samples: list[str] = []
    with psam_path.open() as fh:
        header = fh.readline().strip().split()
        try:
            iid_idx = header.index("IID")
        except ValueError as exc:
            raise ValueError(f"PSAM missing IID column: {psam_path}") from exc
        for line in fh:
            if not line.strip():
                continue
            parts = line.strip().split()
            samples.append(parts[iid_idx])
    return samples


def read_numeric_table(path: Path) -> tuple[list[str], dict[str, list[float]]]:
    with path.open(newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        if len(header) < 2:
            raise ValueError(f"Expected ID plus numeric columns in {path}")
        out: dict[str, list[float]] = {}
        for row in reader:
            if not row:
                continue
            out[row[0]] = [float(x) for x in row[1:]]
    return header[1:], out


def write_party_sample_files(
    party_dir: Path,
    samples: list[str],
    pheno_by_sample: dict[str, list[float]],
    cov_by_sample: dict[str, list[float]],
) -> None:
    party_dir.mkdir(parents=True, exist_ok=True)
    with (party_dir / "sample_keep.txt").open("w") as fh:
        for sample_id in samples:
            fh.write(f"{sample_id}\t{sample_id}\n")
    with (party_dir / "pheno.txt").open("w") as fh:
        for sample_id in samples:
            fh.write(f"{pheno_by_sample[sample_id][0]:.12g}\n")
    with (party_dir / "cov.txt").open("w") as fh:
        for sample_id in samples:
            fh.write("\t".join(f"{x:.12g}" for x in cov_by_sample[sample_id]) + "\n")


def iter_pvar_records(pvar_path: Path) -> Iterable[Variant]:
    if pvar_path.suffix == ".zst":
        proc = subprocess.Popen(["zstdcat", str(pvar_path)], stdout=subprocess.PIPE, text=True)
        assert proc.stdout is not None
        fh: Iterable[str] = proc.stdout
    else:
        proc = None
        fh = pvar_path.open()

    try:
        header: list[str] | None = None
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("##"):
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "#CHROM":
                header = parts
                continue
            if header is None:
                raise ValueError(f"PVAR header not found in {pvar_path}")
            idx = {name: i for i, name in enumerate(header)}
            yield Variant(
                chrom=parts[idx["#CHROM"]],
                pos=int(parts[idx["POS"]]),
                vid=parts[idx["ID"]],
                ref=parts[idx["REF"]],
                alt=parts[idx["ALT"]],
            )
    finally:
        if proc is not None:
            rc = proc.wait()
            if rc != 0:
                raise RuntimeError(f"zstdcat failed for {pvar_path}")
        elif hasattr(fh, "close"):
            fh.close()  # type: ignore[attr-defined]


def pfile_base_args(plink2: str, source_prefix: Path) -> list[str]:
    cmd = [plink2, "--pfile", str(source_prefix)]
    if pfile_path(source_prefix, ".pvar.zst").exists():
        cmd.append("vzs")
    return cmd


def run_freq_counts(plink2: str, source_prefix: Path, keep_path: Path, out_prefix: Path) -> Path:
    cmd = pfile_base_args(plink2, source_prefix)
    cmd.extend(["--keep", str(keep_path), "--freq", "counts", "--out", str(out_prefix)])
    run(cmd)
    acount = out_prefix.with_suffix(".acount")
    if not acount.exists():
        raise FileNotFoundError(f"PLINK did not create {acount}")
    return acount


def read_alt_counts(acount_path: Path) -> dict[str, int]:
    with acount_path.open() as fh:
        header = fh.readline().strip().split()
        id_idx = header.index("ID")
        candidates = ["ALT_CTS", "ALT1_CTS", "ALT_CT", "ALT1_CT", "A1_CT"]
        alt_idx = next((header.index(c) for c in candidates if c in header), None)
        if alt_idx is None:
            raise ValueError(f"Could not find ALT count column in {acount_path}: {header}")
        out: dict[str, int] = {}
        for line in fh:
            if not line.strip():
                continue
            parts = line.strip().split()
            value = parts[alt_idx]
            if value in {".", "NA"}:
                continue
            out[parts[id_idx]] = int(round(float(value)))
    return out


def choose_blocks(
    variants: dict[str, Variant],
    ac_a: dict[str, int],
    ac_m: dict[str, int],
    rare_ac_max: int,
    window_bp: int,
    num_blocks: int,
    max_overlap: int,
    max_mvp_only: int,
    max_aou_hidden: int,
) -> list[SelectedBlock]:
    windows: dict[int, dict[str, list[Variant]]] = {}
    for vid, variant in variants.items():
        a_rare = 1 <= ac_a.get(vid, 0) <= rare_ac_max
        m_rare = 1 <= ac_m.get(vid, 0) <= rare_ac_max
        if not (a_rare or m_rare):
            continue
        window_id = variant.pos // window_bp
        bucket = windows.setdefault(window_id, {"overlap": [], "mvp_only": [], "aou_hidden": []})
        if a_rare and m_rare:
            bucket["overlap"].append(variant)
        elif m_rare:
            bucket["mvp_only"].append(variant)
        else:
            bucket["aou_hidden"].append(variant)

    candidates: list[tuple[int, int, int, int, int]] = []
    for window_id, bucket in windows.items():
        no = len(bucket["overlap"])
        nm = len(bucket["mvp_only"])
        nh = len(bucket["aou_hidden"])
        if no == 0:
            continue
        richness = min(no, max_overlap) + min(nm, max_mvp_only) + min(nh, max_aou_hidden)
        balance = min(no, nm, nh)
        candidates.append((balance, richness, no, nm + nh, window_id))

    candidates.sort(reverse=True)
    selected: list[SelectedBlock] = []
    for _, _, _, _, window_id in candidates:
        bucket = windows[window_id]
        if not bucket["mvp_only"] or not bucket["aou_hidden"]:
            continue
        selected.append(
            SelectedBlock(
                block_id=len(selected) + 1,
                window_id=window_id,
                overlap=sorted(bucket["overlap"], key=lambda v: (v.pos, v.vid))[:max_overlap],
                mvp_only=sorted(bucket["mvp_only"], key=lambda v: (v.pos, v.vid))[:max_mvp_only],
                aou_hidden=sorted(bucket["aou_hidden"], key=lambda v: (v.pos, v.vid))[:max_aou_hidden],
            )
        )
        if len(selected) == num_blocks:
            return selected

    raise RuntimeError(
        "Could not find enough windows with overlap, MVP-only public, and AoU-hidden variants. "
        f"Found {len(selected)} / {num_blocks}."
    )


def write_extract(path: Path, variants: list[Variant]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{v.vid}\n" for v in variants))


def export_vcf(plink2: str, source_prefix: Path, keep_path: Path, extract_path: Path, out_prefix: Path) -> Path:
    cmd = pfile_base_args(plink2, source_prefix)
    cmd.extend(
        [
            "--keep",
            str(keep_path),
            "--extract",
            str(extract_path),
            "--export",
            "vcf",
            "--out",
            str(out_prefix),
        ]
    )
    run(cmd)
    vcf_path = out_prefix.with_suffix(".vcf")
    if not vcf_path.exists():
        raise FileNotFoundError(f"PLINK did not create {vcf_path}")
    return vcf_path


def parse_gt_to_dosage(field: str, gt_index: int) -> float:
    pieces = field.split(":")
    if gt_index >= len(pieces):
        return np.nan
    gt = pieces[gt_index]
    if "." in gt:
        return np.nan
    alleles = gt.replace("|", "/").split("/")
    return float(sum(1 for allele in alleles if allele == "1"))


def normalize_vcf_sample_name(sample: str) -> str:
    # PLINK2 VCF export combines FID/IID as FID_IID.  Our 1000G PSAM uses
    # FID == IID, while secure-skat sample_keep/pheno/cov files use IID.
    if "_" in sample:
        left, right = sample.split("_", 1)
        if left == right:
            return left
    return sample


def read_vcf_alt_matrix(vcf_path: Path, expected_variants: list[Variant]) -> tuple[list[str], np.ndarray]:
    rows_by_id: dict[str, tuple[list[str], list[float]]] = {}
    samples: list[str] | None = None
    with vcf_path.open() as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "#CHROM":
                samples = [normalize_vcf_sample_name(sample) for sample in parts[9:]]
                continue
            if samples is None:
                raise ValueError(f"VCF header not found in {vcf_path}")
            fmt = parts[8].split(":")
            gt_index = fmt.index("GT")
            dosages = [parse_gt_to_dosage(field, gt_index) for field in parts[9:]]
            rows_by_id[parts[2]] = (parts[:8], dosages)
    if samples is None:
        raise ValueError(f"VCF samples not found in {vcf_path}")

    matrix = np.zeros((len(samples), len(expected_variants)), dtype=float)
    for j, variant in enumerate(expected_variants):
        if variant.vid not in rows_by_id:
            raise ValueError(f"Variant {variant.vid} missing from {vcf_path}")
        matrix[:, j] = np.asarray(rows_by_id[variant.vid][1], dtype=float)
    return samples, matrix


def dosage_to_gt(value: float) -> str:
    if not np.isfinite(value):
        return "./."
    dosage = int(round(float(value)))
    if dosage == 0:
        return "0/0"
    if dosage == 1:
        return "0/1"
    if dosage == 2:
        return "1/1"
    raise ValueError(f"Unexpected ALT dosage {value}")


def write_synthetic_vcf(path: Path, samples: list[str], variants: list[Variant], alt_matrix: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##source=secure-skat-1000g-mvp-public-fixture\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        fh.write("\t".join(samples) + "\n")
        for j, variant in enumerate(variants):
            gts = [dosage_to_gt(alt_matrix[i, j]) for i in range(len(samples))]
            fh.write(
                f"{variant.chrom}\t{variant.pos}\t{variant.vid}\t{variant.ref}\t{variant.alt}"
                f"\t.\tPASS\t.\tGT\t" + "\t".join(gts) + "\n"
            )


def vcf_to_pgen(plink2: str, vcf_path: Path, out_prefix: Path) -> None:
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        plink2,
        "--vcf",
        str(vcf_path),
        "--double-id",
        "--indiv-sort",
        "none",
        "--make-pgen",
        "--out",
        str(out_prefix),
    ]
    run(cmd)


def write_block_metadata_files(party_dir: Path, blocks: list[list[Variant]]) -> None:
    party_dir.mkdir(parents=True, exist_ok=True)
    (party_dir / "chrom_sizes.txt").write_text("".join(f"{len(block)}\n" for block in blocks))
    with (party_dir / "snp_ids.txt").open("w") as ids_fh, (party_dir / "snp_pos.txt").open("w") as pos_fh:
        for block in blocks:
            for variant in block:
                ids_fh.write(f"{variant.vid}\n")
                pos_fh.write(f"{variant.chrom}\t{variant.pos}\n")


def secure_oriented(alt_matrix: np.ndarray) -> np.ndarray:
    return np.where(np.isfinite(alt_matrix), 2.0 - alt_matrix, 0.0)


def compute_weight(secure_matrix: np.ndarray) -> np.ndarray:
    n_total = secure_matrix.shape[0]
    p_bar = secure_matrix.sum(axis=0) / (2.0 * n_total)
    p = 1.0 - p_bar
    beta_base = np.maximum(p, p_bar)
    return 25.0 * np.power(beta_base, 24)


def compute_block_stats(secure_parts: list[np.ndarray], residual: np.ndarray) -> tuple[float, float]:
    G = np.vstack(secure_parts)
    w = compute_weight(G)
    score = G.T @ residual
    q_skat = float(np.sum((w * w) * (score * score)))
    b = float(np.sum(w * score))
    return q_skat, b


def fit_residual(pheno_a: np.ndarray, cov_a: np.ndarray, pheno_m: np.ndarray, cov_m: np.ndarray) -> np.ndarray:
    y = np.concatenate([pheno_a, pheno_m])
    X = np.vstack([cov_a, cov_m])
    design = np.column_stack([np.ones(y.shape[0]), X])
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    return y - design @ beta


def write_reference_json(
    path: Path,
    blocks: list[SelectedBlock],
    public_a_secure: list[np.ndarray],
    public_m_secure: list[np.ndarray],
    hidden_a_secure: list[np.ndarray],
    pheno_a: np.ndarray,
    cov_a: np.ndarray,
    pheno_m: np.ndarray,
    cov_m: np.ndarray,
) -> None:
    residual = fit_residual(pheno_a, cov_a, pheno_m, cov_m)
    n_a = len(pheno_a)
    n_m = len(pheno_m)
    r_a = residual[:n_a]
    r_m = residual[n_a:]

    total_full_q = 0.0
    total_full_b = 0.0
    total_decomp_q = 0.0
    total_decomp_b = 0.0
    block_results = []

    for idx, block in enumerate(blocks):
        q_public, b_public = compute_block_stats([public_a_secure[idx], public_m_secure[idx]], residual)

        hidden_m_secure = np.full((n_m, hidden_a_secure[idx].shape[1]), 2.0, dtype=float)
        q_hidden_full, b_hidden_full = compute_block_stats([hidden_a_secure[idx], hidden_m_secure], residual)

        hidden_score_corr = hidden_a_secure[idx].T @ r_a - 2.0 * float(np.sum(r_a))
        hidden_score_full = np.vstack([hidden_a_secure[idx], hidden_m_secure]).T @ residual
        hidden_score_diff = float(np.max(np.abs(hidden_score_corr - hidden_score_full))) if hidden_score_full.size else 0.0

        full_secure = [
            np.concatenate([public_a_secure[idx], hidden_a_secure[idx]], axis=1),
            np.concatenate([public_m_secure[idx], hidden_m_secure], axis=1),
        ]
        q_full, b_full = compute_block_stats(full_secure, residual)
        q_decomp = q_public + q_hidden_full
        b_decomp = b_public + b_hidden_full

        total_full_q += q_full
        total_full_b += b_full
        total_decomp_q += q_decomp
        total_decomp_b += b_decomp
        block_results.append(
            {
                "block_id": block.block_id,
                "public_variant_count": len(block.public),
                "hidden_variant_count": len(block.hidden),
                "q_full": q_full,
                "q_decomposed": q_decomp,
                "burden_linear_full": b_full,
                "burden_linear_decomposed": b_decomp,
                "max_hidden_score_correction_abs_diff": hidden_score_diff,
            }
        )

    q_diff = abs(total_full_q - total_decomp_q)
    b_diff = abs(total_full_b - total_decomp_b)
    if q_diff > 1e-7 * max(1.0, abs(total_full_q)):
        raise AssertionError(f"Full/decomposed SKAT mismatch: {total_full_q} vs {total_decomp_q}")
    if b_diff > 1e-7 * max(1.0, abs(total_full_b)):
        raise AssertionError(f"Full/decomposed burden-linear mismatch: {total_full_b} vs {total_decomp_b}")

    payload = {
        "n_aou": n_a,
        "n_mvp": n_m,
        "residual_sum_total": float(np.sum(residual)),
        "residual_sum_aou": float(np.sum(r_a)),
        "residual_sum_mvp": float(np.sum(r_m)),
        "q_full_total": total_full_q,
        "q_decomposed_total": total_decomp_q,
        "burden_linear_full_total": total_full_b,
        "burden_linear_decomposed_total": total_decomp_b,
        "burden_q_full": float(total_full_b * total_full_b),
        "burden_q_decomposed": float(total_decomp_b * total_decomp_b),
        "blocks": block_results,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def replace_toml_assignment(text: str, key: str, value: str) -> str:
    lines = text.splitlines()
    out: list[str] = []
    replaced = False
    for line in lines:
        if line.strip().startswith(f"{key} "):
            prefix = line.split("=", 1)[0]
            out.append(f"{prefix}= {value}")
            replaced = True
        else:
            out.append(line)
    if not replaced:
        raise ValueError(f"Could not patch TOML key {key}")
    return "\n".join(out) + "\n"


def upsert_root_toml_assignment(text: str, key: str, value: str) -> str:
    lines = text.splitlines()
    first_table = next((i for i, line in enumerate(lines) if line.lstrip().startswith("[")), len(lines))
    for i in range(first_table):
        if lines[i].strip().startswith(f"{key} "):
            prefix = lines[i].split("=", 1)[0]
            lines[i] = f"{prefix}= {value}"
            return "\n".join(lines) + "\n"
    lines.insert(first_table, f"{key} = {value}")
    return "\n".join(lines) + "\n"


def emit_config(template_dir: Path, config_dir: Path, dataset_dir: Path, keys_dir: Path, n_a: int, n_m: int, num_covs: int, num_public_variants: int, num_blocks: int) -> None:
    if config_dir.exists():
        shutil.rmtree(config_dir)
    shutil.copytree(template_dir, config_dir)
    global_path = config_dir / "configGlobal.toml"
    text = global_path.read_text()
    updates = {
        "num_inds": f"[0, {n_a}, {n_m}]",
        "num_snps": str(num_public_variants),
        "num_covs": str(num_covs),
        "geno_file_format": '"pgen"',
        "geno_num_blocks": str(num_blocks),
        "use_precomputed_geno_count": "false",
        "skip_qc": "true",
        "skip_pca": "true",
        "gmiss": "1.0",
        "maf_lb": "0.0",
        "hwe_ub": "1.0e18",
        "snp_dist_thres": "0",
        "blocks_for_assoc_test": "[]",
    }
    for key, value in updates.items():
        text = replace_toml_assignment(text, key, value)
    global_path.write_text(text)

    for party_id in (1, 2):
        party_dir = dataset_dir / f"party{party_id}"
        local_path = config_dir / f"configLocal.Party{party_id}.toml"
        text = local_path.read_text()
        updates = {
            "shared_keys_path": f'"{keys_dir}"',
            "geno_binary_file_prefix": f'"{party_dir / "geno" / "chr%d"}"',
            "geno_num_blocks": str(num_blocks),
            "geno_block_size_file": f'"{party_dir / "chrom_sizes.txt"}"',
            "pheno_file": f'"{party_dir / "pheno.txt"}"',
            "covar_file": f'"{party_dir / "cov.txt"}"',
            "snp_position_file": f'"{party_dir / "snp_pos.txt"}"',
            "sample_keep_file": f'"{party_dir / "sample_keep.txt"}"',
            "snp_ids_file": f'"{party_dir / "snp_ids.txt"}"',
            "geno_count_file": '""',
        }
        for key, value in updates.items():
            text = replace_toml_assignment(text, key, value)
        local_path.write_text(text)

    party0_path = config_dir / "configLocal.Party0.toml"
    text = party0_path.read_text()
    text = replace_toml_assignment(text, "shared_keys_path", f'"{keys_dir}"')
    party0_path.write_text(text)


def emit_config_2party(config_dir: Path, config_2party_dir: Path, dataset_dir: Path, num_hidden_blocks: int) -> None:
    if config_2party_dir.exists():
        shutil.rmtree(config_2party_dir)
    shutil.copytree(config_dir, config_2party_dir)

    hidden_dir = dataset_dir / "party1_hidden"
    global_path = config_2party_dir / "configGlobal.toml"
    text = global_path.read_text()
    updates = {
        "rare_variant_set_mode": '"mvp_public_union"',
        "public_variant_party_id": "2",
        "private_variant_party_id": "1",
        "hidden_geno_file_format": '"pgen"',
        "hidden_geno_num_blocks": str(num_hidden_blocks),
        "hidden_geno_binary_file_prefix": f'"{hidden_dir / "geno" / "chr%d"}"',
        "hidden_geno_block_size_file": f'"{hidden_dir / "chrom_sizes.txt"}"',
        "hidden_snp_ids_file": f'"{hidden_dir / "snp_ids.txt"}"',
        "hidden_snp_position_file": f'"{hidden_dir / "snp_pos.txt"}"',
        "hidden_sample_keep_file": f'"{hidden_dir / "sample_keep.txt"}"',
    }
    for key, value in updates.items():
        text = upsert_root_toml_assignment(text, key, value)
    global_path.write_text(text)


def ensure_shared_keys(keys_dir: Path) -> None:
    keys_dir.mkdir(parents=True, exist_ok=True)
    for name in ["shared_key_global.bin", "shared_key_0_1.bin", "shared_key_0_2.bin", "shared_key_1_2.bin"]:
        path = keys_dir / name
        if not path.exists():
            path.write_bytes(os.urandom(32))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    root = repo_root_from_script()
    parser.add_argument("--source-prefix", default=str(root / "datasets/1000g_source/pgen/1000g.chr22"))
    parser.add_argument("--pheno-file", default=str(root / "datasets/1000g_source/pheno/phenotype_data.csv"))
    parser.add_argument("--cov-file", default=str(root / "datasets/1000g_source/pheno/covariate_data.csv"))
    parser.add_argument("--out-root", default=str(root / "datasets/1000g_all_chr22_mvp_public_2party"))
    parser.add_argument("--config-template-dir", default=str(root / "config"))
    parser.add_argument("--plink2", default="plink2")
    parser.add_argument("--n-per-party", type=int, default=64)
    parser.add_argument("--seed", type=int, default=260530)
    parser.add_argument("--num-blocks", type=int, default=2)
    parser.add_argument("--window-bp", type=int, default=50_000)
    parser.add_argument("--rare-ac-max", type=int, default=6)
    parser.add_argument("--max-overlap", type=int, default=16)
    parser.add_argument("--max-mvp-only", type=int, default=16)
    parser.add_argument("--max-aou-hidden", type=int, default=16)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source_prefix = resolve_path(args.source_prefix)
    out_root = resolve_path(args.out_root)
    dataset_dir = out_root / "dataset"
    metadata_dir = out_root / "metadata"
    reference_dir = out_root / "reference"
    work_dir = out_root / "work"
    config_dir = out_root / "config"
    config_2party_dir = out_root / "config_2party"
    keys_dir = out_root / "keys"

    ensure_clean_dir(out_root, args.force)
    dataset_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)

    pheno_cols, pheno_by_sample = read_numeric_table(resolve_path(args.pheno_file))
    cov_cols, cov_by_sample = read_numeric_table(resolve_path(args.cov_file))
    if len(pheno_cols) != 1:
        raise ValueError("This fixture expects one phenotype column")

    psam_samples = read_psam_samples(pfile_path(source_prefix, ".psam"))
    eligible = [sample for sample in psam_samples if sample in pheno_by_sample and sample in cov_by_sample]
    rng = random.Random(args.seed)
    shuffled = eligible[:]
    rng.shuffle(shuffled)
    if len(shuffled) < 2 * args.n_per_party:
        raise RuntimeError(f"Need {2 * args.n_per_party} eligible samples, found {len(shuffled)}")
    aou_set = set(shuffled[: args.n_per_party])
    mvp_set = set(shuffled[args.n_per_party : 2 * args.n_per_party])
    aou_samples = [sample for sample in psam_samples if sample in aou_set]
    mvp_samples = [sample for sample in psam_samples if sample in mvp_set]

    for party_name, samples in [("party1", aou_samples), ("party2", mvp_samples), ("party1_hidden", aou_samples)]:
        write_party_sample_files(dataset_dir / party_name, samples, pheno_by_sample, cov_by_sample)

    aou_keep = dataset_dir / "party1" / "sample_keep.txt"
    mvp_keep = dataset_dir / "party2" / "sample_keep.txt"
    ac_a = read_alt_counts(run_freq_counts(args.plink2, source_prefix, aou_keep, work_dir / "aou_freq"))
    ac_m = read_alt_counts(run_freq_counts(args.plink2, source_prefix, mvp_keep, work_dir / "mvp_freq"))

    pvar_path = pfile_path(source_prefix, ".pvar.zst")
    if not pvar_path.exists():
        pvar_path = pfile_path(source_prefix, ".pvar")
    variants = {variant.vid: variant for variant in iter_pvar_records(pvar_path)}
    blocks = choose_blocks(
        variants,
        ac_a,
        ac_m,
        args.rare_ac_max,
        args.window_bp,
        args.num_blocks,
        args.max_overlap,
        args.max_mvp_only,
        args.max_aou_hidden,
    )
    nontrivial_mask_variants = [
        variant.vid
        for block in blocks
        for variant in block.mvp_only
        if ac_a.get(variant.vid, 0) > 0
    ]
    if not nontrivial_mask_variants:
        raise AssertionError("Selected blocks have no MVP-only public variant with AoU raw AC > 0; masking override is not tested")

    public_blocks = [block.public for block in blocks]
    hidden_blocks = [block.hidden for block in blocks]
    write_block_metadata_files(dataset_dir / "party1", public_blocks)
    write_block_metadata_files(dataset_dir / "party2", public_blocks)
    write_block_metadata_files(dataset_dir / "party1_hidden", hidden_blocks)

    public_a_secure: list[np.ndarray] = []
    public_m_secure: list[np.ndarray] = []
    hidden_a_secure: list[np.ndarray] = []

    with (metadata_dir / "two_party_variant_sets.tsv").open("w") as meta_fh:
        meta_fh.write("block_id\tchrom\tpos\tid\tref\talt\tcategory\tac_aou\tac_mvp\n")
        for block in blocks:
            for category, variants_for_category in [
                ("overlap", block.overlap),
                ("mvp_only_public", block.mvp_only),
                ("aou_hidden", block.aou_hidden),
            ]:
                for variant in variants_for_category:
                    meta_fh.write(
                        f"{block.block_id}\t{variant.chrom}\t{variant.pos}\t{variant.vid}\t"
                        f"{variant.ref}\t{variant.alt}\t{category}\t"
                        f"{ac_a.get(variant.vid, 0)}\t{ac_m.get(variant.vid, 0)}\n"
                    )

    with (metadata_dir / "block_summary.tsv").open("w") as summary_fh:
        summary_fh.write("block_id\twindow_id\toverlap\tmvp_only_public\taou_hidden\tpublic_total\tunion_total\n")
        for block in blocks:
            summary_fh.write(
                f"{block.block_id}\t{block.window_id}\t{len(block.overlap)}\t{len(block.mvp_only)}\t"
                f"{len(block.aou_hidden)}\t{len(block.public)}\t{len(block.public) + len(block.hidden)}\n"
            )

    for block in blocks:
        public_extract = work_dir / "extract" / f"block{block.block_id}_public.snps.txt"
        hidden_extract = work_dir / "extract" / f"block{block.block_id}_hidden.snps.txt"
        write_extract(public_extract, block.public)
        write_extract(hidden_extract, block.hidden)

        aou_public_vcf = export_vcf(args.plink2, source_prefix, aou_keep, public_extract, work_dir / f"aou_public_block{block.block_id}")
        mvp_public_vcf = export_vcf(args.plink2, source_prefix, mvp_keep, public_extract, work_dir / f"mvp_public_block{block.block_id}")
        aou_hidden_vcf = export_vcf(args.plink2, source_prefix, aou_keep, hidden_extract, work_dir / f"aou_hidden_block{block.block_id}")

        aou_samples_vcf, aou_public_alt = read_vcf_alt_matrix(aou_public_vcf, block.public)
        mvp_samples_vcf, mvp_public_alt = read_vcf_alt_matrix(mvp_public_vcf, block.public)
        hidden_samples_vcf, aou_hidden_alt = read_vcf_alt_matrix(aou_hidden_vcf, block.hidden)

        if aou_samples_vcf != aou_samples:
            raise AssertionError("AoU VCF sample order mismatch")
        if mvp_samples_vcf != mvp_samples:
            raise AssertionError("MVP VCF sample order mismatch")
        if hidden_samples_vcf != aou_samples:
            raise AssertionError("AoU hidden VCF sample order mismatch")

        masked_aou_public_alt = aou_public_alt.copy()
        mvp_only_ids = {variant.vid for variant in block.mvp_only}
        for j, variant in enumerate(block.public):
            if variant.vid in mvp_only_ids:
                masked_aou_public_alt[:, j] = 0.0

        masked_secure = secure_oriented(masked_aou_public_alt)
        for j, variant in enumerate(block.public):
            if variant.vid in mvp_only_ids and not np.all(masked_secure[:, j] == 2.0):
                raise AssertionError(f"AoU masked public variant is not secure value 2: {variant.vid}")

        aou_synth_vcf = work_dir / "synthetic_vcf" / f"party1_public_chr{block.block_id}.vcf"
        mvp_synth_vcf = work_dir / "synthetic_vcf" / f"party2_public_chr{block.block_id}.vcf"
        hidden_synth_vcf = work_dir / "synthetic_vcf" / f"party1_hidden_chr{block.block_id}.vcf"
        write_synthetic_vcf(aou_synth_vcf, aou_samples, block.public, masked_aou_public_alt)
        write_synthetic_vcf(mvp_synth_vcf, mvp_samples, block.public, mvp_public_alt)
        write_synthetic_vcf(hidden_synth_vcf, aou_samples, block.hidden, aou_hidden_alt)

        vcf_to_pgen(args.plink2, aou_synth_vcf, dataset_dir / "party1" / "geno" / f"chr{block.block_id}")
        vcf_to_pgen(args.plink2, mvp_synth_vcf, dataset_dir / "party2" / "geno" / f"chr{block.block_id}")
        vcf_to_pgen(args.plink2, hidden_synth_vcf, dataset_dir / "party1_hidden" / "geno" / f"chr{block.block_id}")

        public_a_secure.append(masked_secure)
        public_m_secure.append(secure_oriented(mvp_public_alt))
        hidden_a_secure.append(secure_oriented(aou_hidden_alt))

    for block in blocks:
        if not block.overlap:
            raise AssertionError(f"Block {block.block_id} has no intersection")
        if {v.vid for v in block.public} & {v.vid for v in block.hidden}:
            raise AssertionError(f"Block {block.block_id} has V_M/H overlap")

    pheno_a = np.asarray([pheno_by_sample[s][0] for s in aou_samples], dtype=float)
    cov_a = np.asarray([cov_by_sample[s] for s in aou_samples], dtype=float)
    pheno_m = np.asarray([pheno_by_sample[s][0] for s in mvp_samples], dtype=float)
    cov_m = np.asarray([cov_by_sample[s] for s in mvp_samples], dtype=float)
    write_reference_json(
        reference_dir / "full_union_reference.json",
        blocks,
        public_a_secure,
        public_m_secure,
        hidden_a_secure,
        pheno_a,
        cov_a,
        pheno_m,
        cov_m,
    )

    ensure_shared_keys(keys_dir)
    emit_config(
        resolve_path(args.config_template_dir),
        config_dir,
        dataset_dir,
        keys_dir,
        len(aou_samples),
        len(mvp_samples),
        len(cov_cols),
        sum(len(block.public) for block in blocks),
        len(blocks),
    )
    emit_config_2party(config_dir, config_2party_dir, dataset_dir, len(blocks))

    two_party_manifest = {
        "roles": {
            "party1": "AoU",
            "party2": "MVP",
            "party0": "non-data hub",
        },
        "variant_set_mode": "mvp_public_union",
        "public_variant_party_id": 2,
        "private_variant_party_id": 1,
        "paths": {
            "dataset_root": str(dataset_dir),
            "public_aou_masked": str(dataset_dir / "party1"),
            "public_mvp_actual": str(dataset_dir / "party2"),
            "hidden_aou": str(dataset_dir / "party1_hidden"),
            "shared_config": str(config_dir),
            "two_party_config": str(config_2party_dir),
            "reference": str(reference_dir / "full_union_reference.json"),
        },
        "block_counts": {
            "public_blocks": len(blocks),
            "hidden_blocks": len(blocks),
            "public_variants": sum(len(block.public) for block in blocks),
            "hidden_variants": sum(len(block.hidden) for block in blocks),
            "union_variants": sum(len(block.public) + len(block.hidden) for block in blocks),
        },
        "validation": {
            "mvp_only_public_variants_with_aou_raw_ac_gt_0": len(nontrivial_mask_variants),
            "raw_reference_note": "full_union_reference.json stores raw unscaled SKAT and burden-linear statistics plus raw burden_q; secure outputs are alpha-scaled.",
        },
    }
    (metadata_dir / "two_party_manifest.json").write_text(json.dumps(two_party_manifest, indent=2, sort_keys=True) + "\n")

    manifest = {
        "source_prefix": str(source_prefix),
        "out_root": str(out_root),
        "seed": args.seed,
        "n_per_party": args.n_per_party,
        "rare_ac_max": args.rare_ac_max,
        "window_bp": args.window_bp,
        "num_blocks": args.num_blocks,
        "max_overlap": args.max_overlap,
        "max_mvp_only": args.max_mvp_only,
        "max_aou_hidden": args.max_aou_hidden,
        "aou_samples": len(aou_samples),
        "mvp_samples": len(mvp_samples),
        "public_variants": sum(len(block.public) for block in blocks),
        "hidden_variants": sum(len(block.hidden) for block in blocks),
        "mvp_only_public_variants_with_aou_raw_ac_gt_0": len(nontrivial_mask_variants),
    }
    (out_root / "fixture_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"Fixture written to {out_root}")
    print(f"Public dataset root: {dataset_dir}")
    print(f"Config directory: {config_dir}")
    print(f"2-party config directory: {config_2party_dir}")


if __name__ == "__main__":
    main()
