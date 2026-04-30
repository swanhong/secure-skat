#!/usr/bin/env python3

"""Build a reduced chr22 dataset around rare coding/splice anchor windows.

This script takes a single-block chr22 dataset plus a VEP annotation table,
selects rare coding/splice anchor variants, chooses a manageable set of
fixed-width windows around those anchors, and materializes a new multi-block
dataset where each selected window/locus becomes one block.

The output dataset reuses the original phenotype/covariate/sample files and
subsets the original PGEN, SNP metadata, and precomputed genotype-count binary
to the selected windows. It can also emit a matching config directory by
copying and patching an existing template config.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import shutil
import subprocess
import sys
from bisect import bisect_left, bisect_right
from pathlib import Path

import numpy as np


ALLOWED_CONSEQUENCES = {
    "missense_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "inframe_insertion",
    "inframe_deletion",
    "splice_donor_variant",
    "splice_acceptor_variant",
}

CONSEQUENCE_SEVERITY = {
    "transcript_ablation": 0,
    "splice_acceptor_variant": 1,
    "splice_donor_variant": 1,
    "stop_gained": 2,
    "frameshift_variant": 3,
    "stop_lost": 4,
    "start_lost": 5,
    "inframe_insertion": 6,
    "inframe_deletion": 6,
    "missense_variant": 7,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Select rare coding/splice anchor windows from a single-block chr22 "
            "dataset and materialize a reduced multi-block dataset."
        )
    )
    parser.add_argument(
        "--source-dataset",
        required=True,
        help="Root of the existing chr22 dataset (for example .local/datasets/1000g_all_chr22).",
    )
    parser.add_argument(
        "--vep",
        required=True,
        help="Path to the bgzipped/tab-delimited VEP output table for the same source dataset.",
    )
    parser.add_argument(
        "--out-dataset",
        required=True,
        help="Root directory for the generated reduced dataset.",
    )
    parser.add_argument(
        "--margin-bp",
        type=int,
        default=50000,
        help="Half-width of each anchor-centered window in base pairs (default: 50000).",
    )
    parser.add_argument(
        "--maf-threshold",
        type=float,
        default=0.01,
        help="Maximum cohort MAF used to define rare variants (default: 0.01).",
    )
    parser.add_argument(
        "--min-anchor-count",
        type=int,
        default=20,
        help="Minimum number of anchor variants inside a candidate window (default: 20).",
    )
    parser.add_argument(
        "--min-window-variants",
        type=int,
        default=500,
        help="Minimum total number of variants inside a candidate window (default: 500).",
    )
    parser.add_argument(
        "--max-windows",
        type=int,
        default=16,
        help=(
            "Maximum number of selected windows/loci to keep "
            "(default: 16, use 0 to keep all after selection)."
        ),
    )
    parser.add_argument(
        "--selection-mode",
        choices=("greedy", "merged"),
        default="greedy",
        help=(
            "Window selection strategy: 'greedy' keeps the highest-scoring "
            "non-overlapping candidates, while 'merged' unions overlapping "
            "candidates into merged loci (default: greedy)."
        ),
    )
    parser.add_argument(
        "--plink2",
        default="plink2",
        help="Path to the plink2 executable (default: plink2).",
    )
    parser.add_argument(
        "--config-template-dir",
        default=".local/config_1kg_all_chr22",
        help="Template config directory to copy and patch (default: .local/config_1kg_all_chr22).",
    )
    parser.add_argument(
        "--config-out-dir",
        default=None,
        help=(
            "Optional output config directory. When omitted, a sibling .local/config_<dataset-name> "
            "directory is created."
        ),
    )
    parser.add_argument(
        "--skip-config",
        action="store_true",
        help="Skip config directory generation.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output dataset/config directories.",
    )
    args = parser.parse_args()

    if args.margin_bp <= 0:
        parser.error("--margin-bp must be positive")
    if not (0.0 < args.maf_threshold < 0.5):
        parser.error("--maf-threshold must be in (0, 0.5)")
    if args.min_anchor_count <= 0:
        parser.error("--min-anchor-count must be positive")
    if args.min_window_variants <= 0:
        parser.error("--min-window-variants must be positive")
    if args.max_windows < 0:
        parser.error("--max-windows must be non-negative")

    return args


def read_lines(path: Path) -> list[str]:
    with path.open() as fh:
        return [line.rstrip("\n") for line in fh]


def read_snp_metadata(dataset_root: Path) -> tuple[list[str], np.ndarray]:
    ids_path = dataset_root / "party1" / "snp_ids.txt"
    pos_path = dataset_root / "party1" / "snp_pos.txt"

    snp_ids = [line for line in read_lines(ids_path) if line]
    positions = []
    with pos_path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            toks = line.split("\t")
            pos = toks[1] if len(toks) >= 2 else toks[0]
            positions.append(int(pos))

    if len(snp_ids) != len(positions):
        raise ValueError(
            f"Variant metadata length mismatch: {ids_path} has {len(snp_ids)} IDs "
            f"but {pos_path} has {len(positions)} positions"
        )

    return snp_ids, np.asarray(positions, dtype=np.int64)


def load_party_gcount(path: Path, expected_variants: int) -> np.ndarray:
    raw = np.fromfile(path, dtype=np.uint32)
    if raw.size != 6 * expected_variants:
        raise ValueError(
            f"{path} has {raw.size} uint32 values, expected {6 * expected_variants} "
            f"for a 6 x {expected_variants} matrix"
        )
    return raw.reshape(6, expected_variants)


def compute_combined_maf(dataset_root: Path, expected_variants: int) -> np.ndarray:
    party_counts = []
    for party_name in ("party1", "party2"):
        gcount_path = dataset_root / party_name / "all.gcount.transpose.bin"
        party_counts.append(load_party_gcount(gcount_path, expected_variants))

    combined = party_counts[0][:3, :].astype(np.uint64) + party_counts[1][:3, :].astype(np.uint64)
    hom_ref = combined[0, :]
    het = combined[1, :]
    hom_alt = combined[2, :]

    alt_alleles = het + 2 * hom_alt
    called_genotypes = hom_ref + het + hom_alt
    denom = 2 * called_genotypes
    af = np.divide(
        alt_alleles,
        denom,
        out=np.zeros_like(alt_alleles, dtype=np.float64),
        where=denom > 0,
    )
    return np.minimum(af, 1.0 - af)


def consequence_terms(raw_value: str) -> set[str]:
    out = set()
    for piece in re.split(r"[,&]", raw_value or ""):
        piece = piece.strip()
        if piece:
            out.add(piece)
    return out


def normalize_uploaded_variation(uploaded: str) -> str | None:
    toks = uploaded.split("_", 2)
    if len(toks) != 3 or "/" not in toks[2]:
        return None
    ref, alt = toks[2].split("/", 1)
    if not ref or not alt:
        return None
    chrom = toks[0].removeprefix("chr")
    pos = toks[1]
    return f"{chrom}:{pos}:{ref}:{alt}"


def annotation_rank(row: dict[str, str]) -> tuple[int, int, int, int]:
    terms = consequence_terms(row["Consequence"])
    protein_coding = int(row["BIOTYPE"] == "protein_coding")
    canonical = int(row["CANONICAL"] == "YES")
    allowed = int(bool(terms & ALLOWED_CONSEQUENCES))
    best_severity = min((CONSEQUENCE_SEVERITY.get(term, 999) for term in terms), default=999)
    return protein_coding, canonical, allowed, -best_severity


def load_best_vep_annotations(vep_path: Path, allowed_variant_ids: set[str]) -> dict[str, dict[str, str]]:
    annotations: dict[str, dict[str, str]] = {}

    with gzip.open(vep_path, "rt") as fh:
        reader = csv.DictReader((line for line in fh if not line.startswith("##")), delimiter="\t")
        uploaded_key = reader.fieldnames[0]

        for row in reader:
            variant_id = normalize_uploaded_variation(row[uploaded_key])
            if variant_id is None or variant_id not in allowed_variant_ids:
                continue

            row = dict(row)
            row["variant_id"] = variant_id

            prev = annotations.get(variant_id)
            if prev is None or annotation_rank(row) > annotation_rank(prev):
                annotations[variant_id] = row

    return annotations


def build_anchor_rows(
    snp_ids: list[str],
    positions: np.ndarray,
    cohort_maf: np.ndarray,
    annotations: dict[str, dict[str, str]],
    maf_threshold: float,
) -> list[dict[str, object]]:
    anchors = []

    for idx, variant_id in enumerate(snp_ids):
        annot = annotations.get(variant_id)
        if annot is None:
            continue
        if cohort_maf[idx] >= maf_threshold:
            continue

        terms = consequence_terms(annot["Consequence"])
        if annot["BIOTYPE"] != "protein_coding" or not (terms & ALLOWED_CONSEQUENCES):
            continue

        anchors.append(
            {
                "variant_index": idx,
                "variant_id": variant_id,
                "pos": int(positions[idx]),
                "gene": annot["SYMBOL"] if annot["SYMBOL"] not in {"", "-"} else annot["Gene"],
                "consequence": annot["Consequence"],
                "canonical": annot["CANONICAL"],
                "maf": float(cohort_maf[idx]),
            }
        )

    return anchors


def candidate_sort_key(row: dict[str, object]) -> tuple[int, int, float, int]:
    return (
        -int(row["anchor_count"]),
        -int(row["n_variants"]),
        float(row["center_maf"]),
        int(row["center_pos"]),
    )


def candidate_preference_key(row: dict[str, object]) -> tuple[int, int, float, int]:
    return (
        int(row["anchor_count"]),
        int(row["n_variants"]),
        -float(row["center_maf"]),
        -int(row["center_pos"]),
    )


def build_candidate_windows(
    anchors: list[dict[str, object]],
    positions: np.ndarray,
    margin_bp: int,
    min_anchor_count: int,
    min_window_variants: int,
) -> list[dict[str, object]]:
    anchor_variant_indices = [int(anchor["variant_index"]) for anchor in anchors]
    dedup: dict[tuple[int, int], dict[str, object]] = {}

    for anchor in anchors:
        center_pos = int(anchor["pos"])
        window_start_bp = max(1, center_pos - margin_bp)
        window_end_bp = center_pos + margin_bp
        start_index = bisect_left(positions, window_start_bp)
        end_index = bisect_right(positions, window_end_bp)
        n_variants = end_index - start_index
        if n_variants <= 0:
            continue

        anchor_left = bisect_left(anchor_variant_indices, start_index)
        anchor_right = bisect_left(anchor_variant_indices, end_index)
        anchor_count = anchor_right - anchor_left
        if anchor_count < min_anchor_count or n_variants < min_window_variants:
            continue

        key = (start_index, end_index)
        candidate = {
            "center_variant_id": anchor["variant_id"],
            "center_variant_index": int(anchor["variant_index"]),
            "center_gene": anchor["gene"],
            "center_pos": center_pos,
            "center_consequence": anchor["consequence"],
            "center_maf": float(anchor["maf"]),
            "window_start_bp": window_start_bp,
            "window_end_bp": window_end_bp,
            "start_index": start_index,
            "end_index": end_index,
            "n_variants": n_variants,
            "anchor_count": anchor_count,
            "retained_start_bp": int(positions[start_index]),
            "retained_end_bp": int(positions[end_index - 1]),
        }

        prev = dedup.get(key)
        if prev is None:
            dedup[key] = candidate
            continue

        prev_sort = candidate_preference_key(prev)
        cand_sort = candidate_preference_key(candidate)
        if cand_sort > prev_sort:
            dedup[key] = candidate

    candidates = list(dedup.values())
    candidates.sort(key=candidate_sort_key)
    return candidates


def finalize_windows(selected: list[dict[str, object]]) -> list[dict[str, object]]:
    selected = [dict(row) for row in selected]
    selected.sort(key=lambda row: int(row["start_index"]))
    for block_id, row in enumerate(selected, start=1):
        row["block_id"] = block_id
    return selected


def select_greedy_windows(
    candidates: list[dict[str, object]],
    max_windows: int,
) -> list[dict[str, object]]:
    selected = []
    for candidate in candidates:
        overlaps_existing = any(
            not (
                int(candidate["end_index"]) <= int(chosen["start_index"])
                or int(candidate["start_index"]) >= int(chosen["end_index"])
            )
            for chosen in selected
        )
        if overlaps_existing:
            continue

        selected.append(candidate)
        if max_windows > 0 and len(selected) >= max_windows:
            break

    return finalize_windows(selected)


def merge_overlapping_windows(
    candidates: list[dict[str, object]],
    anchors: list[dict[str, object]],
    positions: np.ndarray,
    max_windows: int,
) -> list[dict[str, object]]:
    if not candidates:
        return []

    anchor_variant_indices = sorted(int(anchor["variant_index"]) for anchor in anchors)
    by_start = sorted(candidates, key=lambda row: (int(row["start_index"]), int(row["end_index"])))
    merged = []

    current_cluster = [by_start[0]]
    current_start = int(by_start[0]["start_index"])
    current_end = int(by_start[0]["end_index"])

    for candidate in by_start[1:]:
        start_index = int(candidate["start_index"])
        end_index = int(candidate["end_index"])
        if start_index < current_end:
            current_cluster.append(candidate)
            current_end = max(current_end, end_index)
            continue

        merged.append(
            build_merged_window(current_cluster, anchor_variant_indices, positions, current_start, current_end)
        )
        current_cluster = [candidate]
        current_start = start_index
        current_end = end_index

    merged.append(
        build_merged_window(current_cluster, anchor_variant_indices, positions, current_start, current_end)
    )

    merged.sort(key=candidate_sort_key)
    if max_windows > 0:
        merged = merged[:max_windows]

    return finalize_windows(merged)


def build_merged_window(
    cluster: list[dict[str, object]],
    anchor_variant_indices: list[int],
    positions: np.ndarray,
    start_index: int,
    end_index: int,
) -> dict[str, object]:
    representative = max(cluster, key=candidate_preference_key)
    anchor_left = bisect_left(anchor_variant_indices, start_index)
    anchor_right = bisect_left(anchor_variant_indices, end_index)

    return {
        "center_variant_id": representative["center_variant_id"],
        "center_variant_index": representative["center_variant_index"],
        "center_gene": representative["center_gene"],
        "center_pos": representative["center_pos"],
        "center_consequence": representative["center_consequence"],
        "center_maf": representative["center_maf"],
        "window_start_bp": min(int(row["window_start_bp"]) for row in cluster),
        "window_end_bp": max(int(row["window_end_bp"]) for row in cluster),
        "start_index": start_index,
        "end_index": end_index,
        "n_variants": end_index - start_index,
        "anchor_count": anchor_right - anchor_left,
        "retained_start_bp": int(positions[start_index]),
        "retained_end_bp": int(positions[end_index - 1]),
        "merged_candidate_count": len(cluster),
        "selection_mode": "merged",
    }


def ensure_clean_dir(path: Path, force: bool) -> None:
    if path.exists():
        if not force:
            raise FileExistsError(f"{path} already exists; rerun with --force to overwrite it")
        if path.is_symlink() or path.is_file():
            path.unlink()
        else:
            shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def write_tsv(path: Path, header: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=header, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def materialize_keys_link(source_dataset: Path, out_dataset: Path) -> None:
    src = source_dataset / "keys"
    dst = out_dataset / "keys"
    if dst.exists() or dst.is_symlink():
        if dst.is_dir() and not dst.is_symlink():
            shutil.rmtree(dst)
        else:
            dst.unlink()
    try:
        dst.symlink_to(src.resolve())
    except OSError:
        shutil.copytree(src, dst)


def write_window_extract_lists(
    out_dataset: Path,
    windows: list[dict[str, object]],
    snp_ids: list[str],
) -> list[Path]:
    selection_dir = out_dataset / "selection"
    selection_dir.mkdir(parents=True, exist_ok=True)
    extract_paths = []

    for window in windows:
        block_id = int(window["block_id"])
        start_index = int(window["start_index"])
        end_index = int(window["end_index"])
        extract_path = selection_dir / f"window{block_id:02d}.snps.txt"
        with extract_path.open("w") as fh:
            fh.write("\n".join(snp_ids[start_index:end_index]))
            fh.write("\n")
        extract_paths.append(extract_path)

    return extract_paths


def subset_counts_for_windows(
    counts: np.ndarray,
    windows: list[dict[str, object]],
) -> np.ndarray:
    selected_indices = np.concatenate(
        [
            np.arange(int(window["start_index"]), int(window["end_index"]), dtype=np.int64)
            for window in windows
        ]
    )
    return counts[:, selected_indices]


def write_party_metadata(
    source_dataset: Path,
    out_dataset: Path,
    party_name: str,
    windows: list[dict[str, object]],
    snp_ids: list[str],
    positions: np.ndarray,
) -> None:
    source_party = source_dataset / party_name
    out_party = out_dataset / party_name
    (out_party / "geno").mkdir(parents=True, exist_ok=True)

    for filename in ("pheno.txt", "cov.txt", "sample_keep.txt"):
        shutil.copy2(source_party / filename, out_party / filename)

    chrom_sizes = []
    selected_ids = []
    selected_pos = []
    for window in windows:
        start_index = int(window["start_index"])
        end_index = int(window["end_index"])
        chrom_sizes.append(str(end_index - start_index))
        selected_ids.extend(snp_ids[start_index:end_index])
        selected_pos.extend(f"22\t{pos}\n" for pos in positions[start_index:end_index])

    (out_party / "chrom_sizes.txt").write_text("\n".join(chrom_sizes) + "\n")
    (out_party / "snp_ids.txt").write_text("\n".join(selected_ids) + "\n")
    (out_party / "snp_pos.txt").write_text("".join(selected_pos))

    original_counts = load_party_gcount(
        source_party / "all.gcount.transpose.bin",
        len(snp_ids),
    )
    subset_counts = subset_counts_for_windows(original_counts, windows)
    subset_counts.astype(np.uint32, copy=False).tofile(out_party / "all.gcount.transpose.bin")


def pfile_has_zst(prefix: Path) -> bool:
    return prefix.with_suffix(".pvar.zst").exists()


def run_plink_extract(
    plink2: str,
    source_prefix: Path,
    extract_path: Path,
    out_prefix: Path,
) -> None:
    cmd = [plink2, "--pfile", str(source_prefix)]
    if pfile_has_zst(source_prefix):
        cmd.append("vzs")
    cmd.extend(
        [
            "--extract",
            str(extract_path),
            "--make-pgen",
            "vzs",
            "--out",
            str(out_prefix),
        ]
    )
    subprocess.run(cmd, check=True)


def materialize_pgen_blocks(
    plink2: str,
    source_dataset: Path,
    out_dataset: Path,
    extract_paths: list[Path],
) -> None:
    for party_name in ("party1", "party2"):
        source_prefix = source_dataset / party_name / "geno" / "chr1"
        out_geno_dir = out_dataset / party_name / "geno"
        for block_id, extract_path in enumerate(extract_paths, start=1):
            out_prefix = out_geno_dir / f"chr{block_id}"
            print(f"[{party_name}] building block {block_id:02d} with plink2")
            run_plink_extract(plink2, source_prefix, extract_path, out_prefix)


def replace_toml_assignment(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"(?m)^({re.escape(key)}\s*=\s*).*$")
    updated, count = pattern.subn(lambda match: f"{match.group(1)}{value}", text, count=1)
    if count != 1:
        raise ValueError(f"Could not patch TOML key '{key}'")
    return updated


def emit_config_dir(
    template_dir: Path,
    out_dir: Path,
    windows: list[dict[str, object]],
    force: bool,
) -> None:
    if out_dir.exists():
        if not force:
            raise FileExistsError(f"{out_dir} already exists; rerun with --force to overwrite it")
        shutil.rmtree(out_dir)

    shutil.copytree(template_dir, out_dir)

    total_variants = sum(int(window["n_variants"]) for window in windows)
    num_blocks = len(windows)

    global_path = out_dir / "configGlobal.toml"
    global_text = global_path.read_text()
    global_text = replace_toml_assignment(global_text, "num_snps", str(total_variants))
    global_text = replace_toml_assignment(global_text, "geno_num_blocks", str(num_blocks))
    global_path.write_text(global_text)

    for party_id in (1, 2):
        local_path = out_dir / f"configLocal.Party{party_id}.toml"
        local_text = local_path.read_text()
        local_text = replace_toml_assignment(local_text, "geno_num_blocks", str(num_blocks))
        local_path.write_text(local_text)


def build_selection_outputs(
    out_dataset: Path,
    anchors: list[dict[str, object]],
    windows: list[dict[str, object]],
    snp_ids: list[str],
) -> None:
    selection_dir = out_dataset / "selection"
    selection_dir.mkdir(parents=True, exist_ok=True)

    selected_block_by_center = {row["center_variant_id"]: int(row["block_id"]) for row in windows}
    anchor_rows = []
    for anchor in anchors:
        anchor_rows.append(
            {
                "variant_id": anchor["variant_id"],
                "pos": int(anchor["pos"]),
                "gene": anchor["gene"],
                "consequence": anchor["consequence"],
                "canonical": anchor["canonical"],
                "cohort_maf": f"{float(anchor['maf']):.8f}",
                "selected_center_block": selected_block_by_center.get(anchor["variant_id"], ""),
            }
        )

    write_tsv(
        selection_dir / "anchor_variants.tsv",
        ["variant_id", "pos", "gene", "consequence", "canonical", "cohort_maf", "selected_center_block"],
        anchor_rows,
    )

    window_rows = []
    for window in windows:
        start_index = int(window["start_index"])
        end_index = int(window["end_index"])
        window_rows.append(
            {
                "block_id": int(window["block_id"]),
                "center_variant_id": window["center_variant_id"],
                "center_gene": window["center_gene"],
                "center_consequence": window["center_consequence"],
                "center_maf": f"{float(window['center_maf']):.8f}",
                "center_pos": int(window["center_pos"]),
                "window_start_bp": int(window["window_start_bp"]),
                "window_end_bp": int(window["window_end_bp"]),
                "retained_start_bp": int(window["retained_start_bp"]),
                "retained_end_bp": int(window["retained_end_bp"]),
                "anchor_count": int(window["anchor_count"]),
                "n_variants": int(window["n_variants"]),
                "selection_mode": window.get("selection_mode", "greedy"),
                "merged_candidate_count": int(window.get("merged_candidate_count", 1)),
                "start_variant_id": snp_ids[start_index],
                "end_variant_id": snp_ids[end_index - 1],
            }
        )

    write_tsv(
        selection_dir / "selected_windows.tsv",
        [
            "block_id",
            "center_variant_id",
            "center_gene",
            "center_consequence",
            "center_maf",
            "center_pos",
            "window_start_bp",
            "window_end_bp",
            "retained_start_bp",
            "retained_end_bp",
            "anchor_count",
            "n_variants",
            "selection_mode",
            "merged_candidate_count",
            "start_variant_id",
            "end_variant_id",
        ],
        window_rows,
    )


def default_config_out_dir(out_dataset: Path) -> Path:
    return out_dataset.parent.parent / f"config_{out_dataset.name}" if out_dataset.parent.name == "datasets" else out_dataset.parent / f"config_{out_dataset.name}"


def main() -> int:
    args = parse_args()

    source_dataset = Path(args.source_dataset).resolve()
    vep_path = Path(args.vep).resolve()
    out_dataset = Path(args.out_dataset).resolve()

    if not source_dataset.is_dir():
        raise FileNotFoundError(f"Missing source dataset: {source_dataset}")
    if not vep_path.is_file():
        raise FileNotFoundError(f"Missing VEP table: {vep_path}")

    print(f"Loading source metadata from {source_dataset}")
    snp_ids, positions = read_snp_metadata(source_dataset)
    cohort_maf = compute_combined_maf(source_dataset, len(snp_ids))

    print(f"Loading VEP annotations from {vep_path}")
    annotations = load_best_vep_annotations(vep_path, set(snp_ids))

    print("Selecting rare coding/splice anchor variants")
    anchors = build_anchor_rows(
        snp_ids=snp_ids,
        positions=positions,
        cohort_maf=cohort_maf,
        annotations=annotations,
        maf_threshold=args.maf_threshold,
    )
    if not anchors:
        raise RuntimeError("No rare coding/splice anchors were found with the requested filters")

    print(
        "Anchor variants retained:",
        len(anchors),
        f"(MAF < {args.maf_threshold}, protein_coding, allowed consequence)",
    )

    candidates = build_candidate_windows(
        anchors=anchors,
        positions=positions,
        margin_bp=args.margin_bp,
        min_anchor_count=args.min_anchor_count,
        min_window_variants=args.min_window_variants,
    )
    if not candidates:
        raise RuntimeError(
            "No windows met the requested anchor-count/variant-count thresholds; "
            "try lowering --min-anchor-count or --min-window-variants"
        )

    if args.selection_mode == "greedy":
        windows = select_greedy_windows(candidates, args.max_windows)
    else:
        windows = merge_overlapping_windows(candidates, anchors, positions, args.max_windows)

    if not windows:
        raise RuntimeError("Selection produced zero windows/loci after applying the requested mode")

    selection_label = "non-overlapping windows" if args.selection_mode == "greedy" else "merged loci"
    print(
        f"Selected {len(windows)} {selection_label} "
        f"from {len(candidates)} candidate windows (mode={args.selection_mode})"
    )
    for window in windows:
        print(
            f"  block {int(window['block_id']):02d}: "
            f"{int(window['retained_start_bp'])}-{int(window['retained_end_bp'])} "
            f"n_variants={int(window['n_variants'])} "
            f"anchors={int(window['anchor_count'])} "
            f"center={window['center_gene']}:{window['center_variant_id']}"
            f" merged_candidates={int(window.get('merged_candidate_count', 1))}"
        )

    ensure_clean_dir(out_dataset, args.force)
    materialize_keys_link(source_dataset, out_dataset)
    build_selection_outputs(out_dataset, anchors, windows, snp_ids)
    extract_paths = write_window_extract_lists(out_dataset, windows, snp_ids)

    for party_name in ("party1", "party2"):
        write_party_metadata(
            source_dataset=source_dataset,
            out_dataset=out_dataset,
            party_name=party_name,
            windows=windows,
            snp_ids=snp_ids,
            positions=positions,
        )

    materialize_pgen_blocks(args.plink2, source_dataset, out_dataset, extract_paths)

    if not args.skip_config:
        template_dir = Path(args.config_template_dir).resolve()
        if not template_dir.is_dir():
            raise FileNotFoundError(f"Missing config template directory: {template_dir}")
        config_out_dir = Path(args.config_out_dir).resolve() if args.config_out_dir else default_config_out_dir(out_dataset)
        emit_config_dir(template_dir, config_out_dir, windows, args.force)
        print(f"Config directory written to {config_out_dir}")

    total_variants = sum(int(window["n_variants"]) for window in windows)
    print(f"Reduced dataset written to {out_dataset}")
    print(f"Total selected variants: {total_variants}")
    print(f"Selected windows metadata: {out_dataset / 'selection' / 'selected_windows.tsv'}")
    print(f"Anchor table: {out_dataset / 'selection' / 'anchor_variants.tsv'}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
