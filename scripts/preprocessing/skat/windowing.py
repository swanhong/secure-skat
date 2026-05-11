from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from utils import parse_float, pfile_path


def normalize_chrom(value: object) -> int:
    text = str(value)
    text = text[3:] if text.lower().startswith("chr") else text
    return int(text)


def read_source_variants(source_prefix: Path) -> list[dict[str, object]]:
    pvar_path = pfile_path(source_prefix, ".pvar")
    freq_path = pfile_path(source_prefix, ".afreq")

    variants: list[dict[str, object]] = []
    with pvar_path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            toks = line.split()
            if len(toks) < 5:
                continue
            variants.append(
                {
                    "chrom": normalize_chrom(toks[0]),
                    "pos": int(toks[1]),
                    "id": toks[2],
                }
            )
    if not variants:
        raise ValueError(f"No variants found in {pvar_path}")

    maf_by_id: dict[str, float] = {}
    with freq_path.open() as fh:
        header = None
        for line in fh:
            if line.strip():
                header = line.split()
                break
        if header is None:
            raise ValueError(f"Empty AFREQ file: {freq_path}")
        if "ID" not in header or "ALT_FREQS" not in header:
            raise ValueError(f"Expected ID and ALT_FREQS columns in {freq_path}")
        id_idx = header.index("ID")
        af_idx = header.index("ALT_FREQS")
        for line in fh:
            if not line.strip():
                continue
            toks = line.split()
            if len(toks) <= max(id_idx, af_idx):
                continue
            af = parse_float(toks[af_idx])
            if math.isfinite(af):
                af = min(max(af, 0.0), 1.0)
                maf_by_id[toks[id_idx]] = min(af, 1.0 - af)

    missing = 0
    for row in variants:
        maf = maf_by_id.get(str(row["id"]))
        if maf is None:
            missing += 1
            row["maf"] = float("nan")
        else:
            row["maf"] = maf
    if missing:
        raise ValueError(f"PVAR/AFREQ variant ID mismatch for {missing} variants")
    return variants


def write_selected_windows(path: Path, args: argparse.Namespace, windows: list[dict[str, object]]) -> None:
    fields = [
        "block_id",
        "chromosome",
        "window_start_bp",
        "window_end_bp",
        "rare_variant_count",
        "materialized_variant_count",
        "block_mode",
        "maf_threshold",
        "extract_path",
    ]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in windows:
            writer.writerow(
                {
                    "block_id": row["block_id"],
                    "chromosome": args.chromosome,
                    "window_start_bp": row["bin_start"],
                    "window_end_bp": row["bin_end"],
                    "rare_variant_count": row["rare_count"],
                    "materialized_variant_count": row["variant_count"],
                    "block_mode": args.block_mode,
                    "maf_threshold": args.maf_threshold,
                    "extract_path": str(row["extract_path"]),
                }
            )


def select_windows(args: argparse.Namespace, variants: list[dict[str, object]], work_dir: Path, out_dataset: Path) -> list[dict[str, object]]:
    rare = [
        row
        for row in variants
        if math.isfinite(float(row["maf"])) and 0.0 < float(row["maf"]) <= args.maf_threshold
    ]
    if not rare:
        raise ValueError("No rare variants found; increase --n-samples or --maf-threshold")

    all_by_bin: dict[int, list[dict[str, object]]] = {}
    rare_by_bin: dict[int, list[dict[str, object]]] = {}
    for row in variants:
        bin_start = (int(row["pos"]) // args.window_bp) * args.window_bp
        row["bin_start"] = bin_start
        all_by_bin.setdefault(bin_start, []).append(row)
    for row in rare:
        rare_by_bin.setdefault(int(row["bin_start"]), []).append(row)

    candidate_rows = []
    for bin_start in sorted(rare_by_bin):
        rare_group = sorted(rare_by_bin[bin_start], key=lambda row: int(row["pos"]))
        if len(rare_group) < args.min_rare_per_window:
            continue
        if args.block_mode == "rare-only":
            block_variants = rare_group
        else:
            block_variants = sorted(all_by_bin[bin_start], key=lambda row: int(row["pos"]))
        candidate_rows.append(
            {
                "bin_start": int(bin_start),
                "bin_end": int(bin_start + args.window_bp - 1),
                "rare_count": int(len(rare_group)),
                "variant_count": int(len(block_variants)),
                "variants": block_variants,
            }
        )

    if not candidate_rows:
        raise ValueError(
            "No windows met --min-rare-per-window; lower the threshold or increase the sample count"
        )

    candidate_rows.sort(key=lambda row: (-int(row["rare_count"]), int(row["bin_start"])))
    if args.max_windows > 0:
        candidate_rows = candidate_rows[: args.max_windows]
    candidate_rows.sort(key=lambda row: int(row["bin_start"]))

    extract_dir = out_dataset / "selection"
    extract_dir.mkdir(parents=True, exist_ok=True)

    selected_ids: list[str] = []
    selected_pos: list[str] = []
    chrom_sizes: list[int] = []
    windows: list[dict[str, object]] = []

    for block_id, row in enumerate(candidate_rows, start=1):
        block_variants = row["variants"]
        extract_path = extract_dir / f"window{block_id:02d}.snps.txt"
        with extract_path.open("w") as fh:
            for variant in block_variants:
                fh.write(f"{variant['id']}\n")

        chrom_sizes.append(int(row["variant_count"]))
        selected_ids.extend(str(variant["id"]) for variant in block_variants)
        selected_pos.extend(f"{int(variant['chrom'])}\t{int(variant['pos'])}\n" for variant in block_variants)

        windows.append(
            {
                "block_id": block_id,
                "bin_start": row["bin_start"],
                "bin_end": row["bin_end"],
                "rare_count": row["rare_count"],
                "variant_count": row["variant_count"],
                "extract_path": extract_path,
            }
        )

    for party_name in ("party1", "party2"):
        party_dir = out_dataset / party_name
        (party_dir / "chrom_sizes.txt").write_text("\n".join(map(str, chrom_sizes)) + "\n")
        (party_dir / "snp_ids.txt").write_text("\n".join(selected_ids) + "\n")
        (party_dir / "snp_pos.txt").write_text("".join(selected_pos))

    selected_windows_path = extract_dir / "selected_windows.tsv"
    write_selected_windows(selected_windows_path, args, windows)

    print(f"Selected {len(windows)} windows from {len(candidate_rows)} retained candidates")
    print(f"Total materialized variants: {sum(chrom_sizes)}")
    print(f"Window metadata: {selected_windows_path}")
    return windows
