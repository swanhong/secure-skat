from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

from utils import parse_float, pfile_path, resolve_path, split_table_line


def split_cov_cols(cov_cols: str) -> list[str]:
    cols = [col.strip() for col in cov_cols.split(",") if col.strip()]
    if not cols:
        raise ValueError("--cov-cols must contain at least one column")
    return cols


def read_psam_samples(psam_path: Path) -> list[dict[str, object]]:
    with psam_path.open() as fh:
        header = None
        for line in fh:
            if line.strip():
                header = line.split()
                break
        if header is None:
            raise ValueError(f"Empty PSAM file: {psam_path}")

        iid_col = "#IID" if "#IID" in header else "IID"
        fid_col = "#FID" if "#FID" in header else ("FID" if "FID" in header else None)
        if iid_col not in header:
            raise ValueError(f"Could not find IID column in {psam_path}")
        iid_idx = header.index(iid_col)
        fid_idx = header.index(fid_col) if fid_col else None

        samples = []
        for order, line in enumerate(fh):
            if not line.strip():
                continue
            toks = line.split()
            iid = toks[iid_idx]
            fid = toks[fid_idx] if fid_idx is not None else iid
            samples.append(
                {
                    "_psam_order": order,
                    "FID_OUT": fid,
                    "IID_OUT": iid,
                }
            )
    if not samples:
        raise ValueError(f"No samples found in {psam_path}")
    return samples


def load_table_records(path: Path, sep: str | None, required_cols: list[str], id_col: str) -> dict[str, dict[str, str]]:
    with path.open(newline="") as fh:
        header = None
        for line in fh:
            if line.strip():
                header = split_table_line(line, sep)
                break
        if header is None:
            raise ValueError(f"Empty phenotype table: {path}")

        missing = [col for col in required_cols if col not in header]
        if missing:
            raise ValueError(f"Missing phenotype/covariate columns in {path}: {missing}")

        id_idx = header.index(id_col)
        required_indices = {col: header.index(col) for col in required_cols}
        records: dict[str, dict[str, str]] = {}
        for line in fh:
            if not line.strip():
                continue
            toks = split_table_line(line, sep)
            if len(toks) <= id_idx:
                continue
            sample_id = toks[id_idx]
            row = {}
            for col, idx in required_indices.items():
                row[col] = toks[idx] if idx < len(toks) else ""
            records[sample_id] = row
    return records


def write_keep_file(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w") as fh:
        for row in rows:
            fh.write(f"{row['FID_OUT']}\t{row['IID_OUT']}\n")


def write_pheno_file(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w") as fh:
        for row in rows:
            fh.write(f"{float(row['phenotype']):.12g}\n")


def write_cov_file(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w") as fh:
        for row in rows:
            cov = row["covariates"]
            fh.write("\t".join(f"{float(value):.12g}" for value in cov) + "\n")


def build_sample_files(
    args: argparse.Namespace,
    raw_prefix: Path,
    work_dir: Path,
    out_dataset: Path,
) -> tuple[int, int, int]:
    psam_samples = read_psam_samples(pfile_path(raw_prefix, ".psam"))

    if args.pheno_file:
        cov_input_cols = split_cov_cols(args.cov_cols)
        num_covs = len(cov_input_cols)

        pheno_path = resolve_path(args.pheno_file)
        required = [args.id_col, args.pheno_col, *cov_input_cols]
        pheno_by_id = load_table_records(pheno_path, args.pheno_sep, required, args.id_col)
        merged = []
        for sample in psam_samples:
            record = pheno_by_id.get(str(sample["IID_OUT"]))
            if record is None:
                continue
            phenotype = parse_float(record[args.pheno_col])
            covariates = [parse_float(record[col]) for col in cov_input_cols]
            if math.isfinite(phenotype) and all(math.isfinite(value) for value in covariates):
                merged.append({**sample, "phenotype": phenotype, "covariates": covariates})
    else:
        pheno_values = np.atleast_1d(np.loadtxt(resolve_path(args.pheno_vector_file), dtype=float))
        cov_values = np.loadtxt(resolve_path(args.cov_matrix_file), dtype=float, ndmin=2)
        if cov_values.shape[0] != len(psam_samples) and cov_values.shape[1] == len(psam_samples):
            cov_values = cov_values.T
        if pheno_values.shape[0] != len(psam_samples):
            raise ValueError(
                f"Phenotype vector has {pheno_values.shape[0]} rows but PSAM has {len(psam_samples)} samples"
            )
        if cov_values.shape[0] != len(psam_samples):
            raise ValueError(
                f"Covariate matrix has {cov_values.shape[0]} rows but PSAM has {len(psam_samples)} samples"
            )
        num_covs = int(cov_values.shape[1])
        merged = []
        for idx, sample in enumerate(psam_samples):
            phenotype = float(pheno_values[idx])
            covariates = [float(value) for value in cov_values[idx, :]]
            if math.isfinite(phenotype) and all(math.isfinite(value) for value in covariates):
                merged.append({**sample, "phenotype": phenotype, "covariates": covariates})

    if not merged:
        raise ValueError("No samples remain after applying phenotype/covariate filters")

    rng = np.random.default_rng(args.seed)
    if args.n_samples and len(merged) > args.n_samples:
        chosen = np.sort(rng.choice(np.arange(len(merged)), size=args.n_samples, replace=False))
        merged = [merged[int(idx)] for idx in chosen]

    party_labels = np.array(["party2"] * len(merged), dtype=object)
    n_party1 = max(1, int(round(len(merged) * args.party1_frac)))
    n_party1 = min(n_party1, len(merged) - 1)
    shuffled = np.arange(len(merged))
    rng.shuffle(shuffled)
    party_labels[shuffled[:n_party1]] = "party1"
    for row, party_label in zip(merged, party_labels):
        row["_party"] = party_label

    work_dir.mkdir(parents=True, exist_ok=True)
    all_keep = work_dir / "sample_keep_all.txt"
    write_keep_file(all_keep, merged)

    counts = {}
    for party_name in ("party1", "party2"):
        party_rows = sorted(
            [row for row in merged if row["_party"] == party_name],
            key=lambda row: int(row["_psam_order"]),
        )
        party_dir = out_dataset / party_name
        (party_dir / "geno").mkdir(parents=True, exist_ok=True)
        write_keep_file(party_dir / "sample_keep.txt", party_rows)
        write_pheno_file(party_dir / "pheno.txt", party_rows)
        write_cov_file(party_dir / "cov.txt", party_rows)
        counts[party_name] = len(party_rows)

    print(f"Samples retained: {len(merged)}")
    print(f"Party sample counts: party1={counts['party1']} party2={counts['party2']}")
    print(f"Covariate count: {num_covs}")
    return counts["party1"], counts["party2"], num_covs
