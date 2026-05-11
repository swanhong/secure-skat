#!/usr/bin/env python3

from __future__ import annotations

import argparse
import hashlib
import math
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create a global synthetic 1000G phenotype/covariate TSV.")
    parser.add_argument("--panel-file", required=True, help="1000G integrated_call_samples panel file.")
    parser.add_argument("--sample-file", help="Optional sample list/keep file. Uses all panel samples if omitted.")
    parser.add_argument("--out", required=True, help="Output TSV path.")
    parser.add_argument("--seed", default="270318_pheno", help="Deterministic phenotype seed.")
    parser.add_argument(
        "--cov-mode",
        choices=("superpop", "afr-pop"),
        default="superpop",
        help="Covariates to emit. superpop is for all 1000G; afr-pop matches the old AFR subset.",
    )
    return parser.parse_args()


def read_panel(path: Path) -> dict[str, dict[str, str]]:
    panel = {}
    with path.open() as fh:
        next(fh)
        for line in fh:
            if not line.strip():
                continue
            sample, pop, super_pop, gender, *_ = line.rstrip("\n").split("\t")
            panel[sample] = {"pop": pop, "super_pop": super_pop, "gender": gender}
    return panel


def read_sample_ids(path: Path | None, panel: dict[str, dict[str, str]]) -> list[str]:
    if path is None:
        return list(panel)

    sample_ids = []
    with path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            toks = line.split()
            sample_ids.append(toks[-1])
    return sample_ids


def hash_uniform(sample_id: str, seed: str) -> float:
    digest = hashlib.sha256(f"{sample_id}|{seed}".encode()).digest()
    return int.from_bytes(digest[:8], "big") / float(2**64)


def standardize(values: list[float]) -> list[float]:
    mean = sum(values) / len(values)
    var = sum((value - mean) ** 2 for value in values) / (len(values) - 1)
    sd = math.sqrt(var)
    return [(value - mean) / sd for value in values]


def cov_header(cov_mode: str) -> list[str]:
    if cov_mode == "superpop":
        return ["superpop_AFR", "superpop_AMR", "superpop_EAS", "superpop_EUR"]
    return ["sex_male", "pop_GWD", "pop_YRI", "pop_ESN"]


def cov_values(meta: dict[str, str], cov_mode: str) -> list[float]:
    if cov_mode == "superpop":
        return [1.0 if meta["super_pop"] == super_pop else 0.0 for super_pop in ("AFR", "AMR", "EAS", "EUR")]
    return [
        1.0 if meta["gender"] == "male" else 0.0,
        1.0 if meta["pop"] == "GWD" else 0.0,
        1.0 if meta["pop"] == "YRI" else 0.0,
        1.0 if meta["pop"] == "ESN" else 0.0,
    ]


def main() -> int:
    args = parse_args()
    panel = read_panel(Path(args.panel_file).expanduser().resolve())
    sample_ids = read_sample_ids(
        Path(args.sample_file).expanduser().resolve() if args.sample_file else None,
        panel,
    )

    missing = [sample_id for sample_id in sample_ids if sample_id not in panel]
    if missing:
        raise ValueError(f"{len(missing)} samples are missing from the panel; first missing sample: {missing[0]}")
    if len(sample_ids) < 2:
        raise ValueError("Need at least two samples to standardize the synthetic phenotype")

    pheno = standardize([hash_uniform(sample_id, args.seed) for sample_id in sample_ids])
    out = Path(args.out).expanduser().resolve()
    out.parent.mkdir(parents=True, exist_ok=True)

    with out.open("w") as fh:
        fh.write("\t".join(["sample_id", "pheno", *cov_header(args.cov_mode)]) + "\n")
        for sample_id, y in zip(sample_ids, pheno):
            cov = cov_values(panel[sample_id], args.cov_mode)
            row = [sample_id, f"{y:.12f}", *[f"{value:.1f}" for value in cov]]
            fh.write("\t".join(row) + "\n")

    print(f"Wrote {len(sample_ids)} samples to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
