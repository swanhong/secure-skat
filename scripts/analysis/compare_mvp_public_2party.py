#!/usr/bin/env python3
"""Compare AoU/MVP MVP-public secure outputs against a full-union plaintext reference."""

from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def resolve_path(path_text: str | Path) -> Path:
    path = Path(path_text)
    if path.is_absolute():
        return path
    return (repo_root() / path).resolve()


def require_executable(name: str) -> str:
    found = shutil.which(name)
    if found is None:
        raise RuntimeError(f"Required executable not found on PATH: {name}")
    return found


def read_int_lines(path: Path) -> list[int]:
    return [int(line.strip()) for line in path.read_text().splitlines() if line.strip()]


def read_vector(path: Path) -> np.ndarray:
    values = np.loadtxt(path, dtype=float)
    return np.asarray(values, dtype=float).reshape(-1)


def read_matrix(path: Path) -> np.ndarray:
    values = np.loadtxt(path, dtype=float)
    arr = np.asarray(values, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    return arr


def parse_gt_to_alt_dosage(field: str, gt_index: int) -> float:
    pieces = field.split(":")
    if gt_index >= len(pieces):
        return float("nan")
    gt = pieces[gt_index]
    if "." in gt:
        return float("nan")
    return float(sum(1 for allele in gt.replace("|", "/").split("/") if allele == "1"))


def export_alt_matrix(plink2: str, party_dir: Path, block_index: int, scratch_dir: Path) -> tuple[list[str], np.ndarray]:
    out_prefix = scratch_dir / f"{party_dir.name}_block{block_index:02d}"
    vcf_path = out_prefix.with_suffix(".vcf")
    if not vcf_path.exists():
        pfile_prefix = party_dir / "geno" / f"chr{block_index}"
        cmd = [plink2, "--pfile", str(pfile_prefix)]
        if Path(str(pfile_prefix) + ".pvar.zst").exists():
            cmd.append("vzs")
        cmd.extend(["--keep", str(party_dir / "sample_keep.txt"), "--export", "vcf", "--out", str(out_prefix)])
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"plink2 export failed for {party_dir.name} block {block_index}\n{proc.stderr.strip()}")

    variant_ids: list[str] = []
    rows: list[list[float]] = []
    with vcf_path.open() as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "#CHROM":
                continue
            if len(parts) < 10:
                raise RuntimeError(f"Malformed VCF row in {vcf_path}")
            fmt = parts[8].split(":")
            gt_index = fmt.index("GT")
            variant_ids.append(parts[2])
            rows.append([parse_gt_to_alt_dosage(field, gt_index) for field in parts[9:]])
    if not rows:
        raise RuntimeError(f"No genotype rows found in {vcf_path}")
    return variant_ids, np.asarray(rows, dtype=float).T


def resolve_fixture_root(dataset_root: Path) -> Path:
    if (dataset_root / "reference" / "full_union_reference.json").exists():
        return dataset_root
    if dataset_root.name == "dataset" and (dataset_root.parent / "reference" / "full_union_reference.json").exists():
        return dataset_root.parent
    return dataset_root


def secure_oriented(alt_matrix: np.ndarray) -> np.ndarray:
    return np.where(np.isfinite(alt_matrix), 2.0 - alt_matrix, 0.0)


def beta_weight(secure_matrix: np.ndarray) -> np.ndarray:
    n_total = secure_matrix.shape[0]
    p_bar = secure_matrix.sum(axis=0) / (2.0 * n_total)
    p = 1.0 - p_bar
    return 25.0 * np.power(np.maximum(p, p_bar), 24)


def block_stats(secure_parts: list[np.ndarray], residual: np.ndarray) -> dict[str, np.ndarray | float]:
    G = np.vstack(secure_parts)
    weights = beta_weight(G)
    score = G.T @ residual
    return {
        "weights": weights,
        "score": score,
        "q_skat": float(np.sum((weights * weights) * (score * score))),
        "burden_linear": float(np.sum(weights * score)),
    }


def fit_null(pheno_a: np.ndarray, cov_a: np.ndarray, pheno_m: np.ndarray, cov_m: np.ndarray) -> dict[str, float | np.ndarray]:
    y = np.concatenate([pheno_a, pheno_m])
    X = np.vstack([cov_a, cov_m])
    design = np.column_stack([np.ones(y.shape[0]), X])
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    residual = y - design @ beta
    rss = float(np.sum(residual * residual))
    dof = int(y.shape[0] - design.shape[1])
    if dof <= 0 or rss <= 0.0:
        raise RuntimeError("Invalid null model scale")
    return {
        "residual": residual,
        "rss": rss,
        "dof": dof,
        "alpha": float(dof / (2.0 * rss)),
    }


def read_secure_scalar(path: Path) -> float | None:
    if not path.exists():
        return None
    return float(np.asarray(np.loadtxt(path, dtype=float, max_rows=1)).reshape(-1)[0])


def rel_abs_diff(observed: float | None, expected: float) -> dict[str, float | None]:
    if observed is None:
        return {"abs": None, "rel": None}
    abs_diff = abs(observed - expected)
    rel_diff = abs_diff / max(1.0, abs(expected))
    return {"abs": float(abs_diff), "rel": float(rel_diff)}


def is_close(observed: float | None, expected: float, rtol: float, atol: float) -> bool:
    if observed is None:
        return False
    return math.isclose(observed, expected, rel_tol=rtol, abs_tol=atol)


def build_compare(dataset_root: Path, run_root: Path, plink2: str, mode: str, skato_rho: float, rtol: float, atol: float) -> dict:
    fixture_root = resolve_fixture_root(dataset_root)
    dataset_dir = dataset_root / "dataset" if (dataset_root / "dataset").exists() else dataset_root
    party1_dir = dataset_dir / "party1"
    party2_dir = dataset_dir / "party2"
    hidden_dir = dataset_dir / "party1_hidden"
    if not party1_dir.exists() or not party2_dir.exists() or not hidden_dir.exists():
        raise RuntimeError(f"Expected party1, party2, and party1_hidden under {dataset_dir}")

    public_sizes = read_int_lines(party1_dir / "chrom_sizes.txt")
    hidden_sizes = read_int_lines(hidden_dir / "chrom_sizes.txt")
    if len(public_sizes) != len(hidden_sizes):
        raise RuntimeError("This fixture compare expects public and hidden block-count metadata to be paired")

    analysis_dir = run_root / "analysis_2party"
    scratch_dir = analysis_dir / "scratch"
    scratch_dir.mkdir(parents=True, exist_ok=True)

    pheno_a = read_vector(party1_dir / "pheno.txt")
    pheno_m = read_vector(party2_dir / "pheno.txt")
    cov_a = read_matrix(party1_dir / "cov.txt")
    cov_m = read_matrix(party2_dir / "cov.txt")
    model = fit_null(pheno_a, cov_a, pheno_m, cov_m)
    residual = np.asarray(model["residual"], dtype=float)
    n_a = pheno_a.shape[0]
    n_m = pheno_m.shape[0]
    r_a = residual[:n_a]

    total_full_q = 0.0
    total_full_b = 0.0
    total_decomp_q = 0.0
    total_decomp_b = 0.0
    max_hidden_score_diff = 0.0
    block_rows = []

    for block_index in range(1, len(public_sizes) + 1):
        public_a_ids, public_a_alt = export_alt_matrix(plink2, party1_dir, block_index, scratch_dir)
        public_m_ids, public_m_alt = export_alt_matrix(plink2, party2_dir, block_index, scratch_dir)
        hidden_ids, hidden_a_alt = export_alt_matrix(plink2, hidden_dir, block_index, scratch_dir)
        if public_a_ids != public_m_ids:
            raise RuntimeError(f"Public variant order mismatch in block {block_index}")

        public_a_secure = secure_oriented(public_a_alt)
        public_m_secure = secure_oriented(public_m_alt)
        hidden_a_secure = secure_oriented(hidden_a_alt)
        hidden_m_secure = np.full((n_m, hidden_a_secure.shape[1]), 2.0, dtype=float)

        public = block_stats([public_a_secure, public_m_secure], residual)
        hidden = block_stats([hidden_a_secure, hidden_m_secure], residual)
        full = block_stats(
            [
                np.concatenate([public_a_secure, hidden_a_secure], axis=1),
                np.concatenate([public_m_secure, hidden_m_secure], axis=1),
            ],
            residual,
        )

        hidden_score_corr = hidden_a_secure.T @ r_a - 2.0 * float(np.sum(r_a))
        hidden_score_full = np.vstack([hidden_a_secure, hidden_m_secure]).T @ residual
        hidden_score_diff = float(np.max(np.abs(hidden_score_corr - hidden_score_full))) if hidden_score_full.size else 0.0
        max_hidden_score_diff = max(max_hidden_score_diff, hidden_score_diff)

        q_decomp = float(public["q_skat"]) + float(hidden["q_skat"])
        b_decomp = float(public["burden_linear"]) + float(hidden["burden_linear"])
        total_full_q += float(full["q_skat"])
        total_full_b += float(full["burden_linear"])
        total_decomp_q += q_decomp
        total_decomp_b += b_decomp

        block_rows.append(
            {
                "block_id": block_index,
                "public_variants": len(public_a_ids),
                "hidden_variants": len(hidden_ids),
                "q_public": float(public["q_skat"]),
                "q_hidden": float(hidden["q_skat"]),
                "q_full": float(full["q_skat"]),
                "q_decomposed": q_decomp,
                "burden_linear_public": float(public["burden_linear"]),
                "burden_linear_hidden": float(hidden["burden_linear"]),
                "burden_linear_full": float(full["burden_linear"]),
                "burden_linear_decomposed": b_decomp,
                "max_hidden_score_correction_abs_diff": hidden_score_diff,
            }
        )

    alpha = float(model["alpha"])
    expected_skat_scaled = total_full_q * alpha
    expected_burden_scaled = (total_full_b * total_full_b) * alpha
    expected_skato_scaled = (1.0 - skato_rho) * expected_skat_scaled + skato_rho * expected_burden_scaled

    secure_party_dir = run_root / "party1"
    secure_skat = read_secure_scalar(secure_party_dir / "skat_out.txt")
    secure_burden = read_secure_scalar(secure_party_dir / "burden_out.txt")
    secure_skato = read_secure_scalar(secure_party_dir / "skato_out.txt")

    raw_q_abs = abs(total_full_q - total_decomp_q)
    raw_b_abs = abs(total_full_b - total_decomp_b)
    hidden_q_total = float(sum(float(row["q_hidden"]) for row in block_rows))
    hidden_nonzero_ok = all(row["hidden_variants"] == 0 or abs(float(row["q_hidden"])) > 1e-8 for row in block_rows)
    raw_ok = (
        raw_q_abs <= 1e-7 * max(1.0, abs(total_full_q))
        and raw_b_abs <= 1e-7 * max(1.0, abs(total_full_b))
        and max_hidden_score_diff <= 1e-7
        and hidden_nonzero_ok
    )
    reference = None
    reference_path = fixture_root / "reference" / "full_union_reference.json"
    if reference_path.exists():
        reference_payload = json.loads(reference_path.read_text())
        reference_q = float(reference_payload["q_full_total"])
        reference_b = float(reference_payload["burden_linear_full_total"])
        reference_bq = float(reference_payload["burden_q_full"])
        reference_q_diff = abs(total_full_q - reference_q)
        reference_b_diff = abs(total_full_b - reference_b)
        reference_bq_diff = abs((total_full_b * total_full_b) - reference_bq)
        reference_ok = (
            reference_q_diff <= 1e-7 * max(1.0, abs(reference_q))
            and reference_b_diff <= 1e-7 * max(1.0, abs(reference_b))
            and reference_bq_diff <= 1e-7 * max(1.0, abs(reference_bq))
        )
        raw_ok = raw_ok and reference_ok
        reference = {
            "path": str(reference_path),
            "q_full": reference_q,
            "q_abs_diff": reference_q_diff,
            "burden_linear_full": reference_b,
            "burden_linear_abs_diff": reference_b_diff,
            "burden_q_full": reference_bq,
            "burden_q_abs_diff": reference_bq_diff,
            "ok": reference_ok,
        }

    secure_skat_ok = is_close(secure_skat, expected_skat_scaled, rtol, atol)
    secure_burden_ok = is_close(secure_burden, expected_burden_scaled, rtol, atol)
    secure_skato_ok = is_close(secure_skato, expected_skato_scaled, rtol, atol)

    return {
        "dataset": str(dataset_root),
        "run_root": str(run_root),
        "mode": mode,
        "n_aou": int(n_a),
        "n_mvp": int(n_m),
        "alpha_exact_plain": alpha,
        "raw": {
            "q_full": total_full_q,
            "q_decomposed": total_decomp_q,
            "q_abs_diff": raw_q_abs,
            "burden_linear_full": total_full_b,
            "burden_linear_decomposed": total_decomp_b,
            "burden_linear_abs_diff": raw_b_abs,
            "hidden_q_total": hidden_q_total,
            "hidden_nonzero_ok": hidden_nonzero_ok,
            "max_hidden_score_correction_abs_diff": max_hidden_score_diff,
            "ok": raw_ok,
        },
        "reference": reference,
        "scaled": {
            "expected_skat": expected_skat_scaled,
            "secure_skat": secure_skat,
            "skat_diff": rel_abs_diff(secure_skat, expected_skat_scaled),
            "skat_ok": secure_skat_ok,
            "expected_burden": expected_burden_scaled,
            "secure_burden": secure_burden,
            "burden_diff": rel_abs_diff(secure_burden, expected_burden_scaled),
            "burden_ok": secure_burden_ok,
            "expected_skato": expected_skato_scaled,
            "secure_skato": secure_skato,
            "skato_diff": rel_abs_diff(secure_skato, expected_skato_scaled),
            "skato_ok": secure_skato_ok,
            "skato_rho": skato_rho,
            "rtol": rtol,
            "atol": atol,
        },
        "blocks": block_rows,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", required=True, help="Fixture root, e.g. datasets/1000g_all_chr22_mvp_public_2party")
    parser.add_argument("--run-root", required=True, help="Run output root containing party1/skat_out.txt")
    parser.add_argument("--plink2", default="plink2")
    parser.add_argument("--mode", choices=["skat", "burden", "skato"], default="skat")
    parser.add_argument("--skato-rho", type=float, default=0.5)
    parser.add_argument("--rtol", type=float, default=0.01)
    parser.add_argument("--atol", type=float, default=1e-3)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    plink2 = require_executable(args.plink2)
    result = build_compare(resolve_path(args.dataset), resolve_path(args.run_root), plink2, args.mode, args.skato_rho, args.rtol, args.atol)

    out_dir = Path(result["run_root"]) / "analysis_2party"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "compare_mvp_public_2party.json"
    out_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

    print("2-party MVP-public comparison")
    print(f"  raw full-vs-decomposed: {'PASS' if result['raw']['ok'] else 'FAIL'}")
    print(f"  raw hidden q total: {result['raw']['hidden_q_total']:.6e}")
    print(f"  hidden score correction max abs diff: {result['raw']['max_hidden_score_correction_abs_diff']:.6e}")
    if result["reference"] is not None:
        print(f"  reference JSON match: {'PASS' if result['reference']['ok'] else 'FAIL'}")
    if args.mode in {"skat", "skato"}:
        print(
            "  SKAT scaled secure/plain: "
            f"{result['scaled']['secure_skat']} vs {result['scaled']['expected_skat']:.6e} "
            f"(rel diff {result['scaled']['skat_diff']['rel']})"
        )
    print(
        "  Burden scaled secure/plain: "
        f"{result['scaled']['secure_burden']} vs {result['scaled']['expected_burden']:.6e} "
        f"(rel diff {result['scaled']['burden_diff']['rel']})"
    )
    if args.mode == "skato":
        print(
            "  SKAT-O scaled secure/plain: "
            f"{result['scaled']['secure_skato']} vs {result['scaled']['expected_skato']:.6e} "
            f"(rel diff {result['scaled']['skato_diff']['rel']})"
        )
    print(f"  JSON: {out_path}")

    if not result["raw"]["ok"]:
        return 1
    if args.mode in {"skat", "skato"} and not result["scaled"]["skat_ok"]:
        return 1
    if not result["scaled"]["burden_ok"]:
        return 1
    if args.mode == "skato" and not result["scaled"]["skato_ok"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
