#!/usr/bin/env python3
"""Join R v9 and All-by-All v8; retain missing results and separate diagnostics."""

import argparse
from collections import Counter
import json
import math
from pathlib import Path

from common import DEFAULT_OUTPUT, gene_id, index_rows, p_text, read_rows, settings, valid_p, write_rows

FIELDS = ["gene_id", "gene_symbol", "phenotype", "axa_phenotype_id",
          "R_burden", "R_skat", "axa_burden", "axa_skat"]
QC_FIELDS = ("chromosome n_samples n_variants_input n_variants_burden n_variants_skat "
             "R_burden_status R_skat_status R_skat_converged R_burden_message "
             "R_skat_message R_version SKAT_version").split()


def external_index(rows, config):
    phenotype_ids = {p["axa_id"]: p["name"] for p in config["phenotypes"]}
    for row in rows:
        if (row["ancestry"] != config["ancestry"] or row["release"] != "v8" or
                row["annotation"] != config["annotation"] or
                float(row["max_MAF"]) != config["max_maf"]):
            raise ValueError("All-by-All export does not match ancestry/release/mask/MAF")
        if phenotype_ids.get(row["axa_phenotype_id"]) != row["phenotype"]:
            raise ValueError("All-by-All phenotype ID/name mismatch")
    if {r["axa_phenotype_id"] for r in rows} != phenotype_ids.keys():
        raise ValueError("All-by-All export must contain all five phenotypes")
    return index_rows(rows, lambda r: (gene_id(r["gene_id"]), r["axa_phenotype_id"]))


def difference(first, second):
    a, b = valid_p(first), valid_p(second)
    if a is None or b is None:
        return "NA", "NA"
    log_difference = abs(math.log10(a) - math.log10(b)) if a > 0 and b > 0 else "NA"
    return abs(a - b), log_difference


def compare(r_rows, axa_rows, python_rows, config):
    axa = external_index(axa_rows, config)
    python = index_rows(python_rows, lambda r: (gene_id(r["gene_id"]), r["phenotype_name"]))
    r_index = index_rows(r_rows, lambda r: (gene_id(r["gene_id"]), r["axa_phenotype_id"]))
    if not r_rows:
        raise ValueError("R result is empty")
    phenotype_columns = {p["axa_id"]: (p["name"], p["column"]) for p in config["phenotypes"]}
    comparison, qc = [], []
    for (stable, phenotype_id), local in r_index.items():
        if phenotype_columns.get(phenotype_id) != (local["phenotype"], local["phenotype_column"]):
            raise ValueError("R phenotype name/column/ID mismatch")
        external = axa.get((stable, phenotype_id))
        row = {**{key: local[key] for key in FIELDS[:4]}, "gene_id": stable}
        if external and external["gene_symbol"]:
            row["gene_symbol"] = external["gene_symbol"]
        diagnostic = {**row, **{key: local[key] for key in QC_FIELDS}}
        original = python.get((stable, local["phenotype_column"]), {})
        diagnostic["python_match"] = "found" if original else "missing_or_not_supplied"
        for test, python_column in (("burden", "r_burden_p"), ("skat", "r_skat_davies_p")):
            row[f"R_{test}"] = p_text(local[f"R_{test}"])
            row[f"axa_{test}"] = p_text(external[f"axa_{test}"] if external else None)
            status = "ok"
            if external is None:
                status = "missing_row"
            elif valid_p(external[f"axa_{test}"]) is None:
                status = "missing_or_invalid_p"
            diagnostic[f"axa_{test}_status"] = status
            py_p = original.get(python_column)
            diagnostic[f"python_{test}"] = p_text(py_p)
            absolute, log = difference(local[f"R_{test}"], py_p)
            diagnostic[f"python_R_{test}_abs_p_diff"] = absolute
            diagnostic[f"python_R_{test}_abs_log10_diff"] = log
        diagnostic["python_skat_converged"] = original.get("r_skat_davies_converged", "NA")
        comparison.append(row)
        qc.append(diagnostic)
    for table in (comparison, qc):
        table.sort(key=lambda r: (r["gene_id"], r["phenotype"]))
    summary = {"rows": len(comparison), "python_rows_supplied": len(python_rows),
               "note": "Cross-release concordance; R/Python equivalence is checked, not assumed.",
               "phenotypes": {}}
    for pheno in config["phenotypes"]:
        selected = [r for r in qc if r["axa_phenotype_id"] == pheno["axa_id"]]
        details = {"rows": len(selected), "allxall_genes_outside_R": sum(
            key[1] == pheno["axa_id"] and key not in r_index for key in axa)}
        for test in ("burden", "skat"):
            for source in ("R", "axa"):
                key = f"{source}_{test}_status"
                details[key] = dict(Counter(r[key] for r in selected))
            differences = [r[f"python_R_{test}_abs_log10_diff"] for r in selected
                           if r[f"python_R_{test}_abs_log10_diff"] != "NA"]
            details[f"max_python_R_{test}_abs_log10_diff"] = max(differences, default=None)
        summary["phenotypes"][pheno["name"]] = details
    return comparison, qc, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r-dir", type=Path, default=DEFAULT_OUTPUT / "pilot")
    parser.add_argument("--allxall", type=Path, default=DEFAULT_OUTPUT / "allxall_v8.tsv")
    parser.add_argument("--python-results", type=Path,
                        default=Path.home() / "allxall-comparison" / "all_r_results.csv")
    args = parser.parse_args()
    directory = args.r_dir.expanduser().resolve()
    if not (directory / "_SUCCESS").is_file():
        raise ValueError("R run is not complete: missing _SUCCESS")
    config = settings()
    manifest = json.loads((directory / "run_info.json").read_text())
    if manifest["settings"] != config:
        raise ValueError("R run settings differ from current settings.json")
    python_path = args.python_results.expanduser()
    python_rows = read_rows(python_path) if python_path.is_file() else []
    rows, qc, summary = compare(read_rows(directory / "r_v9.csv"),
                                read_rows(args.allxall.expanduser(), "\t"), python_rows, config)
    summary["allxall_file"] = str(args.allxall.expanduser().resolve())
    summary["python_file"] = str(python_path.resolve()) if python_rows else None
    write_rows(directory / "comparison.csv", rows, FIELDS)
    write_rows(directory / "comparison_qc.csv", qc)
    (directory / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"Saved: {directory / 'comparison.csv'}")


if __name__ == "__main__":
    main()
