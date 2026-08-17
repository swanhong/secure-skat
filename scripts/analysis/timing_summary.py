#!/usr/bin/env python3
"""Merge top-level run timing with structured substep records from run logs.

Input logs may contain lines such as::

    [timing] scope=secure.null_model parent=secure party=2 kind=phase \
        status=done milliseconds=2039.1 count=1

The base CSV may be the original three-column ``step,status,milliseconds``
file or a file previously expanded by this script.  Re-running the parser
keeps only top-level rows from the base, so detailed rows from an older run
cannot leak into a new summary.
"""
import argparse
import csv
import json
import os
import re
import shlex
import sys


FIELDS = [
    "step", "status", "milliseconds", "scope", "parent", "party", "kind",
    "mode", "phenotype", "phenotype_name", "count", "min_ms", "q1_ms", "mean_ms",
    "q3_ms", "max_ms",
]
TIMING_PREFIX = "[timing]"


def _number(value, default=""):
    """Return a finite-looking numeric string without changing its precision."""
    if value in (None, ""):
        return default
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    if number != number or number in (float("inf"), float("-inf")):
        return default
    return str(value)


def read_base_rows(path):
    """Read only root step rows; discard any detailed rows from a prior merge."""
    if not path or not os.path.exists(path):
        return []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or not {"step", "status", "milliseconds"}.issubset(reader.fieldnames):
            raise ValueError(f"{path}: expected step,status,milliseconds columns")
        expanded = "kind" in reader.fieldnames
        rows = []
        for source in reader:
            # Expanded files label top-level rows as kind=step.  For backwards
            # compatibility, a blank kind is also root only when parent is blank.
            if expanded and source.get("kind") not in ("", "step"):
                continue
            if expanded and source.get("parent") not in (None, ""):
                continue
            step = (source.get("step") or "").strip()
            if not step:
                continue
            row = {field: "" for field in FIELDS}
            row.update({
                "step": step,
                "status": (source.get("status") or "").strip(),
                "milliseconds": _number(source.get("milliseconds"), "0"),
                "scope": (source.get("scope") or step).strip(),
                "kind": "step",
                "count": _number(source.get("count"), "1"),
            })
            rows.append(row)
    return rows


def parse_timing_line(line, default_step, default_party=""):
    marker = line.find(TIMING_PREFIX)
    if marker < 0:
        return None
    try:
        tokens = shlex.split(line[marker + len(TIMING_PREFIX):].strip())
    except ValueError:
        return None
    values = {}
    for token in tokens:
        if "=" not in token:
            continue
        key, value = token.split("=", 1)
        values[key.strip()] = value.strip()

    # Go records use scope=secure plus stage=<mathematical leaf>, whereas the
    # Python prep/compare records put the leaf directly in scope.  The CSV has
    # no separate stage column, so normalize both forms into its scope column.
    record_scope = values.get("scope", "").strip()
    stage = values.get("stage", "").strip()
    scope = stage or record_scope
    milliseconds = _number(values.get("milliseconds", values.get("ms")))
    if not scope or milliseconds == "":
        return None
    parent = values.get("parent", default_step).strip()
    step = values.get("step", default_step or record_scope).strip()
    row = {field: "" for field in FIELDS}
    row.update({
        "step": step,
        "status": values.get("status", "done").strip(),
        "milliseconds": milliseconds,
        "scope": scope,
        "parent": parent,
        "party": values.get("party", default_party).strip(),
        "kind": values.get("kind", "phase").strip(),
        "mode": values.get("mode", "").strip(),
        "phenotype": _number(values.get("phenotype")),
        "count": _number(values.get("count"), "1"),
        "min_ms": _number(values.get("min_ms")),
        "q1_ms": _number(values.get("q1_ms")),
        "mean_ms": _number(values.get("mean_ms")),
        "q3_ms": _number(values.get("q3_ms")),
        "max_ms": _number(values.get("max_ms")),
    })
    return row


def add_phenotype_names(rows, manifest_path):
    """Map zero-based secure timing indices to manifest phenotype names."""
    if not manifest_path:
        return
    with open(manifest_path) as handle:
        phenotypes = json.load(handle).get("phenotypes", [])
    for row in rows:
        if row["phenotype"] == "":
            continue
        phenotype = int(float(row["phenotype"]))
        if 0 <= phenotype < len(phenotypes):
            row["phenotype_name"] = phenotypes[phenotype]


def print_phenotype_summary(rows):
    """Print party-2 phenotype subspans after all top-level timing output."""
    detail = [row for row in rows if row["kind"] == "phenotype" and row["party"] == "2"]
    if not detail:
        return
    values = {}
    for row in detail:
        key = (row["phenotype"], row["phenotype_name"] or f"phenotype_{row['phenotype']}")
        values.setdefault(key, {})[row["scope"]] = float(row["milliseconds"])
    print("=== PHENOTYPE TIMING (party 2; exact sequential subspans) ===")
    print(f"  {'phenotype':<62} {'null RSS':>12} {'score/Q/L':>12} {'total':>12}")
    for phenotype, name in sorted(values, key=lambda item: int(float(item[0]))):
        stages = values[(phenotype, name)]
        null_rss = stages.get("phenotype_null_rss", 0.0)
        score_ql = stages.get("phenotype_score_ql", 0.0)
        print(f"  {name:<62} {null_rss / 1000:>11.3f}s {score_ql / 1000:>11.3f}s "
              f"{(null_rss + score_ql) / 1000:>11.3f}s")
    print("  Shared/batched null factorization, genotype variance/moments, and finalization are not split by phenotype.")


def infer_party(path):
    match = re.search(r"party[_-]?(\d+)", os.path.basename(path), re.IGNORECASE)
    return match.group(1) if match else ""


def read_log(path, default_step, default_party=""):
    if not path or not os.path.exists(path):
        print(f"warning: timing log not found, skipped: {path}", file=sys.stderr)
        return []
    rows = []
    with open(path, errors="replace") as handle:
        for line in handle:
            row = parse_timing_line(line, default_step, default_party)
            if row is not None:
                rows.append(row)
    return rows


def write_rows(path, rows):
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    temporary = f"{path}.tmp"
    with open(temporary, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base_csv", help="top-level timing CSV (may also be --output)")
    parser.add_argument("--prep-log")
    parser.add_argument("--party-log", action="append", default=[], help="repeat for each party log")
    parser.add_argument("--compare-log")
    parser.add_argument("--manifest", help="manifest.json used to label phenotype timing rows")
    parser.add_argument("--log", action="append", default=[], metavar="STEP:PATH",
                        help="additional structured log, e.g. build:/tmp/build.log")
    parser.add_argument("--output", help="output CSV; default overwrites base_csv")
    args = parser.parse_args(argv)

    rows = read_base_rows(args.base_csv)
    if args.prep_log:
        rows.extend(read_log(args.prep_log, "prep", "driver"))
    for path in args.party_log:
        rows.extend(read_log(path, "secure", infer_party(path)))
    if args.compare_log:
        rows.extend(read_log(args.compare_log, "compare", "driver"))
    for spec in args.log:
        if ":" not in spec:
            parser.error(f"--log expects STEP:PATH, got {spec!r}")
        step, path = spec.split(":", 1)
        rows.extend(read_log(path, step))

    add_phenotype_names(rows, args.manifest)

    output = args.output or args.base_csv
    write_rows(output, rows)
    print_phenotype_summary(rows)
    print(f"timing CSV -> {output} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
