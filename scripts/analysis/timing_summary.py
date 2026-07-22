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
import os
import re
import shlex
import sys


FIELDS = [
    "step", "status", "milliseconds", "scope", "parent", "party", "kind",
    "mode", "count", "min_ms", "q1_ms", "mean_ms", "q3_ms", "max_ms",
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
        "count": _number(values.get("count"), "1"),
        "min_ms": _number(values.get("min_ms")),
        "q1_ms": _number(values.get("q1_ms")),
        "mean_ms": _number(values.get("mean_ms")),
        "q3_ms": _number(values.get("q3_ms")),
        "max_ms": _number(values.get("max_ms")),
    })
    return row


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

    output = args.output or args.base_csv
    write_rows(output, rows)
    print(f"timing CSV -> {output} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
