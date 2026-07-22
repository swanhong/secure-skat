#!/usr/bin/env python3
"""Aggregate total and mathematical-stage secure-SKAT communication.

The Go runtime reports application-layer serialized bytes. Network-wide traffic
is the sum of sent bytes over parties; adding sent and received would count each
transfer twice. ``G^T G``, ``G^T X``, and ``G^T y`` local contractions are
reported explicitly as zero-byte leaves. Secure operations that consume those
local values have separate rows and carry the actual communication cost.
"""

import argparse
import csv
import re
import shlex
from pathlib import Path
from typing import Optional


TOTAL_LINE = re.compile(
    r"\[communication\]\s+"
    r"scope=(?P<scope>\S+)\s+"
    r"(?:mode=(?P<mode>\S+)\s+)?"
    r"party=(?P<party>\d+)\s+"
    r"sent_bytes=(?P<sent_bytes>\d+)\s+"
    r"received_bytes=(?P<received_bytes>\d+)\s+"
    r"sent_messages=(?P<sent_messages>\d+)\s+"
    r"received_messages=(?P<received_messages>\d+)"
)
STAGE_PREFIX = "[communication-stage]"
COUNTERS = ("sent_bytes", "received_bytes", "sent_messages", "received_messages")
FIELDS = ["scope", "stage", "math_objects", "parent", "kind", "mode", "party", *COUNTERS]
MATH_OBJECTS = {
    "protocol_total": "all secure-SKAT traffic",
    "collective_key_setup": "MHE public/relinearization/rotation keys",
    "initialization_other": "post-key-setup initialization",
    "load_private_only": "private genotype file loading (local)",
    "assoc_init": "association-test initialization",
    "null_local_xtx_xty_yty": "X^T X; X^T y; y^T y (local)",
    "null_aggregate_xtx": "sum_i X_i^T X_i",
    "null_aggregate_xty": "sum_i X_i^T y_i",
    "null_aggregate_yty": "sum_i y_i^T y_i",
    "null_solve": "beta_hat; (X^T X)^-1",
    "null_beta_pack": "beta_hat shares -> MHE",
    "null_rss": "residual sum of squares",
    "null_other": "unclassified null-model overhead",
    "pre_block_setup": "per-run gene-output allocation and setup",
    "gene_shape_sync": "m_pub shape metadata",
    "gene_local_public_gtg_gtx_gty": "G_pub^T G_pub; G_pub^T X; G_pub^T y; G_pub^T 1 (local)",
    "gene_public_score_gtx_gty": "G_pub^T y - G_pub^T X beta_hat",
    "gene_public_weight_dosage_maf": "G_pub^T 1 -> MAF -> w",
    "gene_public_stat_share": "Q_pub; burden linear score",
    "gene_local_private_gtg_gtx_gty": "G_priv^T X; G_priv^T y; local burden/moment contractions",
    "gene_private_stat_share": "Q_priv; private burden score",
    "gene_burden_public_gtg": "w_pub^T G_pub^T G_pub w_pub",
    "gene_burden_public_gtx": "X^T G_pub w_pub",
    "gene_burden_private_cross": "G_pub^T G_priv; G_priv^T G_priv cross terms",
    "gene_burden_projection": "(X^T z)^T (X^T X)^-1 (X^T z)",
    "gene_moments_setup_gtx": "G^T X moment setup",
    "gene_moments_public_trace_gtg_gtx": "G^T G action + G^T X projection traces",
    "gene_moments_private_cross": "private/cross kernel-moment corrections",
    "gene_other": "unclassified per-gene overhead",
    "finalize_scale": "Q scaling by sigma_hat^2",
    "finalize_burden_pvalue": "burden sqrt(T/2) pivot",
    "finalize_skat_pvalue": "SKAT moments + Wilson-Hilferty z",
    "finalize_output_pack": "result shares -> MHE",
    "finalize_other": "unclassified post-block finalization overhead",
    "decrypt_outputs": "collective output decryption",
    "run_other": "unclassified secure-run traffic",
}
LOCAL_ZERO_STAGES = {
    "load_private_only",
    "null_local_xtx_xty_yty",
    "gene_local_public_gtg_gtx_gty",
    "gene_local_private_gtg_gtx_gty",
}


def selected_record(path: Path, scope: str, mode: str, occurrence: str) -> dict:
    matches = []
    for line in path.read_text(errors="replace").splitlines():
        match = TOTAL_LINE.search(line)
        if match and match.group("scope") == scope:
            record = {
                "scope": scope,
                "stage": "protocol_total",
                "parent": "",
                "kind": "total",
                "mode": match.group("mode") or "",
                "log": str(path),
            }
            for key in ("party", *COUNTERS):
                record[key] = int(match.group(key))
            matches.append(record)
    if not matches:
        raise SystemExit(f"no communication record with scope={scope!r} in {path}")
    mode_matches = [record for record in matches if record["mode"] == mode]
    if not mode_matches:
        # Backward compatibility for one-run logs emitted before totals carried
        # a mode tag. Multiple untagged totals are inherently ambiguous.
        if len(matches) == 1 and matches[0]["mode"] == "":
            mode_matches = matches
        else:
            modes = [record["mode"] or "<untagged>" for record in matches]
            raise SystemExit(f"no unambiguous total for mode={mode!r} in {path}; found {modes}")
    return mode_matches[0] if occurrence == "first" else mode_matches[-1]


def parse_stage_line(line: str) -> Optional[dict]:
    marker = line.find(STAGE_PREFIX)
    if marker < 0:
        return None
    try:
        tokens = shlex.split(line[marker + len(STAGE_PREFIX):].strip())
    except ValueError:
        return None
    values = {}
    for token in tokens:
        if "=" in token:
            key, value = token.split("=", 1)
            values[key] = value
    required = {"scope", "mode", "party", "stage", "parent", "kind", *COUNTERS}
    if not required.issubset(values):
        return None
    try:
        record = {key: values[key] for key in ("scope", "stage", "parent", "kind", "mode")}
        for key in ("party", *COUNTERS):
            record[key] = int(values[key])
    except ValueError:
        return None
    return record


def all_stage_records(path: Path, scope: str) -> list[dict]:
    records = []
    for line in path.read_text(errors="replace").splitlines():
        record = parse_stage_line(line)
        if record is not None and record["scope"] == scope:
            records.append(record)
    return records


def choose_mode(stage_records_by_log: list[list[dict]], requested: Optional[str]) -> str:
    if requested:
        return requested
    last_modes = [records[-1]["mode"] for records in stage_records_by_log if records]
    if not last_modes:
        raise SystemExit("no [communication-stage] records found")
    if len(set(last_modes)) != 1:
        raise SystemExit(f"latest communication-stage mode differs across party logs: {last_modes}")
    return last_modes[0]


def selected_stage_records(records: list[dict], mode: str, occurrence: str) -> list[dict]:
    """Keep one complete ordered stage block for the requested mode."""
    runs = []
    current = []
    for record in (row for row in records if row["mode"] == mode):
        if record["stage"] == "collective_key_setup" and current:
            runs.append(current)
            current = []
        current.append(record)
    if current:
        runs.append(current)
    if not runs:
        return []
    return runs[0] if occurrence == "first" else runs[-1]


def totals(records: list[dict]) -> dict:
    return {key: sum(int(row[key]) for row in records) for key in COUNTERS}


def validate_balance(label: str, values: dict) -> None:
    if values["sent_bytes"] != values["received_bytes"]:
        raise SystemExit(
            f"{label}: communication counter mismatch: "
            f"sum(sent)={values['sent_bytes']} != sum(received)={values['received_bytes']}"
        )
    if values["sent_messages"] != values["received_messages"]:
        raise SystemExit(
            f"{label}: message-count mismatch: "
            f"sum(sent)={values['sent_messages']} != sum(received)={values['received_messages']}"
        )


def network_row(template: dict, records: list[dict]) -> dict:
    row = {key: template.get(key, "") for key in FIELDS}
    row["math_objects"] = MATH_OBJECTS.get(str(row["stage"]), "")
    row["party"] = "network_total"
    row.update(totals(records))
    return row


def csv_row(record: dict) -> dict:
    row = {key: record.get(key, "") for key in FIELDS}
    row["math_objects"] = MATH_OBJECTS.get(str(row["stage"]), "")
    return row


def write_csv(path: Path, total_records: list[dict], stage_records: list[dict],
              stage_order: list[str], mode: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for record in total_records:
        row = csv_row(record)
        row["mode"] = mode
        rows.append(row)
    rows.append(network_row(rows[0], rows))

    for stage in stage_order:
        group = sorted(
            (row for row in stage_records if row["stage"] == stage),
            key=lambda row: int(row["party"]),
        )
        rows.extend(csv_row(row) for row in group)
        rows.append(network_row(group[0], group))

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logs", nargs="+", type=Path)
    parser.add_argument("--scope", default="skat_fed_total")
    parser.add_argument("--occurrence", choices=("first", "last"), default="last")
    parser.add_argument("--mode", choices=("raw_q", "skat_p"), help="default: latest mode in all logs")
    parser.add_argument(
        "--expected-parties",
        nargs="+",
        type=int,
        default=(0, 1, 2),
        metavar="PARTY",
        help="exact party IDs required in the selected records (default: 0 1 2)",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    all_stages = [all_stage_records(path, args.scope) for path in args.logs]
    mode = choose_mode(all_stages, args.mode)
    total_records = [selected_record(path, args.scope, mode, args.occurrence) for path in args.logs]
    total_records.sort(key=lambda row: int(row["party"]))
    parties = [int(row["party"]) for row in total_records]
    if len(parties) != len(set(parties)):
        raise SystemExit(f"duplicate party records: {parties}")
    expected_parties = sorted(args.expected_parties)
    if len(expected_parties) != len(set(expected_parties)):
        raise SystemExit(f"duplicate expected party IDs: {args.expected_parties}")
    if parties != expected_parties:
        raise SystemExit(f"party set mismatch: expected {expected_parties}, got {parties}")

    per_log_stages = [selected_stage_records(records, mode, args.occurrence) for records in all_stages]
    for path, records in zip(args.logs, per_log_stages):
        if not records:
            raise SystemExit(f"no communication-stage records for mode={mode!r} in {path}")

    expected_stage_names = {row["stage"] for row in per_log_stages[0]}
    stage_order = [row["stage"] for row in per_log_stages[0]]
    for path, records in zip(args.logs[1:], per_log_stages[1:]):
        names = {row["stage"] for row in records}
        if names != expected_stage_names:
            raise SystemExit(
                f"stage set mismatch in {path}: missing={sorted(expected_stage_names - names)}, "
                f"extra={sorted(names - expected_stage_names)}"
            )
    stage_records = [row for records in per_log_stages for row in records]

    for row in stage_records:
        if row["stage"] in LOCAL_ZERO_STAGES and any(int(row[key]) != 0 for key in COUNTERS):
            raise SystemExit(
                f"local-only stage unexpectedly communicated: party={row['party']} "
                f"stage={row['stage']}"
            )

    total_network = totals(total_records)
    validate_balance(args.scope, total_network)
    stage_network = totals(stage_records)
    validate_balance(f"{args.scope} stage leaves", stage_network)
    for stage in stage_order:
        validate_balance(
            f"{args.scope} stage={stage}",
            totals([row for row in stage_records if row["stage"] == stage]),
        )
    for counter in COUNTERS:
        if stage_network[counter] != total_network[counter]:
            raise SystemExit(
                f"stage partition mismatch for {counter}: "
                f"leaves={stage_network[counter]} total={total_network[counter]}"
            )

    # A per-party partition catches swapped traffic between parties that a
    # network-wide equality alone would miss.
    for total_record in total_records:
        party = int(total_record["party"])
        leaves = totals([row for row in stage_records if int(row["party"]) == party])
        for counter in COUNTERS:
            if leaves[counter] != int(total_record[counter]):
                raise SystemExit(
                    f"party {party} stage partition mismatch for {counter}: "
                    f"leaves={leaves[counter]} total={total_record[counter]}"
                )

    print("=== COMMUNICATION (application-layer serialized bytes) ===")
    print(f"  scope={args.scope} mode={mode}")
    for row in total_records:
        print(
            f"  party {row['party']}: sent={row['sent_bytes']} B  "
            f"received={row['received_bytes']} B  "
            f"messages={row['sent_messages']}/{row['received_messages']} (sent/received)"
        )
    print(
        f"  network transmitted: {total_network['sent_bytes']} B "
        f"({total_network['sent_bytes'] / 2**20:.3f} MiB)"
    )
    print("  stage breakdown (network transmitted; local contractions should be zero):")
    for stage in stage_order:
        value = sum(int(row["sent_bytes"]) for row in stage_records if row["stage"] == stage)
        print(f"    {stage:<43} {value:>12} B")

    if args.output:
        write_csv(args.output, total_records, stage_records, stage_order, mode)
        print(f"  wrote {args.output}")


if __name__ == "__main__":
    main()
