#!/usr/bin/env python3
"""Summarize wall-clock and protocol timings from a secure-skat run directory."""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path


TIMESTAMP_RE = re.compile(
    r"(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:[+-]\d{2}:\d{2}|Z))"
)
PARTY_LOG_RE = re.compile(r"stdout_party(\d+)\.txt$")
STEP2_RE = re.compile(r"SKAT Step 2/4: block (\d+)/(\d+) score loading \((\d+) SNPs\)")
STEP3_RE = re.compile(r"SKAT Step 3/4: block weight calculation \((\d+) SNPs\)")
STEP4_RE = re.compile(r"SKAT Step 4/4: block statistic aggregation \((\d+) SNPs\)")
MATMULT_RE = re.compile("MatMult: block\\s+(\\d+)\\s*/\\s*(\\d+)\\s+elapsed time\\s+([0-9a-zA-Z.\\u00b5]+)")
GO_DURATION_RE = re.compile("([0-9]*\\.?[0-9]+)(ns|us|\\u00b5s|ms|s|m|h)")


@dataclass
class BlockTiming:
    index: int
    total_blocks: int | None = None
    snps: int | None = None
    step2_at: datetime | None = None
    step3_at: datetime | None = None
    step4_at: datetime | None = None
    matmult_seconds: float | None = None


@dataclass
class PartyTiming:
    pid: int
    path: Path
    collective_init_seen: bool = False
    protocol_start: datetime | None = None
    qc_start: datetime | None = None
    qc_finished: datetime | None = None
    step1_at: datetime | None = None
    finished_compute: datetime | None = None
    output_saved: datetime | None = None
    blocks: dict[int, BlockTiming] = field(default_factory=dict)

    def ordered_blocks(self) -> list[BlockTiming]:
        return [self.blocks[i] for i in sorted(self.blocks)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_root", type=Path, help="Run output directory produced by run_example.sh")
    return parser.parse_args()


def read_metadata(path: Path) -> dict[str, str]:
    metadata_file = path / "run_metadata.txt"
    metadata: dict[str, str] = {}
    if not metadata_file.exists():
        return metadata
    for raw_line in metadata_file.read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" not in raw_line:
            continue
        key, value = raw_line.split("=", 1)
        metadata[key.strip()] = value.strip()
    return metadata


def parse_datetime(value: str | None) -> datetime | None:
    if not value:
        return None

    value = value.strip()
    if re.search(r"[+-]\d{4}$", value):
        value = f"{value[:-5]}{value[-5:-2]}:{value[-2:]}"
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
        return parsed.replace(tzinfo=None)
    except ValueError:
        pass

    # Older metadata used "YYYY-mm-dd HH:MM:SS EDT". For durations from the
    # same run, dropping the timezone name is enough and avoids platform-specific
    # %Z parsing behavior.
    match = re.match(r"(.+?)\s+[A-Z]{2,5}$", value)
    if match:
        value = match.group(1)
    try:
        return datetime.strptime(value, "%Y-%m-%d %H:%M:%S")
    except ValueError:
        return None


def parse_metadata_time(metadata: dict[str, str], prefix: str) -> datetime | None:
    epoch = metadata.get(f"{prefix}_epoch")
    if epoch:
        try:
            return datetime.fromtimestamp(float(epoch))
        except ValueError:
            pass
    return parse_datetime(metadata.get(f"{prefix}_at_iso")) or parse_datetime(metadata.get(f"{prefix}_at"))


def parse_go_duration_seconds(value: str) -> float | None:
    total = 0.0
    matched = False
    scale = {
        "ns": 1e-9,
        "us": 1e-6,
        "\u00b5s": 1e-6,
        "ms": 1e-3,
        "s": 1.0,
        "m": 60.0,
        "h": 3600.0,
    }
    for number, unit in GO_DURATION_RE.findall(value):
        matched = True
        total += float(number) * scale[unit]
    return total if matched else None


def parse_log_timestamp(line: str) -> datetime | None:
    match = TIMESTAMP_RE.search(line)
    if not match:
        return None
    return parse_datetime(match.group(1))


def block_for(party: PartyTiming, index: int) -> BlockTiming:
    if index not in party.blocks:
        party.blocks[index] = BlockTiming(index=index)
    return party.blocks[index]


def parse_party_log(path: Path) -> PartyTiming | None:
    match = PARTY_LOG_RE.search(path.name)
    if not match:
        return None

    party = PartyTiming(pid=int(match.group(1)), path=path)
    current_block: int | None = None

    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        timestamp = parse_log_timestamp(line)

        if "CollectiveInit started" in line:
            party.collective_init_seen = True
        if timestamp and "Running rare-variant protocol in mode:" in line:
            party.protocol_start = party.protocol_start or timestamp
        elif timestamp and "Starting QC" in line:
            party.qc_start = party.qc_start or timestamp
        elif timestamp and "Finished QC" in line:
            party.qc_finished = party.qc_finished or timestamp
        elif timestamp and "SKAT Step 1/4: Null model residuals" in line:
            party.step1_at = party.step1_at or timestamp
        elif timestamp and "Finished rare-variant statistic computation" in line:
            party.finished_compute = timestamp
        elif timestamp and "Output collectively decrypted and saved to:" in line:
            party.output_saved = timestamp

        if timestamp:
            if step2_match := STEP2_RE.search(line):
                current_block = int(step2_match.group(1))
                block = block_for(party, current_block)
                block.total_blocks = int(step2_match.group(2))
                block.snps = int(step2_match.group(3))
                block.step2_at = timestamp
                continue

            if step3_match := STEP3_RE.search(line):
                if current_block is not None:
                    block = block_for(party, current_block)
                    block.snps = block.snps or int(step3_match.group(1))
                    block.step3_at = timestamp
                continue

            if step4_match := STEP4_RE.search(line):
                if current_block is not None:
                    block = block_for(party, current_block)
                    block.snps = block.snps or int(step4_match.group(1))
                    block.step4_at = timestamp
                continue

            if matmult_match := MATMULT_RE.search(line):
                block_index = int(matmult_match.group(1))
                block = block_for(party, block_index)
                block.total_blocks = int(matmult_match.group(2))
                block.matmult_seconds = parse_go_duration_seconds(matmult_match.group(3))

    return party


def seconds_between(start: datetime | None, end: datetime | None) -> float | None:
    if start is None or end is None:
        return None
    seconds = (end - start).total_seconds()
    return seconds if seconds >= 0 else None


def format_seconds(seconds: float | None) -> str:
    if seconds is None:
        return "n/a"
    if seconds < 0.001:
        return "0s"
    if seconds < 1:
        milliseconds = seconds * 1000
        value = f"{milliseconds:.1f}".rstrip("0").rstrip(".")
        return f"{value}ms"
    if seconds < 60:
        value = f"{seconds:.2f}".rstrip("0").rstrip(".")
        return f"{value}s"

    total = int(round(seconds))
    hours, remainder = divmod(total, 3600)
    minutes, secs = divmod(remainder, 60)
    parts: list[str] = []
    if hours:
        parts.append(f"{hours}h")
    if minutes:
        parts.append(f"{minutes}m")
    if secs or not parts:
        parts.append(f"{secs}s")
    return " ".join(parts)


def max_duration(values: list[float | None]) -> float | None:
    present = [value for value in values if value is not None]
    return max(present) if present else None


def min_duration(values: list[float | None]) -> float | None:
    present = [value for value in values if value is not None]
    return min(present) if present else None


def mean_duration(values: list[float | None]) -> float | None:
    present = [value for value in values if value is not None]
    return sum(present) / len(present) if present else None


def pid_breakdown(values_by_pid: list[tuple[int, float | None]]) -> str:
    present = [(pid, value) for pid, value in values_by_pid if value is not None]
    if not present:
        return ""
    return "(" + ", ".join(f"pid{pid} {format_seconds(value)}" for pid, value in present) + ")"


def first_step2_at(party: PartyTiming) -> datetime | None:
    step2_times = [block.step2_at for block in party.ordered_blocks() if block.step2_at is not None]
    return min(step2_times) if step2_times else None


def next_block_step2_at(party: PartyTiming, block_index: int) -> datetime | None:
    later_step2_times = [
        block.step2_at
        for block in party.ordered_blocks()
        if block.index > block_index and block.step2_at is not None
    ]
    return min(later_step2_times) if later_step2_times else None


def block_step2_to_step3_seconds(party: PartyTiming, block_index: int) -> float | None:
    block = party.blocks.get(block_index)
    if block is None:
        return None
    timestamp_seconds = seconds_between(block.step2_at, block.step3_at)
    if timestamp_seconds is None:
        return block.matmult_seconds
    if block.matmult_seconds is None:
        return timestamp_seconds
    return max(timestamp_seconds, block.matmult_seconds)


def metric_line(label: str, values_by_pid: list[tuple[int, float | None]], indent: str = "    ") -> str:
    max_value = max_duration([value for _, value in values_by_pid])
    breakdown = pid_breakdown(values_by_pid)
    if breakdown:
        return f"{indent}{label}: max {format_seconds(max_value)} {breakdown}"
    return f"{indent}{label}: n/a"


def block_summary_line(
    label: str,
    values_by_block: dict[int, list[tuple[int, float | None]]],
    indent: str = "    ",
) -> str:
    block_maxes = [
        max_duration([value for _, value in values_by_pid])
        for _, values_by_pid in sorted(values_by_block.items())
    ]
    present_block_maxes = [value for value in block_maxes if value is not None]
    if not present_block_maxes:
        return f"{indent}{label}: n/a"

    pids = sorted({pid for values_by_pid in values_by_block.values() for pid, _ in values_by_pid})
    per_pid_parts = []
    for pid in pids:
        pid_values = [
            value
            for values_by_pid in values_by_block.values()
            for value_pid, value in values_by_pid
            if value_pid == pid and value is not None
        ]
        if pid_values:
            per_pid_parts.append(f"pid{pid} avg {format_seconds(mean_duration(pid_values))}")

    per_pid_text = f"; per_pid_avg ({', '.join(per_pid_parts)})" if per_pid_parts else ""
    return (
        f"{indent}{label}: block_max avg {format_seconds(mean_duration(present_block_maxes))}, "
        f"min {format_seconds(min_duration(present_block_maxes))}, "
        f"max {format_seconds(max_duration(present_block_maxes))} "
        f"over {len(present_block_maxes)} blocks{per_pid_text}"
    )


def summarize(run_root: Path) -> list[str]:
    metadata = read_metadata(run_root)
    started_at = parse_metadata_time(metadata, "started")
    finished_at = parse_metadata_time(metadata, "finished")
    parties = [
        party
        for path in sorted(run_root.glob("stdout_party*.txt"))
        if (party := parse_party_log(path)) is not None
    ]

    lines = ["Collapsed timing summary:"]
    total = seconds_between(started_at, finished_at)
    if total is not None:
        lines.append(f"  total_wall_time: {format_seconds(total)}")
    else:
        lines.append("  total_wall_time: n/a")
    if exit_status := metadata.get("exit_status"):
        lines.append(f"  exit_status: {exit_status}")

    if not parties:
        lines.append("  party_logs: none found")
        return lines

    setup_values = [
        (party.pid, seconds_between(started_at, party.protocol_start)) for party in parties
    ]
    protocol_values = [
        (
            party.pid,
            seconds_between(party.protocol_start, party.output_saved or party.finished_compute),
        )
        for party in parties
    ]
    qc_values = [
        (party.pid, seconds_between(party.qc_start, party.qc_finished)) for party in parties
    ]
    step1_values = [
        (party.pid, seconds_between(party.step1_at, first_step2_at(party))) for party in parties
    ]

    lines.append(metric_line("setup_network_keygen_until_protocol_start", setup_values, indent="  "))
    if any(party.collective_init_seen for party in parties):
        lines.append("  keygen_note: exact PubKey/Relin/RotKey substeps are not timestamped in current logs")
    lines.append(metric_line("protocol_total_until_output", protocol_values, indent="  "))
    lines.append("  timing_note: block step2 timings use exact MatMult elapsed times when available")
    lines.append("  ordered_steps:")
    lines.append(metric_line("qc", qc_values))
    lines.append(metric_line("step1_null_model_residuals_until_first_block", step1_values))

    block_indices = sorted({block.index for party in parties for block in party.blocks.values()})
    step2_values_by_block: dict[int, list[tuple[int, float | None]]] = {}
    step3_values_by_block: dict[int, list[tuple[int, float | None]]] = {}
    after_step4_values_by_block: dict[int, list[tuple[int, float | None]]] = {}
    for block_index in block_indices:
        step2_values_by_block[block_index] = [
            (
                party.pid,
                block_step2_to_step3_seconds(party, block_index),
            )
            for party in parties
        ]
        step3_values_by_block[block_index] = [
            (
                party.pid,
                seconds_between(
                    party.blocks.get(block_index).step3_at if party.blocks.get(block_index) else None,
                    party.blocks.get(block_index).step4_at if party.blocks.get(block_index) else None,
                ),
            )
            for party in parties
        ]
        after_step4_values_by_block[block_index] = [
            (
                party.pid,
                seconds_between(
                    party.blocks.get(block_index).step4_at if party.blocks.get(block_index) else None,
                    next_block_step2_at(party, block_index) or party.finished_compute,
                ),
            )
            for party in parties
        ]
    if block_indices:
        lines.append(f"    block_count: {len(block_indices)}")
        lines.append(block_summary_line("block_step2_score_loading_to_step3_weights", step2_values_by_block))
        lines.append(block_summary_line("block_step3_weights_to_step4_aggregation", step3_values_by_block))
        lines.append(block_summary_line("block_step4_aggregation_to_next_or_finish", after_step4_values_by_block))

    output_values = [
        (party.pid, seconds_between(party.finished_compute, party.output_saved)) for party in parties
    ]
    lines.append(metric_line("decrypt_and_write_after_compute", output_values))
    return lines


def main() -> int:
    args = parse_args()
    run_root = args.run_root.expanduser().resolve()
    if not run_root.exists():
        raise SystemExit(f"Run directory does not exist: {run_root}")
    print("\n".join(summarize(run_root)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
