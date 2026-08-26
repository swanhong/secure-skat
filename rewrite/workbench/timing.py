from __future__ import annotations

import csv
import time
from collections import defaultdict
from contextlib import contextmanager
from pathlib import Path
from typing import Callable, Iterator, Sequence


TIMING_FIELDS = [
    "component",
    "scope",
    "chromosome",
    "party",
    "lane",
    "phase",
    "parent_phase",
    "event_index",
    "batch_index",
    "batch_width",
    "batch_gene_count",
    "gene_index",
    "gene_id",
    "phenotype_index",
    "phenotype_name",
    "trace_mode",
    "measurement_kind",
    "elapsed_seconds",
    "sample_count_a",
    "sample_count_b",
    "public_variant_count",
    "private_variant_count",
    "ckks",
    "data_bits",
    "frac_bits",
    "probes",
    "status",
]


def empty_timing_row(**values: object) -> dict[str, object]:
    row = {field: "" for field in TIMING_FIELDS}
    row.update({
        "measurement_kind": "actual",
        "status": "success",
        **values,
    })
    return row


class TimingRecorder:
    def __init__(
        self,
        component: str,
        chromosome: str = "",
        party: str = "",
        lane: str = "",
        clock: Callable[[], float] = time.perf_counter,
    ) -> None:
        self.component = component
        self.chromosome = chromosome
        self.party = party
        self.lane = lane
        self.clock = clock
        self.rows: list[dict[str, object]] = []
        self.defaults: dict[str, object] = {}

    def set_defaults(self, **values: object) -> None:
        self.defaults.update(values)

    @contextmanager
    def measure(self, **values: object) -> Iterator[None]:
        started = self.clock()
        status = "success"
        try:
            yield
        except BaseException:
            status = "failure"
            raise
        finally:
            self.record(self.clock() - started, status=status, **values)

    def record(
        self,
        elapsed_seconds: float,
        status: str = "success",
        **values: object,
    ) -> None:
        self.rows.append(empty_timing_row(
            component=self.component,
            chromosome=self.chromosome,
            party=self.party,
            lane=self.lane,
            event_index=len(self.rows),
            elapsed_seconds=f"{elapsed_seconds:.9f}",
            status=status,
            **values,
        ))

    def write(self, path: Path | None) -> None:
        if path is None:
            return
        for row in self.rows:
            for field, value in self.defaults.items():
                if row[field] == "":
                    row[field] = value
        write_timing_rows(path, self.rows)


def read_timing_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as file:
        reader = csv.DictReader(file)
        if reader.fieldnames != TIMING_FIELDS:
            raise ValueError(f"unexpected timing header in {path}")
        return list(reader)


def write_timing_rows(path: Path, rows: Sequence[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=TIMING_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def combine_timing_files(paths: Sequence[Path], output: Path) -> list[dict[str, str]]:
    rows = []
    for path in paths:
        rows.extend(read_timing_rows(path))
    write_timing_rows(output, rows)
    return rows


def append_timing_event(path: Path, **values: object) -> dict[str, object]:
    rows: list[dict[str, object]] = []
    if path.exists():
        rows.extend(read_timing_rows(path))
    row = empty_timing_row(event_index=len(rows), **values)
    rows.append(row)
    write_timing_rows(path, rows)
    return row


def aggregate_timing_events(
    input_root: Path,
    output_dir: Path,
    chromosomes: Sequence[int],
    coordinator_timing: Path | None = None,
) -> list[dict[str, str]]:
    rows = []
    for chromosome in chromosomes:
        path = input_root / f"chr{chromosome}" / "timing_events.csv"
        if not path.exists():
            raise ValueError(f"missing {path}")
        rows.extend(read_timing_rows(path))
    if coordinator_timing is not None:
        rows.extend(read_timing_rows(coordinator_timing))

    rows.sort(key=timing_sort_key)
    write_timing_rows(output_dir / "all_timing_events.csv", rows)
    summary = build_timing_summary(rows)
    write_timing_rows(output_dir / "timing_summary.csv", summary)
    return summary


def build_timing_summary(
    rows: Sequence[dict[str, str]],
) -> list[dict[str, object]]:
    groups: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row["scope"] in {"gene", "gene_phenotype"}:
            scope = row["scope"]
            gene_index = row["gene_index"]
            gene_id = row["gene_id"]
            phenotype_index = row["phenotype_index"]
            phenotype_name = row["phenotype_name"]
        elif row["scope"] in {"batch_phenotype", "phenotype"}:
            scope = "phenotype"
            gene_index = gene_id = ""
            phenotype_index = row["phenotype_index"]
            phenotype_name = row["phenotype_name"]
        else:
            scope = "run" if row["component"] == "coordinator" else "chromosome"
            gene_index = gene_id = phenotype_index = phenotype_name = ""

        key = (
            row["component"],
            scope,
            row["chromosome"],
            row["party"],
            row["lane"],
            row["phase"],
            row["parent_phase"],
            gene_index,
            gene_id,
            phenotype_index,
            phenotype_name,
            row["trace_mode"],
            row["measurement_kind"],
        )
        groups[key].append(row)

    summary = []
    for key, grouped in groups.items():
        first = grouped[0]
        (
            component,
            scope,
            chromosome,
            party,
            lane,
            phase,
            parent_phase,
            gene_index,
            gene_id,
            phenotype_index,
            phenotype_name,
            trace_mode,
            measurement_kind,
        ) = key
        status = "success" if all(row["status"] == "success" for row in grouped) else "failure"
        public_variant_count = sum_optional_integers(
            row["public_variant_count"] for row in grouped
        )
        private_variant_count = sum_optional_integers(
            row["private_variant_count"] for row in grouped
        )
        summary.append(empty_timing_row(
            component=component,
            scope=scope,
            chromosome=chromosome,
            party=party,
            lane=lane,
            phase=phase,
            parent_phase=parent_phase,
            gene_index=gene_index,
            gene_id=gene_id,
            phenotype_index=phenotype_index,
            phenotype_name=phenotype_name,
            trace_mode=trace_mode,
            measurement_kind=measurement_kind,
            elapsed_seconds=f"{sum(float(row['elapsed_seconds']) for row in grouped):.9f}",
            sample_count_a=first["sample_count_a"],
            sample_count_b=first["sample_count_b"],
            public_variant_count=public_variant_count,
            private_variant_count=private_variant_count,
            ckks=first["ckks"],
            data_bits=first["data_bits"],
            frac_bits=first["frac_bits"],
            probes=first["probes"],
            status=status,
        ))
    summary.extend(build_phenotype_independent_rollups(summary))
    summary.sort(key=timing_sort_key)
    return summary


def build_phenotype_independent_rollups(
    rows: Sequence[dict[str, object]],
) -> list[dict[str, object]]:
    groups: dict[tuple[str, ...], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        if row["component"] != "secure" or row["phase"] not in {
            "prepare_gene_batch",
            "kernel_statistics",
        }:
            continue
        key = (
            str(row["scope"]),
            str(row["chromosome"]),
            str(row["party"]),
            str(row["lane"]),
            str(row["gene_index"]),
            str(row["gene_id"]),
            str(row["trace_mode"]),
            str(row["measurement_kind"]),
        )
        groups[key].append(row)

    rollups = []
    for grouped in groups.values():
        phases = {str(row["phase"]) for row in grouped}
        if phases != {"prepare_gene_batch", "kernel_statistics"}:
            continue
        row = dict(grouped[0])
        row["phase"] = "phenotype_independent_gene_work"
        row["parent_phase"] = "packed_statistics"
        row["elapsed_seconds"] = f"{sum(float(item['elapsed_seconds']) for item in grouped):.9f}"
        row["status"] = (
            "success"
            if all(item["status"] == "success" for item in grouped)
            else "failure"
        )
        rollups.append(row)
    return rollups


def print_timing_summary(rows: Sequence[dict[str, object]]) -> None:
    workflow_phases = [
        "genotype_localization",
        "preprocessing",
        "shared_key_generation",
        "secure_protocol",
        "r_reference",
        "join_results",
        "chromosome_total",
    ]
    by_chromosome: dict[str, dict[str, float]] = defaultdict(dict)
    for row in rows:
        if row["component"] != "workflow" or row["scope"] != "chromosome":
            continue
        by_chromosome[str(row["chromosome"])][str(row["phase"])] = float(
            row["elapsed_seconds"]
        )

    print("\n=== BENCHMARK TIMING BY CHROMOSOME (seconds) ===")
    print("chromosome " + " ".join(f"{phase:>22}" for phase in workflow_phases))
    for chromosome in sorted(by_chromosome, key=chromosome_number):
        values = " ".join(
            f"{by_chromosome[chromosome].get(phase, 0.0):22.3f}"
            for phase in workflow_phases
        )
        print(f"{chromosome:>10} {values}")

    print("\n=== BENCHMARK PHASE DETAIL (seconds) ===")
    print(
        "component chromosome party phase parent phenotype trace_mode seconds"
    )
    phase_rows: dict[tuple[str, ...], float] = {}
    for row in rows:
        if row["scope"] in {"gene", "gene_phenotype"}:
            continue
        party = "max" if row["component"] == "secure" else str(row["party"])
        key = (
            str(row["component"]),
            str(row["chromosome"]),
            party,
            str(row["phase"]),
            str(row["parent_phase"]),
            str(row["phenotype_name"]),
            str(row["trace_mode"]),
        )
        elapsed = float(row["elapsed_seconds"])
        if key not in phase_rows or elapsed > phase_rows[key]:
            phase_rows[key] = elapsed
    for key in sorted(
        phase_rows,
        key=lambda item: (
            component_order(item[0]),
            chromosome_number(item[1]),
            item[3],
            item[5],
            item[6],
        ),
    ):
        print(" ".join((*key, f"{phase_rows[key]:.6f}")))

    print("\n=== BENCHMARK TIMING BY GENE (seconds) ===")
    print("chromosome gene_id secure_amortized_max r_gene_actual")
    secure: dict[tuple[str, str], list[float]] = defaultdict(list)
    reference: dict[tuple[str, str], float] = {}
    for row in rows:
        key = str(row["chromosome"]), str(row["gene_id"])
        if (
            row["component"] == "secure"
            and row["scope"] == "gene"
            and row["phase"] == "batch_total"
            and row["measurement_kind"] == "amortized"
        ):
            secure[key].append(float(row["elapsed_seconds"]))
        if (
            row["component"] == "r_reference"
            and row["scope"] == "gene"
            and row["phase"] == "r_gene_total"
        ):
            reference[key] = float(row["elapsed_seconds"])

    for chromosome, gene_id in sorted(
        set(secure) | set(reference),
        key=lambda key: (chromosome_number(key[0]), key[1]),
    ):
        secure_seconds = max(secure.get((chromosome, gene_id), [0.0]))
        reference_seconds = reference.get((chromosome, gene_id), 0.0)
        print(
            f"{chromosome:>10} {gene_id} "
            f"{secure_seconds:.6f} {reference_seconds:.6f}"
        )

    run_total = next(
        (
            float(row["elapsed_seconds"])
            for row in rows
            if row["component"] == "coordinator"
            and row["phase"] == "run_total"
        ),
        None,
    )
    if run_total is not None:
        print(f"\nRUN WALL TIME: {run_total:.3f} seconds")


def timing_sort_key(row: dict[str, object]) -> tuple[object, ...]:
    return (
        component_order(str(row["component"])),
        chromosome_number(str(row["chromosome"])),
        integer_or_default(row["party"]),
        integer_or_default(row["gene_index"]),
        integer_or_default(row["phenotype_index"]),
        str(row["phase"]),
        str(row["measurement_kind"]),
        integer_or_default(row["event_index"]),
    )


def component_order(component: str) -> int:
    order = {
        "workflow": 0,
        "preprocessing": 1,
        "secure": 2,
        "r_reference": 3,
        "coordinator": 4,
    }
    return order.get(component, len(order))


def chromosome_number(value: str) -> int:
    normalized = value.lower()
    if normalized.startswith("chr"):
        normalized = normalized[3:]
    return int(normalized) if normalized else 10_000


def integer_or_default(value: object) -> int:
    return int(value) if str(value) else -1


def sum_optional_integers(values: Iterator[str]) -> str:
    present = [int(value) for value in values if value]
    return str(sum(present)) if present else ""
