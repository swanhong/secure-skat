import argparse
import time
from pathlib import Path

from .results import (
    SUMMARY_FIELDS,
    aggregate_results,
    build_error_summary,
    join_results,
    write_rows,
)
from .timing import (
    aggregate_timing_events,
    append_timing_event,
    combine_timing_files,
    print_timing_summary,
)


def join_command(args: argparse.Namespace) -> None:
    rows = join_results(args.secure, args.r_results, args.output)
    write_rows(args.summary, build_error_summary(rows, include_all=False), SUMMARY_FIELDS)


def aggregate_command(args: argparse.Namespace) -> None:
    from .plots import write_all_plots

    aggregate_results(args.input_root, args.output_dir, args.chromosomes)
    write_all_plots(
        args.output_dir / "all_gene_results.csv",
        args.output_dir / "plots",
    )


def keys_command(args: argparse.Namespace) -> None:
    from .keys import write_shared_keys

    write_shared_keys(args.output_dir)


def timing_event_command(args: argparse.Namespace) -> None:
    elapsed_seconds = (time.perf_counter_ns() - args.started_ns) / 1e9
    append_timing_event(
        args.output,
        component=args.component,
        scope=args.scope,
        chromosome=args.chromosome,
        party=args.party,
        lane=args.lane,
        phase=args.phase,
        parent_phase=args.parent_phase,
        measurement_kind=args.measurement_kind,
        elapsed_seconds=f"{elapsed_seconds:.9f}",
        status=args.status,
    )


def timing_combine_command(args: argparse.Namespace) -> None:
    combine_timing_files(args.inputs, args.output)


def timing_summary_command(args: argparse.Namespace) -> None:
    summary = aggregate_timing_events(
        args.input_root,
        args.output_dir,
        args.chromosomes,
        args.coordinator_timing,
    )
    print_timing_summary(summary)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m rewrite.workbench")
    subparsers = parser.add_subparsers(dest="command", required=True)

    join = subparsers.add_parser("join")
    join.add_argument("--secure", type=Path, required=True)
    join.add_argument("--r-results", type=Path, required=True)
    join.add_argument("--output", type=Path, required=True)
    join.add_argument("--summary", type=Path, required=True)
    join.set_defaults(handler=join_command)

    aggregate = subparsers.add_parser("aggregate")
    aggregate.add_argument("--input-root", type=Path, required=True)
    aggregate.add_argument("--output-dir", type=Path, required=True)
    aggregate.add_argument(
        "--chromosomes",
        nargs="+",
        type=int,
        default=tuple(range(1, 23)),
    )
    aggregate.set_defaults(handler=aggregate_command)

    keys = subparsers.add_parser("keys")
    keys.add_argument("--output-dir", type=Path, required=True)
    keys.set_defaults(handler=keys_command)

    timing_event = subparsers.add_parser("timing-event")
    timing_event.add_argument("--output", type=Path, required=True)
    timing_event.add_argument("--started-ns", type=int, required=True)
    timing_event.add_argument("--component", required=True)
    timing_event.add_argument("--scope", required=True)
    timing_event.add_argument("--chromosome", default="")
    timing_event.add_argument("--party", default="")
    timing_event.add_argument("--lane", default="")
    timing_event.add_argument("--phase", required=True)
    timing_event.add_argument("--parent-phase", default="")
    timing_event.add_argument("--measurement-kind", default="actual")
    timing_event.add_argument("--status", default="success")
    timing_event.set_defaults(handler=timing_event_command)

    timing_combine = subparsers.add_parser("timing-combine")
    timing_combine.add_argument("--output", type=Path, required=True)
    timing_combine.add_argument("--inputs", type=Path, nargs="+", required=True)
    timing_combine.set_defaults(handler=timing_combine_command)

    timing_summary = subparsers.add_parser("timing-summary")
    timing_summary.add_argument("--input-root", type=Path, required=True)
    timing_summary.add_argument("--output-dir", type=Path, required=True)
    timing_summary.add_argument("--coordinator-timing", type=Path, required=True)
    timing_summary.add_argument(
        "--chromosomes",
        nargs="+",
        type=int,
        default=tuple(range(1, 23)),
    )
    timing_summary.set_defaults(handler=timing_summary_command)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.handler(args)


if __name__ == "__main__":
    main()
