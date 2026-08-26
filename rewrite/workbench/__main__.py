import argparse
from pathlib import Path

from .results import (
    SUMMARY_FIELDS,
    aggregate_results,
    build_error_summary,
    join_results,
    write_rows,
)


def join_command(args: argparse.Namespace) -> None:
    rows = join_results(args.secure, args.r_results, args.output)
    write_rows(args.summary, build_error_summary(rows, include_all=False), SUMMARY_FIELDS)


def aggregate_command(args: argparse.Namespace) -> None:
    from .plots import write_all_plots

    aggregate_results(args.input_root, args.output_dir)
    write_all_plots(
        args.output_dir / "all_gene_results.csv",
        args.output_dir / "plots",
    )


def keys_command(args: argparse.Namespace) -> None:
    from .keys import write_shared_keys

    write_shared_keys(args.output_dir)


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
    aggregate.set_defaults(handler=aggregate_command)

    keys = subparsers.add_parser("keys")
    keys.add_argument("--output-dir", type=Path, required=True)
    keys.set_defaults(handler=keys_command)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.handler(args)


if __name__ == "__main__":
    main()
