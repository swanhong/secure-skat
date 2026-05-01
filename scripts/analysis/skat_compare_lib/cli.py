"""CLI entrypoint for the modular SKAT comparison pipeline."""

from __future__ import annotations

import argparse
import sys
from typing import Sequence

from .context import build_context
from .workflow import run_compare, run_manual, run_reference_only, run_secure_only


def add_shared_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--repo-root", default=".", help="Repository root (default: current directory).")
    parser.add_argument("--run-id", help="Secure run id suffix, for example ca92.")
    parser.add_argument("--dataset", help="Dataset root. If omitted, infer from run metadata or local datasets.")
    parser.add_argument(
        "--blocks",
        help="Analysis block specification, for example 'all', '1,last', or '1-4,8'. Defaults to all blocks.",
    )
    parser.add_argument(
        "--detail-blocks",
        default="1,last",
        help="Blocks to print detailed diagnostics for (default: 1,last).",
    )
    parser.add_argument("--debug", action="store_true", help="Write per-variant debug CSVs and print extra diagnostics.")
    parser.add_argument("--window-bp", type=int, help="Optional window size in base pairs for sliding-window summaries.")
    parser.add_argument("--step-bp", type=int, help="Optional window step in base pairs. Defaults to --window-bp.")
    parser.add_argument("--min-window-variants", type=int, default=1, help="Minimum variants per window (default: 1).")
    parser.add_argument("--window-limit", type=int, help="Maximum number of windows to retain across all blocks.")
    parser.add_argument("--window-output-tag", help="Optional output tag for window CSV/plot filenames.")
    parser.add_argument(
        "--skip-reference",
        action="store_true",
        help="Skip the R SKAT package reference step.",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Python-first secure-vs-plain SKAT comparison pipeline.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    compare_parser = subparsers.add_parser("compare", help="Run reference, manual, and secure comparisons together.")
    add_shared_arguments(compare_parser)

    manual_parser = subparsers.add_parser("manual", help="Run only the manual plain recomputation.")
    add_shared_arguments(manual_parser)

    secure_parser = subparsers.add_parser("secure", help="Load only secure outputs for the selected analysis blocks.")
    add_shared_arguments(secure_parser)

    reference_parser = subparsers.add_parser("reference", help="Run only the R SKAT package reference.")
    add_shared_arguments(reference_parser)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    try:
        parser = build_parser()
        args = parser.parse_args(argv)
        if args.command == "secure" and args.run_id is None:
            raise RuntimeError("secure mode requires --run-id")

        ctx = build_context(args)

        if args.command == "compare":
            return run_compare(ctx)
        if args.command == "manual":
            return run_manual(ctx)
        if args.command == "secure":
            return run_secure_only(ctx)
        if args.command == "reference":
            return run_reference_only(ctx)
        raise RuntimeError(f"Unknown command: {args.command}")
    except RuntimeError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
