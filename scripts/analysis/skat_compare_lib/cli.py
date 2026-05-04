"""CLI entrypoint for the simplified SKAT comparison pipeline."""

from __future__ import annotations

import argparse
import sys
from typing import Sequence

from .pipeline import build_context, run_compare, run_reference_only


def add_shared_arguments(arg_parser: argparse.ArgumentParser) -> None:
    arg_parser.add_argument("--repo-root", default=".", help="Repository root (default: current directory).")
    arg_parser.add_argument("--run-id", required=True, help="Secure run id suffix, for example ca92.")
    arg_parser.add_argument("--dataset", help="Dataset root. If omitted, infer from run metadata or local datasets.")
    arg_parser.add_argument(
        "--blocks",
        help="Analysis block specification, for example 'all', '1,last', or '1-4,8'. Defaults to all blocks.",
    )


def build_parser() -> argparse.ArgumentParser:
    arg_parser = argparse.ArgumentParser(description="Plain-vs-secure SKAT block comparison pipeline.")
    arg_subparsers = arg_parser.add_subparsers(dest="command", required=True)

    arg_compare_parser = arg_subparsers.add_parser("compare", help="Run plain-vs-secure block comparison.")
    add_shared_arguments(arg_compare_parser)
    arg_compare_parser.add_argument(
        "--skip-reference",
        action="store_true",
        help="Skip the R SKAT package reference step.",
    )
    arg_compare_parser.add_argument(
        "--plain-mode",
        choices=["standard", "local-weight-burden"],
        default="standard",
        help=(
            "Manual/plain calculation mode. "
            "'standard' uses pooled weights for SKAT and burden; "
            "'local-weight-burden' keeps SKAT standard but computes burden from "
            "party-local weights and party-local partial sums."
        ),
    )
    arg_compare_parser.add_argument(
        "--local-weight-mode",
        choices=["direct-total", "product-approx"],
        default="direct-total",
        help=(
            "Experimental submode used only when --plain-mode local-weight-burden. "
            "'direct-total' applies the beta weight separately per party using the "
            "party-local numerator over the global 2N denominator. "
            "'product-approx' builds a shared approximate weight "
            "25 * product_p (1 - x_p)^24 from each party's local alt-frequency contribution x_p."
        ),
    )

    arg_reference_parser = arg_subparsers.add_parser("reference", help="Run only the R SKAT package reference.")
    add_shared_arguments(arg_reference_parser)

    return arg_parser


def main(arg_argv: Sequence[str] | None = None) -> int:
    try:
        arg_parser = build_parser()
        arg_ns = arg_parser.parse_args(arg_argv)
        ctx = build_context(arg_ns)

        if arg_ns.command == "compare":
            return run_compare(ctx)
        if arg_ns.command == "reference":
            return run_reference_only(ctx)
        raise RuntimeError(f"Unknown command: {arg_ns.command}")
    except RuntimeError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
