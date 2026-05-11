from __future__ import annotations

import argparse
from pathlib import Path

from utils import run


def cloud_billing_args(args: argparse.Namespace) -> list[str]:
    if not args.billing_project:
        return []
    if args.cloud_cli == "gcloud":
        return [f"--billing-project={args.billing_project}"]
    return ["-u", args.billing_project]


def cloud_ls(uri: str, args: argparse.Namespace) -> bool:
    if args.cloud_cli == "gcloud":
        cmd = ["gcloud", "storage", "ls", uri, *cloud_billing_args(args)]
    else:
        cmd = ["gsutil", *cloud_billing_args(args), "ls", uri]
    proc = run(cmd, check=False, capture=True)
    return proc.returncode == 0


def cloud_cp(src: str, dst: Path | str, args: argparse.Namespace) -> None:
    if args.cloud_cli == "gcloud":
        cmd = ["gcloud", "storage", "cp", src, str(dst), *cloud_billing_args(args)]
    else:
        cmd = ["gsutil", *cloud_billing_args(args), "cp", src, str(dst)]
    run(cmd)


def cloud_cp_recursive(src: Path | str, dst: str, args: argparse.Namespace) -> None:
    if args.cloud_cli == "gcloud":
        cmd = ["gcloud", "storage", "cp", "-r", str(src), dst, *cloud_billing_args(args)]
    else:
        cmd = ["gsutil", *cloud_billing_args(args), "-m", "cp", "-r", str(src), dst]
    run(cmd)
