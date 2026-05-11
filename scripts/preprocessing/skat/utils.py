from __future__ import annotations

import csv
import math
import shutil
import subprocess
from pathlib import Path


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def resolve_path(path_text: str) -> Path:
    return Path(path_text).expanduser().resolve()


def pfile_path(prefix: Path, ext: str) -> Path:
    """Return a PLINK component path without treating dots in the prefix as suffixes."""
    return Path(str(prefix) + ext)


def is_gcs_uri(path_text: str) -> bool:
    return path_text.startswith("gs://")


def run(cmd: list[str], *, check: bool = True, capture: bool = False) -> subprocess.CompletedProcess[str]:
    printable = " ".join(cmd)
    print(f"+ {printable}", flush=True)
    return subprocess.run(
        cmd,
        check=check,
        text=True,
        capture_output=capture,
    )


def ensure_clean_dir(path: Path, force: bool) -> None:
    if path.exists():
        if not force:
            raise FileExistsError(f"{path} already exists; rerun with --force to overwrite")
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def copy_local_file(src: Path, dst: Path) -> None:
    if src.resolve() == dst.resolve():
        return
    shutil.copy2(src, dst)


def basename_from_path_text(path_text: str, default: str) -> str:
    name = Path(path_text.rstrip("/")).name
    return name if name else default


def split_table_line(line: str, sep: str | None) -> list[str]:
    line = line.rstrip("\n")
    if sep is not None:
        return next(csv.reader([line], delimiter=sep))
    if "\t" in line:
        return line.split("\t")
    if "," in line:
        return next(csv.reader([line]))
    return line.split()


def parse_float(value: object) -> float:
    try:
        out = float(str(value).strip())
    except ValueError:
        return float("nan")
    return out if math.isfinite(out) else float("nan")
