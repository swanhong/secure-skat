"""Shared helpers for the SKAT comparison pipeline."""

from __future__ import annotations

import math
import re
import shutil
from pathlib import Path
from typing import Sequence

import numpy as np


def require_executable(name: str) -> str:
    path = shutil.which(name)
    if not path:
        raise RuntimeError(f"{name} not found in PATH")
    return path


def safe_rel_diff(x: float, y: float) -> float:
    if not (math.isfinite(x) and math.isfinite(y)):
        return float("nan")
    return abs(x - y) / max(abs(y), 1e-12)


def safe_corr(x: Sequence[float], y: Sequence[float]) -> float:
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    keep = np.isfinite(x_arr) & np.isfinite(y_arr)
    if int(np.sum(keep)) < 2:
        return float("nan")
    corr = np.corrcoef(x_arr[keep], y_arr[keep])[0, 1]
    return float(corr)


def format_float(value: float) -> str:
    return "NA" if not math.isfinite(value) else f"{value:.10e}"


def sanitize_path_tag(path: Path | str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(path))


def count_nonblank_lines(path: Path) -> int:
    return sum(1 for line in path.read_text().splitlines() if line.strip())


def read_kv_file(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    if not path.exists():
        return out
    for line in path.read_text().splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        out[key.strip()] = value.strip()
    return out


def trim_or_none(vec: np.ndarray | None, target_len: int) -> np.ndarray | None:
    if vec is None or vec.size < target_len:
        return None
    return np.asarray(vec[:target_len], dtype=float)

