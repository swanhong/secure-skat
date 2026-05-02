"""Scatter-plot rendering for block-level comparison outputs."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np


def draw_scatter_png(
    arg_plain_values: np.ndarray,
    arg_secure_values: np.ndarray,
    arg_out_path: Path,
    *,
    title: str,
    subtitle: str,
) -> bool:
    keep = np.isfinite(arg_plain_values) & np.isfinite(arg_secure_values)
    if int(np.sum(keep)) == 0:
        return False

    cache_root = arg_out_path.parent / ".plot-cache"
    cache_root.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_root / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_root))

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        raise RuntimeError(
            "matplotlib is required for scatter-plot output. "
            "Install scripts/analysis/requirements-skat-compare.txt first."
        ) from exc

    eps = 1e-6
    x = np.log10(np.maximum(np.asarray(arg_plain_values, dtype=float)[keep], eps))
    y = np.log10(np.maximum(np.asarray(arg_secure_values, dtype=float)[keep], eps))
    min_v = float(min(np.min(x), np.min(y)))
    max_v = float(max(np.max(x), np.max(y)))
    if min_v == max_v:
        min_v -= 1.0
        max_v += 1.0
    pad = 0.05 * (max_v - min_v)
    min_v -= pad
    max_v += pad

    fig, ax = plt.subplots(figsize=(9, 8), dpi=200)
    fig.subplots_adjust(left=0.16, right=0.96, bottom=0.14, top=0.86)
    ax.scatter(x, y, s=24, c=[(0.10, 0.35, 0.70, 0.72)], edgecolors="none")
    ax.grid(True, color="0.9", linewidth=0.8)
    ax.plot([min_v, max_v], [min_v, max_v], color="firebrick", linewidth=2, linestyle="--")
    ax.set_xlim(min_v, max_v)
    ax.set_ylim(min_v, max_v)
    ax.set_xlabel(f"log10(plain + {eps:.2e})", fontsize=12)
    ax.set_ylabel(f"log10(secure + {eps:.2e})", fontsize=12)
    ax.tick_params(axis="both", labelsize=11)
    fig.suptitle(title, fontsize=16, y=0.97)
    ax.set_title(subtitle, fontsize=11, pad=10)
    fig.savefig(arg_out_path, dpi=200)
    plt.close(fig)
    return True
