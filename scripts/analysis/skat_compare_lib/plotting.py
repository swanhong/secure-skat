"""PNG plotting helpers for comparison scatter plots."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont


def _load_truetype_font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Arial Unicode.ttf",
        "/System/Library/Fonts/SFNS.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ]
    for candidate in candidates:
        path = Path(candidate)
        if path.exists():
            try:
                return ImageFont.truetype(str(path), size=size)
            except OSError:
                continue
    return ImageFont.load_default()


def _compute_plot_range(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    min_v = float(min(np.min(x), np.min(y)))
    max_v = float(max(np.max(x), np.max(y)))
    if math.isclose(min_v, max_v):
        min_v -= 1.0
        max_v += 1.0
    pad = 0.05 * (max_v - min_v)
    return min_v - pad, max_v + pad


def _draw_scatter_with_matplotlib(
    x: np.ndarray,
    y: np.ndarray,
    out_path: Path,
    *,
    title: str,
    subtitle: str,
    eps: float,
) -> bool:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return False

    min_v, max_v = _compute_plot_range(x, y)
    fig, ax = plt.subplots(figsize=(9, 8), dpi=200)
    fig.subplots_adjust(left=0.16, right=0.96, bottom=0.14, top=0.86)

    ax.scatter(
        x,
        y,
        s=24,
        c=[(0.10, 0.35, 0.70, 0.72)],
        edgecolors="none",
    )
    ax.grid(True, color="0.9", linewidth=0.8)
    ax.plot(
        [min_v, max_v],
        [min_v, max_v],
        color="firebrick",
        linewidth=2,
        linestyle="--",
    )

    ax.set_xlim(min_v, max_v)
    ax.set_ylim(min_v, max_v)
    ax.set_xlabel(f"log10(plain + {eps:.2e})", fontsize=12)
    ax.set_ylabel(f"log10(secure + {eps:.2e})", fontsize=12)
    ax.tick_params(axis="both", labelsize=11)
    fig.suptitle(title, fontsize=16, y=0.97)
    ax.set_title(subtitle, fontsize=11, pad=10)

    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return True


def _draw_scatter_with_pillow(
    x: np.ndarray,
    y: np.ndarray,
    out_path: Path,
    *,
    title: str,
    subtitle: str,
    eps: float,
) -> bool:
    width, height = 1800, 1600
    margin_left, margin_right = 210, 70
    margin_top, margin_bottom = 150, 190

    image = Image.new("RGBA", (width, height), (255, 255, 255, 255))
    draw = ImageDraw.Draw(image, "RGBA")

    title_font = _load_truetype_font(34)
    subtitle_font = _load_truetype_font(24)
    axis_font = _load_truetype_font(24)
    tick_font = _load_truetype_font(20)

    min_v, max_v = _compute_plot_range(x, y)
    plot_left = margin_left
    plot_right = width - margin_right
    plot_top = margin_top
    plot_bottom = height - margin_bottom

    def project(px: float, py: float) -> tuple[int, int]:
        nx = (px - min_v) / (max_v - min_v)
        ny = (py - min_v) / (max_v - min_v)
        x_pix = int(plot_left + nx * (plot_right - plot_left))
        y_pix = int(plot_bottom - ny * (plot_bottom - plot_top))
        return x_pix, y_pix

    for i in range(6):
        frac = i / 5
        x_pix = int(plot_left + frac * (plot_right - plot_left))
        y_pix = int(plot_top + frac * (plot_bottom - plot_top))
        draw.line([(x_pix, plot_top), (x_pix, plot_bottom)], fill=(230, 230, 230, 255), width=2)
        draw.line([(plot_left, y_pix), (plot_right, y_pix)], fill=(230, 230, 230, 255), width=2)
        tick_value = min_v + frac * (max_v - min_v)
        tick_label = f"{tick_value:.2f}"
        bbox = draw.textbbox((0, 0), tick_label, font=tick_font)
        tick_w = bbox[2] - bbox[0]
        tick_h = bbox[3] - bbox[1]
        draw.text((x_pix - tick_w / 2, plot_bottom + 22), tick_label, fill="black", font=tick_font)
        draw.text((plot_left - tick_w - 24, y_pix - tick_h / 2), tick_label, fill="black", font=tick_font)

    draw.rectangle([(plot_left, plot_top), (plot_right, plot_bottom)], outline="black", width=3)

    diag_start = project(min_v, min_v)
    diag_end = project(max_v, max_v)
    dash_len = 16
    gap_len = 10
    dx = diag_end[0] - diag_start[0]
    dy = diag_end[1] - diag_start[1]
    dist = math.hypot(dx, dy)
    if dist > 0:
        ux = dx / dist
        uy = dy / dist
        pos = 0.0
        while pos < dist:
            seg_end = min(dist, pos + dash_len)
            start = (diag_start[0] + ux * pos, diag_start[1] + uy * pos)
            end = (diag_start[0] + ux * seg_end, diag_start[1] + uy * seg_end)
            draw.line([start, end], fill=(178, 34, 34, 255), width=4)
            pos += dash_len + gap_len

    for px, py in zip(x, y):
        point = project(float(px), float(py))
        draw.ellipse(
            [(point[0] - 6, point[1] - 6), (point[0] + 6, point[1] + 6)],
            fill=(26, 89, 179, 190),
            outline=None,
        )

    draw.text((plot_left, 24), title, fill="black", font=title_font)
    draw.text((plot_left, 68), subtitle, fill="black", font=subtitle_font)

    x_label = f"log10(plain + {eps:.2e})"
    x_bbox = draw.textbbox((0, 0), x_label, font=axis_font)
    x_w = x_bbox[2] - x_bbox[0]
    draw.text(((plot_left + plot_right - x_w) / 2, height - 60), x_label, fill="black", font=axis_font)

    y_label = f"log10(secure + {eps:.2e})"
    y_bbox = draw.textbbox((0, 0), y_label, font=axis_font)
    y_w = y_bbox[2] - y_bbox[0]
    y_h = y_bbox[3] - y_bbox[1]
    y_image = Image.new("RGBA", (y_w + 10, y_h + 10), (255, 255, 255, 0))
    y_draw = ImageDraw.Draw(y_image)
    y_draw.text((5, 5), y_label, fill="black", font=axis_font)
    y_rotated = y_image.rotate(90, expand=True)
    image.alpha_composite(y_rotated, (40, int((plot_top + plot_bottom - y_rotated.height) / 2)))

    image.convert("RGB").save(out_path, dpi=(200, 200))
    return True


def draw_scatter_png(
    plain_values: np.ndarray,
    secure_values: np.ndarray,
    out_path: Path,
    *,
    title: str,
    subtitle: str,
) -> bool:
    keep = np.isfinite(plain_values) & np.isfinite(secure_values)
    if int(np.sum(keep)) == 0:
        return False

    eps = 1e-6
    x = np.log10(np.maximum(plain_values[keep], eps))
    y = np.log10(np.maximum(secure_values[keep], eps))

    if _draw_scatter_with_matplotlib(x, y, out_path, title=title, subtitle=subtitle, eps=eps):
        return True
    return _draw_scatter_with_pillow(x, y, out_path, title=title, subtitle=subtitle, eps=eps)
