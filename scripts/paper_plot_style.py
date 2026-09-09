"""Shared matplotlib style for paper figures."""

from __future__ import annotations

import os


def apply_paper_figure_style(scale: float | None = None) -> float:
    """Apply large, print-readable matplotlib fonts.

    The default scale is 3.0 because the paper is now in a two-column layout and
    the raster figures are often reduced to column or page width.
    """

    import matplotlib as mpl

    if scale is None:
        scale = float(os.environ.get("PICOS_FIGURE_FONT_SCALE", "3.0"))

    mpl.rcParams.update(
        {
            "font.size": 10.0 * scale,
            "axes.titlesize": 11.0 * scale,
            "axes.labelsize": 10.0 * scale,
            "xtick.labelsize": 9.0 * scale,
            "ytick.labelsize": 9.0 * scale,
            "legend.fontsize": 8.0 * scale,
            "figure.titlesize": 12.0 * scale,
            "lines.linewidth": 1.5,
            "axes.linewidth": 1.2,
            "xtick.major.width": 1.2,
            "ytick.major.width": 1.2,
        }
    )
    return scale
