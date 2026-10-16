# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
"""Shared matplotlib style for the project's figures.

Static PNGs for the wiki: a light surface, recessive grid and axes, thin
lines, and a fixed categorical order (never cycled) from a palette that
passes the colour-vision checks for adjacent series. Three of the hues
sit below 3:1 contrast on the surface, so every series also gets its own
line style and marker, and every figure with more than one series has a
legend.
"""

from __future__ import annotations

from typing import Dict

SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#e4e3df"
NEUTRAL = "#8a8984"

# Categorical slots in fixed order: blue, orange, aqua, yellow, magenta,
# green, violet, red.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
LINESTYLES = ["-", "--", "-.", ":", (0, (5, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1)), (0, (1, 1)), (0, (8, 2))]
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]

# Diverging: blue (below the midpoint) to red (above), neutral grey midpoint.
DIVERGING = ["#104281", "#2a78d6", "#86b6ef", "#f0efec", "#f2a3a2", "#e34948", "#8f1f1f"]


def apply_style() -> None:
    """Set matplotlib rcParams for the project's figures."""
    import matplotlib

    matplotlib.use("Agg")
    matplotlib.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": NEUTRAL,
        "axes.labelcolor": TEXT_PRIMARY,
        "axes.titlecolor": TEXT_PRIMARY,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "axes.grid": True,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "xtick.color": TEXT_SECONDARY,
        "ytick.color": TEXT_SECONDARY,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.frameon": False,
        "legend.fontsize": 8.5,
        "lines.linewidth": 1.6,
        "lines.markersize": 6,
        "font.size": 10,
        "savefig.dpi": 130,
        "savefig.bbox": "tight",
    })


def series_style(index: int) -> Dict[str, object]:
    """Colour, line style and marker of categorical slot ``index`` (no cycling)."""
    if index >= len(SERIES):
        raise ValueError("more than {} series: fold into small multiples instead".format(len(SERIES)))
    return {"color": SERIES[index], "linestyle": LINESTYLES[index], "marker": MARKERS[index]}


def diverging_cmap():
    """Blue-grey-red colormap for a quantity with a meaningful midpoint."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list("project_diverging", DIVERGING)
