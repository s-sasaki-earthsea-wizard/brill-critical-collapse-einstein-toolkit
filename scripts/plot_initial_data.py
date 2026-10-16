#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
"""Figures of the Brill wave initial data.

``mass``  ADM mass against amplitude, with the published values, and the
          small-amplitude behaviour M ~ A^2 on log axes.
``psi``   The conformal factor psi(rho, z) of one amplitude, and its
          profiles along the symmetry axis and the equator.

Usage:
    plot_initial_data.py mass --out initial_data_mass.png
    plot_initial_data.py psi -A 5.0 --out initial_data_psi_A5.png
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import List, Optional

import numpy as np

import brill_initial_data as bid
import plot_style as ps

# Published ADM masses of the centred Holz family (sigma = 1, e^{2q}).
PUBLISHED = [
    (4.6, 0.597, "Hilditch et al. 2017"),
    (4.69, 0.619, "Hilditch et al. 2017"),
    (4.697, 0.621, "Hilditch et al. 2017"),
    (4.7, 0.622, "Hilditch et al. 2017"),
    (6.073, 1.02, "Hilditch et al. 2017"),
    (5.0, 0.70, "Khirnov & Ledvinka 2019"),
    (5.5, 0.84, "Khirnov & Ledvinka 2019"),
]
A_STAR = 4.6967


def adm_mass(amplitude: float, grid: bid.Grid) -> float:
    return bid.adm_mass(bid.solve_brill(bid.BrillParams(amplitude), grid)).flux


def plot_mass(out: Path, n_s: int, n_theta: int) -> None:
    """ADM mass against amplitude, and M ~ A^2 at small amplitude."""
    import matplotlib.pyplot as plt

    grid = bid.Grid(n_s, n_theta, 2.0)
    positive = np.linspace(0.0, 6.5, 66)
    negative = np.linspace(-6.0, 0.0, 61)
    m_pos = [adm_mass(a, grid) for a in positive]
    m_neg = [adm_mass(a, grid) for a in negative]

    fig, (left, right) = plt.subplots(1, 2, figsize=(10.5, 4.2))
    style = ps.series_style(0)
    left.plot(np.concatenate([negative, positive[1:]]), np.concatenate([m_neg, m_pos[1:]]),
              color=style["color"], linestyle="-", label="this project (6th-order solver)")
    for k, source in enumerate(("Hilditch et al. 2017", "Khirnov & Ledvinka 2019")):
        st = ps.series_style(k + 1)
        points = [(a, m) for a, m, s in PUBLISHED if s == source]
        left.plot([a for a, _ in points], [m for _, m in points], linestyle="none",
                  marker=st["marker"], markersize=7, markerfacecolor="none", markeredgewidth=1.6,
                  color=st["color"], label=source)
    left.axvline(A_STAR, color=ps.NEUTRAL, linewidth=1.0, linestyle=":")
    left.annotate("A* = 4.6967", xy=(A_STAR, 2.6), xytext=(A_STAR - 0.15, 2.6), ha="right",
                  color=ps.TEXT_SECONDARY, fontsize=8.5)
    left.set_ylim(0.0, 3.4)
    left.set_xlabel("amplitude A")
    left.set_ylabel("ADM mass M")
    left.set_title("ADM mass of the Holz-type Brill wave (sigma = 1)")
    left.legend(loc="upper center")

    small = np.logspace(-2, math.log10(6.0), 30)
    m_small_pos = [adm_mass(a, grid) for a in small]
    m_small_neg = [adm_mass(-a, grid) for a in small]
    for k, (values, label) in enumerate(((m_small_pos, "A > 0 (prolate)"), (m_small_neg, "A < 0 (oblate)"))):
        st = ps.series_style(k)
        right.loglog(small, values, color=st["color"], linestyle=st["linestyle"], label=label)
    right.loglog(small, 0.0372 * small ** 2, color=ps.NEUTRAL, linewidth=1.0, linestyle=":",
                 label="0.0372 A^2")
    right.set_xlabel("|A|")
    right.set_ylabel("ADM mass M")
    right.set_title("Small amplitude: M ~ A^2")
    right.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(str(out))
    plt.close(fig)


def plot_psi(out: Path, amplitude: float, n_s: int, n_theta: int, extent: float, profile_extent: float) -> None:
    """psi(rho, z) map and its axis and equator profiles."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    grid = bid.Grid(n_s, n_theta, 2.0)
    solution = bid.solve_brill(bid.BrillParams(amplitude), grid)
    u, scale = solution.u, grid.scale

    coords = np.linspace(0.0, extent, 241)
    rho, z = np.meshgrid(coords, coords, indexing="xy")
    r = np.hypot(rho, z)
    theta = np.arctan2(rho, z)
    psi = 1.0 + bid.evaluate_table(u, scale, r, theta)

    fig, (left, right) = plt.subplots(1, 2, figsize=(10.5, 4.4), gridspec_kw={"width_ratios": [1.0, 1.15]})
    norm = TwoSlopeNorm(vcenter=1.0, vmin=min(float(psi.min()), 0.999), vmax=max(float(psi.max()), 1.001))
    image = left.pcolormesh(rho, z, psi, cmap=ps.diverging_cmap(), norm=norm, shading="auto")
    contours = left.contour(rho, z, psi, levels=8, colors=ps.TEXT_SECONDARY, linewidths=0.6)
    left.clabel(contours, fontsize=7, fmt="%.2f")
    left.set_aspect("equal")
    left.grid(False)
    left.set_xlabel("rho")
    left.set_ylabel("z (symmetry axis)")
    left.set_title("Conformal factor psi, A = {:g}".format(amplitude))
    fig.colorbar(image, ax=left, fraction=0.046, pad=0.04, label="psi (1 = flat)")

    line = np.linspace(1e-6, profile_extent, 600)
    profiles = (
        (bid.evaluate_table(u, scale, line, np.zeros_like(line)), "along the axis (rho = 0)"),
        (bid.evaluate_table(u, scale, line, np.full_like(line, 0.5 * math.pi)), "along the equator (z = 0)"),
    )
    for k, (values, label) in enumerate(profiles):
        st = ps.series_style(k)
        right.plot(line, 1.0 + values, color=st["color"], linestyle=st["linestyle"], label=label)
    masses = bid.adm_mass(solution)
    far = line > 3.0
    right.plot(line[far], 1.0 + masses.flux / (2.0 * line[far]), color=ps.NEUTRAL,
               linewidth=1.0, linestyle=":", label="1 + M/(2r), M = {:.4f}".format(masses.flux))
    right.axhline(1.0, color=ps.NEUTRAL, linewidth=0.8)
    right.set_xlabel("distance from the origin")
    right.set_ylabel("psi")
    right.axvspan(0.0, extent, color=ps.GRID, alpha=0.5, linewidth=0)
    right.text(0.5 * extent, 0.97, "map range", ha="center", va="top", color=ps.TEXT_SECONDARY, fontsize=8,
               transform=right.get_xaxis_transform())
    right.set_title("Profiles: below 1 inside the wave, 1 + M/2r far away")
    right.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(str(out))
    plt.close(fig)


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("figure", choices=("mass", "psi"))
    parser.add_argument("--out", type=Path, required=True, help="output PNG")
    parser.add_argument("--amplitude", "-A", type=float, default=5.0, help="amplitude for psi (default: 5)")
    parser.add_argument("--extent", type=float, default=3.0, help="plotted range in rho and z for the psi map")
    parser.add_argument("--profile-extent", type=float, default=12.0, help="range of the psi profiles")
    parser.add_argument("--n-s", type=int, default=256)
    parser.add_argument("--n-theta", type=int, default=32)
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    ps.apply_style()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if args.figure == "mass":
        plot_mass(args.out, args.n_s, args.n_theta)
    else:
        plot_psi(args.out, args.amplitude, args.n_s, args.n_theta, args.extent, args.profile_extent)
    print("wrote {}".format(args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
