#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
"""Critical exponent gamma from both sides of the threshold.

Reads one or more ``bisection_history.json`` files and fits

    supercritical  ln M_AH  = gamma ln(A - A*) + c       (collapse runs)
    subcritical    ln K_max = -4 gamma ln(A* - A) + c'   (dispersion runs)

by least squares, where K_max is the largest value of the curvature
invariant during the run (a Kretschmann-type scalar, length^-4). Two
independent estimates of gamma that agree are the strongest evidence for
type II behaviour this project can produce (goal G3). Prints gamma, its
standard error and R^2 for each side and writes the data with the fitted
lines to ``critical_scaling.png``.

How A* is chosen matters more than the fit itself, because the runs
closest to the threshold have |A - A*| comparable to the final bracket
width:

* ``lower`` (default): the last subcritical bound a_min. Keeps every
  supercritical Delta A positive, but biases the points nearest A*.
* ``mid``: the midpoint of the final bracket.
* ``value``: an explicit ``--a-star``.
* ``fit``: fit A*, gamma and c together on each side (nonlinear least
  squares, needs at least four points). The recommended mode once a
  dedicated scan exists. ``--min-delta-factor`` does not apply.

Points with |A - A*| below ``--min-delta-factor`` times the final bracket
width are excluded, since their Delta A is dominated by the uncertainty
in A*. Measured exponents belong to this family and this numerical setup;
there is no universal value to compare against (Baumgarte et al. 2023).
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

A_STAR_MODES = ("lower", "mid", "value", "fit")
MASS_CHOICES = ("settled", "final", "first")
SIDES = ("super", "sub", "both")

# Exponent of the observable in units of gamma: M ~ dA^gamma, K ~ dA^(-4 gamma).
GAMMA_FACTOR = {"super": 1.0, "sub": -4.0}


@dataclass
class ScalingPoint:
    """One run used in a fit.

    Attributes:
        amplitude: Brill wave amplitude A.
        value: Horizon mass (supercritical) or peak curvature (subcritical).
        run_dir: Run directory, for traceability.
        side: ``super`` or ``sub``.
    """

    amplitude: float
    value: float
    run_dir: str
    side: str = "super"


@dataclass
class FitResult:
    """Result of a power-law fit on one side.

    Attributes:
        gamma: Critical exponent.
        gamma_stderr: Standard error of gamma.
        intercept: c in ln(value) = slope ln|A - A*| + c.
        r_squared: Coefficient of determination in log-log space.
        a_star: A* used (fitted in ``fit`` mode).
        n_points: Number of points in the fit.
        side: ``super`` or ``sub``.
        a_star_stderr: Standard error of A* (``fit`` mode only).
    """

    gamma: float
    gamma_stderr: float
    intercept: float
    r_squared: float
    a_star: float
    n_points: int
    side: str = "super"
    a_star_stderr: Optional[float] = None

    @property
    def slope(self) -> float:
        return GAMMA_FACTOR[self.side] * self.gamma


def _load(path: Path) -> dict:
    return json.loads(path.read_text())


def load_points(history_paths: List[Path], mass: str) -> List[ScalingPoint]:
    """Collect the supercritical runs from one or more history files.

    Args:
        history_paths: ``bisection_history.json`` files.
        mass: ``settled``, ``final`` or ``first`` horizon mass.

    Returns:
        One point per collapse run with a positive mass of that kind.
    """
    key = {"settled": "m_ah_settled", "final": "m_ah_final", "first": "m_ah_first"}[mass]
    points = []
    for path in history_paths:
        for run in _load(path)["runs"]:
            value = run.get(key)
            if run.get("outcome") == "collapse" and value is not None and value > 0.0:
                points.append(ScalingPoint(run["amplitude"], value, run.get("run_dir", ""), "super"))
    return points


def load_subcritical_points(history_paths: List[Path]) -> List[ScalingPoint]:
    """Collect the dispersion runs with their peak curvature."""
    points = []
    for path in history_paths:
        for run in _load(path)["runs"]:
            value = run.get("curvature_max")
            if run.get("outcome") == "dispersion" and value is not None and value > 0.0:
                points.append(ScalingPoint(run["amplitude"], value, run.get("run_dir", ""), "sub"))
    return points


def final_bracket(history_paths: List[Path]) -> Tuple[float, float]:
    """Final bracket of the first history file."""
    bracket = _load(history_paths[0])["bracket"]
    return bracket["a_min"], bracket["a_max"]


def resolve_a_star(history_paths: List[Path], mode: str, explicit: Optional[float]) -> float:
    """Determine A* for the ``lower``, ``mid`` and ``value`` modes.

    Args:
        history_paths: History files; the final bracket comes from the first.
        mode: One of ``A_STAR_MODES`` other than ``fit``.
        explicit: Value for ``value`` mode.

    Returns:
        A*.
    """
    if mode == "value":
        if explicit is None:
            raise ValueError("value mode needs an explicit A*")
        return explicit
    a_min, a_max = final_bracket(history_paths)
    if mode == "lower":
        return a_min
    if mode == "mid":
        return 0.5 * (a_min + a_max)
    raise ValueError("resolve_a_star does not handle mode '{}'".format(mode))


def select_points(points: List[ScalingPoint], a_star: float, min_delta: float) -> List[ScalingPoint]:
    """Keep the points on the right side of A* and at least ``min_delta`` from it."""
    kept = []
    for p in points:
        delta = p.amplitude - a_star if p.side == "super" else a_star - p.amplitude
        if delta > 0.0 and delta >= min_delta:
            kept.append(p)
    return kept


def _deltas(points: List[ScalingPoint], a_star: float) -> np.ndarray:
    a = np.array([p.amplitude for p in points])
    return a - a_star if points[0].side == "super" else a_star - a


def fit_power_law(points: List[ScalingPoint], a_star: float) -> FitResult:
    """Linear fit of ln(value) against ln|A - A*| on one side.

    Args:
        points: Points of one side, all on the correct side of A*.
        a_star: Critical amplitude.

    Returns:
        The fit; gamma is the slope (supercritical) or -slope/4 (subcritical).

    Raises:
        ValueError: Fewer than two points, mixed sides, or a point on the
            wrong side of A*.
    """
    if len(points) < 2:
        raise ValueError("need at least two points, got {}".format(len(points)))
    side = points[0].side
    if any(p.side != side for p in points):
        raise ValueError("points from both sides of the threshold")
    delta = _deltas(points, a_star)
    if np.any(delta <= 0.0):
        raise ValueError("a point lies on the wrong side of A* = {}".format(a_star))
    x = np.log(delta)
    y = np.log([p.value for p in points])
    slope, intercept = np.polyfit(x, y, 1)
    residuals = y - (slope * x + intercept)
    ss_res = float(np.sum(residuals ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    n = len(points)
    if n > 2:
        slope_stderr = math.sqrt(ss_res / (n - 2) / float(np.sum((x - x.mean()) ** 2)))
    else:
        slope_stderr = float("nan")
    factor = GAMMA_FACTOR[side]
    return FitResult(
        gamma=float(slope) / factor,
        gamma_stderr=slope_stderr / abs(factor),
        intercept=float(intercept),
        r_squared=1.0 - ss_res / ss_tot if ss_tot > 0.0 else 1.0,
        a_star=a_star,
        n_points=n,
        side=side,
    )


def fit_power_law_free_a_star(points: List[ScalingPoint], a_star_guess: float) -> FitResult:
    """Nonlinear fit of value = C |A - A*|^slope with A* as a free parameter.

    Args:
        points: Points of one side.
        a_star_guess: Starting value, e.g. the final bracket's midpoint.

    Returns:
        The fit, with ``a_star`` and its standard error fitted.

    Raises:
        ValueError: Fewer than four points.
    """
    from scipy.optimize import least_squares

    if len(points) < 4:
        raise ValueError("a free A* needs at least four points, got {}".format(len(points)))
    side = points[0].side
    amplitudes = np.array([p.amplitude for p in points])
    y = np.log([p.value for p in points])
    span = float(np.ptp(amplitudes)) or 1.0
    eps = 1e-12 * span
    if side == "super":
        bound = (-np.inf, float(amplitudes.min()) - eps)
        guess = min(a_star_guess, bound[1] - 1e-3 * span)
    else:
        bound = (float(amplitudes.max()) + eps, np.inf)
        guess = max(a_star_guess, bound[0] + 1e-3 * span)

    def model(params: np.ndarray) -> np.ndarray:
        a_star, slope, intercept = params
        delta = amplitudes - a_star if side == "super" else a_star - amplitudes
        return slope * np.log(delta) + intercept

    start = fit_power_law(points, guess)
    solution = least_squares(
        lambda params: model(params) - y,
        x0=[guess, start.slope, start.intercept],
        bounds=([bound[0], -np.inf, -np.inf], [bound[1], np.inf, np.inf]),
        x_scale=[span, 1.0, 1.0],
    )
    a_star, slope, intercept = solution.x
    residuals = solution.fun
    n = len(points)
    dof = max(n - 3, 1)
    sigma2 = float(np.sum(residuals ** 2)) / dof
    try:
        covariance = np.linalg.inv(solution.jac.T @ solution.jac) * sigma2
        a_star_stderr = math.sqrt(max(covariance[0, 0], 0.0))
        slope_stderr = math.sqrt(max(covariance[1, 1], 0.0))
    except np.linalg.LinAlgError:
        a_star_stderr = slope_stderr = float("nan")
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    factor = GAMMA_FACTOR[side]
    return FitResult(
        gamma=float(slope) / factor,
        gamma_stderr=slope_stderr / abs(factor),
        intercept=float(intercept),
        r_squared=1.0 - float(np.sum(residuals ** 2)) / ss_tot if ss_tot > 0.0 else 1.0,
        a_star=float(a_star),
        n_points=n,
        side=side,
        a_star_stderr=a_star_stderr,
    )


def plot_scaling(points: Dict[str, List[ScalingPoint]], fits: Dict[str, FitResult], output: Path) -> None:
    """Plot ln(value) against ln|A - A*| with the fitted lines, one panel per side.

    Args:
        points: Points used in each fit, keyed by side.
        fits: Fit results, keyed by side.
        output: PNG file to write.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sides = [s for s in ("super", "sub") if s in fits]
    fig, axes = plt.subplots(1, len(sides), figsize=(5.5 * len(sides), 4.5), squeeze=False)
    labels = {"super": ("A - A*", "M_AH"), "sub": ("A* - A", "K_max")}
    for ax, side in zip(axes[0], sides):
        fit = fits[side]
        delta = _deltas(points[side], fit.a_star)
        values = np.array([p.value for p in points[side]])
        ax.loglog(delta, values, "o", label="runs")
        grid = np.logspace(np.log10(delta.min()), np.log10(delta.max()), 50)
        ax.loglog(grid, np.exp(fit.intercept) * grid ** fit.slope, "-",
                  label="gamma = {:.3f} +/- {:.3f}".format(fit.gamma, fit.gamma_stderr))
        ax.set_xlabel(labels[side][0])
        ax.set_ylabel(labels[side][1])
        ax.set_title("{} (A* = {:.6g}, n = {})".format(
            "supercritical" if side == "super" else "subcritical", fit.a_star, fit.n_points))
        ax.grid(True, which="both", alpha=0.3)
        ax.legend()
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output), dpi=120)
    plt.close(fig)


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Parse the command line."""
    parser = argparse.ArgumentParser(
        description="Fit the critical exponent gamma from horizon masses and peak curvature."
    )
    parser.add_argument(
        "--history", type=Path, action="append", required=True,
        help="bisection_history.json; repeat to combine a bisection with later scans",
    )
    parser.add_argument(
        "--a-star-mode", choices=A_STAR_MODES, default="lower",
        help="how to choose A* (default: lower)",
    )
    parser.add_argument("--a-star", type=float, help="A* for --a-star-mode value")
    parser.add_argument(
        "--mass", choices=MASS_CHOICES, default="settled",
        help="horizon mass: settled, at the last find, or at the first find (default: settled)",
    )
    parser.add_argument(
        "--side", choices=SIDES, default="both",
        help="fit the supercritical masses, the subcritical curvature, or both (default: both)",
    )
    parser.add_argument(
        "--min-delta-factor", type=float, default=2.0,
        help="drop points with |A - A*| below this many final bracket widths (default: 2)",
    )
    parser.add_argument(
        "--output", type=Path, default=Path("critical_scaling.png"),
        help="output figure (default: critical_scaling.png)",
    )
    args = parser.parse_args(argv)
    if args.a_star_mode == "value" and args.a_star is None:
        parser.error("--a-star-mode value requires --a-star")
    return args


def main(argv: Optional[List[str]] = None) -> int:
    """Entry point."""
    args = parse_args(argv)
    a_min, a_max = final_bracket(args.history)
    width = a_max - a_min
    candidates = {}
    if args.side in ("super", "both"):
        candidates["super"] = load_points(args.history, args.mass)
    if args.side in ("sub", "both"):
        candidates["sub"] = load_subcritical_points(args.history)

    print("final bracket: [{}, {}] (width {:.3g})".format(a_min, a_max, width))
    fits: Dict[str, FitResult] = {}
    used: Dict[str, List[ScalingPoint]] = {}
    for side, points in candidates.items():
        try:
            if args.a_star_mode == "fit":
                selected = points
                fit = fit_power_law_free_a_star(selected, 0.5 * (a_min + a_max))
            else:
                a_star = resolve_a_star(args.history, args.a_star_mode, args.a_star)
                selected = select_points(points, a_star, args.min_delta_factor * width)
                fit = fit_power_law(selected, a_star)
        except ValueError as exc:
            print("{:<5}: no fit ({}; {} candidate runs)".format(side, exc, len(points)))
            continue
        fits[side], used[side] = fit, selected
        print("{:<5}: gamma = {:.4f} +/- {:.4f}, R^2 = {:.4f}, A* = {:.8g}{}, n = {}".format(
            side, fit.gamma, fit.gamma_stderr, fit.r_squared, fit.a_star,
            "" if fit.a_star_stderr is None else " +/- {:.2g}".format(fit.a_star_stderr), fit.n_points))
    if not fits:
        return 1
    if len(fits) == 2:
        difference = fits["super"].gamma - fits["sub"].gamma
        combined = math.hypot(fits["super"].gamma_stderr, fits["sub"].gamma_stderr)
        print("both : gamma_super - gamma_sub = {:.4f} +/- {:.4f}".format(difference, combined))
    plot_scaling(used, fits, args.output)
    print("wrote {}".format(args.output))
    return 0


if __name__ == "__main__":
    sys.exit(main())
