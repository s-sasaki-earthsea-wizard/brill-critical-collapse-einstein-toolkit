#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
"""Critical exponent gamma from the mass scaling M_AH ~ (A - A*)^gamma.

Reads one or more ``bisection_history.json`` files, keeps the runs whose
outcome is ``collapse``, computes Delta A = A - A* and fits

    ln M_AH = gamma ln(Delta A) + c

by linear least squares. Prints gamma, its standard error and R^2, and
writes the data points with the fitted line to ``critical_scaling.png``.

How A* is chosen matters more than the fit itself, because the
supercritical points closest to the threshold have Delta A comparable to
the final bracket width:

* ``lower`` (default): the last subcritical bound a_min. Keeps every
  Delta A positive, but biases the points nearest A* upwards.
* ``mid``: the midpoint of the final bracket. Drops points with
  Delta A <= 0.
* ``value``: an explicit ``--a-star``.
* ``fit``: fit A*, gamma and c together (nonlinear least squares, needs
  at least four points). The recommended mode once a dedicated
  supercritical scan exists.

Points whose Delta A is below ``--min-delta-factor`` times the final
bracket width are excluded, since their Delta A is dominated by the
uncertainty in A*.

Status: skeleton. The CLI is fixed here; loading, fitting and plotting
are not implemented yet.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

A_STAR_MODES = ("lower", "mid", "value", "fit")
MASS_CHOICES = ("final", "first")


@dataclass
class ScalingPoint:
    """One supercritical run used in the fit.

    Attributes:
        amplitude: Brill wave amplitude A.
        mass: Apparent horizon mass (first or final find, see ``--mass``).
        run_dir: Run directory, for traceability.
    """

    amplitude: float
    mass: float
    run_dir: str


@dataclass
class FitResult:
    """Result of the power-law fit.

    Attributes:
        gamma: Critical exponent.
        gamma_stderr: Standard error of gamma.
        intercept: c in ln M = gamma ln(Delta A) + c.
        r_squared: Coefficient of determination of the linear fit.
        a_star: A* used for Delta A (fitted in ``fit`` mode).
        n_points: Number of points in the fit.
    """

    gamma: float
    gamma_stderr: float
    intercept: float
    r_squared: float
    a_star: float
    n_points: int


def load_points(history_paths: List[Path], mass: str) -> List[ScalingPoint]:
    """Collect the supercritical runs from one or more history files.

    Args:
        history_paths: ``bisection_history.json`` files.
        mass: ``final`` or ``first`` horizon mass.

    Returns:
        One point per collapse run with a positive mass.
    """
    raise NotImplementedError("load_points is not implemented yet")


def resolve_a_star(history_paths: List[Path], mode: str, explicit: Optional[float]) -> float:
    """Determine A* for the ``lower``, ``mid`` and ``value`` modes.

    Args:
        history_paths: History files; the final bracket comes from the first.
        mode: One of ``A_STAR_MODES`` other than ``fit``.
        explicit: Value for ``value`` mode.

    Returns:
        A*.
    """
    raise NotImplementedError("resolve_a_star is not implemented yet")


def fit_power_law(points: List[ScalingPoint], a_star: float) -> FitResult:
    """Linear fit of ln M against ln(A - A*).

    Args:
        points: Supercritical runs with A > A*.
        a_star: Critical amplitude.

    Returns:
        The fit.
    """
    raise NotImplementedError("fit_power_law is not implemented yet")


def fit_power_law_free_a_star(points: List[ScalingPoint], a_star_guess: float) -> FitResult:
    """Nonlinear fit of M = C (A - A*)^gamma with A* as a free parameter.

    Args:
        points: Supercritical runs.
        a_star_guess: Starting value, e.g. the final bracket's lower bound.

    Returns:
        The fit, with ``a_star`` fitted.
    """
    raise NotImplementedError("fit_power_law_free_a_star is not implemented yet")


def plot_scaling(points: List[ScalingPoint], fit: FitResult, output: Path) -> None:
    """Plot ln M against ln(Delta A) with the fitted line.

    Args:
        points: Points used in the fit.
        fit: Fit result.
        output: PNG file to write.
    """
    raise NotImplementedError("plot_scaling is not implemented yet")


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Parse the command line."""
    parser = argparse.ArgumentParser(
        description="Fit the critical exponent gamma from apparent horizon masses."
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
        "--mass", choices=MASS_CHOICES, default="final",
        help="horizon mass at the last or the first find (default: final)",
    )
    parser.add_argument(
        "--min-delta-factor", type=float, default=2.0,
        help="drop points with Delta A below this many final bracket widths (default: 2)",
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
    raise NotImplementedError("the analysis is not implemented yet")


if __name__ == "__main__":
    sys.exit(main())
