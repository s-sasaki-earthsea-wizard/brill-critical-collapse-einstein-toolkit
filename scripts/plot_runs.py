#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
"""Compare finished runs in figures.

``horizon``    Irreducible mass of each run's horizon against time. A cross
               marks a crash, a dot the time the mass counts as settled
               under the current classifier.
``histories``  Global minimum of the lapse and global maximum of the
               curvature invariant |Re I| against time, one panel each
               (log scale), a cross at the last value before a crash.
``profiles``   One variable of one run along a coordinate axis at several
               times, from the CarpetIOASCII 1D output (finest level at
               each point).

Runs are given as LABEL=RUN_DIR, where RUN_DIR holds run_summary.json,
the parfile and output/.

Usage:
    plot_runs.py horizon --run "kappa=1=/path/run" ... --out horizon.png
    plot_runs.py histories --run ... --out histories.png
    plot_runs.py profiles --run "A=4=/path/run" --file admbase-shift.x.asc \
        --var betax --times 4.2,4.6,5.0 --out shift.png
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import cactus_outputs as co
import outcome_classifier as oc
import plot_style as ps


@dataclass
class Run:
    """A finished run to plot.

    Attributes:
        label: Legend label.
        path: Run directory.
        summary: Contents of run_summary.json.
    """

    label: str
    path: Path
    summary: dict

    @property
    def output(self) -> Path:
        return self.path / "output"

    @property
    def crashed(self) -> bool:
        reasons = " ".join(self.summary.get("reasons", []))
        return self.summary.get("outcome") == oc.CRASH or "crashed afterwards" in reasons


def parse_run(spec: str) -> Run:
    """Parse LABEL=RUN_DIR (the label may contain '=')."""
    label, _, path = spec.rpartition("=")
    if not label or not path:
        raise argparse.ArgumentTypeError("expected LABEL=RUN_DIR, got {!r}".format(spec))
    run_dir = Path(path)
    summary_path = run_dir / "run_summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    return Run(label, run_dir, summary)


def finite(series: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
    return [(t, v) for t, v in series if math.isfinite(v)]


def horizon_track(run: Run) -> Optional[co.HorizonSeries]:
    """The confirmed horizon, else the finder with the largest last mass."""
    number = run.summary.get("horizon_number")
    if number is None:
        candidates = [h for h in run.summary.get("horizons", []) if h.get("n_finds")]
        if not candidates:
            return None
        number = max(candidates, key=lambda h: h["m_last"])["number"]
    return co.read_bh_diagnostics(run.output, number)


def settle_point(series: co.HorizonSeries, config: oc.ClassifierConfig) -> Optional[Tuple[float, float]]:
    """First find at which the mass counts as settled under ``config``."""
    for k in range(2, len(series.times)):
        mass = oc.settled_mass(series.times[: k + 1], series.irreducible_masses[: k + 1], config)
        if mass is not None:
            return series.times[k], mass
    return None


def plot_horizon(runs: List[Run], out: Path, title: str, band: Optional[Tuple[float, float, str]]) -> None:
    import matplotlib.pyplot as plt

    config = oc.ClassifierConfig()
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    if band is not None:
        ax.axhspan(band[0], band[1], color=ps.GRID, alpha=0.7, linewidth=0, label=band[2])
    for k, run in enumerate(runs):
        st = ps.series_style(k)
        series = horizon_track(run)
        if series is None:
            continue
        ax.plot(series.times, series.irreducible_masses, color=st["color"], linestyle=st["linestyle"],
                label=run.label)
        point = settle_point(series, config)
        if point is not None:
            ax.plot([point[0]], [point[1]], linestyle="none", marker="o", markersize=7, color=st["color"],
                    markeredgecolor=ps.SURFACE, markeredgewidth=1.2)
        if run.crashed:
            ax.plot([series.times[-1]], [series.irreducible_masses[-1]], linestyle="none", marker="x",
                    markersize=9, markeredgewidth=2.0, color=st["color"])
    ax.plot([], [], linestyle="none", marker="o", color=ps.TEXT_SECONDARY, label="mass settled (rate criterion)")
    ax.plot([], [], linestyle="none", marker="x", markeredgewidth=2.0, color=ps.TEXT_SECONDARY, label="crash")
    ax.set_xlabel("time t")
    ax.set_ylabel("irreducible mass M_irr of the apparent horizon")
    ax.set_title(title)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(str(out))
    plt.close(fig)


def plot_histories(runs: List[Run], out: Path, title: str, curvature_cap: float) -> None:
    import matplotlib.pyplot as plt

    fig, (top, bottom) = plt.subplots(2, 1, figsize=(7.5, 7.2), sharex=True)
    for k, run in enumerate(runs):
        st = ps.series_style(k)
        panels = (
            (top, finite(co.read_scalar_series(run.output, co.MIN_LAPSE_FILENAME))),
            (bottom, [(t, abs(v)) for t, v in finite(co.read_scalar_series(run.output, co.CURVATURE_NORM_INF_FILENAME))]),
        )
        for ax, series in panels:
            if not series:
                continue
            ax.semilogy([t for t, _ in series], [v for _, v in series], color=st["color"],
                        linestyle=st["linestyle"], label=run.label if ax is top else None)
            if run.crashed:
                last = series[-1][1] if ax is top else min(series[-1][1], 0.7 * curvature_cap)
                ax.plot([series[-1][0]], [last], linestyle="none", marker="x", markersize=9,
                        markeredgewidth=2.0, color=st["color"])
    top.plot([], [], linestyle="none", marker="x", markeredgewidth=2.0, color=ps.TEXT_SECONDARY, label="crash")
    top.set_ylabel("minimum lapse")
    top.set_title(title)
    bottom.set_ylim(top=curvature_cap)
    bottom.set_ylabel("maximum |Re I|  (curvature)")
    bottom.set_xlabel("time t")
    bottom.set_title("curvature above {:.0e} clipped (it diverges at a crash)".format(curvature_cap),
                     fontsize=9, color=ps.TEXT_SECONDARY)
    handles, labels = top.get_legend_handles_labels()
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.legend(handles, labels, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.0))
    fig.savefig(str(out))
    plt.close(fig)


DATA_COLUMNS_PATTERN = re.compile(r"^#\s*data columns:\s*(.+)$")
# Coordinate columns of CarpetIOASCII 1D output (0-based): x, y, z.
AXIS_COLUMN = {"x": 9, "y": 10, "z": 11}
TIME_COLUMN = 8
LEVEL_COLUMN = 2


def read_profiles(path: Path, var: str, times: List[float]) -> Dict[float, List[Tuple[float, float]]]:
    """Profiles of ``var`` along the file's axis at ``times``.

    Where refinement levels overlap, the finest level's value is kept.

    Args:
        path: CarpetIOASCII 1D file, e.g. ``admbase-shift.x.asc``.
        var: Data column name from the header, e.g. ``betax``.
        times: Output times to extract (matched to 1e-6).

    Returns:
        (coordinate, value) pairs sorted by coordinate, per requested time.
    """
    axis = path.name.split(".")[-2]
    if axis not in AXIS_COLUMN:
        raise ValueError("cannot tell the axis of {}".format(path.name))
    column = None
    finest: Dict[float, Dict[float, Tuple[int, float]]] = {t: {} for t in times}
    with path.open() as stream:
        for line in stream:
            if line.startswith("#"):
                match = DATA_COLUMNS_PATTERN.match(line)
                if match and column is None:
                    names = dict(item.split(":")[::-1] for item in match.group(1).split())
                    if var not in names:
                        raise ValueError("{} has no column {!r} (has {})".format(path.name, var, sorted(names)))
                    column = int(names[var]) - 1
                continue
            fields = line.split()
            if column is None or len(fields) <= column:
                continue
            time = float(fields[TIME_COLUMN])
            target = next((t for t in times if abs(t - time) < 1e-6), None)
            if target is None:
                continue
            coordinate = float(fields[AXIS_COLUMN[axis]])
            level = int(fields[LEVEL_COLUMN])
            kept = finest[target].get(coordinate)
            if kept is None or level > kept[0]:
                finest[target][coordinate] = (level, float(fields[column]))
    return {t: sorted((c, v) for c, (_, v) in points.items()) for t, points in finest.items()}


def plot_profiles(run: Run, path: Path, var: str, times: List[float], out: Path, title: str,
                  max_coordinate: Optional[float], log: bool, ylabel: str) -> None:
    import matplotlib.pyplot as plt

    profiles = read_profiles(path, var, times)
    axis = path.name.split(".")[-2]
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    for k, t in enumerate(times):
        points = [(c, v) for c, v in profiles[t] if max_coordinate is None or c <= max_coordinate + 1e-9]
        if not points:
            raise ValueError("no {} output at t = {} in {}".format(var, t, path))
        values = [abs(v) if log else v for _, v in points]
        st = ps.sequential_style(k, len(times))
        draw = ax.semilogy if log else ax.plot
        draw([c for c, _ in points], values, color=st["color"], linestyle=st["linestyle"], label="t = {:g}".format(t))
    if not log:
        ax.axhline(0.0, color=ps.NEUTRAL, linewidth=0.8)
    ax.set_xlabel({"x": "x on the equator (z = 0)", "z": "z on the symmetry axis (x = 0)"}.get(axis, axis))
    ax.set_ylabel(ylabel or ("|{}|".format(var) if log else var))
    ax.set_title(title or "{}: {}".format(run.label, var))
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(str(out))
    plt.close(fig)


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("figure", choices=("horizon", "histories", "profiles"))
    parser.add_argument("--run", type=parse_run, action="append", required=True, metavar="LABEL=RUN_DIR")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--title", default="")
    parser.add_argument("--curvature-cap", type=float, default=1e6,
                        help="upper limit of the curvature axis in the histories figure (default: 1e6)")
    parser.add_argument("--band", metavar="LOW,HIGH,LABEL",
                        help="shade a reference mass range in the horizon figure")
    profiles = parser.add_argument_group("profiles")
    profiles.add_argument("--file", help="1D output file under the run's output/, e.g. admbase-shift.x.asc")
    profiles.add_argument("--var", help="data column, e.g. betax")
    profiles.add_argument("--times", help="comma-separated output times, earliest first")
    profiles.add_argument("--max-coordinate", type=float, help="right end of the plotted range")
    profiles.add_argument("--log", action="store_true", help="plot |value| on a log axis")
    profiles.add_argument("--ylabel", default="")
    args = parser.parse_args(argv)
    if args.figure == "profiles":
        if len(args.run) != 1 or not (args.file and args.var and args.times):
            parser.error("profiles needs one --run, --file, --var and --times")
        args.times = [float(t) for t in args.times.split(",")]
    if args.band:
        low, high, label = args.band.split(",", 2)
        args.band = (float(low), float(high), label)
    return args


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    ps.apply_style()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if args.figure == "horizon":
        plot_horizon(args.run, args.out, args.title, args.band)
    elif args.figure == "profiles":
        run = args.run[0]
        plot_profiles(run, run.output / args.file, args.var, args.times, args.out, args.title,
                      args.max_coordinate, args.log, args.ylabel)
    else:
        plot_histories(args.run, args.out, args.title, args.curvature_cap)
    print("wrote {}".format(args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
