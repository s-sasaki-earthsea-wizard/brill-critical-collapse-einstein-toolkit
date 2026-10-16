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

Runs are given as LABEL=RUN_DIR, where RUN_DIR holds run_summary.json,
the parfile and output/.

Usage:
    plot_runs.py horizon --run "kappa=1=/path/run" ... --out horizon.png
    plot_runs.py histories --run ... --out histories.png
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

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


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("figure", choices=("horizon", "histories"))
    parser.add_argument("--run", type=parse_run, action="append", required=True, metavar="LABEL=RUN_DIR")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--title", default="")
    parser.add_argument("--curvature-cap", type=float, default=1e6,
                        help="upper limit of the curvature axis in the histories figure (default: 1e6)")
    parser.add_argument("--band", metavar="LOW,HIGH,LABEL",
                        help="shade a reference mass range in the horizon figure")
    args = parser.parse_args(argv)
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
    else:
        plot_histories(args.run, args.out, args.title, args.curvature_cap)
    print("wrote {}".format(args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
