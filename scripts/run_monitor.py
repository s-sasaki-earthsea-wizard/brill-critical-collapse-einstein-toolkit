#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
"""Resource monitoring and ETA estimation for Cactus runs.

Both are mandatory for every run in this project.

Memory (OOM protection)
    A background thread samples the container's cgroup v2 accounting
    (``memory.current``, ``memory.peak``, ``memory.max`` and the
    ``oom_kill`` counter in ``memory.events``) at a fixed interval and
    appends each sample to ``<run_dir>/resources.csv``. When
    ``memory.current`` exceeds ``guard_fraction * memory.max`` the guard
    trips: the runner asks Cactus to stop through the TerminationTrigger
    file, escalates to SIGTERM and then SIGKILL on mpirun, and classifies
    the run as ``oom``. The container itself has swap disabled
    (``memswap_limit = mem_limit``), so the kernel OOM killer is the
    backstop when the guard is too slow; an increment of ``oom_kill``
    also marks the run as ``oom``.

ETA
    Per run: the runner feeds every stdout line to ``EtaEstimator``,
    which picks up the ``(iteration, time)`` pairs from the CarpetIOBasic
    table, timestamps them with the wall clock, and estimates the
    simulation-time rate over a sliding window. Run ETA is
    ``(t_final - t_now) / rate``.

    Per campaign: remaining bisection steps times the median wall time of
    the completed runs, plus the current run's ETA. Collapse runs may be
    stopped early once the horizon mass settles, so the median is a
    conservative estimate.

Both estimates go into ``<campaign>/status.json``, which ``make status``
prints through the ``status`` subcommand of this script.

Status: skeleton. The interfaces are fixed here; the bodies are not
implemented yet.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

CGROUP_ROOT = Path("/sys/fs/cgroup")
STATUS_FILENAME = "status.json"
RESOURCES_FILENAME = "resources.csv"

# One row of the CarpetIOBasic screen table, e.g.
#   "     4800    30.000 |  0.0123456 |  1.2345e-05 | ..."
INFO_LINE_PATTERN = re.compile(r"^\s*(\d+)\s+([-+0-9.eE]+)\s*\|")


@dataclass
class MemorySample:
    """One sample of the container's memory accounting.

    Attributes:
        timestamp: Wall-clock time (Unix epoch seconds).
        current_bytes: ``memory.current``.
        peak_bytes: ``memory.peak`` (high-water mark since container start).
        limit_bytes: ``memory.max``; ``None`` when unlimited.
        oom_kill_count: ``oom_kill`` counter from ``memory.events``.
        load_1min: 1-minute load average from ``/proc/loadavg``.
    """

    timestamp: float
    current_bytes: int
    peak_bytes: Optional[int]
    limit_bytes: Optional[int]
    oom_kill_count: int
    load_1min: float


def read_memory_sample(cgroup_root: Path = CGROUP_ROOT) -> MemorySample:
    """Take one sample of the cgroup v2 memory accounting.

    Args:
        cgroup_root: Mount point of the container's cgroup v2 hierarchy.

    Returns:
        The sample.
    """
    raise NotImplementedError("read_memory_sample is not implemented yet")


class ResourceMonitor:
    """Background memory sampler with an OOM guard.

    Args:
        csv_path: File that receives one CSV row per sample.
        interval_s: Sampling interval in seconds.
        guard_fraction: Fraction of ``memory.max`` at which the guard trips.
    """

    def __init__(self, csv_path: Path, interval_s: float, guard_fraction: float) -> None:
        self.csv_path = csv_path
        self.interval_s = interval_s
        self.guard_fraction = guard_fraction

    def start(self) -> None:
        """Start the sampling thread."""
        raise NotImplementedError("ResourceMonitor.start is not implemented yet")

    def stop(self) -> None:
        """Stop the sampling thread and flush the CSV file."""
        raise NotImplementedError("ResourceMonitor.stop is not implemented yet")

    @property
    def guard_tripped(self) -> bool:
        """Whether memory use has crossed the guard threshold."""
        raise NotImplementedError("ResourceMonitor.guard_tripped is not implemented yet")

    @property
    def peak_bytes(self) -> Optional[int]:
        """Highest ``memory.current`` seen during this run."""
        raise NotImplementedError("ResourceMonitor.peak_bytes is not implemented yet")


def parse_info_line(line: str) -> Optional[Tuple[int, float]]:
    """Extract ``(iteration, time)`` from a CarpetIOBasic table row.

    Args:
        line: One line of Cactus stdout.

    Returns:
        The pair, or ``None`` for any other line.
    """
    raise NotImplementedError("parse_info_line is not implemented yet")


class EtaEstimator:
    """Sliding-window estimate of a run's simulation-time rate and ETA.

    Args:
        t_final: Physical end time of the run.
        window_s: Width of the wall-clock window used for the rate.
    """

    def __init__(self, t_final: float, window_s: float = 600.0) -> None:
        self.t_final = t_final
        self.window_s = window_s

    def observe(self, line: str, wall_time: float) -> bool:
        """Feed one stdout line.

        Args:
            line: One line of Cactus stdout.
            wall_time: Wall-clock time the line was read.

        Returns:
            Whether the line was a table row and updated the estimate.
        """
        raise NotImplementedError("EtaEstimator.observe is not implemented yet")

    @property
    def rate_per_hour(self) -> Optional[float]:
        """Simulation time advanced per wall-clock hour."""
        raise NotImplementedError("EtaEstimator.rate_per_hour is not implemented yet")

    @property
    def eta_seconds(self) -> Optional[float]:
        """Wall-clock seconds until ``t_final`` at the current rate."""
        raise NotImplementedError("EtaEstimator.eta_seconds is not implemented yet")


def campaign_eta_seconds(
    completed_wall_seconds: List[float],
    remaining_runs: int,
    current_run_eta: Optional[float],
) -> Optional[float]:
    """Estimate the wall-clock time left in a bisection campaign.

    Args:
        completed_wall_seconds: Wall time of every finished run.
        remaining_runs: Runs still to start after the current one.
        current_run_eta: ETA of the run in progress, if any.

    Returns:
        Seconds left, or ``None`` before any rate is known.
    """
    raise NotImplementedError("campaign_eta_seconds is not implemented yet")


def write_status(runs_dir: Path, status: dict) -> None:
    """Atomically replace ``<runs_dir>/status.json``.

    Args:
        runs_dir: Campaign directory.
        status: JSON-serialisable status snapshot.
    """
    raise NotImplementedError("write_status is not implemented yet")


def print_status(runs_dir: Path) -> int:
    """Pretty-print ``<runs_dir>/status.json`` for ``make status``.

    Args:
        runs_dir: Campaign directory.

    Returns:
        Process exit code.
    """
    raise NotImplementedError("print_status is not implemented yet")


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Parse the command line."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    subparsers = parser.add_subparsers(dest="command", required=True)
    status = subparsers.add_parser("status", help="print the campaign status")
    status.add_argument("--runs-dir", type=Path, required=True, help="campaign directory")
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    """Entry point."""
    args = parse_args(argv)
    if args.command == "status":
        return print_status(args.runs_dir)
    return 2


if __name__ == "__main__":
    sys.exit(main())
