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
"""

from __future__ import annotations

import argparse
import collections
import csv
import datetime
import json
import os
import re
import statistics
import sys
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Deque, List, Optional, Tuple

CGROUP_ROOT = Path("/sys/fs/cgroup")
STATUS_FILENAME = "status.json"
RESOURCES_FILENAME = "resources.csv"

GIB = 1024 ** 3

# One row of the CarpetIOBasic screen table, e.g.
#   "     4800    30.000 |  0.0123456 |  1.2345e-05 | ..."
INFO_LINE_PATTERN = re.compile(r"^\s*(\d+)\s+([-+0-9.eE]+)\s*\|")

RESOURCES_COLUMNS = (
    "timestamp",
    "current_bytes",
    "peak_bytes",
    "limit_bytes",
    "oom_kill_count",
    "load_1min",
)


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


def _read_int_file(path: Path) -> Optional[int]:
    """Read a cgroup file holding one integer, or ``max`` for unlimited."""
    try:
        text = path.read_text().strip()
    except FileNotFoundError:
        return None
    if text == "max":
        return None
    return int(text)


def _read_oom_kill_count(events_path: Path) -> int:
    """Read the ``oom_kill`` counter from ``memory.events``."""
    try:
        for line in events_path.read_text().splitlines():
            key, _, value = line.partition(" ")
            if key == "oom_kill":
                return int(value)
    except FileNotFoundError:
        pass
    return 0


def _read_load_1min(loadavg_path: Path = Path("/proc/loadavg")) -> float:
    """Read the 1-minute load average."""
    try:
        return float(loadavg_path.read_text().split()[0])
    except (FileNotFoundError, IndexError, ValueError):
        return float("nan")


def read_memory_sample(cgroup_root: Path = CGROUP_ROOT) -> MemorySample:
    """Take one sample of the cgroup v2 memory accounting.

    Args:
        cgroup_root: Mount point of the container's cgroup v2 hierarchy.

    Returns:
        The sample.

    Raises:
        FileNotFoundError: ``memory.current`` is missing, i.e. this is not
            a cgroup v2 hierarchy.
    """
    current = _read_int_file(cgroup_root / "memory.current")
    if current is None:
        raise FileNotFoundError("no memory.current under {}".format(cgroup_root))
    return MemorySample(
        timestamp=time.time(),
        current_bytes=current,
        peak_bytes=_read_int_file(cgroup_root / "memory.peak"),
        limit_bytes=_read_int_file(cgroup_root / "memory.max"),
        oom_kill_count=_read_oom_kill_count(cgroup_root / "memory.events"),
        load_1min=_read_load_1min(),
    )


class ResourceMonitor:
    """Background memory sampler with an OOM guard.

    Args:
        csv_path: File that receives one CSV row per sample.
        interval_s: Sampling interval in seconds.
        guard_fraction: Fraction of ``memory.max`` at which the guard trips.
        cgroup_root: Mount point of the cgroup v2 hierarchy.
    """

    def __init__(
        self,
        csv_path: Path,
        interval_s: float,
        guard_fraction: float,
        cgroup_root: Path = CGROUP_ROOT,
    ) -> None:
        self.csv_path = csv_path
        self.interval_s = interval_s
        self.guard_fraction = guard_fraction
        self.cgroup_root = cgroup_root
        self._stop_event = threading.Event()
        self._thread: Optional[threading.Thread] = None
        self._lock = threading.Lock()
        self._guard_tripped = False
        self._peak_bytes: Optional[int] = None
        self._baseline_oom_kills: Optional[int] = None
        self._oom_killed = False
        self._latest: Optional[MemorySample] = None

    def start(self) -> None:
        """Take a first sample synchronously, then start the sampling thread."""
        self.csv_path.parent.mkdir(parents=True, exist_ok=True)
        with self.csv_path.open("w", newline="") as handle:
            csv.writer(handle).writerow(RESOURCES_COLUMNS)
        self.sample_once()
        self._thread = threading.Thread(target=self._run, name="resource-monitor", daemon=True)
        self._thread.start()

    def stop(self) -> None:
        """Stop the sampling thread after one last sample."""
        self._stop_event.set()
        if self._thread is not None:
            self._thread.join(timeout=self.interval_s + 5.0)
        self.sample_once()

    def sample_once(self) -> MemorySample:
        """Take one sample, record it and update the guard state."""
        sample = read_memory_sample(self.cgroup_root)
        with self._lock:
            self._latest = sample
            if self._peak_bytes is None or sample.current_bytes > self._peak_bytes:
                self._peak_bytes = sample.current_bytes
            if self._baseline_oom_kills is None:
                self._baseline_oom_kills = sample.oom_kill_count
            elif sample.oom_kill_count > self._baseline_oom_kills:
                self._oom_killed = True
            if (
                sample.limit_bytes is not None
                and sample.current_bytes > self.guard_fraction * sample.limit_bytes
            ):
                self._guard_tripped = True
        with self.csv_path.open("a", newline="") as handle:
            csv.writer(handle).writerow(
                [
                    "{:.3f}".format(sample.timestamp),
                    sample.current_bytes,
                    "" if sample.peak_bytes is None else sample.peak_bytes,
                    "" if sample.limit_bytes is None else sample.limit_bytes,
                    sample.oom_kill_count,
                    sample.load_1min,
                ]
            )
        return sample

    def _run(self) -> None:
        while not self._stop_event.wait(self.interval_s):
            try:
                self.sample_once()
            except OSError as exc:  # keep sampling through transient read errors
                print("[monitor] memory sample failed: {}".format(exc), file=sys.stderr)

    @property
    def guard_tripped(self) -> bool:
        """Whether memory use has crossed the guard threshold."""
        with self._lock:
            return self._guard_tripped

    @property
    def oom_killed(self) -> bool:
        """Whether the kernel OOM killer fired in this cgroup during the run."""
        with self._lock:
            return self._oom_killed

    @property
    def peak_bytes(self) -> Optional[int]:
        """Highest ``memory.current`` seen during this run."""
        with self._lock:
            return self._peak_bytes

    @property
    def latest(self) -> Optional[MemorySample]:
        """Most recent sample."""
        with self._lock:
            return self._latest


def parse_info_line(line: str) -> Optional[Tuple[int, float]]:
    """Extract ``(iteration, time)`` from a CarpetIOBasic table row.

    Args:
        line: One line of Cactus stdout.

    Returns:
        The pair, or ``None`` for any other line.
    """
    match = INFO_LINE_PATTERN.match(line)
    if match is None:
        return None
    try:
        return int(match.group(1)), float(match.group(2))
    except ValueError:
        return None


class EtaEstimator:
    """Sliding-window estimate of a run's simulation-time rate and ETA.

    Args:
        t_final: Physical end time of the run.
        window_s: Width of the wall-clock window used for the rate.
    """

    def __init__(self, t_final: float, window_s: float = 600.0) -> None:
        self.t_final = t_final
        self.window_s = window_s
        self._points: Deque[Tuple[float, float]] = collections.deque()
        self.iteration: Optional[int] = None
        self.sim_time: Optional[float] = None

    def observe(self, line: str, wall_time: float) -> bool:
        """Feed one stdout line.

        Args:
            line: One line of Cactus stdout.
            wall_time: Wall-clock time the line was read.

        Returns:
            Whether the line was a table row and updated the estimate.
        """
        parsed = parse_info_line(line)
        if parsed is None:
            return False
        self.iteration, self.sim_time = parsed
        self._points.append((wall_time, self.sim_time))
        # Keep one point older than the window so the rate always spans
        # the full window once the run is long enough.
        while len(self._points) > 2 and self._points[1][0] < wall_time - self.window_s:
            self._points.popleft()
        return True

    @property
    def rate_per_hour(self) -> Optional[float]:
        """Simulation time advanced per wall-clock hour."""
        if len(self._points) < 2:
            return None
        (w0, t0), (w1, t1) = self._points[0], self._points[-1]
        if w1 <= w0 or t1 <= t0:
            return None
        return (t1 - t0) / (w1 - w0) * 3600.0

    @property
    def eta_seconds(self) -> Optional[float]:
        """Wall-clock seconds until ``t_final`` at the current rate."""
        rate = self.rate_per_hour
        if rate is None or self.sim_time is None:
            return None
        return max(self.t_final - self.sim_time, 0.0) / rate * 3600.0


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
    if remaining_runs > 0 and not completed_wall_seconds:
        return None
    future = remaining_runs * statistics.median(completed_wall_seconds) if remaining_runs > 0 else 0.0
    if current_run_eta is None:
        return future if remaining_runs > 0 else None
    return future + current_run_eta


def format_duration(seconds: Optional[float]) -> str:
    """Format seconds as ``1d 02:03:04``, ``02:03:04`` or ``-``."""
    if seconds is None:
        return "-"
    total = int(round(seconds))
    days, rest = divmod(total, 86400)
    hours, rest = divmod(rest, 3600)
    minutes, secs = divmod(rest, 60)
    clock = "{:02d}:{:02d}:{:02d}".format(hours, minutes, secs)
    return "{}d {}".format(days, clock) if days else clock


def now_iso() -> str:
    """Current local time in ISO 8601 with offset, to the second."""
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def write_status(runs_dir: Path, status: dict) -> None:
    """Atomically replace ``<runs_dir>/status.json``.

    Args:
        runs_dir: Campaign directory.
        status: JSON-serialisable status snapshot.
    """
    runs_dir.mkdir(parents=True, exist_ok=True)
    target = runs_dir / STATUS_FILENAME
    tmp = target.with_name(target.name + ".tmp")
    with tmp.open("w") as handle:
        json.dump(status, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, target)


def _gib(value: Optional[int]) -> str:
    return "-" if value is None else "{:.2f} GiB".format(value / GIB)


def print_status(runs_dir: Path) -> int:
    """Pretty-print ``<runs_dir>/status.json`` for ``make status``.

    Args:
        runs_dir: Campaign directory.

    Returns:
        Process exit code.
    """
    path = runs_dir / STATUS_FILENAME
    if not path.exists():
        print("no status file at {}".format(path))
        return 1
    status = json.loads(path.read_text())
    memory = status.get("memory", {})
    rate = status.get("rate_per_hour")
    rows = [
        ("run", status.get("run_name", "-")),
        ("state", status.get("state", "-")),
        ("updated", status.get("updated_at", "-")),
        ("iteration", status.get("iteration", "-")),
        ("time", "{} / {}".format(status.get("sim_time", "-"), status.get("t_final", "-"))),
        ("rate", "-" if rate is None else "{:.3f} per hour".format(rate)),
        ("run ETA", "{} (at {})".format(
            format_duration(status.get("eta_seconds")), status.get("eta_at") or "-")),
        ("campaign ETA", format_duration(status.get("campaign_eta_seconds"))),
        ("memory", "{} now, {} peak this run, limit {}, guard at {}".format(
            _gib(memory.get("current_bytes")),
            _gib(memory.get("run_peak_bytes")),
            _gib(memory.get("limit_bytes")),
            _gib(memory.get("guard_bytes")))),
        ("elapsed", format_duration(status.get("elapsed_seconds"))),
    ]
    if status.get("outcome"):
        rows.append(("outcome", status["outcome"]))
    width = max(len(name) for name, _ in rows)
    for name, value in rows:
        print("{:<{w}} : {}".format(name, value, w=width))
    return 0


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
