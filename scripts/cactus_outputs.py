# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
"""Readers for the Cactus output files that the bisection relies on.

Shared by ``bisection_runner.py`` (outcome classification) and
``analyze_critical_exponent.py``. All paths are relative to a run's
``output/`` directory (``IO::out_dir`` in the parfile template).

Files of interest:

* ``BH_diagnostics.ah1.gp`` (AHFinderDirect): one row per successful
  horizon find. The column layout is read from the file's ``# column N =
  name`` header rather than hard-coded, since it differs between
  AHFinderDirect versions.
* ``h.t<iteration>.ah1.gp`` (AHFinderDirect): horizon shapes, only used
  for inspection.
* ``admbase-lapse..asc`` (CarpetIOASCII 0D): central lapse alpha(0,0,0)
  per refinement level.
* ``admbase-lapse.minimum.asc`` (CarpetIOScalar): global minimum lapse.
* ``ml_admconstraints-ml_ham.norm2.asc`` (CarpetIOScalar): Hamiltonian
  constraint violation.
* the run log (stdout and stderr of mpirun): NaNChecker reports, MPI
  aborts and the kernel's "Killed" message.

``horizon`` selects the AHFinderDirect horizon number (``ah<N>``) for
runs with several finders.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

BH_DIAGNOSTICS_FILENAME = "BH_diagnostics.ah1.gp"
BH_DIAGNOSTICS_PATTERN = "BH_diagnostics.ah{}.gp"
CENTRAL_LAPSE_FILENAME = "admbase-lapse..asc"
MIN_LAPSE_FILENAME = "admbase-lapse.minimum.asc"
HAMILTONIAN_NORM2_FILENAME = "ml_admconstraints-ml_ham.norm2.asc"


@dataclass
class HorizonSeries:
    """Apparent horizon history of one run, as found by AHFinderDirect.

    Attributes:
        iterations: Cactus iterations at which the horizon was found.
        times: Simulation times of those finds.
        irreducible_masses: ``m_irreducible`` column, sqrt(A / 16 pi).
        areal_radii: ``areal_radius`` column.
        mean_radii: ``mean radius`` column (coordinate radius).
        min_radii: ``min radius`` column (coordinate radius).
        centroids_z: ``centroid_z`` column.
    """

    iterations: List[int]
    times: List[float]
    irreducible_masses: List[float]
    areal_radii: List[float]
    mean_radii: List[float]
    min_radii: List[float]
    centroids_z: List[float]

    @property
    def found(self) -> bool:
        """Whether a horizon with positive mass was found at least once."""
        return any(m > 0.0 for m in self.irreducible_masses)


@dataclass
class LogFailures:
    """Failure signatures extracted from a run log.

    Attributes:
        nan_reports: Variables NaNChecker reported, in order of first report.
        mpi_abort: Whether MPI reported an abort or a dead rank.
        killed: Whether the kernel's "Killed" message (OOM) appears.
        last_iteration: Last iteration seen in the IOBasic table.
        last_time: Simulation time of that iteration.
    """

    nan_reports: List[str]
    mpi_abort: bool
    killed: bool
    last_iteration: Optional[int]
    last_time: Optional[float]


COLUMN_HEADER_PATTERN = re.compile(r"^#\s*column\s+(\d+)\s*=\s*(.+?)\s*$")
NAN_REPORT_PATTERN = re.compile(r"There were \d+ NaN/Inf value\(s\) found in variable '([^']+)'")
INFO_ROW_PATTERN = re.compile(r"^\s*(\d+)\s+([-+0-9.eE]+)\s*\|")
MPI_ABORT_PATTERN = re.compile(r"MPI_ABORT|BAD TERMINATION|application called MPI_Abort", re.I)
KILLED_PATTERN = re.compile(r"\bKilled\b|signal 9|Out of memory", re.I)


def read_bh_diagnostics(output_dir: Path, horizon: int = 1) -> Optional[HorizonSeries]:
    """Read ``BH_diagnostics.ah<N>.gp`` from a run's output directory.

    Args:
        output_dir: The run's ``output/`` directory.
        horizon: AHFinderDirect horizon number.

    Returns:
        The horizon history, or ``None`` when the file does not exist
        (AHFinderDirect never found that horizon).
    """
    path = output_dir / BH_DIAGNOSTICS_PATTERN.format(horizon)
    if not path.exists():
        return None
    lines = path.read_text().splitlines()
    columns = parse_column_header([line for line in lines if line.startswith("#")])
    required = (
        "cctk_iteration", "cctk_time", "m_irreducible", "areal radius",
        "mean radius", "min radius", "centroid_z",
    )
    missing = [name for name in required if name not in columns]
    if missing:
        raise ValueError("{}: missing columns {}".format(path, ", ".join(missing)))
    series = HorizonSeries([], [], [], [], [], [], [])
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split()
        series.iterations.append(int(float(fields[columns["cctk_iteration"]])))
        series.times.append(float(fields[columns["cctk_time"]]))
        series.irreducible_masses.append(float(fields[columns["m_irreducible"]]))
        series.areal_radii.append(float(fields[columns["areal radius"]]))
        series.mean_radii.append(float(fields[columns["mean radius"]]))
        series.min_radii.append(float(fields[columns["min radius"]]))
        series.centroids_z.append(float(fields[columns["centroid_z"]]))
    return series


def parse_column_header(lines: List[str]) -> Dict[str, int]:
    """Map column names to zero-based indices from ``# column N = name`` lines.

    Args:
        lines: Header lines of an AHFinderDirect ``.gp`` file.

    Returns:
        Column name to zero-based column index.
    """
    columns = {}
    for line in lines:
        match = COLUMN_HEADER_PATTERN.match(line)
        if match:
            columns[match.group(2)] = int(match.group(1)) - 1
    return columns


def read_central_lapse(output_dir: Path) -> List[Tuple[float, float]]:
    """Read the central lapse alpha(0,0,0) on the finest level containing the origin.

    The 0D file has one row per refinement level and output iteration;
    at each time the row of the highest refinement level is kept.

    Args:
        output_dir: The run's ``output/`` directory.

    Returns:
        ``(time, alpha)`` pairs in time order.
    """
    finest: Dict[float, Tuple[int, float]] = {}
    for line in (output_dir / CENTRAL_LAPSE_FILENAME).read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split()
        level, time, value = int(fields[2]), float(fields[8]), float(fields[12])
        if time not in finest or level > finest[time][0]:
            finest[time] = (level, value)
    return [(time, finest[time][1]) for time in sorted(finest)]


def read_scalar_series(output_dir: Path, filename: str) -> List[Tuple[float, float]]:
    """Read a CarpetIOScalar reduction file.

    Args:
        output_dir: The run's ``output/`` directory.
        filename: e.g. ``MIN_LAPSE_FILENAME``.

    Returns:
        ``(time, value)`` pairs in time order; for a restarted or
        overlapping file the last value at a given time wins.
    """
    values: Dict[float, float] = {}
    for line in (output_dir / filename).read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split()
        values[float(fields[1])] = float(fields[2])
    return sorted(values.items())


def scan_log(log_path: Path) -> LogFailures:
    """Extract failure signatures and the last reached iteration from a run log.

    Args:
        log_path: The run's combined stdout/stderr log.

    Returns:
        The failure signatures found in the log.
    """
    failures = LogFailures([], False, False, None, None)
    with log_path.open(errors="replace") as handle:
        for line in handle:
            nan = NAN_REPORT_PATTERN.search(line)
            if nan and nan.group(1) not in failures.nan_reports:
                failures.nan_reports.append(nan.group(1))
            if MPI_ABORT_PATTERN.search(line):
                failures.mpi_abort = True
            if KILLED_PATTERN.search(line):
                failures.killed = True
            row = INFO_ROW_PATTERN.match(line)
            if row:
                failures.last_iteration = int(row.group(1))
                failures.last_time = float(row.group(2))
    return failures
