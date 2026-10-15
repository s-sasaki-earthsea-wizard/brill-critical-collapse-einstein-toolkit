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

Status: skeleton. The interfaces are fixed here; the bodies are not
implemented yet.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

BH_DIAGNOSTICS_FILENAME = "BH_diagnostics.ah1.gp"
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
    """

    iterations: List[int]
    times: List[float]
    irreducible_masses: List[float]
    areal_radii: List[float]

    @property
    def found(self) -> bool:
        """Whether a horizon with positive mass was found at least once."""
        return any(m > 0.0 for m in self.irreducible_masses)


@dataclass
class LogFailures:
    """Failure signatures extracted from a run log.

    Attributes:
        nan_reports: NaNChecker lines naming the offending variables.
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


def read_bh_diagnostics(output_dir: Path) -> Optional[HorizonSeries]:
    """Read ``BH_diagnostics.ah1.gp`` from a run's output directory.

    Args:
        output_dir: The run's ``output/`` directory.

    Returns:
        The horizon history, or ``None`` when the file does not exist
        (AHFinderDirect never found a horizon).
    """
    raise NotImplementedError("read_bh_diagnostics is not implemented yet")


def parse_column_header(lines: List[str]) -> Dict[str, int]:
    """Map column names to zero-based indices from ``# column N = name`` lines.

    Args:
        lines: Header lines of an AHFinderDirect ``.gp`` file.

    Returns:
        Column name to zero-based column index.
    """
    raise NotImplementedError("parse_column_header is not implemented yet")


def read_central_lapse(output_dir: Path) -> List[tuple]:
    """Read the central lapse alpha(0,0,0) on the finest level containing the origin.

    Args:
        output_dir: The run's ``output/`` directory.

    Returns:
        ``(time, alpha)`` pairs in time order.
    """
    raise NotImplementedError("read_central_lapse is not implemented yet")


def read_scalar_series(output_dir: Path, filename: str) -> List[tuple]:
    """Read a CarpetIOScalar reduction file.

    Args:
        output_dir: The run's ``output/`` directory.
        filename: e.g. ``MIN_LAPSE_FILENAME``.

    Returns:
        ``(time, value)`` pairs in time order.
    """
    raise NotImplementedError("read_scalar_series is not implemented yet")


def scan_log(log_path: Path) -> LogFailures:
    """Extract failure signatures and the last reached iteration from a run log.

    Args:
        log_path: The run's combined stdout/stderr log.

    Returns:
        The failure signatures found in the log.
    """
    raise NotImplementedError("scan_log is not implemented yet")
