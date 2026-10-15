#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
"""Bisection driver for the critical amplitude A* of Brill wave collapse.

Each iteration renders ``templates/brill_collapse.par.template`` for the
midpoint of the current bracket ``[a_min, a_max]``, runs Cactus under
mpirun, classifies the outcome and narrows the bracket:

* ``collapse``   -- an apparent horizon with M_AH > 0 was found
  (supercritical): ``a_max = A``.
* ``dispersion`` -- no horizon, and the central lapse recovered towards 1
  as the wave dispersed (subcritical): ``a_min = A``.
* anything else (``crash``, ``oom``, ``timeout``, ``undecided``) leaves
  the bracket untouched and stops the campaign for a human decision. A
  crash close to A* is expected with moving-puncture gauge
  (arXiv:1309.5008); guessing a side would silently corrupt A*.

Every run is appended to ``<runs_dir>/bisection_history.json`` with its
amplitude, outcome and horizon masses, so a campaign can be resumed and
analysed by ``analyze_critical_exponent.py``.

Each run is monitored for memory (OOM guard) and reports its ETA; see
``run_monitor.py``.

Meant to run inside the container (Python 3.8), where mpirun.mpich and
cactus_sim live; ``make bisection`` and friends wrap it.

Status: skeleton. The CLI, template rendering and ``--dry-run`` work.
Running Cactus, monitoring, classification and the bisection loop are
not implemented yet.
"""

from __future__ import annotations

import argparse
import enum
import os
import re
import shlex
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_TEMPLATE = REPO_ROOT / "templates" / "brill_collapse.par.template"
DEFAULT_CACTUS_BIN = "/home/etuser/Cactus/exe/cactus_sim"
DEFAULT_RUNS_DIR = Path("/home/etuser/simulations/bisection")
DEFAULT_MPIRUN = "mpirun.mpich"

HISTORY_FILENAME = "bisection_history.json"
HISTORY_SCHEMA_VERSION = 1
PARFILE_NAME = "brill_collapse.par"
LOG_NAME = "cactus.log"

# Template placeholders look like __AMPLITUDE__. Lower-case double
# underscores (AHFinderDirect's initial_guess__coord_sphere__radius)
# deliberately do not match.
PLACEHOLDER_PATTERN = re.compile(r"__[A-Z][A-Z0-9_]*__")


class Outcome(str, enum.Enum):
    """Verdict on a single run."""

    COLLAPSE = "collapse"
    DISPERSION = "dispersion"
    CRASH = "crash"
    OOM = "oom"
    TIMEOUT = "timeout"
    UNDECIDED = "undecided"


@dataclass
class RunResult:
    """Everything recorded about one run in ``bisection_history.json``.

    Attributes:
        iteration: Bisection step, starting at 0.
        amplitude: Brill wave amplitude A.
        run_dir: Run directory (parfile, log, output/).
        outcome: One of ``Outcome``.
        exit_code: mpirun exit code; ``None`` if it never started.
        m_ah_first: Irreducible mass at the first horizon find.
        t_ah_first: Simulation time of the first horizon find.
        m_ah_final: Irreducible mass at the last horizon find.
        t_ah_final: Simulation time of the last horizon find.
        min_lapse: Global minimum of the lapse over the whole run.
        final_central_lapse: Central lapse at the end of the run.
        t_reached: Last simulation time reached.
        wall_seconds: Wall-clock duration of the run.
        peak_memory_bytes: Highest container memory seen during the run.
        started_at: ISO 8601 start time.
        finished_at: ISO 8601 end time.
        notes: Free-form diagnosis, e.g. the variable NaNChecker flagged.
    """

    iteration: int
    amplitude: float
    run_dir: str
    outcome: str
    exit_code: Optional[int] = None
    m_ah_first: Optional[float] = None
    t_ah_first: Optional[float] = None
    m_ah_final: Optional[float] = None
    t_ah_final: Optional[float] = None
    min_lapse: Optional[float] = None
    final_central_lapse: Optional[float] = None
    t_reached: Optional[float] = None
    wall_seconds: Optional[float] = None
    peak_memory_bytes: Optional[int] = None
    started_at: Optional[str] = None
    finished_at: Optional[str] = None
    notes: str = ""


# ---------------------------------------------------------------------
# Parfile rendering (implemented: needed by --dry-run and the tests)
# ---------------------------------------------------------------------


def format_amplitude(amplitude: float) -> str:
    """Format an amplitude losslessly for parfiles and directory names.

    Args:
        amplitude: Brill wave amplitude.

    Returns:
        The shortest string that round-trips to the same float.
    """
    return repr(float(amplitude))


def render_parfile(template_text: str, values: Dict[str, str]) -> str:
    """Substitute ``__NAME__`` placeholders in a parfile template.

    Args:
        template_text: Template contents.
        values: Placeholder name (without underscores) to replacement text.

    Returns:
        The rendered parfile.

    Raises:
        KeyError: A value was given for a placeholder the template lacks.
        ValueError: The rendered text still contains a placeholder.
    """
    rendered = template_text
    for name, value in values.items():
        token = "__{}__".format(name)
        if token not in rendered:
            raise KeyError("placeholder {} not found in the template".format(token))
        rendered = rendered.replace(token, value)
    leftover = sorted(set(PLACEHOLDER_PATTERN.findall(rendered)))
    if leftover:
        raise ValueError("unrendered placeholders: {}".format(", ".join(leftover)))
    return rendered


def run_dir_name(iteration: int, amplitude: float) -> str:
    """Name of the directory for one bisection step.

    Args:
        iteration: Bisection step, starting at 0.
        amplitude: Brill wave amplitude.

    Returns:
        e.g. ``iter00_A5.0``.
    """
    return "iter{:02d}_A{}".format(iteration, format_amplitude(amplitude))


def prepare_run_dir(
    runs_dir: Path,
    iteration: int,
    amplitude: float,
    template_path: Path,
    t_final: float,
) -> Path:
    """Create a run directory and write the rendered parfile into it.

    Args:
        runs_dir: Campaign directory.
        iteration: Bisection step, starting at 0.
        amplitude: Brill wave amplitude.
        template_path: Parfile template.
        t_final: Physical end time.

    Returns:
        Path of the rendered parfile.

    Raises:
        FileExistsError: The run directory already exists. Runs are never
            overwritten; resuming is handled by ``--resume``.
    """
    run_dir = runs_dir / run_dir_name(iteration, amplitude)
    run_dir.mkdir(parents=True, exist_ok=False)
    rendered = render_parfile(
        template_path.read_text(),
        {
            "AMPLITUDE": format_amplitude(amplitude),
            "FINAL_TIME": repr(float(t_final)),
        },
    )
    parfile = run_dir / PARFILE_NAME
    parfile.write_text(rendered)
    return parfile


def build_mpirun_command(
    mpirun: str, cactus_bin: str, parfile: Path, n_procs: int, omp_threads: int
) -> List[str]:
    """Build the mpirun command line for one run.

    OMP_NUM_THREADS is always passed to every rank explicitly: left unset,
    each rank starts one thread per visible core and the run slows down
    by an order of magnitude without failing (GW230529 finding 7).

    Args:
        mpirun: MPICH launcher.
        cactus_bin: Path to cactus_sim.
        parfile: Rendered parfile.
        n_procs: Number of MPI ranks.
        omp_threads: OpenMP threads per rank.

    Returns:
        The command as an argument list.
    """
    return [
        mpirun,
        "-np", str(n_procs),
        "-genv", "OMP_NUM_THREADS", str(omp_threads),
        "-genv", "HDF5_USE_FILE_LOCKING", "FALSE",
        cactus_bin,
        str(parfile),
    ]


# ---------------------------------------------------------------------
# Running, classification and the bisection loop (not implemented yet)
# ---------------------------------------------------------------------


def execute_run(
    command: List[str],
    run_dir: Path,
    t_final: float,
    timeout_s: float,
    mem_guard_fraction: float,
    monitor_interval_s: float,
) -> Tuple[Optional[int], float, Optional[int], bool]:
    """Run Cactus in ``run_dir`` with memory monitoring and ETA reporting.

    Streams stdout/stderr into ``run_dir/cactus.log``, feeds every line to
    an ``EtaEstimator``, samples memory with a ``ResourceMonitor``, keeps
    ``status.json`` current, and stops the run (TerminationTrigger file,
    then SIGTERM, then SIGKILL) on timeout or when the memory guard trips.

    Args:
        command: mpirun argument list.
        run_dir: Working directory of the run.
        t_final: Physical end time, for the ETA.
        timeout_s: Wall-clock limit for this run.
        mem_guard_fraction: Fraction of the container memory limit at
            which the run is aborted.
        monitor_interval_s: Memory sampling interval.

    Returns:
        ``(exit_code, wall_seconds, peak_memory_bytes, guard_tripped)``.
    """
    raise NotImplementedError("execute_run is not implemented yet")


def classify_run(
    iteration: int,
    amplitude: float,
    run_dir: Path,
    exit_code: Optional[int],
    timed_out: bool,
    guard_tripped: bool,
) -> RunResult:
    """Decide the outcome of a finished run from its output and log.

    Order of precedence:

    1. memory guard tripped or kernel OOM kill        -> ``oom``
    2. ``BH_diagnostics.ah1.gp`` has a find with M > 0 -> ``collapse``
       (even if the run crashed afterwards: the horizon settles the verdict)
    3. NaNChecker report or abnormal exit              -> ``crash``
    4. wall-clock limit reached                        -> ``timeout``
    5. reached t_final, no horizon, central lapse back
       above a recovery threshold                      -> ``dispersion``
    6. otherwise (e.g. lapse collapsed but no horizon) -> ``undecided``

    Args:
        iteration: Bisection step.
        amplitude: Brill wave amplitude.
        run_dir: Run directory.
        exit_code: mpirun exit code.
        timed_out: Whether the wall-clock limit stopped the run.
        guard_tripped: Whether the memory guard stopped the run.

    Returns:
        The populated run record.
    """
    raise NotImplementedError("classify_run is not implemented yet")


def update_bracket(a_min: float, a_max: float, result: RunResult) -> Tuple[float, float]:
    """Narrow the bracket after a decisive run.

    Args:
        a_min: Current subcritical bound.
        a_max: Current supercritical bound.
        result: The run just classified.

    Returns:
        The new ``(a_min, a_max)``.

    Raises:
        RuntimeError: The outcome is not decisive; the campaign stops.
    """
    raise NotImplementedError("update_bracket is not implemented yet")


def load_history(runs_dir: Path) -> dict:
    """Load ``bisection_history.json``, or start a new history."""
    raise NotImplementedError("load_history is not implemented yet")


def save_history(runs_dir: Path, history: dict) -> None:
    """Atomically replace ``bisection_history.json`` (write, fsync, rename)."""
    raise NotImplementedError("save_history is not implemented yet")


def run_campaign(args: argparse.Namespace) -> int:
    """Run (or resume) the bisection campaign.

    Args:
        args: Parsed command line.

    Returns:
        Process exit code.
    """
    raise NotImplementedError("run_campaign is not implemented yet")


def run_single(args: argparse.Namespace) -> int:
    """Run one amplitude outside the bisection (smoke tests, benchmarks).

    Args:
        args: Parsed command line; ``args.single`` is the amplitude.

    Returns:
        Process exit code.
    """
    raise NotImplementedError("run_single is not implemented yet")


# ---------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------


def dry_run(args: argparse.Namespace) -> int:
    """Render the first parfile and print the plan without running Cactus.

    Args:
        args: Parsed command line.

    Returns:
        Process exit code.
    """
    amplitude = args.single if args.single is not None else 0.5 * (args.a_min + args.a_max)
    dry_dir = args.runs_dir / "dry-run"
    target = dry_dir / run_dir_name(0, amplitude)
    if target.exists():
        print("error: {} already exists; remove it first".format(target), file=sys.stderr)
        return 1
    parfile = prepare_run_dir(dry_dir, 0, amplitude, args.template, args.t_final)
    command = build_mpirun_command(
        args.mpirun, args.cactus_bin, parfile, args.np, args.omp
    )
    print("bracket            : [{}, {}]".format(args.a_min, args.a_max))
    print("max iterations     : {}".format(args.max_iter))
    print("final bracket width: {:.3g}".format((args.a_max - args.a_min) / 2 ** args.max_iter))
    print("first amplitude    : {}".format(format_amplitude(amplitude)))
    print("t_final            : {}".format(args.t_final))
    print("parallelism        : np={} x OMP={}".format(args.np, args.omp))
    print("memory guard       : {:.0%} of the container limit".format(args.mem_guard_fraction))
    print("timeout per run    : {} h".format(args.timeout_hours))
    print("rendered parfile   : {}".format(parfile))
    print("command            : {}".format(" ".join(shlex.quote(c) for c in command)))
    print("ETA                : unknown until the first run measures a rate")
    return 0


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Parse the command line."""
    parser = argparse.ArgumentParser(
        description="Bisection search for the critical amplitude of Brill wave collapse."
    )
    parser.add_argument("--a-min", type=float, default=4.0, help="subcritical bound (default: 4.0)")
    parser.add_argument("--a-max", type=float, default=6.0, help="supercritical bound (default: 6.0)")
    parser.add_argument("--max-iter", type=int, default=7, help="bisection steps (default: 7)")
    parser.add_argument(
        "--cactus-bin", default=DEFAULT_CACTUS_BIN,
        help="Cactus executable (default: %(default)s)",
    )
    parser.add_argument(
        "--template", type=Path, default=DEFAULT_TEMPLATE,
        help="parfile template (default: templates/brill_collapse.par.template)",
    )
    parser.add_argument(
        "--runs-dir", type=Path, default=DEFAULT_RUNS_DIR,
        help="campaign directory (default: %(default)s)",
    )
    parser.add_argument("--t-final", type=float, default=30.0, help="physical end time (default: 30.0)")
    parser.add_argument(
        "--np", type=int, default=int(os.environ.get("SIM_MPI_PROCS", "4")),
        help="MPI ranks (default: $SIM_MPI_PROCS or 4)",
    )
    parser.add_argument(
        "--omp", type=int, default=int(os.environ.get("SIM_OMP_THREADS", "4")),
        help="OpenMP threads per rank (default: $SIM_OMP_THREADS or 4)",
    )
    parser.add_argument("--mpirun", default=DEFAULT_MPIRUN, help="MPI launcher (default: %(default)s)")
    parser.add_argument(
        "--timeout-hours", type=float, default=6.0,
        help="wall-clock limit per run in hours (default: 6.0)",
    )
    parser.add_argument(
        "--mem-guard-fraction", type=float, default=0.9,
        help="abort a run above this fraction of the container memory limit (default: 0.9)",
    )
    parser.add_argument(
        "--monitor-interval", type=float, default=10.0,
        help="memory sampling interval in seconds (default: 10)",
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--resume", action="store_true", help="continue from bisection_history.json")
    mode.add_argument("--single", type=float, metavar="A", help="run one amplitude, no bisection")
    parser.add_argument(
        "--dry-run", action="store_true",
        help="render the first parfile and print the plan without running Cactus",
    )
    args = parser.parse_args(argv)
    if not args.a_min < args.a_max:
        parser.error("--a-min must be smaller than --a-max")
    if not 0.0 < args.mem_guard_fraction < 1.0:
        parser.error("--mem-guard-fraction must lie in (0, 1)")
    return args


def main(argv: Optional[List[str]] = None) -> int:
    """Entry point."""
    args = parse_args(argv)
    if args.dry_run:
        return dry_run(args)
    if args.single is not None:
        return run_single(args)
    return run_campaign(args)


if __name__ == "__main__":
    sys.exit(main())
