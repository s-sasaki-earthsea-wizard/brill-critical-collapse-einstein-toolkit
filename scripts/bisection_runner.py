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

``--parfile`` runs any parfile under the same monitor, for feasibility
tests and benchmarks that do not use the template.

Status: the CLI, template rendering, ``--dry-run``, and monitored single
runs (``--single``, ``--parfile``) work. Outcome classification and the
bisection loop are not implemented yet (Phase 3).
"""

from __future__ import annotations

import argparse
import datetime
import enum
import json
import os
import queue
import re
import shlex
import shutil
import signal
import subprocess
import sys
import threading
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from run_monitor import (
    CGROUP_ROOT,
    RESOURCES_FILENAME,
    EtaEstimator,
    ResourceMonitor,
    format_duration,
    now_iso,
    write_status,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_TEMPLATE = REPO_ROOT / "templates" / "brill_collapse.par.template"
DEFAULT_CACTUS_BIN = "/home/etuser/Cactus/exe/cactus_sim"
DEFAULT_RUNS_DIR = Path("/home/etuser/simulations/bisection")
DEFAULT_MPIRUN = "mpirun.mpich"

HISTORY_FILENAME = "bisection_history.json"
HISTORY_SCHEMA_VERSION = 1
PARFILE_NAME = "brill_collapse.par"
LOG_NAME = "cactus.log"
RUN_SUMMARY_NAME = "run_summary.json"

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
# Running Cactus under the resource monitor
# ---------------------------------------------------------------------


@dataclass
class ExecutionResult:
    """How one cactus_sim invocation ended.

    Attributes:
        exit_code: mpirun exit code (negative: killed by that signal).
        wall_seconds: Wall-clock duration.
        peak_memory_bytes: Highest container ``memory.current`` seen.
        guard_tripped: Whether the memory guard stopped the run.
        oom_killed: Whether the kernel OOM killer fired during the run.
        timed_out: Whether the wall-clock limit stopped the run.
        last_iteration: Last iteration in the IOBasic table.
        t_reached: Simulation time of that iteration.
        rate_per_hour: Final simulation-time rate estimate.
    """

    exit_code: Optional[int]
    wall_seconds: float
    peak_memory_bytes: Optional[int]
    guard_tripped: bool
    oom_killed: bool
    timed_out: bool
    last_iteration: Optional[int]
    t_reached: Optional[float]
    rate_per_hour: Optional[float]


# Grace periods when stopping a run: first ask Cactus through the
# TerminationTrigger file (checked every few iterations), then SIGTERM,
# then SIGKILL on the whole process group.
TERMINATE_GRACE_S = 120.0
SIGTERM_GRACE_S = 30.0
STATUS_INTERVAL_S = 15.0


def request_termination(run_dir: Path) -> None:
    """Ask Cactus to stop cleanly through TerminationTrigger.

    The parfiles create ``IO::out_dir/TERMINATE`` (``output/TERMINATE``)
    and poll it; writing a non-zero integer ends the run after the
    current iteration.

    Args:
        run_dir: Working directory of the run.
    """
    trigger = run_dir / "output" / "TERMINATE"
    trigger.parent.mkdir(parents=True, exist_ok=True)
    trigger.write_text("1\n")


def _signal_group(process: subprocess.Popen, sig: int) -> None:
    try:
        os.killpg(process.pid, sig)
    except ProcessLookupError:
        pass


def _stream_lines(stream, sink: "queue.Queue[Optional[str]]") -> None:
    for line in iter(stream.readline, ""):
        sink.put(line)
    sink.put(None)


def execute_run(
    command: List[str],
    run_dir: Path,
    run_name: str,
    status_dir: Path,
    t_final: float,
    timeout_s: float,
    mem_guard_fraction: float,
    monitor_interval_s: float,
    cgroup_root: Path = CGROUP_ROOT,
) -> ExecutionResult:
    """Run Cactus in ``run_dir`` with memory monitoring and ETA reporting.

    Streams stdout/stderr into ``run_dir/cactus.log``, feeds every line to
    an ``EtaEstimator``, samples memory with a ``ResourceMonitor``, keeps
    ``status.json`` current, and stops the run (TerminationTrigger file,
    then SIGTERM, then SIGKILL) on timeout or when the memory guard trips.

    Args:
        command: mpirun argument list.
        run_dir: Working directory of the run.
        run_name: Label shown by ``make status``.
        status_dir: Directory that receives ``status.json``.
        t_final: Physical end time, for the ETA.
        timeout_s: Wall-clock limit for this run.
        mem_guard_fraction: Fraction of the container memory limit at
            which the run is aborted.
        monitor_interval_s: Memory sampling interval.
        cgroup_root: cgroup v2 hierarchy to sample (overridable for tests).

    Returns:
        How the run ended.
    """
    monitor = ResourceMonitor(
        run_dir / RESOURCES_FILENAME, monitor_interval_s, mem_guard_fraction, cgroup_root
    )
    eta = EtaEstimator(t_final)
    started = time.time()
    started_at = now_iso()
    monitor.start()

    def status_snapshot(state: str) -> dict:
        sample = monitor.latest
        limit = sample.limit_bytes if sample else None
        eta_s = eta.eta_seconds
        eta_at = None
        if eta_s is not None:
            eta_at = (
                datetime.datetime.now().astimezone() + datetime.timedelta(seconds=eta_s)
            ).isoformat(timespec="seconds")
        return {
            "run_name": run_name,
            "run_dir": str(run_dir),
            "state": state,
            "started_at": started_at,
            "updated_at": now_iso(),
            "elapsed_seconds": time.time() - started,
            "iteration": eta.iteration,
            "sim_time": eta.sim_time,
            "t_final": t_final,
            "rate_per_hour": eta.rate_per_hour,
            "eta_seconds": eta_s,
            "eta_at": eta_at,
            "campaign_eta_seconds": None,
            "memory": {
                "current_bytes": sample.current_bytes if sample else None,
                "run_peak_bytes": monitor.peak_bytes,
                "limit_bytes": limit,
                "guard_bytes": None if limit is None else int(mem_guard_fraction * limit),
                "load_1min": sample.load_1min if sample else None,
            },
        }

    timed_out = False
    stop_reason: Optional[str] = None
    stop_requested_at: Optional[float] = None
    sigterm_sent_at: Optional[float] = None
    last_status = 0.0
    lines: "queue.Queue[Optional[str]]" = queue.Queue()

    # Line-buffered so that `tail -f cactus.log` follows the run.
    with (run_dir / LOG_NAME).open("w", buffering=1) as log:
        process = subprocess.Popen(
            command,
            cwd=str(run_dir),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            bufsize=1,
            start_new_session=True,
        )
        reader = threading.Thread(target=_stream_lines, args=(process.stdout, lines), daemon=True)
        reader.start()
        write_status(status_dir, status_snapshot("running"))
        stream_open = True
        while stream_open or process.poll() is None:
            try:
                line = lines.get(timeout=1.0)
                if line is None:
                    stream_open = False
                else:
                    log.write(line)
                    eta.observe(line, time.time())
            except queue.Empty:
                pass

            now = time.time()
            if stop_reason is None:
                if monitor.guard_tripped or monitor.oom_killed:
                    stop_reason = "memory guard" if monitor.guard_tripped else "kernel OOM kill"
                elif now - started > timeout_s:
                    stop_reason = "wall-clock limit"
                    timed_out = True
                if stop_reason is not None:
                    message = "[runner] stopping the run: {}\n".format(stop_reason)
                    log.write(message)
                    print(message, end="", file=sys.stderr)
                    request_termination(run_dir)
                    stop_requested_at = now
            elif sigterm_sent_at is None and stop_requested_at is not None:
                if now - stop_requested_at > TERMINATE_GRACE_S and process.poll() is None:
                    log.write("[runner] no clean stop, sending SIGTERM\n")
                    _signal_group(process, signal.SIGTERM)
                    sigterm_sent_at = now
            elif sigterm_sent_at is not None and now - sigterm_sent_at > SIGTERM_GRACE_S:
                if process.poll() is None:
                    log.write("[runner] sending SIGKILL\n")
                    _signal_group(process, signal.SIGKILL)
                    sigterm_sent_at = float("inf")

            if now - last_status > STATUS_INTERVAL_S:
                write_status(status_dir, status_snapshot("running"))
                last_status = now
        exit_code = process.wait()
        reader.join(timeout=5.0)

    monitor.stop()
    result = ExecutionResult(
        exit_code=exit_code,
        wall_seconds=time.time() - started,
        peak_memory_bytes=monitor.peak_bytes,
        guard_tripped=monitor.guard_tripped,
        oom_killed=monitor.oom_killed,
        timed_out=timed_out,
        last_iteration=eta.iteration,
        t_reached=eta.sim_time,
        rate_per_hour=eta.rate_per_hour,
    )
    final = status_snapshot("finished")
    final["outcome"] = describe_execution(result)
    write_status(status_dir, final)
    return result


def describe_execution(result: ExecutionResult) -> str:
    """One-line summary of how a run ended (not a physics verdict)."""
    if result.oom_killed:
        return "oom (kernel OOM kill)"
    if result.guard_tripped:
        return "oom (memory guard)"
    if result.timed_out:
        return "timeout"
    if result.exit_code == 0:
        return "completed"
    return "failed (exit code {})".format(result.exit_code)


FINAL_TIME_PATTERN = re.compile(
    r"^\s*Cactus::cctk_final_time\s*=\s*([-+0-9.eE]+)", re.IGNORECASE | re.MULTILINE
)


def parse_final_time(parfile_text: str) -> Optional[float]:
    """Read ``Cactus::cctk_final_time`` from a parfile, for the ETA.

    Args:
        parfile_text: Parfile contents.

    Returns:
        The final time, or ``None`` when the parfile does not set it.
    """
    match = FINAL_TIME_PATTERN.search(parfile_text)
    return float(match.group(1)) if match else None


def run_prepared(args: argparse.Namespace, parfile: Path, run_name: str, t_final: float) -> int:
    """Execute an already prepared run directory and record the result.

    Args:
        args: Parsed command line.
        parfile: Parfile inside its run directory.
        run_name: Label for ``status.json`` and the summary.
        t_final: Physical end time, for the ETA.

    Returns:
        Process exit code: 0 when Cactus exited cleanly.
    """
    run_dir = parfile.parent
    command = build_mpirun_command(args.mpirun, args.cactus_bin, parfile, args.np, args.omp)
    print("[runner] {}: np={} x OMP={}, t_final={}, log {}".format(
        run_name, args.np, args.omp, t_final, run_dir / LOG_NAME))
    result = execute_run(
        command,
        run_dir,
        run_name,
        args.runs_dir,
        t_final,
        args.timeout_hours * 3600.0,
        args.mem_guard_fraction,
        args.monitor_interval,
    )
    summary = asdict(result)
    summary.update(
        {
            "run_name": run_name,
            "parfile": str(parfile),
            "command": command,
            "outcome": describe_execution(result),
            "finished_at": now_iso(),
        }
    )
    (run_dir / RUN_SUMMARY_NAME).write_text(json.dumps(summary, indent=2) + "\n")
    print("[runner] {}: {} after {}, t={} (it {}), peak memory {}".format(
        run_name,
        summary["outcome"],
        format_duration(result.wall_seconds),
        result.t_reached,
        result.last_iteration,
        "-" if result.peak_memory_bytes is None
        else "{:.2f} GiB".format(result.peak_memory_bytes / 1024 ** 3),
    ))
    return 0 if result.exit_code == 0 and not (result.guard_tripped or result.timed_out) else 1


def run_parfile(args: argparse.Namespace) -> int:
    """Run an arbitrary parfile under the monitor (feasibility tests, benchmarks).

    The parfile is copied into ``<runs_dir>/<run_name>/`` and run there,
    so the run directory is self-contained.

    Args:
        args: Parsed command line; ``args.parfile`` is the source parfile.

    Returns:
        Process exit code.
    """
    source = args.parfile
    text = source.read_text()
    run_name = args.run_name or source.stem
    run_dir = args.runs_dir / run_name
    if run_dir.exists():
        print("error: {} already exists; pick another --run-name".format(run_dir), file=sys.stderr)
        return 1
    run_dir.mkdir(parents=True)
    parfile = run_dir / source.name
    shutil.copyfile(str(source), str(parfile))
    t_final = parse_final_time(text)
    if t_final is None:
        t_final = args.t_final
    return run_prepared(args, parfile, run_name, t_final)


# ---------------------------------------------------------------------
# Classification and the bisection loop (not implemented yet)
# ---------------------------------------------------------------------


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

    The outcome is not classified yet (Phase 3); the run summary records
    how Cactus ended.

    Args:
        args: Parsed command line; ``args.single`` is the amplitude.

    Returns:
        Process exit code.
    """
    single_dir = args.runs_dir / "single"
    target = single_dir / run_dir_name(0, args.single)
    if target.exists():
        print("error: {} already exists; remove it first".format(target), file=sys.stderr)
        return 1
    parfile = prepare_run_dir(single_dir, 0, args.single, args.template, args.t_final)
    return run_prepared(args, parfile, "single A={}".format(format_amplitude(args.single)), args.t_final)


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
    mode.add_argument(
        "--parfile", type=Path, metavar="PAR",
        help="run this parfile as is under the monitor (feasibility tests, benchmarks)",
    )
    parser.add_argument(
        "--run-name", metavar="NAME",
        help="run directory name for --parfile (default: the parfile stem)",
    )
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
    if args.parfile is not None:
        return run_parfile(args)
    if args.dry_run:
        return dry_run(args)
    if args.single is not None:
        return run_single(args)
    return run_campaign(args)


if __name__ == "__main__":
    sys.exit(main())
