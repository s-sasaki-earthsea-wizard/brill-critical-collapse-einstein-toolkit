#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
"""Bisection driver for the critical amplitude A* of Brill wave collapse.

Each iteration renders ``templates/brill_collapse.par.template`` for the
midpoint of the current bracket ``[a_min, a_max]``, runs Cactus under
mpirun, classifies the outcome and narrows the bracket:

* ``collapse``   -- a confirmed apparent horizon (supercritical):
  ``a_max = A``.
* ``dispersion`` -- positive evidence that the wave dispersed
  (subcritical): ``a_min = A``.
* anything else (``crash``, ``oom``, ``timeout``, ``undecided``) leaves
  the bracket untouched and stops the campaign for a human decision. A
  crash close to A* is expected with moving-puncture gauge
  (arXiv:1309.5008); guessing a side would silently corrupt A*.

The criteria live in ``outcome_classifier.py``. The same evidence stops
a run early once the horizon mass has settled or dispersion is
established.

Every run is appended to ``<runs_dir>/bisection_history.json`` with its
amplitude, outcome, reasons and evidence, so a campaign can be resumed
and analysed by ``analyze_critical_exponent.py``. The history records a
fingerprint of the numerical setup (template, classifier thresholds,
end time, Cactus executable, initial data grid); resuming with a
different setup is refused, because runs from two setups do not share
one A*. ``--set-outcome`` records a human verdict for a run the
classifier left undecided, and the bracket is replayed from the runs.

If the template contains ``__TABLE_FILE__``, the runner solves the Brill
initial data for each amplitude (``brill_initial_data.py``) and writes
the table into the run directory.

Each run is monitored for memory (OOM guard) and reports its ETA; see
``run_monitor.py``.

Meant to run inside the container (Python 3.8), where mpirun.mpich and
cactus_sim live; ``make bisection`` and friends wrap it.

``--parfile`` runs any parfile under the same monitor, for feasibility
tests and benchmarks that do not use the template.

"""

from __future__ import annotations

import argparse
import datetime
import enum
import hashlib
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
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import outcome_classifier as oc
from run_monitor import (
    CGROUP_ROOT,
    RESOURCES_FILENAME,
    EtaEstimator,
    ResourceMonitor,
    campaign_eta_seconds,
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
TABLE_NAME = "psi_table.h5"
EARLY_STOP_PREFIX = "decided: "

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
        iteration: Index of the run in the campaign, starting at 0.
        amplitude: Brill wave amplitude A.
        run_dir: Run directory (parfile, log, output/).
        outcome: One of ``Outcome``.
        role: ``endpoint``, ``bisection`` or ``single``.
        reasons: Why the classifier decided as it did.
        exit_code: mpirun exit code; ``None`` if it never started.
        stop_reason: Why the runner stopped Cactus, if it did.
        horizon_number: AHFinderDirect horizon that confirmed collapse.
        horizon_kind: ``origin`` or ``axis``.
        m_ah_first: Irreducible mass at its first find.
        t_ah_first: Simulation time of its first find.
        t_ah_confirmed: Time the collapse criteria were met.
        m_ah_final: Irreducible mass at its last find.
        t_ah_final: Simulation time of its last find.
        m_ah_settled: Settled irreducible mass; ``None`` if still changing.
        r_ah_min: Smallest coordinate radius after confirmation.
        z_ah: Centroid z at the last find.
        horizons: Summary of every finder that found something.
        curvature_initial: Global max of the curvature invariant at t = 0.
        curvature_max: Its maximum over the run.
        t_curvature_max: Time of that maximum.
        curvature_final: Its last value.
        t_dispersed: Time the dispersion window criteria were first met.
        min_lapse: Global minimum of the lapse over the run.
        final_min_lapse: Global minimum of the lapse at the end.
        t_reached: Last simulation time reached.
        wall_seconds: Wall-clock duration of the run.
        peak_memory_bytes: Highest container memory seen during the run.
        started_at: ISO 8601 start time.
        finished_at: ISO 8601 end time.
        human_override: Classifier outcome replaced by ``--set-outcome``.
        notes: Free-form remarks.
    """

    iteration: int
    amplitude: float
    run_dir: str
    outcome: str
    role: str = "bisection"
    reasons: List[str] = field(default_factory=list)
    exit_code: Optional[int] = None
    stop_reason: Optional[str] = None
    horizon_number: Optional[int] = None
    horizon_kind: Optional[str] = None
    m_ah_first: Optional[float] = None
    t_ah_first: Optional[float] = None
    t_ah_confirmed: Optional[float] = None
    m_ah_final: Optional[float] = None
    t_ah_final: Optional[float] = None
    m_ah_settled: Optional[float] = None
    r_ah_min: Optional[float] = None
    z_ah: Optional[float] = None
    horizons: List[dict] = field(default_factory=list)
    curvature_initial: Optional[float] = None
    curvature_max: Optional[float] = None
    t_curvature_max: Optional[float] = None
    curvature_final: Optional[float] = None
    t_dispersed: Optional[float] = None
    min_lapse: Optional[float] = None
    final_min_lapse: Optional[float] = None
    t_reached: Optional[float] = None
    wall_seconds: Optional[float] = None
    peak_memory_bytes: Optional[int] = None
    started_at: Optional[str] = None
    finished_at: Optional[str] = None
    human_override: Optional[str] = None
    notes: str = ""

    @classmethod
    def from_dict(cls, data: dict) -> "RunResult":
        """Rebuild from a history entry, ignoring unknown keys."""
        names = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in data.items() if k in names})


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


@dataclass
class InitialDataSettings:
    """Grid of the initial data solve when the template uses a table.

    Attributes:
        n_s: Cells in the compactified radius.
        n_theta: Cells in the polar angle.
        scale: L in r = L tan(pi s / 2).
    """

    n_s: int = 512
    n_theta: int = 64
    scale: float = 2.0


def template_uses_table(template_text: str) -> bool:
    """Whether the template expects a per-amplitude initial data table."""
    return "__TABLE_FILE__" in template_text


def write_initial_data_table(path: Path, amplitude: float, settings: InitialDataSettings) -> float:
    """Solve the Brill initial data for one amplitude and write its table.

    Args:
        path: Output HDF5 file.
        amplitude: Brill amplitude A (sigma = 1, e^{2q}).
        settings: Solver grid.

    Returns:
        The ADM mass (flux estimate).
    """
    import brill_initial_data as bid  # deferred: numpy and scipy are not needed otherwise

    solution = bid.solve_brill(
        bid.BrillParams(amplitude), bid.Grid(settings.n_s, settings.n_theta, settings.scale)
    )
    masses = bid.adm_mass(solution)
    bid.write_table(path, solution, masses)
    return masses.flux


def prepare_run_dir(
    runs_dir: Path,
    iteration: int,
    amplitude: float,
    template_path: Path,
    t_final: float,
    id_settings: Optional[InitialDataSettings] = None,
) -> Path:
    """Create a run directory and write the rendered parfile into it.

    Args:
        runs_dir: Campaign directory.
        iteration: Index of the run, starting at 0.
        amplitude: Brill wave amplitude.
        template_path: Parfile template.
        t_final: Physical end time.
        id_settings: Initial data grid, used when the template has
            ``__TABLE_FILE__``.

    Returns:
        Path of the rendered parfile.

    Raises:
        FileExistsError: The run directory already exists. Runs are never
            overwritten; resuming is handled by ``--resume``.
    """
    run_dir = runs_dir / run_dir_name(iteration, amplitude)
    run_dir.mkdir(parents=True, exist_ok=False)
    template_text = template_path.read_text()
    values = {
        "AMPLITUDE": format_amplitude(amplitude),
        "FINAL_TIME": repr(float(t_final)),
    }
    if template_uses_table(template_text):
        table = run_dir / TABLE_NAME
        write_initial_data_table(table, amplitude, id_settings or InitialDataSettings())
        values["TABLE_FILE"] = str(table)
    parfile = run_dir / PARFILE_NAME
    parfile.write_text(render_parfile(template_text, values))
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
        stop_reason: Why the runner asked Cactus to stop, if it did.
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
    stop_reason: Optional[str] = None


# Grace periods when stopping a run: first ask Cactus through the
# TerminationTrigger file (checked every few iterations), then SIGTERM,
# then SIGKILL on the whole process group.
TERMINATE_GRACE_S = 120.0
SIGTERM_GRACE_S = 30.0
STATUS_INTERVAL_S = 15.0
PROBE_INTERVAL_S = 60.0


@dataclass
class CampaignContext:
    """What a run's status needs to know about its campaign.

    Attributes:
        completed_wall_seconds: Wall time of the finished runs.
        remaining_runs: Runs still to start after this one.
    """

    completed_wall_seconds: List[float]
    remaining_runs: int


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
    probe: Optional[Callable[[], Optional[str]]] = None,
    probe_interval_s: float = PROBE_INTERVAL_S,
    campaign: Optional[CampaignContext] = None,
) -> ExecutionResult:
    """Run Cactus in ``run_dir`` with memory monitoring and ETA reporting.

    Streams stdout/stderr into ``run_dir/cactus.log``, feeds every line to
    an ``EtaEstimator``, samples memory with a ``ResourceMonitor``, keeps
    ``status.json`` current, and stops the run (TerminationTrigger file,
    then SIGTERM, then SIGKILL) on timeout or when the memory guard trips.
    ``probe`` is called every ``probe_interval_s``; a non-empty answer
    stops the run the same way (early termination once decided).

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
        probe: Returns a reason to stop early, or ``None``.
        probe_interval_s: Seconds between probe calls.
        campaign: Campaign progress, for the campaign ETA.

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
            "campaign_eta_seconds": None if campaign is None else campaign_eta_seconds(
                campaign.completed_wall_seconds, campaign.remaining_runs, eta_s
            ),
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
    last_probe = time.time()
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
                elif probe is not None and now - last_probe > probe_interval_s:
                    last_probe = now
                    decided = probe()
                    if decided:
                        stop_reason = EARLY_STOP_PREFIX + decided
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
        stop_reason=stop_reason,
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
    if result.exit_code == 0 and result.stop_reason:
        return "stopped early ({})".format(result.stop_reason[len(EARLY_STOP_PREFIX):])
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
# Classification
# ---------------------------------------------------------------------


def load_classifier_config(path: Optional[Path]) -> oc.ClassifierConfig:
    """Classifier thresholds: defaults, overridden by a JSON file.

    Args:
        path: JSON object with ``ClassifierConfig`` field names, or ``None``.

    Returns:
        The configuration.

    Raises:
        ValueError: The file names an unknown threshold.
    """
    config = oc.ClassifierConfig()
    if path is None:
        return config
    overrides = json.loads(path.read_text())
    known = {f.name for f in fields(oc.ClassifierConfig)}
    unknown = sorted(set(overrides) - known)
    if unknown:
        raise ValueError("unknown classifier settings: {}".format(", ".join(unknown)))
    for name, value in overrides.items():
        setattr(config, name, value)
    return config


def classify_run(
    iteration: int,
    amplitude: float,
    run_dir: Path,
    execution: ExecutionResult,
    config: oc.ClassifierConfig,
    role: str = "bisection",
) -> RunResult:
    """Decide the outcome of a finished run from its output and log.

    See ``outcome_classifier`` for the criteria. The order of precedence
    is: memory abort -> ``oom``; confirmed horizon -> ``collapse`` (even
    if the run crashed afterwards, unless the dispersion criteria were
    met before the horizon appeared); crash -> ``crash``; dispersion
    evidence -> ``dispersion``; wall-clock stop -> ``timeout``; anything
    else -> ``undecided``.

    Args:
        iteration: Index of the run in the campaign.
        amplitude: Brill wave amplitude.
        run_dir: Run directory.
        execution: How the Cactus process ended.
        config: Classifier thresholds.
        role: ``endpoint``, ``bisection`` or ``single``.

    Returns:
        The populated run record.
    """
    parfile = run_dir / PARFILE_NAME
    inputs = oc.load_run_inputs(run_dir, parfile.read_text() if parfile.exists() else "", LOG_NAME)
    inputs.exit_code = execution.exit_code
    inputs.guard_tripped = execution.guard_tripped
    inputs.oom_killed = execution.oom_killed
    inputs.timed_out = execution.timed_out
    inputs.t_reached = execution.t_reached
    verdict = oc.classify(inputs, config)

    result = RunResult(
        iteration=iteration,
        amplitude=amplitude,
        run_dir=str(run_dir),
        outcome=verdict.outcome,
        role=role,
        reasons=list(verdict.reasons),
        exit_code=execution.exit_code,
        stop_reason=execution.stop_reason,
        horizons=[asdict(e) for e in verdict.horizons],
        t_reached=execution.t_reached,
        wall_seconds=execution.wall_seconds,
        peak_memory_bytes=execution.peak_memory_bytes,
    )
    horizon = verdict.horizon
    if horizon is not None:
        result.horizon_number = horizon.number
        result.horizon_kind = horizon.kind
        result.m_ah_first = horizon.m_first
        result.t_ah_first = horizon.t_first
        result.t_ah_confirmed = horizon.t_confirmed
        result.m_ah_final = horizon.m_last
        result.t_ah_final = horizon.t_last
        result.m_ah_settled = horizon.m_settled
        result.r_ah_min = horizon.r_min
        result.z_ah = horizon.z_centroid
    curvature = verdict.curvature
    if curvature is not None:
        result.curvature_initial = curvature.k_initial
        result.curvature_max = curvature.k_max
        result.t_curvature_max = curvature.t_max
        result.curvature_final = curvature.k_final
        result.t_dispersed = curvature.t_dispersed
    if inputs.min_lapse:
        result.min_lapse = min(a for _, a in inputs.min_lapse)
        result.final_min_lapse = inputs.min_lapse[-1][1]
    return result


# ---------------------------------------------------------------------
# Bracket and history
# ---------------------------------------------------------------------


class BracketError(RuntimeError):
    """A decisive outcome contradicts the bracket (non-monotone verdicts)."""


def update_bracket(a_min: float, a_max: float, result: RunResult) -> Tuple[float, float]:
    """Narrow the bracket after a decisive run.

    A collapse at or below the subcritical bound, or a dispersion at or
    above the supercritical bound, contradicts monotonicity in A.

    Args:
        a_min: Current subcritical bound.
        a_max: Current supercritical bound.
        result: The run just classified.

    Returns:
        The new ``(a_min, a_max)``.

    Raises:
        RuntimeError: The outcome is not decisive; the campaign stops.
        BracketError: The outcome contradicts the bracket.
    """
    a = result.amplitude
    if result.outcome == oc.COLLAPSE:
        if a <= a_min:
            raise BracketError("collapse at A={} at or below the subcritical bound {}".format(a, a_min))
        return a_min, min(a_max, a)
    if result.outcome == oc.DISPERSION:
        if a >= a_max:
            raise BracketError("dispersion at A={} at or above the supercritical bound {}".format(a, a_max))
        return max(a_min, a), a_max
    raise RuntimeError("outcome '{}' at A={} is not decisive".format(result.outcome, a))


def replay_bracket(a_min: float, a_max: float, runs: List[RunResult]) -> Tuple[float, float]:
    """Bracket after applying every decisive run in order.

    Endpoint runs that confirm their own bound leave it unchanged.

    Raises:
        BracketError: Some run contradicts the bracket.
    """
    for run in runs:
        if run.outcome not in oc.DECISIVE:
            continue
        if run.role == "endpoint" and run.amplitude in (a_min, a_max):
            expected = oc.DISPERSION if run.amplitude == a_min else oc.COLLAPSE
            if run.outcome != expected:
                raise BracketError(
                    "endpoint A={} is {}, expected {}".format(run.amplitude, run.outcome, expected)
                )
            continue
        a_min, a_max = update_bracket(a_min, a_max, run)
    return a_min, a_max


def monotonicity_violation(runs: List[RunResult]) -> Optional[str]:
    """Return a description if some dispersion lies above some collapse."""
    collapses = [r.amplitude for r in runs if r.outcome == oc.COLLAPSE]
    dispersions = [r.amplitude for r in runs if r.outcome == oc.DISPERSION]
    if collapses and dispersions and max(dispersions) >= min(collapses):
        return "dispersion at A={} but collapse at A={}".format(max(dispersions), min(collapses))
    return None


def _sha256(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def campaign_fingerprint(
    template_text: str,
    t_final: float,
    config: oc.ClassifierConfig,
    cactus_bin: str,
    id_settings: InitialDataSettings,
) -> Dict[str, object]:
    """Fingerprint of the numerical setup that every run must share.

    Parallelism (np, OMP) is left out: it does not change the numerical
    family beyond round-off.

    Returns:
        JSON-serialisable components.
    """
    binary = Path(cactus_bin)
    stat = binary.stat() if binary.exists() else None
    fingerprint: Dict[str, object] = {
        "template_sha256": _sha256(template_text),
        "t_final": float(t_final),
        "classifier": config.as_dict(),
        "cactus_bin": None if stat is None else {
            "path": str(binary), "size": stat.st_size, "mtime": int(stat.st_mtime),
        },
        "initial_data": asdict(id_settings) if template_uses_table(template_text) else None,
    }
    if template_uses_table(template_text):
        import brill_initial_data as bid

        fingerprint["initial_data_solver"] = bid.SOLVER_VERSION
    return fingerprint


def fingerprint_differences(stored: Dict[str, object], current: Dict[str, object]) -> List[str]:
    """Names of the fingerprint components that differ."""
    return sorted(k for k in set(stored) | set(current) if stored.get(k) != current.get(k))


def new_history(args: argparse.Namespace, fingerprint: Dict[str, object]) -> dict:
    """A fresh campaign history."""
    return {
        "schema_version": HISTORY_SCHEMA_VERSION,
        "created_at": now_iso(),
        "updated_at": now_iso(),
        "campaign": {
            "a_min_initial": args.a_min,
            "a_max_initial": args.a_max,
            "max_iter": args.max_iter,
            "t_final": args.t_final,
            "template": str(args.template),
            "verify_endpoints": bool(args.verify_endpoints),
        },
        "fingerprint": fingerprint,
        "fingerprint_sha256": _sha256(json.dumps(fingerprint, sort_keys=True)),
        "bracket": {"a_min": args.a_min, "a_max": args.a_max},
        "status": "running",
        "stop_reason": None,
        "runs": [],
    }


def load_history(runs_dir: Path) -> Optional[dict]:
    """Load ``bisection_history.json``, or ``None`` if there is none.

    Raises:
        ValueError: Unsupported schema version.
    """
    path = runs_dir / HISTORY_FILENAME
    if not path.exists():
        return None
    history = json.loads(path.read_text())
    if history.get("schema_version") != HISTORY_SCHEMA_VERSION:
        raise ValueError("{}: unsupported schema_version {}".format(path, history.get("schema_version")))
    return history


def save_history(runs_dir: Path, history: dict) -> None:
    """Atomically replace ``bisection_history.json`` (write, fsync, rename)."""
    runs_dir.mkdir(parents=True, exist_ok=True)
    history["updated_at"] = now_iso()
    target = runs_dir / HISTORY_FILENAME
    tmp = target.with_name(target.name + ".tmp")
    with tmp.open("w") as handle:
        json.dump(history, handle, indent=2)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, target)


def history_runs(history: dict) -> List[RunResult]:
    return [RunResult.from_dict(entry) for entry in history["runs"]]


def apply_overrides(history: dict, overrides: List[str]) -> None:
    """Apply ``--set-outcome ITER=OUTCOME`` human verdicts to the history.

    Raises:
        ValueError: Malformed override, unknown run or non-decisive outcome.
    """
    for item in overrides:
        index_text, _, outcome = item.partition("=")
        if outcome not in oc.DECISIVE:
            raise ValueError("--set-outcome {}: outcome must be collapse or dispersion".format(item))
        index = int(index_text)
        matches = [entry for entry in history["runs"] if entry["iteration"] == index]
        if not matches:
            raise ValueError("--set-outcome {}: no run with iteration {}".format(item, index))
        entry = matches[0]
        if entry.get("human_override") is None:
            entry["human_override"] = entry["outcome"]
        entry["outcome"] = outcome
        entry.setdefault("reasons", []).append("outcome set to {} by a human at {}".format(outcome, now_iso()))


# ---------------------------------------------------------------------
# Running one classified run and the campaign
# ---------------------------------------------------------------------


def run_one(
    args: argparse.Namespace,
    config: oc.ClassifierConfig,
    iteration: int,
    amplitude: float,
    role: str,
    campaign: Optional[CampaignContext],
) -> RunResult:
    """Prepare, execute and classify one run.

    Args:
        args: Parsed command line.
        config: Classifier thresholds.
        iteration: Index of the run in the campaign.
        amplitude: Brill wave amplitude.
        role: ``endpoint``, ``bisection`` or ``single``.
        campaign: Campaign progress for the ETA, or ``None``.

    Returns:
        The classified run.
    """
    runs_dir = args.runs_dir if role != "single" else args.runs_dir / "single"
    started_at = now_iso()
    parfile = prepare_run_dir(
        runs_dir, iteration, amplitude, args.template, args.t_final,
        InitialDataSettings(args.id_n_s, args.id_n_theta, args.id_scale),
    )
    run_dir = parfile.parent
    parfile_text = parfile.read_text()
    command = build_mpirun_command(args.mpirun, args.cactus_bin, parfile, args.np, args.omp)
    run_name = "{} A={}".format(role, format_amplitude(amplitude))
    print("[runner] {}: np={} x OMP={}, t_final={}, log {}".format(
        run_name, args.np, args.omp, args.t_final, run_dir / LOG_NAME), flush=True)

    probe = None
    if not args.no_early_stop:
        def probe() -> Optional[str]:
            try:
                return oc.early_stop_reason(oc.load_run_inputs(run_dir, parfile_text, LOG_NAME), config)
            except (OSError, ValueError, IndexError):
                return None  # files caught mid-write; try again next time

    execution = execute_run(
        command, run_dir, run_name, args.runs_dir, args.t_final, args.timeout_hours * 3600.0,
        args.mem_guard_fraction, args.monitor_interval, probe=probe, campaign=campaign,
    )
    result = classify_run(iteration, amplitude, run_dir, execution, config, role)
    result.started_at = started_at
    result.finished_at = now_iso()
    summary = asdict(result)
    summary["command"] = command
    summary["execution"] = describe_execution(execution)
    (run_dir / RUN_SUMMARY_NAME).write_text(json.dumps(summary, indent=2) + "\n")
    print("[runner] {}: {} ({}) after {}, t={}".format(
        run_name, result.outcome, "; ".join(result.reasons) or "-",
        format_duration(execution.wall_seconds), execution.t_reached), flush=True)
    return result


def _stop(args: argparse.Namespace, history: dict, reason: str) -> int:
    history["status"] = "stopped"
    history["stop_reason"] = reason
    save_history(args.runs_dir, history)
    print("[campaign] stopped for a human decision: {}".format(reason), file=sys.stderr)
    print("[campaign] inspect the run, then resume with --resume (and --set-outcome ITER=collapse|dispersion "
          "or --retry)", file=sys.stderr)
    return 2


def run_campaign(args: argparse.Namespace) -> int:
    """Run (or resume) the bisection campaign.

    Exit codes: 0 when the requested number of bisection steps finished,
    2 when the campaign stopped for a human decision, 1 on usage errors.

    Args:
        args: Parsed command line.

    Returns:
        Process exit code.
    """
    config = load_classifier_config(args.classifier_config)
    template_text = args.template.read_text()
    id_settings = InitialDataSettings(args.id_n_s, args.id_n_theta, args.id_scale)
    fingerprint = campaign_fingerprint(template_text, args.t_final, config, args.cactus_bin, id_settings)

    history = load_history(args.runs_dir)
    if history is None:
        if args.resume:
            print("error: nothing to resume in {}".format(args.runs_dir), file=sys.stderr)
            return 1
        history = new_history(args, fingerprint)
    else:
        if not args.resume:
            print("error: {} already holds a campaign; use --resume or another CAMPAIGN".format(
                args.runs_dir), file=sys.stderr)
            return 1
        differences = fingerprint_differences(history["fingerprint"], fingerprint)
        if differences:
            print("error: the numerical setup changed since the campaign started ({}); runs from two "
                  "setups do not share one A*. Start a new campaign.".format(", ".join(differences)),
                  file=sys.stderr)
            return 1

    campaign_setup = history["campaign"]
    a_min0, a_max0 = campaign_setup["a_min_initial"], campaign_setup["a_max_initial"]
    max_iter = campaign_setup["max_iter"]
    if args.set_outcome:
        try:
            apply_overrides(history, args.set_outcome)
        except ValueError as exc:
            print("error: {}".format(exc), file=sys.stderr)
            return 1

    runs = history_runs(history)
    try:
        a_min, a_max = replay_bracket(a_min0, a_max0, runs)
    except BracketError as exc:
        return _stop(args, history, str(exc))
    history["bracket"] = {"a_min": a_min, "a_max": a_max}

    last = runs[-1] if runs else None
    if last is not None and last.outcome not in oc.DECISIVE and not args.retry:
        return _stop(args, history, "run {} at A={} is {}".format(last.iteration, last.amplitude, last.outcome))
    history["status"] = "running"
    history["stop_reason"] = None
    save_history(args.runs_dir, history)

    def context(remaining: int) -> CampaignContext:
        walls = [r.wall_seconds for r in history_runs(history) if r.wall_seconds]
        return CampaignContext(walls, remaining)

    def record(result: RunResult) -> None:
        history["runs"].append(asdict(result))
        save_history(args.runs_dir, history)

    # Endpoints first: the bisection is meaningless unless a_min disperses
    # and a_max collapses in this numerical setup.
    if campaign_setup.get("verify_endpoints"):
        done = {(r.amplitude, r.outcome) for r in history_runs(history) if r.role == "endpoint"}
        for amplitude, expected in ((a_min0, oc.DISPERSION), (a_max0, oc.COLLAPSE)):
            if (amplitude, expected) in done:
                continue
            bisection_left = max_iter - sum(1 for r in history_runs(history) if r.role == "bisection")
            result = run_one(args, config, len(history["runs"]), amplitude, "endpoint",
                             context(bisection_left + (1 if amplitude == a_min0 else 0)))
            record(result)
            if result.outcome != expected:
                return _stop(args, history, "endpoint A={} is {}, expected {}".format(
                    amplitude, result.outcome, expected))

    while True:
        steps_done = sum(1 for r in history_runs(history) if r.role == "bisection")
        if steps_done >= max_iter:
            break
        amplitude = 0.5 * (a_min + a_max)
        result = run_one(args, config, len(history["runs"]), amplitude, "bisection",
                         context(max_iter - steps_done - 1))
        record(result)
        if result.outcome not in oc.DECISIVE:
            return _stop(args, history, "run {} at A={} is {}: {}".format(
                result.iteration, amplitude, result.outcome, "; ".join(result.reasons)))
        try:
            a_min, a_max = update_bracket(a_min, a_max, result)
        except BracketError as exc:
            return _stop(args, history, str(exc))
        violation = monotonicity_violation(history_runs(history))
        if violation:
            return _stop(args, history, "non-monotone verdicts: " + violation)
        history["bracket"] = {"a_min": a_min, "a_max": a_max}
        save_history(args.runs_dir, history)
        print("[campaign] bracket [{}, {}] (width {:.3g})".format(a_min, a_max, a_max - a_min), flush=True)

    history["status"] = "completed"
    save_history(args.runs_dir, history)
    print("[campaign] completed: A* in [{}, {}]".format(a_min, a_max))
    return 0


def run_single(args: argparse.Namespace) -> int:
    """Run and classify one amplitude outside the bisection (smoke tests, benchmarks).

    Args:
        args: Parsed command line; ``args.single`` is the amplitude.

    Returns:
        Process exit code: 0 for a decisive outcome, 2 otherwise.
    """
    target = args.runs_dir / "single" / run_dir_name(0, args.single)
    if target.exists():
        print("error: {} already exists; remove it first".format(target), file=sys.stderr)
        return 1
    config = load_classifier_config(args.classifier_config)
    result = run_one(args, config, 0, args.single, "single", None)
    return 0 if result.outcome in oc.DECISIVE else 2


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
    id_settings = InitialDataSettings(args.id_n_s, args.id_n_theta, args.id_scale)
    parfile = prepare_run_dir(dry_dir, 0, amplitude, args.template, args.t_final, id_settings)
    command = build_mpirun_command(
        args.mpirun, args.cactus_bin, parfile, args.np, args.omp
    )
    config = load_classifier_config(args.classifier_config)
    fingerprint = campaign_fingerprint(
        args.template.read_text(), args.t_final, config, args.cactus_bin, id_settings
    )
    print("bracket            : [{}, {}]".format(args.a_min, args.a_max))
    print("max iterations     : {}".format(args.max_iter))
    print("final bracket width: {:.3g}".format((args.a_max - args.a_min) / 2 ** args.max_iter))
    print("first amplitude    : {}".format(format_amplitude(amplitude)))
    print("t_final            : {}".format(args.t_final))
    print("parallelism        : np={} x OMP={}".format(args.np, args.omp))
    print("memory guard       : {:.0%} of the container limit".format(args.mem_guard_fraction))
    print("timeout per run    : {} h".format(args.timeout_hours))
    print("endpoint runs      : {}".format("yes" if args.verify_endpoints else "no"))
    print("early stop         : {}".format("no" if args.no_early_stop else "yes"))
    print("classifier         : {}".format(", ".join(
        "{}={}".format(k, v) for k, v in config.as_dict().items())))
    print("setup fingerprint  : {}".format(_sha256(json.dumps(fingerprint, sort_keys=True))[:16]))
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
    campaign = parser.add_argument_group("campaign")
    campaign.add_argument(
        "--verify-endpoints", action="store_true",
        help="run --a-min and --a-max first and require dispersion and collapse there",
    )
    campaign.add_argument(
        "--set-outcome", action="append", metavar="ITER=OUTCOME", default=[],
        help="with --resume: record a human verdict (collapse or dispersion) for run ITER",
    )
    campaign.add_argument(
        "--retry", action="store_true",
        help="with --resume: run the next midpoint again after a non-decisive run",
    )
    campaign.add_argument(
        "--no-early-stop", action="store_true",
        help="run every amplitude to t_final instead of stopping once decided",
    )
    campaign.add_argument(
        "--classifier-config", type=Path, metavar="JSON",
        help="override classifier thresholds (see outcome_classifier.ClassifierConfig)",
    )
    initial_data = parser.add_argument_group("initial data (templates with __TABLE_FILE__)")
    initial_data.add_argument("--id-n-s", type=int, default=512, help="cells in s (default: 512)")
    initial_data.add_argument("--id-n-theta", type=int, default=64, help="cells in theta (default: 64)")
    initial_data.add_argument("--id-scale", type=float, default=2.0, help="L in r = L tan(pi s/2) (default: 2)")
    args = parser.parse_args(argv)
    if not args.a_min < args.a_max:
        parser.error("--a-min must be smaller than --a-max")
    if not 0.0 < args.mem_guard_fraction < 1.0:
        parser.error("--mem-guard-fraction must lie in (0, 1)")
    if (args.set_outcome or args.retry) and not args.resume:
        parser.error("--set-outcome and --retry need --resume")
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
