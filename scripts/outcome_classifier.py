#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
"""Three-way outcome classification of a Brill wave collapse run.

A bisection on the amplitude needs a verdict that is monotone in A and
correct for the numerical spacetime. Both ways of getting it wrong bias
A* upwards: missing a horizon calls a collapse a dispersion, and calling
dispersion too early does the same. So dispersion is never inferred from
the absence of a horizon; it needs positive evidence of its own, and
anything short of either kind of evidence is ``undecided``.

Collapse (all of)
    C1  one horizon finder found a horizon at ``min_consecutive_finds``
        consecutive find times;
    C2  its irreducible mass did not drop between those finds by more
        than ``mass_decrease_tolerance`` (relative), as the area theorem
        requires of a real dynamical horizon;
    C3  its minimum coordinate radius was at least ``min_horizon_radius``
        (disabled at 0 until the resolution floor is calibrated).

Dispersion (all of)
    D1  the global maximum of the curvature invariant stayed below
        ``1 / curvature_drop_factor`` of its initial value over a window
        of ``dispersion_window``;
    D2  no horizon finder found anything, ever;
    D3  the window ends at or after ``dispersion_min_time``;
    D4  the curvature did not grow back above that level afterwards;
    D5  the minimum lapse is at least ``lapse_recovered`` at the end and
        not falling over the last window (gauge, so only supporting).

The verdict reports the finder confirmed first. When several are confirmed
at the same time, it reports the one with the largest mass: a horizon
forms as a pair of marginally trapped surfaces, an outer one that grows
and an inner one that shrinks, and the apparent horizon is the outer one.

A horizon that appears after D1-D3 were met is a contradiction and makes
the run ``undecided``. A memory abort is ``oom``. A crash (NaN, MPI
abort, kill, failed exit) is ``crash`` unless the collapse evidence was
already complete. A wall-clock stop without evidence is ``timeout``.

The thresholds are provisional until the pilot campaign (Phase 2.5)
calibrates them; they belong to the frozen numerical setup.

The same evidence drives early termination: ``early_stop_reason`` tells
the runner to stop once the horizon mass has settled or dispersion is
established.
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import cactus_outputs as co

COLLAPSE = "collapse"
DISPERSION = "dispersion"
CRASH = "crash"
OOM = "oom"
TIMEOUT = "timeout"
UNDECIDED = "undecided"
DECISIVE = (COLLAPSE, DISPERSION)

ORIGIN_Z_PATTERN = re.compile(
    r"^\s*AHFinderDirect::origin_z\[(\d+)\]\s*=\s*([-+0-9.eE]+)", re.IGNORECASE | re.MULTILINE
)
N_HORIZONS_PATTERN = re.compile(r"^\s*AHFinderDirect::N_horizons\s*=\s*(\d+)", re.IGNORECASE | re.MULTILINE)


@dataclass
class ClassifierConfig:
    """Thresholds of the classifier (provisional until calibrated).

    Attributes:
        min_consecutive_finds: C1.
        mass_decrease_tolerance: C2, relative.
        min_horizon_radius: C3, coordinate radius; 0 disables.
        find_interval: Iterations between horizon finds; ``None`` infers
            it as the smallest gap between finds.
        settle_rate_window: Time span, ending at the last find, over which
            the rate of change of M_irr is fitted.
        settle_rate_threshold: The mass counts as settled once
            |dM/dt| / M over that span falls below this (per unit time).
            A Brill black hole keeps accreting the rest of the wave slowly
            (A = 5: about 0.3 % over 6 time units), so a criterion on the
            rate rather than on the total change.
        curvature_drop_factor: D1.
        dispersion_window: D1 window.
        dispersion_min_time: D3.
        lapse_recovered: D5.
    """

    min_consecutive_finds: int = 3
    mass_decrease_tolerance: float = 1.0e-3
    min_horizon_radius: float = 0.0
    find_interval: Optional[int] = None
    settle_rate_window: float = 1.0
    settle_rate_threshold: float = 2.0e-3
    curvature_drop_factor: float = 10.0
    dispersion_window: float = 2.0
    dispersion_min_time: float = 5.0
    lapse_recovered: float = 0.5

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


# ---------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------


@dataclass
class HorizonTrack:
    """The finds of one AHFinderDirect horizon.

    Attributes:
        number: Horizon number (``ah<number>``).
        origin_z: z of the finder's origin; 0 is the origin finder.
        series: The finds.
    """

    number: int
    origin_z: float
    series: co.HorizonSeries

    @property
    def kind(self) -> str:
        return "origin" if self.origin_z == 0.0 else "axis"


@dataclass
class RunInputs:
    """Everything the classifier looks at, read from one run.

    Attributes:
        horizons: One track per finder that wrote a diagnostics file.
        curvature: (time, global maximum of |curvature invariant|).
        min_lapse: (time, global minimum of the lapse).
        log: Failure signatures from the run log.
        exit_code: mpirun exit code.
        guard_tripped: Memory guard stopped the run.
        oom_killed: Kernel OOM killer fired.
        timed_out: Wall-clock limit stopped the run.
        t_reached: Last simulation time reached.
    """

    horizons: List[HorizonTrack] = field(default_factory=list)
    curvature: List[Tuple[float, float]] = field(default_factory=list)
    min_lapse: List[Tuple[float, float]] = field(default_factory=list)
    log: Optional[co.LogFailures] = None
    exit_code: Optional[int] = 0
    guard_tripped: bool = False
    oom_killed: bool = False
    timed_out: bool = False
    t_reached: Optional[float] = None


def finder_origins(parfile_text: str) -> Dict[int, float]:
    """z of every AHFinderDirect origin set in a parfile (default 0)."""
    count = N_HORIZONS_PATTERN.search(parfile_text)
    origins = {n: 0.0 for n in range(1, int(count.group(1)) + 1)} if count else {}
    for match in ORIGIN_Z_PATTERN.finditer(parfile_text):
        origins[int(match.group(1))] = float(match.group(2))
    return origins


def load_run_inputs(run_dir: Path, parfile_text: str, log_name: str = "cactus.log") -> RunInputs:
    """Read the classifier inputs from a run directory.

    Missing files give empty series; the classifier decides what that
    means. Execution facts (exit code, guard, timeout) are filled in by
    the caller.

    Args:
        run_dir: Run directory with ``output/`` and the log.
        parfile_text: The run's parfile, for the finder origins.
        log_name: Log file name inside ``run_dir``.

    Returns:
        The inputs.
    """
    output = run_dir / "output"
    origins = finder_origins(parfile_text)
    numbers = set(origins)
    for path in output.glob("BH_diagnostics.ah*.gp"):
        match = re.match(r"BH_diagnostics\.ah(\d+)\.gp$", path.name)
        if match:
            numbers.add(int(match.group(1)))
    horizons = []
    for number in sorted(numbers):
        series = co.read_bh_diagnostics(output, number)
        if series is not None and series.iterations:
            horizons.append(HorizonTrack(number, origins.get(number, 0.0), series))

    def scalar(filename: str) -> List[Tuple[float, float]]:
        return co.read_scalar_series(output, filename) if (output / filename).exists() else []

    log_path = run_dir / log_name
    return RunInputs(
        horizons=horizons,
        curvature=[(t, abs(v)) for t, v in scalar(co.CURVATURE_NORM_INF_FILENAME)],
        min_lapse=scalar(co.MIN_LAPSE_FILENAME),
        log=co.scan_log(log_path) if log_path.exists() else None,
    )


# ---------------------------------------------------------------------
# Evidence
# ---------------------------------------------------------------------


@dataclass
class HorizonEvidence:
    """Collapse evidence from one horizon finder.

    Attributes:
        number: Horizon number.
        kind: ``origin`` or ``axis``.
        n_finds: Number of finds.
        t_first: Time of the first find.
        m_first: M_irr at the first find.
        t_confirmed: Time at which C1-C3 were first met; ``None`` if never.
        t_last: Time of the last find.
        m_last: M_irr at the last find.
        m_settled: Settled M_irr; ``None`` while the mass still changes.
        r_min: Smallest coordinate radius over the finds after confirmation.
        z_centroid: Centroid z at the last find.
        failed: Why C1-C3 were not met, if they were not.
    """

    number: int
    kind: str
    n_finds: int
    t_first: float
    m_first: float
    t_confirmed: Optional[float]
    t_last: float
    m_last: float
    m_settled: Optional[float]
    r_min: float
    z_centroid: float
    failed: str = ""


def _find_interval(horizons: Sequence[HorizonTrack], config: ClassifierConfig) -> Optional[int]:
    if config.find_interval is not None:
        return config.find_interval
    gaps = [
        b - a
        for track in horizons
        for a, b in zip(track.series.iterations, track.series.iterations[1:])
        if b > a
    ]
    return min(gaps) if gaps else None


def settled_mass(times: Sequence[float], masses: Sequence[float], config: ClassifierConfig) -> Optional[float]:
    """Last M_irr once its relative rate of change has dropped below the threshold.

    The rate is the least-squares slope of M_irr over the finds in the
    last ``settle_rate_window``, divided by the last mass. The finds must
    cover the whole window and there must be at least three of them.
    """
    if not masses or masses[-1] <= 0.0:
        return None
    m_last, t_last = masses[-1], times[-1]
    window = config.settle_rate_window
    if t_last - times[0] < window:
        return None
    recent = [(t, m) for t, m in zip(times, masses) if t >= t_last - window - 1e-9]
    if len(recent) < 3:
        return None
    t_mean = sum(t for t, _ in recent) / len(recent)
    m_mean = sum(m for _, m in recent) / len(recent)
    denominator = sum((t - t_mean) ** 2 for t, _ in recent)
    slope = sum((t - t_mean) * (m - m_mean) for t, m in recent) / denominator
    if abs(slope) / m_last <= config.settle_rate_threshold:
        return m_last
    return None


def horizon_evidence(track: HorizonTrack, interval: Optional[int], config: ClassifierConfig) -> HorizonEvidence:
    """Check C1-C3 for one finder and summarise its finds."""
    s = track.series
    evidence = HorizonEvidence(
        number=track.number,
        kind=track.kind,
        n_finds=len(s.iterations),
        t_first=s.times[0],
        m_first=s.irreducible_masses[0],
        t_confirmed=None,
        t_last=s.times[-1],
        m_last=s.irreducible_masses[-1],
        m_settled=None,
        r_min=min(s.min_radii),
        z_centroid=s.centroids_z[-1],
    )
    need = config.min_consecutive_finds
    if interval is None and need > 1:
        evidence.failed = "too few finds to know the find interval"
        return evidence

    # Length of the current streak of good, consecutive finds ending at k.
    streak_start: Optional[int] = None
    for k in range(len(s.iterations)):
        if not (s.irreducible_masses[k] > 0.0 and s.min_radii[k] >= config.min_horizon_radius):
            streak_start = None
            continue
        if streak_start is None:
            streak_start = k
        else:
            consecutive = s.iterations[k] - s.iterations[k - 1] == interval
            drop = (s.irreducible_masses[k - 1] - s.irreducible_masses[k]) / s.irreducible_masses[k - 1]
            if not (consecutive and drop <= config.mass_decrease_tolerance):
                streak_start = k
        if k - streak_start + 1 >= need:
            evidence.t_confirmed = s.times[k]
            confirmed_from = streak_start
            break
    else:
        evidence.failed = "no {} consecutive finds with non-decreasing mass{}".format(
            need, "" if config.min_horizon_radius == 0.0 else " and radius >= {}".format(config.min_horizon_radius)
        )
        return evidence

    times = s.times[confirmed_from:]
    masses = s.irreducible_masses[confirmed_from:]
    evidence.r_min = min(s.min_radii[confirmed_from:])
    evidence.m_settled = settled_mass(times, masses, config)
    return evidence


@dataclass
class CurvatureEvidence:
    """Dispersion evidence from the curvature invariant.

    Attributes:
        k_initial: Value at the first sample.
        k_max: Maximum over the run.
        t_max: Time of the maximum.
        k_final: Value at the last sample.
        threshold: k_initial / curvature_drop_factor.
        t_dispersed: End of the first window that satisfies D1 and D3.
        regrew: Whether the curvature exceeded the threshold after that.
    """

    k_initial: float
    k_max: float
    t_max: float
    k_final: float
    threshold: float
    t_dispersed: Optional[float]
    regrew: bool

    @property
    def length_scale(self) -> float:
        """Curvature length k_max^(-1/4) (the invariant scales as length^-4)."""
        return self.k_max ** -0.25 if self.k_max > 0.0 else float("inf")


def curvature_evidence(curvature: Sequence[Tuple[float, float]], config: ClassifierConfig) -> Optional[CurvatureEvidence]:
    """Check D1, D3 and D4 on the curvature series."""
    if not curvature:
        return None
    times = [t for t, _ in curvature]
    values = [k for _, k in curvature]
    k_initial = values[0]
    i_max = max(range(len(values)), key=values.__getitem__)
    threshold = k_initial / config.curvature_drop_factor
    t_dispersed = None
    for k, t_end in enumerate(times):
        t_start = t_end - config.dispersion_window
        if t_end < config.dispersion_min_time or t_start < times[0]:
            continue
        window = [v for t, v in zip(times, values) if t_start <= t <= t_end]
        if window and max(window) <= threshold:
            t_dispersed = t_end
            break
    regrew = t_dispersed is not None and any(v > threshold for t, v in curvature if t >= t_dispersed)
    return CurvatureEvidence(
        k_initial=k_initial,
        k_max=values[i_max],
        t_max=times[i_max],
        k_final=values[-1],
        threshold=threshold,
        t_dispersed=t_dispersed,
        regrew=regrew,
    )


def lapse_recovered(min_lapse: Sequence[Tuple[float, float]], config: ClassifierConfig) -> Tuple[bool, str]:
    """D5: the minimum lapse is high enough and not falling at the end."""
    if not min_lapse:
        return False, "no minimum-lapse output"
    t_end, a_end = min_lapse[-1]
    if a_end < config.lapse_recovered:
        return False, "minimum lapse {:.3g} below {}".format(a_end, config.lapse_recovered)
    window = [a for t, a in min_lapse if t >= t_end - config.dispersion_window]
    if a_end < window[0]:
        return False, "minimum lapse still falling at the end"
    return True, ""


# ---------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------


@dataclass
class Verdict:
    """The classifier's decision with its evidence.

    Attributes:
        outcome: One of the outcome strings.
        reasons: Human-readable reasons, most important first.
        horizon: Evidence of the confirmed horizon (earliest), if any.
        horizons: Evidence of every finder that found something.
        curvature: Curvature evidence, if the run wrote the invariant.
    """

    outcome: str
    reasons: List[str]
    horizon: Optional[HorizonEvidence] = None
    horizons: List[HorizonEvidence] = field(default_factory=list)
    curvature: Optional[CurvatureEvidence] = None


def _collapse(inputs: RunInputs, config: ClassifierConfig) -> Tuple[Optional[HorizonEvidence], List[HorizonEvidence]]:
    interval = _find_interval(inputs.horizons, config)
    every = [horizon_evidence(track, interval, config) for track in inputs.horizons]
    confirmed = [e for e in every if e.t_confirmed is not None]
    # Ties go to the outermost surface (largest mass), not the inner MOTS.
    first = min(confirmed, key=lambda e: (e.t_confirmed, -e.m_last, e.number)) if confirmed else None
    return first, every


def _crashed(inputs: RunInputs) -> Optional[str]:
    log = inputs.log
    if log is not None and log.nan_reports:
        return "NaN in " + ", ".join(log.nan_reports[:4])
    if log is not None and log.mpi_abort:
        return "MPI abort"
    if log is not None and log.killed:
        return "process killed"
    if inputs.exit_code not in (0, None) and not inputs.timed_out:
        return "exit code {}".format(inputs.exit_code)
    return None


def classify(inputs: RunInputs, config: ClassifierConfig) -> Verdict:
    """Decide collapse, dispersion or one of the non-decisive outcomes.

    Args:
        inputs: Outputs and execution facts of one run.
        config: Thresholds.

    Returns:
        The verdict.
    """
    horizon, horizons = _collapse(inputs, config)
    curvature = curvature_evidence(inputs.curvature, config)
    verdict = Verdict(UNDECIDED, [], horizon, horizons, curvature)

    if inputs.guard_tripped or inputs.oom_killed:
        verdict.outcome = OOM
        verdict.reasons.append("memory guard" if inputs.guard_tripped else "kernel OOM kill")
        return verdict

    crash = _crashed(inputs)
    if horizon is not None:
        if curvature is not None and curvature.t_dispersed is not None and curvature.t_dispersed < horizon.t_first:
            verdict.reasons.append(
                "horizon {} found at t={:.4g} after the dispersion criteria were met at t={:.4g}".format(
                    horizon.number, horizon.t_first, curvature.t_dispersed
                )
            )
            return verdict
        verdict.outcome = COLLAPSE
        verdict.reasons.append(
            "horizon {} ({}) confirmed at t={:.4g}, M_irr={:.6g}{}".format(
                horizon.number, horizon.kind, horizon.t_confirmed, horizon.m_last,
                "" if horizon.m_settled is not None else " (not settled)",
            )
        )
        if crash:
            verdict.reasons.append("run crashed afterwards: " + crash)
        return verdict

    if crash:
        verdict.outcome = CRASH
        verdict.reasons.append(crash)
        return verdict

    dispersion_failures = []
    if horizons:
        dispersion_failures.append(
            "horizon finds without confirmation ({})".format(
                "; ".join("ah{}: {}".format(e.number, e.failed) for e in horizons)
            )
        )
    if curvature is None:
        dispersion_failures.append("no curvature output, dispersion cannot be shown")
    elif curvature.t_dispersed is None:
        dispersion_failures.append(
            "curvature not below {:.3g} (initial / {}) for {} after t={}".format(
                curvature.threshold, config.curvature_drop_factor,
                config.dispersion_window, config.dispersion_min_time,
            )
        )
    elif curvature.regrew:
        dispersion_failures.append("curvature grew back after t={:.4g}".format(curvature.t_dispersed))
    lapse_ok, lapse_reason = lapse_recovered(inputs.min_lapse, config)
    if not lapse_ok:
        dispersion_failures.append(lapse_reason)

    if not dispersion_failures:
        verdict.outcome = DISPERSION
        verdict.reasons.append(
            "curvature below initial / {} from t={:.4g}, peak {:.4g} at t={:.4g}".format(
                config.curvature_drop_factor, curvature.t_dispersed, curvature.k_max, curvature.t_max
            )
        )
        return verdict

    if inputs.timed_out:
        verdict.outcome = TIMEOUT
    verdict.reasons.extend(dispersion_failures)
    return verdict


def early_stop_reason(inputs: RunInputs, config: ClassifierConfig) -> Optional[str]:
    """Why a running simulation may stop now, or ``None`` to keep going.

    Stops once a confirmed horizon's mass has settled, or once dispersion
    is established (D1-D5 hold and nothing has been found).
    """
    verdict = classify(inputs, config)
    if verdict.outcome == COLLAPSE and verdict.horizon is not None and verdict.horizon.m_settled is not None:
        return "collapse, M_irr settled at {:.6g}".format(verdict.horizon.m_settled)
    if verdict.outcome == DISPERSION:
        return "dispersion established"
    return None
