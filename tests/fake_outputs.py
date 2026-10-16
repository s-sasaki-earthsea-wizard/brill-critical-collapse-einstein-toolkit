# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
"""Synthetic Cactus output for the classifier and campaign tests."""

import math

import cactus_outputs as co

BH_COLUMNS = [
    "cctk_iteration", "cctk_time", "centroid_x", "centroid_y", "centroid_z",
    "min radius", "max radius", "mean radius",
] + ["filler{}".format(k) for k in range(18)] + ["m_irreducible", "areal radius"]

# Iterations per unit of simulation time (finest-level step 0.00625).
ITERATIONS_PER_TIME = 160
FIND_EVERY = 16


def horizon_series(finds):
    """HorizonSeries from (iteration, time, m_irr, min_radius, centroid_z) tuples."""
    return co.HorizonSeries(
        iterations=[f[0] for f in finds],
        times=[f[1] for f in finds],
        irreducible_masses=[f[2] for f in finds],
        areal_radii=[2.0 * f[2] for f in finds],
        mean_radii=[f[3] for f in finds],
        min_radii=[f[3] for f in finds],
        centroids_z=[f[4] for f in finds],
    )


def growing_horizon(t_form, t_end, m_final, rise_time=1.0, radius=0.3, z=0.0, dt=0.1):
    """Finds every ``dt`` from ``t_form`` with a mass relaxing to ``m_final``."""
    finds = []
    t = t_form
    while t <= t_end + 1e-9:
        mass = m_final * (1.0 - 0.1 * math.exp(-(t - t_form) / rise_time * 5.0))
        finds.append((int(round(t * ITERATIONS_PER_TIME)), round(t, 6), mass, radius, z))
        t += dt
    return finds


def write_bh_diagnostics(output_dir, number, finds):
    output_dir.mkdir(parents=True, exist_ok=True)
    lines = ["# column {:2d} = {}".format(k + 1, name) for k, name in enumerate(BH_COLUMNS)]
    for iteration, time, mass, radius, z in finds:
        values = [0.0] * len(BH_COLUMNS)
        values[0], values[1], values[4] = iteration, time, z
        values[5] = values[6] = values[7] = radius
        values[26], values[27] = mass, 2.0 * mass
        lines.append("\t".join(repr(v) for v in values))
    (output_dir / co.BH_DIAGNOSTICS_PATTERN.format(number)).write_text("\n".join(lines) + "\n")


def write_scalar(output_dir, filename, series):
    output_dir.mkdir(parents=True, exist_ok=True)
    lines = ["# Scalar ASCII output created by CarpetIOScalar", "# 1:iteration 2:time 3:data"]
    for time, value in series:
        lines.append("{} {} {}".format(int(round(time * ITERATIONS_PER_TIME)), time, value))
    (output_dir / filename).write_text("\n".join(lines) + "\n")


def times(t_end, dt=0.05):
    return [round(k * dt, 6) for k in range(int(round(t_end / dt)) + 1)]


def dispersing_curvature(t_end, k0=1.0, t_peak=2.0, peak=3.0, t_decay=1.0):
    """Curvature that rises to a peak and then decays well below k0."""
    series = []
    for t in times(t_end):
        if t <= t_peak:
            value = k0 + (peak - k0) * t / t_peak
        else:
            value = peak * math.exp(-(t - t_peak) / t_decay * 2.0)
        series.append((t, value))
    return series


def recovering_lapse(t_end, dip=0.3, t_dip=2.5):
    series = []
    for t in times(t_end):
        value = 1.0 - (1.0 - dip) * math.exp(-((t - t_dip) / 1.0) ** 2)
        series.append((t, value))
    return series


def collapsing_lapse(t_end):
    return [(t, max(1e-6, 1.0 - t / 3.0)) for t in times(t_end)]


def write_log(run_dir, t_end, extra=""):
    lines = ["INFO (Carpet): Carpet is running on 4 processes"]
    for t in times(t_end, 0.05):
        lines.append("{:>9d} {:>9.3f} |   1.0   1.0 |".format(int(round(t * ITERATIONS_PER_TIME)), t))
    (run_dir / "cactus.log").write_text("\n".join(lines) + "\n" + extra)
