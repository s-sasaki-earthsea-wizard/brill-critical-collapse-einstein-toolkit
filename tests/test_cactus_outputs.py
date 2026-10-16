# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
"""Tests for the Cactus output readers."""

import pytest

import cactus_outputs as co

BH_HEADER_NAMES = [
    "cctk_iteration", "cctk_time", "centroid_x", "centroid_y", "centroid_z",
    "min radius", "max radius", "mean radius",
] + ["filler{}".format(k) for k in range(18)] + ["m_irreducible", "areal radius"]


def write_bh_diagnostics(path, rows):
    header = ["# column {:2d} = {}".format(k + 1, name) for k, name in enumerate(BH_HEADER_NAMES)]
    body = []
    for iteration, time, m_irr, centroid_z, mean_radius in rows:
        values = [0.0] * len(BH_HEADER_NAMES)
        values[0], values[1], values[4] = iteration, time, centroid_z
        values[5] = values[6] = values[7] = mean_radius
        values[26], values[27] = m_irr, 2.0 * m_irr
        body.append("\t".join(repr(v) for v in values))
    path.write_text("\n".join(header + body) + "\n")


def test_parse_column_header():
    columns = co.parse_column_header(["# column  1 = cctk_iteration", "# column 27 = m_irreducible"])
    assert columns == {"cctk_iteration": 0, "m_irreducible": 26}


def test_read_bh_diagnostics(tmp_path):
    write_bh_diagnostics(
        tmp_path / "BH_diagnostics.ah2.gp",
        [(16, 0.1, 0.25, 0.6125, 0.12), (32, 0.2, 0.2501, 0.61, 0.13)],
    )
    series = co.read_bh_diagnostics(tmp_path, horizon=2)
    assert series.iterations == [16, 32]
    assert series.irreducible_masses == [0.25, 0.2501]
    assert series.areal_radii == [0.5, 0.5002]
    assert series.centroids_z == [0.6125, 0.61]
    assert series.mean_radii == [0.12, 0.13]
    assert series.found


def test_read_bh_diagnostics_missing_file(tmp_path):
    assert co.read_bh_diagnostics(tmp_path) is None


def test_read_scalar_series(tmp_path):
    (tmp_path / "x.asc").write_text(
        "# Scalar ASCII output created by CarpetIOScalar\n"
        "# 1:iteration 2:time 3:data\n"
        "0 0 1\n8 0.05 0.99\n16 0.1 0.95\n"
    )
    assert co.read_scalar_series(tmp_path, "x.asc") == [(0.0, 1.0), (0.05, 0.99), (0.1, 0.95)]


def test_read_central_lapse_keeps_finest_level(tmp_path):
    rows = [
        (0, 0, 0.0, 1.0), (0, 1, 0.0, 1.0),
        (8, 0, 0.05, 0.90), (8, 1, 0.05, 0.80), (8, 3, 0.05, 0.70),
        (12, 3, 0.075, 0.65),
    ]
    lines = ["# 0D ASCII output created by CarpetIOASCII"]
    for iteration, level, time, value in rows:
        lines.append("{}\t0\t{} 0 0\t24 24 24\t{}\t0 0 0\t{}".format(iteration, level, time, value))
        lines.append("")
    (tmp_path / co.CENTRAL_LAPSE_FILENAME).write_text("\n".join(lines) + "\n")
    assert co.read_central_lapse(tmp_path) == [(0.0, 1.0), (0.05, 0.70), (0.075, 0.65)]


def test_scan_log(tmp_path):
    log = tmp_path / "cactus.log"
    log.write_text(
        "        0     0.000 |    1.0000000 |\n"
        "        8     0.050 |    0.9999000 |\n"
        "  -> There were 827 NaN/Inf value(s) found in variable 'ML_BSSN::Xt1'\n"
        "  -> There were 210 NaN/Inf value(s) found in variable 'ML_BSSN::Xt1'\n"
        "  -> There were 3 NaN/Inf value(s) found in variable 'ML_BSSN::Xt2'\n"
        "application called MPI_Abort(MPI_COMM_WORLD, 1) - process 0\n"
    )
    failures = co.scan_log(log)
    assert failures.nan_reports == ["ML_BSSN::Xt1", "ML_BSSN::Xt2"]
    assert failures.mpi_abort
    assert not failures.killed
    assert (failures.last_iteration, failures.last_time) == (8, pytest.approx(0.05))
