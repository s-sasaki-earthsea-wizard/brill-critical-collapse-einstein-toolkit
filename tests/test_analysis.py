# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
"""Tests for the critical exponent analysis."""

import json

import numpy as np
import pytest

import analyze_critical_exponent as ace

A_STAR, GAMMA = 4.7, 0.37


def super_points(deltas, a_star=A_STAR, gamma=GAMMA, scale=0.5):
    return [ace.ScalingPoint(a_star + d, scale * d ** gamma, "r", "super") for d in deltas]


def sub_points(deltas, a_star=A_STAR, gamma=GAMMA, scale=2.0):
    return [ace.ScalingPoint(a_star - d, scale * d ** (-4.0 * gamma), "r", "sub") for d in deltas]


DELTAS = list(np.logspace(-3, -0.5, 6))


def test_supercritical_fit_recovers_gamma():
    fit = ace.fit_power_law(super_points(DELTAS), A_STAR)
    assert fit.gamma == pytest.approx(GAMMA, abs=1e-10)
    assert fit.r_squared == pytest.approx(1.0)
    assert fit.gamma_stderr < 1e-8


def test_subcritical_fit_recovers_gamma_from_curvature():
    fit = ace.fit_power_law(sub_points(DELTAS), A_STAR)
    assert fit.gamma == pytest.approx(GAMMA, abs=1e-10)
    assert fit.slope == pytest.approx(-4.0 * GAMMA)


def test_fit_rejects_bad_input():
    with pytest.raises(ValueError):
        ace.fit_power_law(super_points(DELTAS[:1]), A_STAR)
    with pytest.raises(ValueError):
        ace.fit_power_law(super_points(DELTAS), A_STAR + 0.01)  # closest point below A*
    with pytest.raises(ValueError):
        ace.fit_power_law(super_points(DELTAS[:2]) + sub_points(DELTAS[:2]), A_STAR)


@pytest.mark.parametrize("side", ["super", "sub"])
def test_free_a_star_fit(side):
    points = super_points(DELTAS) if side == "super" else sub_points(DELTAS)
    fit = ace.fit_power_law_free_a_star(points, A_STAR - 0.05 if side == "super" else A_STAR + 0.05)
    assert fit.a_star == pytest.approx(A_STAR, abs=1e-7)
    assert fit.gamma == pytest.approx(GAMMA, abs=1e-6)
    with pytest.raises(ValueError):
        ace.fit_power_law_free_a_star(points[:3], A_STAR)


def write_history(path, runs, bracket=(4.69, 4.71)):
    path.write_text(json.dumps({"bracket": {"a_min": bracket[0], "a_max": bracket[1]}, "runs": runs}))
    return path


def history_runs():
    runs = []
    for d in DELTAS:
        runs.append({"amplitude": A_STAR + d, "outcome": "collapse", "run_dir": "c",
                     "m_ah_settled": 0.5 * d ** GAMMA, "m_ah_final": 0.6 * d ** GAMMA, "m_ah_first": None,
                     "curvature_max": 99.0})
        runs.append({"amplitude": A_STAR - d, "outcome": "dispersion", "run_dir": "d",
                     "curvature_max": 2.0 * d ** (-4 * GAMMA)})
    runs.append({"amplitude": 4.705, "outcome": "undecided", "m_ah_settled": 0.1})
    runs.append({"amplitude": 4.8, "outcome": "collapse", "m_ah_settled": None, "m_ah_final": 0.2})
    return runs


def test_load_points_by_side_and_mass(tmp_path):
    path = write_history(tmp_path / "h.json", history_runs())
    assert len(ace.load_points([path], "settled")) == len(DELTAS)
    assert len(ace.load_points([path], "final")) == len(DELTAS) + 1
    assert ace.load_points([path], "first") == []
    sub = ace.load_subcritical_points([path])
    assert len(sub) == len(DELTAS) and all(p.side == "sub" for p in sub)


def test_resolve_a_star_and_selection(tmp_path):
    path = write_history(tmp_path / "h.json", [], bracket=(4.69, 4.71))
    assert ace.resolve_a_star([path], "lower", None) == 4.69
    assert ace.resolve_a_star([path], "mid", None) == pytest.approx(4.70)
    assert ace.resolve_a_star([path], "value", 4.7) == 4.7
    points = super_points([0.001, 0.01, 0.1])
    assert [p.amplitude for p in ace.select_points(points, A_STAR, 0.005)] == [points[1].amplitude, points[2].amplitude]


def test_main_fits_both_sides_and_plots(tmp_path, capsys):
    path = write_history(tmp_path / "h.json", history_runs(), bracket=(A_STAR - 1e-5, A_STAR + 1e-5))
    output = tmp_path / "scaling.png"
    code = ace.main(["--history", str(path), "--a-star-mode", "value", "--a-star", str(A_STAR),
                     "--output", str(output)])
    assert code == 0
    assert output.stat().st_size > 0
    text = capsys.readouterr().out
    assert "super: gamma = 0.3700" in text
    assert "sub  : gamma = 0.3700" in text
    assert "gamma_super - gamma_sub = 0.0000" in text


def test_main_without_enough_points(tmp_path):
    path = write_history(tmp_path / "h.json", [{"amplitude": 5.0, "outcome": "collapse", "m_ah_settled": 0.3}])
    assert ace.main(["--history", str(path), "--output", str(tmp_path / "x.png")]) == 1
