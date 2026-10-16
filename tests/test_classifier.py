# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
"""Tests for the three-way outcome classifier."""

import pytest

import cactus_outputs as co
import outcome_classifier as oc
from fake_outputs import (
    collapsing_lapse,
    dispersing_curvature,
    growing_horizon,
    horizon_series,
    recovering_lapse,
    write_bh_diagnostics,
    write_log,
    write_scalar,
)

CONFIG = oc.ClassifierConfig()


def track(finds, number=1, origin_z=0.0):
    return oc.HorizonTrack(number, origin_z, horizon_series(finds))


def collapse_inputs(t_form=6.0, t_end=12.0, mass=0.6, **kwargs):
    return oc.RunInputs(
        horizons=[track(growing_horizon(t_form, t_end, mass))],
        curvature=[(t, 1.0 + t) for t, _ in collapsing_lapse(t_end)],
        min_lapse=collapsing_lapse(t_end),
        **kwargs
    )


def dispersion_inputs(t_end=12.0, **kwargs):
    return oc.RunInputs(
        curvature=dispersing_curvature(t_end),
        min_lapse=recovering_lapse(t_end),
        **kwargs
    )


# --- collapse -------------------------------------------------------------


def test_confirmed_horizon_is_collapse():
    verdict = oc.classify(collapse_inputs(), CONFIG)
    assert verdict.outcome == oc.COLLAPSE
    assert verdict.horizon.t_confirmed == pytest.approx(6.2)
    assert verdict.horizon.kind == "origin"


def test_mass_settles_after_the_horizon_relaxes():
    # The settle window is 10 M = 6 back from the last find.
    verdict = oc.classify(collapse_inputs(t_end=20.0), CONFIG)
    assert verdict.horizon.m_settled == pytest.approx(0.6, rel=1e-3)
    early = oc.classify(collapse_inputs(t_end=6.5), CONFIG)
    assert early.outcome == oc.COLLAPSE
    assert early.horizon.m_settled is None


def test_two_finds_are_not_enough():
    finds = growing_horizon(6.0, 6.1, 0.6)
    inputs = oc.RunInputs(horizons=[track(finds)], curvature=dispersing_curvature(12.0),
                          min_lapse=recovering_lapse(12.0))
    verdict = oc.classify(inputs, CONFIG)
    assert verdict.outcome == oc.UNDECIDED
    assert any("without confirmation" in r for r in verdict.reasons)


def test_finds_with_gaps_are_not_consecutive():
    finds = growing_horizon(6.0, 7.0, 0.6)[::2]  # every other find time
    other = growing_horizon(1.0, 1.1, 0.01, z=0.5)  # sets the find interval to 16
    inputs = oc.RunInputs(horizons=[track(finds), track(other, number=2, origin_z=0.5)])
    assert oc.classify(inputs, CONFIG).outcome == oc.UNDECIDED


def test_shrinking_mass_breaks_the_streak():
    finds = [(160 * 6 + 16 * k, 6.0 + 0.1 * k, 0.6 * (1 - 0.01 * k), 0.3, 0.0) for k in range(6)]
    verdict = oc.classify(oc.RunInputs(horizons=[track(finds)]), CONFIG)
    assert verdict.outcome == oc.UNDECIDED


def test_minimum_radius_is_enforced_when_set():
    config = oc.ClassifierConfig(min_horizon_radius=0.5)
    assert oc.classify(collapse_inputs(), config).outcome == oc.UNDECIDED
    assert oc.classify(collapse_inputs(), CONFIG).outcome == oc.COLLAPSE


def test_earliest_confirmed_finder_wins():
    origin = track(growing_horizon(8.0, 12.0, 0.6), number=1, origin_z=0.0)
    axis = track(growing_horizon(6.0, 12.0, 0.2, z=0.7), number=2, origin_z=0.7)
    verdict = oc.classify(oc.RunInputs(horizons=[origin, axis]), CONFIG)
    assert verdict.outcome == oc.COLLAPSE
    assert (verdict.horizon.number, verdict.horizon.kind) == (2, "axis")
    assert len(verdict.horizons) == 2


def test_collapse_survives_a_later_crash():
    inputs = collapse_inputs(log=co.LogFailures(["ML_BSSN::alpha"], False, False, 1900, 11.9), exit_code=0)
    verdict = oc.classify(inputs, CONFIG)
    assert verdict.outcome == oc.COLLAPSE
    assert any("crashed afterwards" in r for r in verdict.reasons)


# --- dispersion -----------------------------------------------------------


def test_dispersion_needs_and_finds_positive_evidence():
    verdict = oc.classify(dispersion_inputs(), CONFIG)
    assert verdict.outcome == oc.DISPERSION
    assert verdict.curvature.t_dispersed >= CONFIG.dispersion_min_time
    assert verdict.curvature.k_max == pytest.approx(3.0)


def test_no_curvature_output_is_undecided():
    inputs = oc.RunInputs(min_lapse=recovering_lapse(12.0))
    verdict = oc.classify(inputs, CONFIG)
    assert verdict.outcome == oc.UNDECIDED
    assert any("no curvature output" in r for r in verdict.reasons)


def test_too_early_to_tell_is_undecided():
    verdict = oc.classify(dispersion_inputs(t_end=4.0), CONFIG)
    assert verdict.outcome == oc.UNDECIDED


def test_regrowing_curvature_is_undecided():
    curvature = dispersing_curvature(12.0)
    curvature = [(t, k if t < 10.0 else 2.0) for t, k in curvature]
    inputs = oc.RunInputs(curvature=curvature, min_lapse=recovering_lapse(12.0))
    verdict = oc.classify(inputs, CONFIG)
    assert verdict.outcome == oc.UNDECIDED
    assert any("grew back" in r for r in verdict.reasons)


def test_unrecovered_lapse_is_undecided():
    inputs = oc.RunInputs(curvature=dispersing_curvature(12.0), min_lapse=collapsing_lapse(12.0))
    verdict = oc.classify(inputs, CONFIG)
    assert verdict.outcome == oc.UNDECIDED
    assert any("minimum lapse" in r for r in verdict.reasons)


def test_horizon_after_dispersion_is_a_contradiction():
    inputs = dispersion_inputs(t_end=14.0)
    inputs.horizons = [track(growing_horizon(11.0, 14.0, 0.05))]
    verdict = oc.classify(inputs, CONFIG)
    assert verdict.outcome == oc.UNDECIDED
    assert any("after the dispersion criteria" in r for r in verdict.reasons)


# --- other outcomes ------------------------------------------------------


def test_memory_abort_takes_precedence():
    assert oc.classify(collapse_inputs(guard_tripped=True), CONFIG).outcome == oc.OOM
    assert oc.classify(dispersion_inputs(oom_killed=True), CONFIG).outcome == oc.OOM


def test_nan_without_horizon_is_a_crash():
    inputs = oc.RunInputs(log=co.LogFailures(["ML_BSSN::Xt1"], False, False, 784, 4.9))
    verdict = oc.classify(inputs, CONFIG)
    assert verdict.outcome == oc.CRASH
    assert "ML_BSSN::Xt1" in verdict.reasons[0]


def test_failed_exit_is_a_crash():
    assert oc.classify(oc.RunInputs(exit_code=134), CONFIG).outcome == oc.CRASH


def test_timeout_without_evidence_and_with_evidence():
    assert oc.classify(oc.RunInputs(timed_out=True, exit_code=-15), CONFIG).outcome == oc.TIMEOUT
    assert oc.classify(dispersion_inputs(timed_out=True, exit_code=-15), CONFIG).outcome == oc.DISPERSION


# --- early termination ---------------------------------------------------


def test_early_stop_reasons():
    assert "settled" in oc.early_stop_reason(collapse_inputs(t_end=20.0), CONFIG)
    assert oc.early_stop_reason(collapse_inputs(t_end=6.5), CONFIG) is None
    assert oc.early_stop_reason(dispersion_inputs(), CONFIG) == "dispersion established"
    assert oc.early_stop_reason(dispersion_inputs(t_end=4.0), CONFIG) is None


# --- reading a run directory --------------------------------------------


def test_finder_origins_from_parfile():
    text = (
        "AHFinderDirect::N_horizons = 3\n"
        "AHFinderDirect::origin_z[2] = 0.6125\n"
        "AHFinderDirect::origin_z[3]                            = 1.2\n"
    )
    assert oc.finder_origins(text) == {1: 0.0, 2: 0.6125, 3: 1.2}


def test_load_run_inputs(tmp_path):
    output = tmp_path / "output"
    write_bh_diagnostics(output, 2, growing_horizon(6.0, 8.0, 0.3, z=0.6))
    write_scalar(output, co.CURVATURE_NORM_INF_FILENAME, [(t, -k) for t, k in dispersing_curvature(5.0)])
    write_scalar(output, co.MIN_LAPSE_FILENAME, recovering_lapse(5.0))
    write_log(tmp_path, 5.0)
    inputs = oc.load_run_inputs(tmp_path, "AHFinderDirect::N_horizons = 2\nAHFinderDirect::origin_z[2] = 0.6\n")
    assert [(h.number, h.kind) for h in inputs.horizons] == [(2, "axis")]
    assert inputs.curvature[0][1] == pytest.approx(1.0)  # absolute value
    assert inputs.min_lapse[0][0] == 0.0
    assert inputs.log.last_time == pytest.approx(5.0)
