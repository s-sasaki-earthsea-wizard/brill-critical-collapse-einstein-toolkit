# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
"""Tests for the bisection campaign: bracket, history, resume, end to end."""

import json
import re

import pytest

import bisection_runner as br
import cactus_outputs as co
import outcome_classifier as oc
from fake_outputs import (
    collapsing_lapse,
    dispersing_curvature,
    growing_horizon,
    recovering_lapse,
    write_bh_diagnostics,
    write_log,
    write_scalar,
)

A_STAR = 4.7
TEMPLATE = (
    'Cactus::cctk_run_title = "test, A = __AMPLITUDE__"\n'
    "Cactus::cctk_final_time = __FINAL_TIME__\n"
    "AHFinderDirect::N_horizons = 1\n"
)


def run(amplitude, outcome, role="bisection", iteration=0):
    return br.RunResult(iteration=iteration, amplitude=amplitude, run_dir="x", outcome=outcome, role=role)


# --- bracket ---------------------------------------------------------------


def test_update_bracket():
    assert br.update_bracket(4.0, 6.0, run(5.0, oc.COLLAPSE)) == (4.0, 5.0)
    assert br.update_bracket(4.0, 6.0, run(5.0, oc.DISPERSION)) == (5.0, 6.0)
    with pytest.raises(RuntimeError):
        br.update_bracket(4.0, 6.0, run(5.0, oc.UNDECIDED))
    with pytest.raises(br.BracketError):
        br.update_bracket(4.5, 6.0, run(4.25, oc.COLLAPSE))
    with pytest.raises(br.BracketError):
        br.update_bracket(4.0, 5.0, run(5.5, oc.DISPERSION))


def test_replay_bracket_with_endpoints():
    runs = [
        run(4.0, oc.DISPERSION, role="endpoint"),
        run(6.0, oc.COLLAPSE, role="endpoint"),
        run(5.0, oc.COLLAPSE),
        run(4.5, oc.DISPERSION),
        run(4.75, oc.UNDECIDED),
    ]
    assert br.replay_bracket(4.0, 6.0, runs) == (4.5, 5.0)
    with pytest.raises(br.BracketError):
        br.replay_bracket(4.0, 6.0, [run(6.0, oc.DISPERSION, role="endpoint")])


def test_monotonicity_violation():
    assert br.monotonicity_violation([run(4.5, oc.DISPERSION), run(5.0, oc.COLLAPSE)]) is None
    assert "4.8" in br.monotonicity_violation([run(4.8, oc.DISPERSION), run(4.6, oc.COLLAPSE)])


# --- history -------------------------------------------------------------


def make_args(tmp_path, *extra):
    template = tmp_path / "template.par"
    if not template.exists():
        template.write_text(TEMPLATE)
    return br.parse_args([
        "--runs-dir", str(tmp_path / "campaign"),
        "--template", str(template),
        "--cactus-bin", str(tmp_path / "no-such-cactus"),
        "--a-min", "4.0", "--a-max", "6.0", "--max-iter", "6", "--t-final", "20.0",
    ] + list(extra))


def test_history_round_trip(tmp_path):
    args = make_args(tmp_path)
    history = br.new_history(args, {"template_sha256": "abc"})
    history["runs"].append(json.loads(json.dumps(br.asdict(run(5.0, oc.COLLAPSE)))))
    br.save_history(args.runs_dir, history)
    assert not list(args.runs_dir.glob("*.tmp"))
    loaded = br.load_history(args.runs_dir)
    assert br.history_runs(loaded)[0].outcome == oc.COLLAPSE
    loaded["schema_version"] = 99
    br.save_history(args.runs_dir, loaded)
    with pytest.raises(ValueError):
        br.load_history(args.runs_dir)


def test_run_result_ignores_unknown_keys():
    data = br.asdict(run(5.0, oc.COLLAPSE))
    data["from_the_future"] = 1
    assert br.RunResult.from_dict(data).amplitude == 5.0


def test_apply_overrides():
    history = {"runs": [br.asdict(run(4.75, oc.UNDECIDED, iteration=2))]}
    br.apply_overrides(history, ["2=collapse"])
    entry = history["runs"][0]
    assert (entry["outcome"], entry["human_override"]) == (oc.COLLAPSE, oc.UNDECIDED)
    with pytest.raises(ValueError):
        br.apply_overrides(history, ["2=crash"])
    with pytest.raises(ValueError):
        br.apply_overrides(history, ["9=collapse"])


def test_fingerprint_tracks_the_numerical_setup():
    config = oc.ClassifierConfig()
    settings = br.InitialDataSettings()
    base = br.campaign_fingerprint(TEMPLATE, 20.0, config, "/no/such", settings)
    assert base["initial_data"] is None
    assert br.fingerprint_differences(base, br.campaign_fingerprint(TEMPLATE, 20.0, config, "/no/such", settings)) == []
    changed = br.campaign_fingerprint(TEMPLATE + "# edit\n", 20.0, oc.ClassifierConfig(curvature_drop_factor=5.0),
                                      "/no/such", settings)
    assert br.fingerprint_differences(base, changed) == ["classifier", "template_sha256"]


def test_classifier_config_from_json(tmp_path):
    path = tmp_path / "classifier.json"
    path.write_text('{"min_consecutive_finds": 4, "lapse_recovered": 0.6}')
    config = br.load_classifier_config(path)
    assert (config.min_consecutive_finds, config.lapse_recovered) == (4, 0.6)
    path.write_text('{"no_such_threshold": 1}')
    with pytest.raises(ValueError):
        br.load_classifier_config(path)


# --- end to end with a fake Cactus --------------------------------------


def physical_outcome(amplitude):
    return oc.COLLAPSE if amplitude > A_STAR else oc.DISPERSION


def fake_executor(behaviour, calls):
    """execute_run stand-in that writes the output of the given behaviour."""

    def execute(command, run_dir, run_name, status_dir, t_final, timeout_s, mem_guard_fraction,
                monitor_interval_s, cgroup_root=None, probe=None, probe_interval_s=60.0, campaign=None):
        text = (run_dir / br.PARFILE_NAME).read_text()
        amplitude = float(re.search(r"A = ([-+0-9.eE]+)", text).group(1))
        kind = behaviour(amplitude)
        output = run_dir / "output"
        output.mkdir(parents=True, exist_ok=True)
        exit_code, t_end, extra = 0, 20.0, ""
        if kind == oc.COLLAPSE:
            mass = 0.5 * (amplitude - A_STAR) ** 0.37
            write_bh_diagnostics(output, 1, growing_horizon(5.0, t_end, mass))
            write_scalar(output, co.CURVATURE_NORM_INF_FILENAME, [(t, 1.0 + t) for t, _ in collapsing_lapse(t_end)])
            write_scalar(output, co.MIN_LAPSE_FILENAME, collapsing_lapse(t_end))
        elif kind == oc.DISPERSION:
            t_end = 12.0
            peak = 1.0 / (A_STAR - amplitude)
            write_scalar(output, co.CURVATURE_NORM_INF_FILENAME, dispersing_curvature(t_end, peak=peak))
            write_scalar(output, co.MIN_LAPSE_FILENAME, recovering_lapse(t_end))
        elif kind == oc.CRASH:
            t_end, exit_code = 4.9, 0
            extra = "  -> There were 12 NaN/Inf value(s) found in variable 'ML_BSSN::alpha'\n"
        else:
            t_end = 3.0
        write_log(run_dir, t_end, extra)
        decided = probe() if probe is not None else None
        calls.append((amplitude, kind, decided, campaign))
        return br.ExecutionResult(
            exit_code=exit_code, wall_seconds=60.0, peak_memory_bytes=1 << 30, guard_tripped=False,
            oom_killed=False, timed_out=False, last_iteration=int(t_end * 160), t_reached=t_end,
            rate_per_hour=10.0, stop_reason=(br.EARLY_STOP_PREFIX + decided) if decided else None,
        )

    return execute


def test_campaign_converges_on_the_threshold(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(br, "execute_run", fake_executor(physical_outcome, calls))
    args = make_args(tmp_path)
    assert br.run_campaign(args) == 0

    history = br.load_history(args.runs_dir)
    assert history["status"] == "completed"
    bracket = history["bracket"]
    assert bracket["a_min"] < A_STAR < bracket["a_max"]
    assert bracket["a_max"] - bracket["a_min"] == pytest.approx(2.0 / 2 ** 6)
    runs = br.history_runs(history)
    assert len(runs) == 6 and all(r.outcome in oc.DECISIVE for r in runs)
    assert [r.amplitude for r in runs][:3] == [5.0, 4.5, 4.75]
    collapse = runs[0]
    assert collapse.m_ah_settled == pytest.approx(0.5 * 0.3 ** 0.37, rel=1e-3)
    assert collapse.horizon_kind == "origin"
    dispersion = runs[1]
    assert dispersion.curvature_max == pytest.approx(1.0 / 0.2)
    # The early-stop probe saw the decisive evidence, and the status had campaign context.
    assert all(decided for _, _, decided, _ in calls)
    assert calls[0][3].remaining_runs == 5
    summary = json.loads((tmp_path / "campaign" / "iter00_A5.0" / br.RUN_SUMMARY_NAME).read_text())
    assert summary["outcome"] == oc.COLLAPSE


def test_campaign_stops_for_a_human_and_resumes(tmp_path, monkeypatch):
    calls = []

    def behaviour(amplitude):
        return oc.UNDECIDED if abs(amplitude - A_STAR) < 0.1 else physical_outcome(amplitude)

    monkeypatch.setattr(br, "execute_run", fake_executor(behaviour, calls))
    args = make_args(tmp_path)
    assert br.run_campaign(args) == 2
    history = br.load_history(args.runs_dir)
    assert history["status"] == "stopped"
    assert "A=4.75" in history["stop_reason"]

    # Resuming without a decision stops again without running anything.
    assert br.run_campaign(make_args(tmp_path, "--resume")) == 2
    assert len(calls) == 3

    # A human decides run 2 collapsed; the bisection continues in [4.5, 4.75]
    # and stops again at 4.625, which is also within 0.1 of A*.
    assert br.run_campaign(make_args(tmp_path, "--resume", "--set-outcome", "2=collapse")) == 2
    runs = br.history_runs(br.load_history(args.runs_dir))
    assert (runs[2].outcome, runs[2].human_override) == (oc.COLLAPSE, oc.UNDECIDED)
    assert (runs[3].amplitude, runs[3].outcome) == (pytest.approx(4.625), oc.UNDECIDED)


def test_resume_refuses_a_changed_setup(tmp_path, monkeypatch):
    monkeypatch.setattr(br, "execute_run", fake_executor(lambda a: oc.UNDECIDED, []))
    args = make_args(tmp_path)
    assert br.run_campaign(args) == 2
    (tmp_path / "template.par").write_text(TEMPLATE + "ML_BSSN::harmonicF = 1.0\n")
    assert br.run_campaign(make_args(tmp_path, "--resume", "--retry")) == 1


def test_existing_campaign_needs_resume(tmp_path, monkeypatch):
    monkeypatch.setattr(br, "execute_run", fake_executor(physical_outcome, []))
    assert br.run_campaign(make_args(tmp_path)) == 0
    assert br.run_campaign(make_args(tmp_path)) == 1


def test_endpoints_are_verified_first(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(br, "execute_run", fake_executor(physical_outcome, calls))
    args = make_args(tmp_path, "--verify-endpoints")
    assert br.run_campaign(args) == 0
    runs = br.history_runs(br.load_history(args.runs_dir))
    assert [(r.role, r.amplitude) for r in runs[:3]] == [("endpoint", 4.0), ("endpoint", 6.0), ("bisection", 5.0)]
    assert sum(r.role == "bisection" for r in runs) == 6


def test_wrong_endpoint_stops(tmp_path, monkeypatch):
    monkeypatch.setattr(br, "execute_run", fake_executor(lambda a: oc.DISPERSION, []))
    args = make_args(tmp_path, "--verify-endpoints")
    assert br.run_campaign(args) == 2
    assert "endpoint A=6.0" in br.load_history(args.runs_dir)["stop_reason"]


def test_crash_is_not_assigned_a_side(tmp_path, monkeypatch):
    monkeypatch.setattr(br, "execute_run", fake_executor(lambda a: oc.CRASH, []))
    args = make_args(tmp_path)
    assert br.run_campaign(args) == 2
    history = br.load_history(args.runs_dir)
    assert history["runs"][0]["outcome"] == oc.CRASH
    assert history["bracket"] == {"a_min": 4.0, "a_max": 6.0}


def test_single_run_is_classified(tmp_path, monkeypatch):
    monkeypatch.setattr(br, "execute_run", fake_executor(physical_outcome, []))
    assert br.run_single(make_args(tmp_path, "--single", "5.5")) == 0
    assert (tmp_path / "campaign" / "single" / "iter00_A5.5" / br.RUN_SUMMARY_NAME).exists()


def test_reclassify_rewrites_the_verdict_but_keeps_the_facts(tmp_path, monkeypatch):
    monkeypatch.setattr(br, "execute_run", fake_executor(physical_outcome, []))
    assert br.run_single(make_args(tmp_path, "--single", "5.5")) == 0
    summary_path = tmp_path / "campaign" / "single" / "iter00_A5.5" / br.RUN_SUMMARY_NAME
    summary = json.loads(summary_path.read_text())
    summary.update({"outcome": oc.UNDECIDED, "horizon_number": None, "m_ah_final": None})
    summary_path.write_text(json.dumps(summary))

    assert br.main(["--reclassify", str(summary_path.parent)]) == 0
    redone = json.loads(summary_path.read_text())
    assert redone["outcome"] == oc.COLLAPSE
    assert redone["horizon_number"] == 1
    assert redone["wall_seconds"] == 60.0
    assert redone["execution"] == summary["execution"]
    assert "reclassified" in redone["notes"] and "was undecided" in redone["notes"]

    strict = tmp_path / "strict.json"
    strict.write_text('{"min_consecutive_finds": 1000}')
    assert br.main(["--reclassify", str(summary_path.parent), "--classifier-config", str(strict)]) == 0
    assert json.loads(summary_path.read_text())["outcome"] == oc.UNDECIDED
    assert br.main(["--reclassify", str(tmp_path / "nowhere")]) == 1


def test_table_template_generates_initial_data(tmp_path):
    template = tmp_path / "table.par"
    template.write_text(TEMPLATE + 'BrillPsiTable::table_file = "__TABLE_FILE__"\n')
    settings = br.InitialDataSettings(n_s=32, n_theta=8, scale=2.0)
    parfile = br.prepare_run_dir(tmp_path / "runs", 0, 5.0, template, 10.0, settings)
    table = parfile.parent / br.TABLE_NAME
    assert table.exists()
    assert 'table_file = "{}"'.format(table) in parfile.read_text()
    fingerprint = br.campaign_fingerprint(template.read_text(), 10.0, oc.ClassifierConfig(), "/x", settings)
    assert fingerprint["initial_data"] == {"n_s": 32, "n_theta": 8, "scale": 2.0}
