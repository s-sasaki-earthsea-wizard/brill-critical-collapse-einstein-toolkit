# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
"""Tests for the resource monitor, the ETA estimator and monitored runs."""

import json
import sys

import pytest

import bisection_runner as br
import run_monitor as rm


def make_cgroup(root, current, limit="max", peak=None, oom_kill=0):
    root.mkdir(parents=True, exist_ok=True)
    (root / "memory.current").write_text("{}\n".format(current))
    (root / "memory.max").write_text("{}\n".format(limit))
    if peak is not None:
        (root / "memory.peak").write_text("{}\n".format(peak))
    (root / "memory.events").write_text(
        "low 0\nhigh 0\nmax 0\noom 0\noom_kill {}\n".format(oom_kill)
    )
    return root


def table_row(iteration, time):
    return "{:>9d} {:>9.3f} |   0.9876543   1.0000000 |\n".format(iteration, time)


# --- parsing -------------------------------------------------------------


def test_parse_info_line_reads_table_rows():
    assert rm.parse_info_line(table_row(4800, 30.0)) == (4800, 30.0)


@pytest.mark.parametrize(
    "line",
    [
        "Iteration      Time |  ADMBASE::alp\n",
        "-------------------------------------------\n",
        "INFO (Carpet): Carpet is running on 4 processes\n",
        "",
    ],
)
def test_parse_info_line_ignores_other_lines(line):
    assert rm.parse_info_line(line) is None


def test_parse_final_time():
    text = 'Cactus::terminate = "time"\nCactus::cctk_final_time = 12.5\n'
    assert br.parse_final_time(text) == 12.5
    assert br.parse_final_time("Cactus::cctk_itlast = 0\n") is None


# --- ETA -----------------------------------------------------------------


def test_eta_needs_two_rows():
    eta = rm.EtaEstimator(t_final=10.0)
    assert eta.observe(table_row(0, 0.0), 1000.0)
    assert eta.rate_per_hour is None
    assert eta.eta_seconds is None


def test_eta_from_constant_rate():
    eta = rm.EtaEstimator(t_final=10.0)
    for k in range(3):
        eta.observe(table_row(8 * k, 1.0 * k), 1000.0 + 100.0 * k)
    assert eta.rate_per_hour == pytest.approx(36.0)
    assert eta.eta_seconds == pytest.approx(800.0)
    assert (eta.iteration, eta.sim_time) == (16, 2.0)


def test_eta_window_tracks_a_slowdown():
    eta = rm.EtaEstimator(t_final=100.0, window_s=250.0)
    for k in range(5):  # fast: 1 per 100 s
        eta.observe(table_row(k, float(k)), 100.0 * k)
    for k in range(1, 6):  # slow: 1 per 400 s
        eta.observe(table_row(4 + k, 4.0 + k), 400.0 + 400.0 * k)
    assert eta.rate_per_hour == pytest.approx(9.0)


def test_campaign_eta():
    assert rm.campaign_eta_seconds([], 3, None) is None
    assert rm.campaign_eta_seconds([100.0, 300.0, 200.0], 2, 50.0) == pytest.approx(450.0)
    assert rm.campaign_eta_seconds([100.0], 0, 50.0) == pytest.approx(50.0)
    assert rm.campaign_eta_seconds([100.0], 0, None) is None


def test_format_duration():
    assert rm.format_duration(None) == "-"
    assert rm.format_duration(3725) == "01:02:05"
    assert rm.format_duration(90000) == "1d 01:00:00"


# --- memory --------------------------------------------------------------


def test_read_memory_sample(tmp_path):
    root = make_cgroup(tmp_path / "cg", current=123, limit=1000, peak=456, oom_kill=2)
    sample = rm.read_memory_sample(root)
    assert sample.current_bytes == 123
    assert sample.limit_bytes == 1000
    assert sample.peak_bytes == 456
    assert sample.oom_kill_count == 2


def test_read_memory_sample_unlimited(tmp_path):
    sample = rm.read_memory_sample(make_cgroup(tmp_path / "cg", current=1))
    assert sample.limit_bytes is None
    assert sample.peak_bytes is None


def test_read_memory_sample_requires_cgroup_v2(tmp_path):
    with pytest.raises(FileNotFoundError):
        rm.read_memory_sample(tmp_path)


def test_monitor_guard_and_oom_detection(tmp_path):
    root = make_cgroup(tmp_path / "cg", current=500, limit=1000)
    monitor = rm.ResourceMonitor(tmp_path / "res.csv", 3600.0, 0.9, root)
    monitor.start()
    assert not monitor.guard_tripped
    assert not monitor.oom_killed

    make_cgroup(root, current=950, limit=1000, oom_kill=1)
    monitor.sample_once()
    assert monitor.guard_tripped
    assert monitor.oom_killed
    assert monitor.peak_bytes == 950
    monitor.stop()

    rows = (tmp_path / "res.csv").read_text().splitlines()
    assert rows[0].split(",") == list(rm.RESOURCES_COLUMNS)
    assert len(rows) == 4  # header, start, explicit sample, stop


def test_preexisting_oom_kills_do_not_count(tmp_path):
    root = make_cgroup(tmp_path / "cg", current=1, limit=1000, oom_kill=5)
    monitor = rm.ResourceMonitor(tmp_path / "res.csv", 3600.0, 0.9, root)
    monitor.start()
    monitor.stop()
    assert not monitor.oom_killed


# --- status file ---------------------------------------------------------


def test_write_and_print_status(tmp_path, capsys):
    rm.write_status(tmp_path, {"run_name": "demo", "state": "running", "eta_seconds": 61})
    assert json.loads((tmp_path / rm.STATUS_FILENAME).read_text())["run_name"] == "demo"
    assert not list(tmp_path.glob("*.tmp"))
    assert rm.print_status(tmp_path) == 0
    assert "00:01:01" in capsys.readouterr().out


def test_print_status_without_file(tmp_path):
    assert rm.print_status(tmp_path) == 1


# --- monitored execution -------------------------------------------------


def run_fake(tmp_path, script, timeout_s=60.0):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    root = make_cgroup(tmp_path / "cg", current=100, limit=1000)
    result = br.execute_run(
        [sys.executable, "-c", script],
        run_dir,
        "fake",
        tmp_path,
        t_final=10.0,
        timeout_s=timeout_s,
        mem_guard_fraction=0.9,
        monitor_interval_s=0.2,
        cgroup_root=root,
    )
    return run_dir, result


def test_execute_run_streams_log_and_status(tmp_path):
    script = (
        "import time\n"
        "for k in range(3):\n"
        "    print('{:>9d} {:>9.3f} |   1.0   1.0 |'.format(8 * k, 0.5 * k), flush=True)\n"
        "    time.sleep(0.05)\n"
    )
    run_dir, result = run_fake(tmp_path, script)
    assert result.exit_code == 0
    assert (result.last_iteration, result.t_reached) == (16, 1.0)
    assert not (result.guard_tripped or result.timed_out or result.oom_killed)
    assert result.peak_memory_bytes == 100
    assert "       16     1.000 |" in (run_dir / br.LOG_NAME).read_text()
    status = json.loads((tmp_path / rm.STATUS_FILENAME).read_text())
    assert status["state"] == "finished"
    assert status["outcome"] == "completed"


def test_execute_run_escalates_on_timeout(tmp_path, monkeypatch):
    monkeypatch.setattr(br, "TERMINATE_GRACE_S", 0.5)
    monkeypatch.setattr(br, "SIGTERM_GRACE_S", 0.5)
    run_dir, result = run_fake(tmp_path, "import time; time.sleep(60)", timeout_s=0.5)
    assert result.timed_out
    assert result.exit_code is not None and result.exit_code < 0
    assert (run_dir / "output" / "TERMINATE").read_text().strip() == "1"
    assert br.describe_execution(result) == "timeout"
