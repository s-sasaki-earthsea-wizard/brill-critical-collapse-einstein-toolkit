# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
"""Tests for the parfile template and its rendering."""

import re

import pytest

import bisection_runner as br


@pytest.fixture
def template_text():
    return br.DEFAULT_TEMPLATE.read_text()


def test_template_declares_required_placeholders(template_text):
    assert "__AMPLITUDE__" in template_text
    assert "__FINAL_TIME__" in template_text


def test_render_leaves_no_placeholders(template_text):
    rendered = br.render_parfile(
        template_text, {"AMPLITUDE": "4.75", "FINAL_TIME": "30.0"}
    )
    assert br.PLACEHOLDER_PATTERN.search(rendered) is None
    assert re.search(
        r'^Cactus::cctk_run_title\s*=\s*"Brill wave collapse, A = 4\.75"$', rendered, re.M
    )
    assert re.search(r"^Cactus::cctk_final_time\s*=\s*30\.0$", rendered, re.M)


def test_render_rejects_missing_value(template_text):
    with pytest.raises(ValueError):
        br.render_parfile(template_text, {"AMPLITUDE": "4.75"})


def test_render_rejects_unknown_placeholder(template_text):
    with pytest.raises(KeyError):
        br.render_parfile(
            template_text,
            {"AMPLITUDE": "4.75", "FINAL_TIME": "30.0", "NOT_THERE": "1"},
        )


def test_ahfinderdirect_parameters_are_not_placeholders(template_text):
    # initial_guess__coord_sphere__radius must survive rendering untouched.
    assert "initial_guess__coord_sphere__radius" in template_text
    assert not any(
        "coord_sphere" in match for match in br.PLACEHOLDER_PATTERN.findall(template_text)
    )


def test_amplitude_round_trips_through_directory_name():
    assert br.run_dir_name(3, 4.6875) == "iter03_A4.6875"
    assert float(br.format_amplitude(4.696703125)) == 4.696703125


def test_mpirun_command_pins_openmp_threads(tmp_path):
    command = br.build_mpirun_command(
        "mpirun.mpich", "/bin/cactus_sim", tmp_path / "x.par", 4, 4
    )
    index = command.index("OMP_NUM_THREADS")
    assert command[index - 1] == "-genv"
    assert command[index + 1] == "4"
