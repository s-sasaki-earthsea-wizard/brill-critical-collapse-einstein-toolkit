# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
"""Tests for adding refinement levels to a template."""

import re
from pathlib import Path

import refine_template as rt

TEMPLATE = Path(__file__).resolve().parent.parent / "templates" / "brill_cartoon.par.template"


def value(text, name):
    return re.search(r"^{}\s*=\s*(.+?)\s*$".format(re.escape(name)), text, re.MULTILINE).group(1)


def test_two_extra_levels():
    source = TEMPLATE.read_text()
    refined = rt.refine(source, 2)
    assert value(refined, "Carpet::max_refinement_levels") == "6"
    assert value(refined, "CarpetRegrid2::num_levels_1") == "6"
    assert value(refined, "Carpet::time_refinement_factors") == '"[1, 2, 4, 8, 16, 32]"'
    assert float(value(refined, "CarpetRegrid2::radius_1[4]")) == 1.0
    assert float(value(refined, "CarpetRegrid2::radius_1[5]")) == 0.5
    for name in ("AHFinderDirect::find_every", "IOBasic::outInfo_every", "IOASCII::out1D_every",
                 "NaNChecker::check_every", "TerminationTrigger::check_file_every",
                 "WeylScal4::WeylScal4_invars_calc_Nth_calc_every"):
        assert int(value(refined, name)) == 4 * int(value(source, name)), name
    # Everything else is untouched, blank lines included.
    assert refined.count("\n\n") == source.count("\n\n")
    assert value(refined, "Time::dtfac") == value(source, "Time::dtfac")
