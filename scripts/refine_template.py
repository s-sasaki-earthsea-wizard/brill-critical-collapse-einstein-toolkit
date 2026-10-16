#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
"""Add fine refinement levels to a Cartoon2D Brill template.

The templates use 4 levels (finest dx = 0.025) and count every cadence in
iterations of the finest level: 8 iterations per coarse step. Each extra
level halves the finest spacing and time step, so every cadence has to
double to keep the same physical interval. This script adds levels with
radii halving from the last one and rescales the cadences:

    out*_every, find_every, output_h_every, check_every,
    check_file_every, update_ah_every and WeylScal4's *_calc_every.

Usage: refine_template.py SOURCE DEST [--extra-levels 2]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import List, Optional

CADENCE_PATTERN = re.compile(
    r"^(\s*[A-Za-z0-9_]+::(?:[A-Za-z0-9_]*_every|check_file_every|update_ah_every)\s*=\s*)(\d+)[ \t]*$",
    re.MULTILINE,
)


def refine(text: str, extra_levels: int) -> str:
    """Return the template with ``extra_levels`` more refinement levels.

    Raises:
        ValueError: The template does not have the expected grid block.
    """
    factor = 2 ** extra_levels
    levels = int(re.search(r"Carpet::max_refinement_levels\s*=\s*(\d+)", text).group(1))
    new_levels = levels + extra_levels

    text = re.sub(r"(Carpet::max_refinement_levels\s*=\s*)\d+", r"\g<1>{}".format(new_levels), text)
    text = re.sub(r"(CarpetRegrid2::num_levels_1\s*=\s*)\d+", r"\g<1>{}".format(new_levels), text)
    factors = ", ".join(str(2 ** k) for k in range(new_levels))
    text, n = re.subn(r'(Carpet::time_refinement_factors\s*=\s*)"\[[^\]]*\]"',
                      r'\g<1>"[{}]"'.format(factors), text)
    if n != 1:
        raise ValueError("no Carpet::time_refinement_factors")

    radius_lines = re.findall(r"^CarpetRegrid2::radius_1\[(\d+)\]\s*=\s*([0-9.]+)\s*$", text, re.MULTILINE)
    last_index, last_radius = int(radius_lines[-1][0]), float(radius_lines[-1][1])
    if last_index != levels - 1:
        raise ValueError("radius_1 entries do not match max_refinement_levels")
    anchor = re.search(r"^CarpetRegrid2::radius_1\[{}\].*$".format(last_index), text, re.MULTILINE)
    extra = "".join(
        "\nCarpetRegrid2::radius_1[{}]  = {:5.2f}".format(last_index + k, last_radius / 2 ** k)
        for k in range(1, extra_levels + 1)
    )
    text = text[:anchor.end()] + extra + text[anchor.end():]

    def scale(match: "re.Match[str]") -> str:
        return "{}{}".format(match.group(1), int(match.group(2)) * factor)

    text = CADENCE_PATTERN.sub(scale, text)
    note = ("# Refined by scripts/refine_template.py: {} levels, finest dx = 0.2 / 2^{}, every\n"
            "# cadence scaled by {} to keep its physical interval.\n").format(new_levels, new_levels - 1, factor)
    return note + text


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("source", type=Path)
    parser.add_argument("dest", type=Path)
    parser.add_argument("--extra-levels", type=int, default=2)
    args = parser.parse_args(argv)
    args.dest.write_text(refine(args.source.read_text(), args.extra_levels))
    print("wrote {}".format(args.dest))
    return 0


if __name__ == "__main__":
    sys.exit(main())
