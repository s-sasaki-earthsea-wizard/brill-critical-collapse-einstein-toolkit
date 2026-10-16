#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
#
# Apply the project's patches (patches/*.patch) to a Cactus tree.
# Idempotent: a patch that is already applied is skipped, and a patch that
# neither applies nor reverses cleanly is an error.
#
# Usage: apply_cactus_patches.sh <patch-dir> <cactus-home>
set -euo pipefail

patch_dir=${1:?patch directory}
cactus_home=${2:?Cactus home}

shopt -s nullglob
for patch_file in "${patch_dir}"/*.patch; do
    name=$(basename "${patch_file}")
    if patch -d "${cactus_home}" -p1 -R --dry-run -s -f < "${patch_file}" > /dev/null 2>&1; then
        echo "  already applied: ${name}"
    elif patch -d "${cactus_home}" -p1 --dry-run -s -f < "${patch_file}" > /dev/null 2>&1; then
        patch -d "${cactus_home}" -p1 -s -f < "${patch_file}"
        echo "  applied: ${name}"
    else
        echo "  ERROR: ${name} neither applies nor is already applied" >&2
        exit 1
    fi
done
