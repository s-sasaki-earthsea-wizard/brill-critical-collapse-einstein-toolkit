# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
# Brill wave critical collapse -- Einstein Toolkit run image
# =================================================================
# Thin overlay on top of the GW230529 project's Einstein Toolkit image
# (../gw230529-einstein-toolkit: ET_2025_05 "Kruskal", MPICH only, the
# full einsteintoolkit.th thornlist). Every thorn this project needs is
# already compiled into that image -- McLachlan/ML_BSSN,
# ML_ADMConstraints, IDBrillData, Dissipation, ReflectionSymmetry,
# AHFinderDirect, NaNChecker, SystemStatistics, TerminationTrigger,
# CarpetRegrid2 and the Carpet I/O thorns; `make docker-check` asserts it.
# Reusing it avoids a 60-120 min from-source build and ~17 GB of
# duplicated layers.
#
# Build the base image first, in the GW230529 repository:
#   cd ../gw230529-einstein-toolkit && make docker-build
#
# The project's own thorns (thorns/, e.g. BrillCollapse/BrillPsiTable,
# which loads the initial data solved by scripts/brill_initial_data.py)
# and its patches to stock thorns (patches/, e.g. WeylScal4 tensor types
# for Cartoon2D) are compiled into the base configuration incrementally:
# the base keeps its object files, so this step takes a few minutes.

ARG BASE_IMAGE=gw230529-et:local
FROM ${BASE_IMAGE}

# The base image predates the GW230529 fix that bakes these in (see
# finding 7 in that project's CLAUDE.md). Without OMP_NUM_THREADS every
# MPI rank starts one OpenMP thread per visible core, which oversubscribes
# the CPU and slows a run down by more than an order of magnitude without
# any error. The runner also passes it to every rank explicitly.
ENV OMP_NUM_THREADS=1
ENV PATH=/home/etuser/Cactus/exe:${PATH}

ARG PROJECT_THORNS="BrillCollapse/BrillPsiTable"
COPY --chown=etuser:etuser thorns/ /home/etuser/Cactus/arrangements/
COPY --chown=etuser:etuser patches/ /tmp/project-patches/
COPY --chown=etuser:etuser scripts/apply_cactus_patches.sh /tmp/apply_cactus_patches.sh
RUN bash /tmp/apply_cactus_patches.sh /tmp/project-patches /home/etuser/Cactus \
    && cd /home/etuser/Cactus \
    && for thorn in ${PROJECT_THORNS}; do \
         grep -qx "${thorn}" configs/sim/ThornList || echo "${thorn}" >> configs/sim/ThornList; \
       done \
    && make sim PROMPT=no -j"$(getconf _NPROCESSORS_ONLN)" > /tmp/build-project-thorns.log 2>&1 \
       || (tail -50 /tmp/build-project-thorns.log; exit 1)

WORKDIR /home/etuser/work

# The default command (Jupyter Lab) is inherited from the base image.
