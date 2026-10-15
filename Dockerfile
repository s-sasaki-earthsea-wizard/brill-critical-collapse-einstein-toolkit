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
# TODO(initial-data): IDBrillData solves the Hamiltonian constraint with
# the PUGH-era EllBase/EllSOR solver, which is unigrid-only and cannot
# solve on a Carpet refinement hierarchy. Whatever replaces it (a small
# custom thorn, CT_MultiLevel coefficients, or an imported 2D solution)
# gets added here as an incremental Cactus rebuild on top of the base.

ARG BASE_IMAGE=gw230529-et:local
FROM ${BASE_IMAGE}

# The base image predates the GW230529 fix that bakes these in (see
# finding 7 in that project's CLAUDE.md). Without OMP_NUM_THREADS every
# MPI rank starts one OpenMP thread per visible core, which oversubscribes
# the CPU and slows a run down by more than an order of magnitude without
# any error. The runner also passes it to every rank explicitly.
ENV OMP_NUM_THREADS=1
ENV PATH=/home/etuser/Cactus/exe:${PATH}

WORKDIR /home/etuser/work

# The default command (Jupyter Lab) is inherited from the base image.
