# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
# Simulation targets
# =======================================================
# Drive scripts/bisection_runner.py inside the container. The runner owns
# parfile rendering, mpirun, memory monitoring, ETA reporting and outcome
# classification; these targets only pass settings through.
#
# Every setting can be overridden on the command line or in .env, e.g.
#   make bisection A_MIN=4.5 A_MAX=5.0 MAX_ITER=10 CAMPAIGN=holz-l4

##@ Simulation

# The recipes below rely on bash (PIPESTATUS, pgrep exit codes).
SHELL := /bin/bash

# Parallelism. Placeholder until the np/OMP/cpuset benchmark (see .env.example).
SIM_MPI_PROCS   ?= 4
SIM_OMP_THREADS ?= 4

# Bisection bracket and length.
A_MIN    ?= 4.0
A_MAX    ?= 6.0
MAX_ITER ?= 7

# Physical end time of every run.
T_FINAL ?= 30.0

# Campaign name: one directory per campaign under SIM_OUTPUT_DIR.
CAMPAIGN ?= bisection
RUNS_DIR := /home/etuser/simulations/$(CAMPAIGN)

RUNNER := python3 scripts/bisection_runner.py
RUNNER_ARGS = \
	--a-min $(A_MIN) --a-max $(A_MAX) --max-iter $(MAX_ITER) \
	--t-final $(T_FINAL) \
	--np $(SIM_MPI_PROCS) --omp $(SIM_OMP_THREADS) \
	--runs-dir $(RUNS_DIR)

# Refuse to start while another cactus_sim is running in the container:
# two concurrent runs halve the speed of both and break the ETA, and the
# GW150914 project lost a run to exactly this.
.PHONY: _precheck-idle
_precheck-idle:
	@running=$$($(COMPOSE) exec -T et pgrep -f cactus_sim 2>/dev/null | wc -l | tr -d '\r '); \
	if [ "$$running" != "0" ]; then \
		echo "[error] $$running cactus_sim process(es) already running in $(CONTAINER_NAME)."; \
		echo "        Refusing to start a second run. Check with 'make status'."; \
		exit 1; \
	fi

.PHONY: bisection-dry-run
bisection-dry-run: ## Render the first bisection parfile and print the run plan (no Cactus)
	$(COMPOSE) exec -T et $(RUNNER) $(RUNNER_ARGS) --dry-run

.PHONY: single-run
single-run: _precheck-idle ## Run one amplitude (AMP=...) for smoke tests and benchmarks
	@test -n "$(AMP)" || (echo "Usage: make single-run AMP=<amplitude>" && exit 1)
	$(COMPOSE) exec -T et $(RUNNER) $(RUNNER_ARGS) --single $(AMP)

.PHONY: bisection
bisection: _precheck-idle ## Run the bisection campaign (A_MIN, A_MAX, MAX_ITER, CAMPAIGN)
	$(COMPOSE) exec -T et $(RUNNER) $(RUNNER_ARGS)

.PHONY: bisection-resume
bisection-resume: _precheck-idle ## Resume an interrupted campaign from its bisection_history.json
	$(COMPOSE) exec -T et $(RUNNER) $(RUNNER_ARGS) --resume
