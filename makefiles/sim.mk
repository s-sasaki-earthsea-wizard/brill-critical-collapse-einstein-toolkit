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

# Extra runner options, e.g.
#   make bisection EXTRA="--verify-endpoints"
#   make bisection-resume EXTRA="--set-outcome 3=collapse"
#   make bisection-resume EXTRA="--retry"
#   make bisection EXTRA="--classifier-config classifier.json"
EXTRA ?=

# Campaign name: one directory per campaign under SIM_OUTPUT_DIR.
CAMPAIGN ?= bisection
RUNS_DIR := /home/etuser/simulations/$(CAMPAIGN)

RUNNER := python3 scripts/bisection_runner.py
RUNNER_ARGS = \
	--a-min $(A_MIN) --a-max $(A_MAX) --max-iter $(MAX_ITER) \
	--t-final $(T_FINAL) \
	--np $(SIM_MPI_PROCS) --omp $(SIM_OMP_THREADS) \
	--runs-dir $(RUNS_DIR) $(EXTRA)

# Refuse to start while another cactus_sim is running in the container:
# two concurrent runs halve the speed of both and break the ETA, and the
# GW150914 project lost a run to exactly this. Zombies (state Z) are
# ignored: they are dead ranks of an aborted run that nobody reaped.
.PHONY: _precheck-idle
_precheck-idle:
	@running=$$($(COMPOSE) exec -T et ps -eo stat=,comm= 2>/dev/null \
		| awk '$$2 == "cactus_sim" && $$1 !~ /^Z/' | wc -l | tr -d '\r '); \
	if [ "$$running" != "0" ]; then \
		echo "[error] $$running cactus_sim process(es) already running in $(CONTAINER_NAME)."; \
		echo "        Refusing to start a second run. Check with 'make status'."; \
		exit 1; \
	fi

.PHONY: bisection-dry-run
bisection-dry-run: ## Render the first bisection parfile and print the run plan (no Cactus)
	$(COMPOSE) exec -T et $(RUNNER) $(RUNNER_ARGS) --dry-run

.PHONY: single-run
single-run: _precheck-idle ## Run and classify one amplitude (AMP=...) for smoke tests
	@test -n "$(AMP)" || (echo "Usage: make single-run AMP=<amplitude>" && exit 1)
	$(COMPOSE) exec -T et $(RUNNER) $(RUNNER_ARGS) --single $(AMP)

# Initial data tables: psi on the compactified polar grid, read by the
# BrillPsiTable thorn. Written to SIM_OUTPUT_DIR/tables.
ID_N_S     ?= 512
ID_N_THETA ?= 64
ID_SCALE   ?= 2.0
TABLES_DIR := /home/etuser/simulations/tables

.PHONY: id-table
id-table: ## Solve the Brill initial data for A=... and write its psi table (ID_N_S, ID_N_THETA, ID_SCALE)
	@test -n "$(A)" || (echo "Usage: make id-table A=<amplitude>" && exit 1)
	$(COMPOSE) exec -T et python3 scripts/brill_initial_data.py -A $(A) \
		--n-s $(ID_N_S) --n-theta $(ID_N_THETA) --scale $(ID_SCALE) \
		--out $(TABLES_DIR)/holz_A$(A)_k2_$(ID_N_S)x$(ID_N_THETA).h5

.PHONY: run-par
run-par: _precheck-idle ## Run a parfile as is under the monitor (PAR=..., RUN_NAME=..., CAMPAIGN=...)
	@test -n "$(PAR)" || (echo "Usage: make run-par PAR=par/<file>.par [RUN_NAME=<name>] [CAMPAIGN=<dir>]" && exit 1)
	$(COMPOSE) exec -T et $(RUNNER) $(RUNNER_ARGS) --parfile $(PAR) $(if $(RUN_NAME),--run-name $(RUN_NAME))

.PHONY: reclassify
reclassify: ## Classify finished runs again and rewrite run_summary.json (RUN="<container dir> ...", EXTRA)
	@test -n "$(RUN)" || (echo "Usage: make reclassify RUN=/home/etuser/simulations/<campaign>/single/iter00_A4.8" && exit 1)
	$(COMPOSE) exec -T et $(RUNNER) --reclassify $(RUN) $(EXTRA)

.PHONY: bisection
bisection: _precheck-idle ## Run the bisection campaign (A_MIN, A_MAX, MAX_ITER, CAMPAIGN, EXTRA)
	$(COMPOSE) exec -T et $(RUNNER) $(RUNNER_ARGS)

.PHONY: bisection-resume
bisection-resume: _precheck-idle ## Resume a stopped campaign (EXTRA="--set-outcome ITER=collapse|dispersion" or "--retry")
	$(COMPOSE) exec -T et $(RUNNER) $(RUNNER_ARGS) --resume
