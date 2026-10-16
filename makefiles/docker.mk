# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
# Docker targets
# =======================================================
# Build, start and manage the overlay image. The heavy lifting (building
# the Einstein Toolkit from source) happens in the GW230529 project; this
# image only adds runtime settings on top of it (see Dockerfile).

##@ Docker

# Load .env when present (BASE_IMAGE, SIM_OUTPUT_DIR, SIM_MPI_PROCS, ...).
ifneq (,$(wildcard .env))
include .env
export
endif

# Defaults used when .env does not set them.
BASE_IMAGE ?= gw230529-et:local
SIM_OUTPUT_DIR ?= $(HOME)/brill-collapse-output

COMPOSE := docker compose
CONTAINER_NAME := brill-collapse-et
CACTUS_HOME := /home/etuser/Cactus

.PHONY: docker-setup
docker-setup: ## First-time setup (create .env and the output directory)
	@test -f .env || (cp .env.example .env && echo "Created .env from .env.example")
	@mkdir -p $(SIM_OUTPUT_DIR)
	@echo "Base image       : $(BASE_IMAGE)"
	@echo "Output directory : $(SIM_OUTPUT_DIR)"

.PHONY: docker-build
docker-build: docker-setup ## Build the overlay image on top of BASE_IMAGE
	@docker image inspect $(BASE_IMAGE) >/dev/null 2>&1 \
		|| (echo "Base image $(BASE_IMAGE) not found."; \
		    echo "Build it first: cd ../gw230529-einstein-toolkit && make docker-build"; \
		    exit 1)
	$(COMPOSE) build

.PHONY: docker-up
docker-up: docker-setup ## Start the container in the background (with Jupyter Lab)
	$(COMPOSE) up -d

.PHONY: docker-down
docker-down: ## Stop and remove the container
	$(COMPOSE) down

.PHONY: docker-restart
docker-restart: docker-down docker-up ## Restart the container

.PHONY: docker-shell
docker-shell: ## Open a bash shell inside the container
	$(COMPOSE) exec et bash

.PHONY: docker-logs
docker-logs: ## Follow the container logs
	$(COMPOSE) logs -f et

.PHONY: docker-ps
docker-ps: ## Show container status
	$(COMPOSE) ps

.PHONY: thorn-build
thorn-build: ## Recompile the project thorns in the running container (dev loop; docker-build makes it permanent)
	$(COMPOSE) exec -T et bash -lc '\
		cp -a /home/etuser/work/thorns/. $(CACTUS_HOME)/arrangements/ && \
		cd $(CACTUS_HOME) && \
		make sim PROMPT=no -j$$(getconf _NPROCESSORS_ONLN) > /home/etuser/simulations/build_project_thorns.log 2>&1 \
		|| (tail -40 /home/etuser/simulations/build_project_thorns.log; exit 1)'
	@echo "Project thorns rebuilt in $(CONTAINER_NAME); run make docker-build to keep them."

.PHONY: docker-check
docker-check: ## Verify thorns, OpenMP pinning and MPI linkage inside the container
	@echo "== Thorns compiled into the binary (CST ThornList) =="
	@$(COMPOSE) exec -T et bash -lc '\
		compiled() { \
			if grep -qx "$$1" $(CACTUS_HOME)/configs/sim/ThornList 2>/dev/null; then \
				echo "  OK      $$1"; \
			else \
				echo "  MISSING $$1"; \
			fi; \
		}; \
		compiled McLachlan/ML_BSSN; \
		compiled McLachlan/ML_BSSN_Helper; \
		compiled McLachlan/ML_ADMConstraints; \
		compiled EinsteinInitialData/IDBrillData; \
		compiled CactusNumerical/Dissipation; \
		compiled CactusNumerical/ReflectionSymmetry; \
		compiled CactusBase/CartGrid3D; \
		compiled CactusBase/CoordBase; \
		compiled EinsteinAnalysis/AHFinderDirect; \
		compiled CactusUtils/NaNChecker; \
		compiled CactusNumerical/Cartoon2D; \
		compiled EinsteinInitialData/IDAnalyticBH; \
		compiled BrillCollapse/BrillPsiTable; \
		compiled CactusUtils/SystemStatistics; \
		compiled CactusUtils/TerminationTrigger; \
		compiled Carpet/CarpetRegrid2; \
		compiled Carpet/CarpetIOBasic; \
		compiled Carpet/CarpetIOScalar; \
		compiled Carpet/CarpetIOASCII'
	@echo ""
	@echo "== OpenMP pinning (must be 1) =="
	@$(COMPOSE) exec -T et bash -lc \
		'[ "$$OMP_NUM_THREADS" = "1" ] \
			&& echo "  OK      OMP_NUM_THREADS=1" \
			|| echo "  FAIL    OMP_NUM_THREADS=$${OMP_NUM_THREADS:-unset}"'
	@echo ""
	@echo "== MPI linkage of the binary (must be MPICH only) =="
	@$(COMPOSE) exec -T et bash -lc \
		'ldd $(CACTUS_HOME)/exe/cactus_sim | grep -q libmpich \
			&& echo "  OK      linked against libmpich" \
			|| echo "  FAIL    not linked against libmpich"; \
		 if ldd $(CACTUS_HOME)/exe/cactus_sim | grep -qE "libopen-pal|libopen-rte|libmpi\.so"; then \
			echo "  FAIL    Open MPI is also linked (mixed MPI stack)"; \
		 else \
			echo "  OK      no Open MPI contamination"; \
		 fi'
	@echo ""
	@echo "== MPI communicator test (mpirun.mpich -np 2) =="
	@$(COMPOSE) exec -T et bash -lc \
		'cd /home/etuser/simulations && rm -rf .mpi_check && mkdir -p .mpi_check && cd .mpi_check && \
		 n=$$(mpirun.mpich -np 2 $(CACTUS_HOME)/exe/cactus_sim \
		        /home/etuser/work/par/mpi_check.par 2>&1 \
		      | grep -m1 -oP "Carpet is running on \K[0-9]+"); \
		 if [ "$$n" = "2" ]; then \
			echo "  OK      2 ranks formed a single communicator"; \
		 else \
			echo "  FAIL    Carpet reported $${n:-unknown} processes (expected 2)"; \
		 fi'
