# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
# Brill Wave Vacuum Critical Collapse (Einstein Toolkit)
# ======================================================
# Run `make help` for the list of targets.
# Functionality is split across sub-makefiles under makefiles/:
#   docker.mk   image and container lifecycle
#   sim.mk      single runs and the bisection campaign
#   monitor.mk  resource monitoring and ETA
#   analyze.mk  critical exponent fit
#   test.mk     syntax checks and the pytest suite

.DEFAULT_GOAL := help

# docker.mk must come first: it loads .env, which the others rely on.
include makefiles/docker.mk
include makefiles/sim.mk
include makefiles/monitor.mk
include makefiles/analyze.mk
include makefiles/test.mk

.PHONY: help
help: ## Show this help
	@echo "Brill Wave Vacuum Critical Collapse -- available targets:"
	@echo ""
	@awk 'BEGIN{FS=":.*?## "} \
		/^##@/ { printf "\n\033[1m%s\033[0m\n", substr($$0, 5); next } \
		/^[a-zA-Z0-9_.-]+:.*?## / { printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2 }' \
		$(MAKEFILE_LIST)
	@echo ""
