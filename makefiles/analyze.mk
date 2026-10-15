# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
# Analysis targets
# =======================================================
# Fit the critical exponent gamma from the supercritical runs of a
# campaign. Runs inside the container (numpy / scipy / matplotlib come
# with the base image). Figures land in reports/, which is gitignored.

##@ Analysis

# How A* is chosen for Delta A = A - A* (see the script's --help).
A_STAR_MODE ?= lower

.PHONY: analyze
analyze: ## Fit gamma from M_AH vs (A - A*) and write reports/critical_scaling.png
	@mkdir -p reports
	$(COMPOSE) exec -T et python3 scripts/analyze_critical_exponent.py \
		--history $(RUNS_DIR)/bisection_history.json \
		--a-star-mode $(A_STAR_MODE) \
		--output reports/critical_scaling.png
