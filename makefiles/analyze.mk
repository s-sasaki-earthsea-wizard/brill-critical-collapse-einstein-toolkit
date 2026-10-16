# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
# Analysis targets
# =======================================================
# Fit the critical exponent gamma on both sides of the threshold:
# horizon masses of the collapse runs and peak curvature of the
# dispersion runs. Runs inside the container (numpy / scipy / matplotlib come
# with the base image). Figures land in reports/, which is gitignored.

##@ Analysis

# How A* is chosen for Delta A = A - A* (see the script's --help).
A_STAR_MODE ?= lower
# super, sub or both; horizon mass settled, final or first.
SIDE ?= both
MASS ?= settled

.PHONY: analyze
analyze: ## Fit gamma from M_AH and peak curvature vs |A - A*| (SIDE, MASS, A_STAR_MODE)
	@mkdir -p reports
	$(COMPOSE) exec -T et python3 scripts/analyze_critical_exponent.py \
		--history $(RUNS_DIR)/bisection_history.json \
		--a-star-mode $(A_STAR_MODE) --side $(SIDE) --mass $(MASS) \
		--output reports/critical_scaling.png
