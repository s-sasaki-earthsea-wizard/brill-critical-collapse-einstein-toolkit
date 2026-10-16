# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
# Figure targets
# =======================================================
# Render the figures shown in the wiki into reports/figures/ (gitignored)
# and copy them into the wiki clone. The run figures read finished runs
# from SIM_OUTPUT_DIR, so they only work where those campaigns exist.

##@ Figures

FIG_DIR  := reports/figures
SIM_DIR  := /home/etuser/simulations
WIKI_DIR ?= ../brill-critical-collapse-einstein-toolkit.wiki
PLOT_ID   = $(COMPOSE) exec -T et python3 scripts/plot_initial_data.py
PLOT_RUNS = $(COMPOSE) exec -T et python3 scripts/plot_runs.py

.PHONY: figures figures-id figures-gauge wiki-figures

figures: figures-id figures-gauge ## Render every wiki figure into reports/figures/

figures-id: ## Render the initial-data figures (ADM mass vs A, psi at A = 5)
	@mkdir -p $(FIG_DIR)
	$(PLOT_ID) mass --out $(FIG_DIR)/initial_data_mass.png
	$(PLOT_ID) psi -A 5.0 --out $(FIG_DIR)/initial_data_psi_A5.png

figures-gauge: ## Render the Phase 2 gauge figures (needs the gauge_*, near*, diag_* campaigns)
	@mkdir -p $(FIG_DIR)
	$(PLOT_RUNS) horizon \
		--run "SA kappa=1, eta=1=$(SIM_DIR)/gauge_s2/single/iter00_A5.0" \
		--run "SA kappa=2/3, eta=1=$(SIM_DIR)/gauge_s2_k23/single/iter00_A5.0" \
		--run "SA kappa=1/2, eta=1=$(SIM_DIR)/gauge_s2_k12/single/iter00_A5.0" \
		--run "SA kappa=2/3, eta=0=$(SIM_DIR)/gauge_s2_k23_eta0/single/iter00_A5.0" \
		--run "SA kappa=1, zero shift=$(SIM_DIR)/gauge_s2z_wide/single/iter00_A5.0" \
		--band "0.56,0.58,Khirnov & Ledvinka 2019 (quasi-maximal slicing)" \
		--title "A = 5: apparent horizon mass with shock-avoiding (SA) slicing" \
		--out $(FIG_DIR)/gauge_A5_horizon_mass.png
	$(PLOT_RUNS) histories \
		--run "1+log, Gamma-driver=$(SIM_DIR)/gauge_s0/single/iter00_A5.0" \
		--run "harmonic, Gamma-driver=$(SIM_DIR)/gauge_s1/single/iter00_A5.0" \
		--run "SA kappa=1, eta=1=$(SIM_DIR)/gauge_s2/single/iter00_A5.0" \
		--run "SA kappa=2/3, eta=1=$(SIM_DIR)/gauge_s2_k23/single/iter00_A5.0" \
		--run "SA kappa=2/3, eta=0=$(SIM_DIR)/gauge_s2_k23_eta0/single/iter00_A5.0" \
		--run "SA kappa=1, zero shift=$(SIM_DIR)/gauge_s2z_wide/single/iter00_A5.0" \
		--title "A = 5: lapse and curvature under different slicings" \
		--out $(FIG_DIR)/gauge_A5_histories.png
	$(PLOT_RUNS) histories \
		--run "SA kappa=2/3, eta=1=$(SIM_DIR)/near_k23/single/iter00_A4.8" \
		--run "SA kappa=2/3, eta=1, 6 levels=$(SIM_DIR)/near_k23_l6/single/iter00_A4.8" \
		--run "SA kappa=1/2, eta=1=$(SIM_DIR)/near48_k12/single/iter00_A4.8" \
		--run "SA kappa=1/3, eta=1=$(SIM_DIR)/near48_k13/single/iter00_A4.8" \
		--run "SA kappa=2/3, eta=0=$(SIM_DIR)/near48_k23_eta0/single/iter00_A4.8" \
		--title "A = 4.8 (near the threshold): every variant crashes" \
		--out $(FIG_DIR)/near_critical_A4.8_histories.png
	$(PLOT_RUNS) profiles \
		--run "A=4, eta=0=$(SIM_DIR)/diag_eta0/single/iter00_A4.0" \
		--file admbase-shift.x.asc --var betax --times 4.2,4.4,4.6,4.8,5.0 --max-coordinate 1.2 \
		--ylabel "shift component beta^x" \
		--title "A = 4, SA kappa=2/3, eta = 0: a shock in the shift moves inwards" \
		--out $(FIG_DIR)/eta0_A4_shift_shock.png
	$(PLOT_RUNS) profiles \
		--run "A=4, eta=1=$(SIM_DIR)/diag_k23_noes/dense_A4.0" \
		--file weylscal4-curvir_group.z.asc --var curvIr --times 6.3,6.4,6.45,6.5,6.55 \
		--max-coordinate 0.5 --log --ylabel "|Re I|  (curvature)" \
		--title "A = 4, SA kappa=2/3, eta = 1: an instability grows from the origin" \
		--out $(FIG_DIR)/eta1_A4_origin_instability.png

wiki-figures: ## Copy reports/figures/*.png into the wiki clone's images/ (WIKI_DIR)
	@test -d $(WIKI_DIR)/images || { echo "no wiki clone at $(WIKI_DIR)"; exit 1; }
	cp $(FIG_DIR)/*.png $(WIKI_DIR)/images/
	@echo "copied; review and commit in $(WIKI_DIR)"
