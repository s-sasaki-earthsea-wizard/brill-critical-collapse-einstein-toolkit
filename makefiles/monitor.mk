# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
# Monitoring targets
# =======================================================
# Resource monitoring and ETA are mandatory for every run in this
# project. The runner samples the container's cgroup memory accounting,
# enforces the memory guard and keeps <campaign>/status.json up to date
# (current amplitude, simulation time, run ETA, campaign ETA, memory).
# These targets read that state from the host side.

##@ Monitoring

.PHONY: status
status: ## Show the campaign status (current run, ETA, memory) from status.json
	$(COMPOSE) exec -T et python3 scripts/run_monitor.py status --runs-dir $(RUNS_DIR)

.PHONY: monitor
monitor: ## One-shot snapshot: docker stats for the container plus the campaign status
	@docker stats --no-stream $(CONTAINER_NAME)
	@echo ""
	@$(MAKE) --no-print-directory status
