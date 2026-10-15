# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Syota Sasaki
# Check and test targets
# =======================================================
# Everything runs inside the container on purpose: the base image ships
# Python 3.8, so a syntax check on the host's newer Python would not
# catch 3.9+ constructs that break at run time.

##@ Checks

.PHONY: check
check: ## Byte-compile the scripts and tests with the container's Python 3.8
	$(COMPOSE) exec -T et python3 -m py_compile scripts/*.py tests/*.py
	@echo "py_compile OK"

.PHONY: test
test: ## Run the pytest suite inside the container
	$(COMPOSE) exec -T et python3 -m pytest -q
