.PHONY: clean submodules test_submodules test verify copy_remote copy_logs test_remote verify_remote status_remote status
REMOTE_DIR := /tmp/$(shell basename $(CURDIR))

# Delete the compiled Magma package files (*.sig), also those in the submodules. Skips
# scratch/, which may contain separate git worktrees with their own *.sig files.
clean:
	find . -path ./scratch -prune -o -type f -name "*.sig" -print -delete

# Initialise the given submodules (all if none are given) that are not checked out yet.
# Submodules that are already checked out are left alone, so this never moves them to a
# different commit.
init_missing_submodules = git submodule status $(1) | awk '/^-/ {print $$2}' | xargs -r git submodule update --init --

submodules:
	@$(call init_missing_submodules,)

# The tests only need magma-unittest, not the large ShimuraCurveALQuotients.
test_submodules:
	@$(call init_missing_submodules,tests/magma-unittest)

# Magma exits 0 even after a failed assert, so check the output as well. Any occurrence of
# "error" (case insensitive) counts as a failure: better a false alarm than a missed error.
test: test_submodules
	@out=$$(magma -n tests/test_all.m < /dev/null 2>&1); echo "$$out"; \
	if echo "$$out" | grep -qi "error"; then echo "TESTS FAILED"; exit 1; fi; \
	echo "$$out" | grep -q "All tests finished!"

verify: submodules
	bash verify_all.sh

copy_remote:
	rsync -avz --delete ./ $(ssh):$(REMOTE_DIR)

copy_logs:
	rsync -avz $(ssh):$(REMOTE_DIR)/logs/ logs/

# Usage: make test_remote ssh="user@hostname"
test_remote: copy_remote
	ssh $(ssh) "cd $(REMOTE_DIR) && make test"

# Usage: make status_remote ssh="user@hostname"
status_remote:
	ssh $(ssh) "cd $(REMOTE_DIR) && comm -23 \
	  <(find computations -name '*.m' | sort) \
	  <(grep -oP 'computations/\S+\.m' logs/magma_verify_joblog.txt | sort) \
	  > /tmp/.pending && cat /tmp/.pending && echo \"\$$(wc -l < /tmp/.pending) jobs pending\""

status: status_remote

# Usage: make verify_remote ssh="user@hostname"
verify_remote: copy_remote
	ssh -f $(ssh) "cd $(REMOTE_DIR) && nohup make verify > $(REMOTE_DIR)/logs/nohup.out 2>&1"
