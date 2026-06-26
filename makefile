.PHONY: clean test verify copy_remote copy_logs test_remote verify_remote status_remote status
REMOTE_DIR := /tmp/$(shell basename $(CURDIR))

clean:
	find . -type f -name "*.sig" -delete

test:
	@echo "No tests configured yet - override this target"

verify:
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
