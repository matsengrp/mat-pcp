# Remote training configuration
REMOTE_HOST ?= ermine
REMOTE_DIR ?= /fh/fast/matsen_e/shared/viral-dasm/mat-pcp

# Fetch only key result files from remote server using rsync
# Fetches: aamut_fitness/, clade_founder_nts.csv, expected_vs_actual_mut_counts.csv, spike_codon_sites_to_mask.txt
remote-fetch-key-files:
	@if [ -z "$(DIR)" ]; then \
		echo "Error: Please specify DIR='path/to/sync'"; \
		echo "Example: make remote-fetch-key-files DIR='results/2024-12-01'"; \
		exit 1; \
	fi
	@echo "Syncing key files from $(DIR) on $(REMOTE_HOST)..."
	rsync -av \
		--include='aamut_fitness/***' \
		--include='clade_founder_nts.csv' \
		--include='expected_vs_actual_mut_counts.csv' \
		--include='spike_codon_sites_to_mask.txt' \
		--exclude='*' \
		--no-links \
		-e "ssh -o RemoteCommand=none -o RequestTTY=no" \
		$(REMOTE_HOST):$(REMOTE_DIR)/$(DIR)/ $(DIR)/
	@echo "Key files synced to $(DIR)/"
