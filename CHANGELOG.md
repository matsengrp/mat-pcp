# Changelog

This file tracks changes to the mat-pcp pipeline outputs and schema that affect downstream consumers.

## Version Policy

- **Simple integer versioning**: v1, v2, v3, etc.
- **Version bumps**: Any change to output format, filtering logic, or schema that affects downstream consumers should bump the version
- **Breaking changes**: In our use case, essentially any change is a breaking change since models trained on one version's outputs may not work with another's
- **Version in outputs**: The pipeline version is automatically appended to the output directory (e.g., `results/public_2022-06-01_v1`)

## How to Bump Version

1. Update the version number in the `VERSION` file in the repository root
2. Document the changes in this CHANGELOG under a new version heading
3. Update `OUTPUT_SCHEMA.md` if the schema has changed
4. Commit the changes and create a git tag matching the version (e.g., `git tag v2`)

## Versions

### v3

**Unify the count table, make the GLM the single neutral expectation, and simplify masking**
— tracks issue [#101](https://github.com/matsengrp/mat-pcp/issues/101).

**Changes**:

- **Single unified count table** (`build_count_table`, `workflow/scripts/build_count_table.py`):
  - New `counts_all.csv` replaces `expected_vs_actual_mut_counts.csv`. It enumerates the full
    grid of all possible mutations per clade (supplying the zero-`actual_count` rows the GLM
    trains on), joins the aggregated actual counts, and computes `exclude` / `masked_in_usher`
    **once**.
  - Removed rules/scripts: `synonymous_mut_rates`, `expected_mut_counts`,
    `merge_expected_and_actual_counts`, `aggregate_mutations_to_exclude`. Removed intermediates:
    `synonymous_mut_rates/rates_by_clade.csv`, `expected_mut_counts/expected_counts.csv`,
    `mutations_to_exclude.csv`, `expected_vs_actual_mut_counts.csv`.

- **Masking computed once** (`workflow/modules/masking.py`):
  - `count_mutations.py` no longer computes masking; `mutation_counts/{clade}.csv` and
    `mutation_counts/aggregated.csv` **lose the `exclude` column** (raw counts only).
  - `exclude` / `masked_in_usher` are computed in `build_count_table` via
    `masking.compute_exclusions`, from `sites_to_exclude`, the ProblematicSites VCF, the UShER
    yaml `sites`, and `site_include_range`.
  - The `exclude_ref_to_founder_muts` config key is removed. The count grid enumerates only
    mutations *away from the clade founder* (`nt_mutation` = founder_nt + site + to_nt), so
    clade-defining ref→founder mutations (`ref_nt` + site + founder_nt) can never appear in it
    and need no explicit exclusion — including in the GLM training set.

- **The GLM is the single neutral-expectation estimator**:
  - The synonymous 4-fold `expected_count` is gone. The GLM output column `predicted_count` is
    **renamed `expected_count`** everywhere, and is used as the neutral baseline for **both** the
    pseudocount (`naive_delta_fitness`) and the Bayesian (`delta_fitness`) estimators.
  - Fitness outputs (`fitness/*`, `aamut_fitness/*`) carry `expected_count` in place of
    `predicted_count`; the separate synonymous `expected_count` column no longer exists.
  - The `predicted_counts_by_clade` rule and its output `fitness/predicted_counts_by_clade.csv`
    are renamed to `expected_counts_by_clade` / `fitness/expected_counts_by_clade.csv`.
  - The `synonymous_spectra_min_counts` config key is removed.

- **Behavior:** `naive_delta_fitness` now compares actual counts against the GLM expectation
  (previously the synonymous-rate expectation); `delta_fitness` shifts slightly because
  ref→founder mutations are now excluded from GLM training.

### v2

**Bundled v2 cut** — tracks issues [#75](https://github.com/matsengrp/mat-pcp/issues/75), [#76](https://github.com/matsengrp/mat-pcp/issues/76), and [#80](https://github.com/matsengrp/mat-pcp/issues/80). [#81](https://github.com/matsengrp/mat-pcp/issues/81) ships the first three production runs against this version.

**Changes**:

- **Fused PCP extraction and curation** ([#75](https://github.com/matsengrp/mat-pcp/issues/75)):
  - The `curate_pcps_by_clade` rule, `workflow/scripts/curate_pcps_by_clade.py`, and `workflow/modules/pcp_validation.py` are removed.
  - `extract_pcps_by_clade` now writes the final curated `pcps/{clade}_pcps.csv` directly (no `mat/clades/{clade}_pcps_raw.csv` intermediate) and emits `pcps/{clade}_filtered_branches.csv` alongside it.
  - The identical-parent/child filter is applied inside `pcp_from_tree.py` after the gene-region slice, replacing the separate curation step.
  - Output schema is unchanged: `pcps/{clade}_pcps.csv`, `pcps/all_pcps.csv`, `pcps/{clade}_filtered_branches.csv`, `pcps/filtered_branch_ids.csv`, and `pcps/pcp_curation_report.txt` all keep their v1 columns and content.
  - The unit-test file `tests/test_pcp_validation.py` is removed along with the module it covered.

- **`min_clade_size` guardrail** ([#75](https://github.com/matsengrp/mat-pcp/issues/75)):
  - The Snakefile prints a `WARNING` to stderr at config-load time if any config matching `config-public-*.yaml` sets `min_clade_size > 0`. v1 production configs shipped with `min_clade_size: 10000` and biased the GLM rate estimation; the guardrail makes a recurrence of that mistake noisy.

- **DASM-ready PCP outputs** ([#80](https://github.com/matsengrp/mat-pcp/issues/80)):
  - New `pcps/{clade}_pcps_dasm_ready.csv` and `pcps/all_pcps_dasm_ready.csv` outputs land alongside the existing `*_pcps.csv` files. These absorb the post-processing that `data-central/scripts/filter_viral_pcps.sh` does today (invalid-PCP filtering, masked-only-mutation filtering, 90/10 train/val split).
  - Schema: all columns from `*_pcps.csv` are preserved verbatim; one new column `in_train` (boolean) is appended. `True` means the row is in the train split, `False` means val.
  - Seed convention: per-clade `seed=1729`, aggregate `seed=42` (matches the guardrail in `filter_viral_pcps.sh`).
  - Byte-equivalence to `dasm2 preprocess` is asserted on the test pipeline: the aggregate has md5 `5e538bd088e676b485a11a04745221e9` (sentinel in `validate_test.py`).
  - The filter is a NumPy-only reimplementation of `dasm2.preprocess_command.preprocess`, vendored into `workflow/scripts/curate_pcps_for_dasm.py` so the new rule stays in `envs/default.yml` (no torch/dasm2 dep).
  - Companion files emitted alongside each output: `{clade}_dasm_filter_report.txt` (human-readable per-filter drop counts) and `{clade}_dasm_filter_stats.json` (machine-readable).
  - Once v2 ships, these outputs make `data-central/scripts/filter_viral_pcps.sh` redundant.

- **Gene-sliced clade founder FASTAs** ([#76](https://github.com/matsengrp/mat-pcp/issues/76)):
  - `rule clade_founder_fasta` now also emits a gene-sliced founder FASTA at `clade_founders/{gene_dir}/{clade}.fa` whenever `pcp_gene_region` is configured. For the canonical Spike slice this is `clade_founders/spike/{clade}.fa` (3819 nt, positions 21563-25381, 1-based inclusive, excludes the trailing stop codon).
  - The slice uses the same `seq[start-1:end]` convention as PCP extraction, so the founder slice byte-aligns with the parent/child sequences in `pcps/{clade}_pcps.csv`.
  - Cross-version integrity sentinel: the 19A spike slice has md5 `b948f2a1f20edbce64cfff1dc540c119`; this is asserted in `validate_test.py`.
  - The `pcp_gene_region` block in all production configs and the slicing helper in `workflow/scripts/pcp_from_tree.py` are now documented as 1-based inclusive on both ends.
  - Once v2 ships, this replaces the post-hoc `data-central/scripts/slice_spike_founders.sh` step.

### v1

**Initial release** - Pipeline combining mat-pcp, SARS2-mut-fitness, and SARS2-mut-fitness-v2

**Output Schema**: See [OUTPUT_SCHEMA.md](OUTPUT_SCHEMA.md) for complete schema documentation.

**Key Features**:
- Parent-child pair (PCP) extraction from mutation-annotated trees
- Expected vs actual mutation count computation
- Neutral mutation rate estimation via GLM
- Nucleotide and amino acid fitness estimation
- Per-clade and aggregated outputs

**Main Output Files**:
- `pcps/all_pcps.csv` - Curated parent-child pairs across all clades
- `aamut_fitness/all_aamut_fitness.csv` - Amino acid fitness estimates
- `fitness/ntmut_fitness_aggregated.csv` - Nucleotide fitness estimates
- `rates_table.csv` - GLM-estimated neutral mutation rates
- `validation/pcp_vs_counts_report.txt` - Internal validation report

**Configuration Parameters**: Pipeline behavior is controlled by configuration files (see `config/config.yaml` for full parameter documentation).
