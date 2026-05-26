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

### v2 (unreleased)

**Bundled v2 cut** — VERSION bump and tag will land in the final PR of the v2 bundle. Tracks issues [#75](https://github.com/matsengrp/mat-pcp/issues/75), [#76](https://github.com/matsengrp/mat-pcp/issues/76), [#80](https://github.com/matsengrp/mat-pcp/issues/80), and [#81](https://github.com/matsengrp/mat-pcp/issues/81).

**Changes landed so far** (this entry will be appended to as the bundle progresses):

- **Fused PCP extraction and curation** ([#75](https://github.com/matsengrp/mat-pcp/issues/75)):
  - The `curate_pcps_by_clade` rule, `workflow/scripts/curate_pcps_by_clade.py`, and `workflow/modules/pcp_validation.py` are removed.
  - `extract_pcps_by_clade` now writes the final curated `pcps/{clade}_pcps.csv` directly (no `mat/clades/{clade}_pcps_raw.csv` intermediate) and emits `pcps/{clade}_filtered_branches.csv` alongside it.
  - The identical-parent/child filter is applied inside `pcp_from_tree.py` after the gene-region slice, replacing the separate curation step.
  - Output schema is unchanged: `pcps/{clade}_pcps.csv`, `pcps/all_pcps.csv`, `pcps/{clade}_filtered_branches.csv`, `pcps/filtered_branch_ids.csv`, and `pcps/pcp_curation_report.txt` all keep their v1 columns and content.
  - The unit-test file `tests/test_pcp_validation.py` is removed along with the module it covered.

- **`min_clade_size` guardrail** ([#75](https://github.com/matsengrp/mat-pcp/issues/75)):
  - The Snakefile prints a `WARNING` to stderr at config-load time if any config matching `config-public-*.yaml` sets `min_clade_size > 0`. v1 production configs shipped with `min_clade_size: 10000` and biased the GLM rate estimation; the guardrail makes a recurrence of that mistake noisy.

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
