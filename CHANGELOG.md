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

### v1 (Current)

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
