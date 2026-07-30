# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`mat-pcp` is a unified Snakemake pipeline for analyzing SARS-CoV-2 genomic evolution from mutation-annotated trees (MATs). The pipeline combines three previously separate analysis workflows into a single coherent pipeline:

1. **mat-pcp** - Parent-child pair extraction from MATs
2. **SARS2-mut-fitness** - Expected vs actual mutation count computation
3. **SARS2-mut-fitness-v2** - Neutral mutation rate estimation via GLM

### Pipeline Stages

| Stage | Description | Key Rules |
|-------|-------------|-----------|
| 1-2 | Input Processing & Reference | `filter_metadata`, `filter_mat_tree`, `export_newick`, `chronumental`, `chronumental_qc`, `window_subtree_counts`, `get_ref_fasta`, `ref_coding_sites` |
| 3 | Clade Processing & PCPs | `get_clades` (checkpoint), `mat_clade_subset`, `extract_pcps_by_clade` |
| 4 | Mutation Counting | `count_mutations_clade`, `aggregate_mutation_counts` |
| 5 | Unified Count Table | `build_count_table` (full mutation grid + actual counts + masking flags, via `modules/masking.py`) |
| 6 | GLM Rate Estimation | `annotate_counts`, `curate_counts`, `master_table`, `get_rates_table`, `window_expected_counts` (windowed expected counts + `delta_fitness_ratio`) |
| 7 | Fitness Estimation | `ntmut_fitness_by_clade`, `aamut_fitness_by_clade`, `aggregate_*_fitness` |
| 8 | Validation | `extract_mutations_from_pcps`, `validate_pcp_vs_counts` |

## Versioning

The pipeline uses simple integer versioning (v1, v2, v3, etc.) for reproducibility and output tracking.

### Version Management

- **Version File**: `VERSION` in the repository root contains the current pipeline version
- **Automatic Versioning**: The Snakefile derives each tree's output directory as `{output_root}/{label}_{version}` (e.g., `results/public_2022-06-01_v5`)
- **Git Tags**: Each version should be tagged in git (e.g., `git tag v1`) when the branch is ready to merge
- **Changelog**: All output schema changes must be documented in `CHANGELOG.md`

### When to Bump Version

Any change that affects downstream consumers requires a version bump:
- Changes to output file formats or schemas (see `OUTPUT_SCHEMA.md`)
- Changes to filtering logic or curation criteria
- Changes to fitness estimation algorithms
- Changes to GLM rate estimation

In general, any change to the pipeline is a breaking change since models trained on one version's outputs may not work with another's.

**Additive-output carve-out.** A change that only *adds* new output files — reading no existing schema and read by no current downstream consumer — may ship without a version bump; the coordinated bump then rides with the next schema-breaking change. (Example: the additive outputs of EPIC [#115](https://github.com/matsengrp/mat-pcp/issues/115) — `mat/node_dates*` (Phase 1), `windows/subtree_counts.csv` + `founder_muts.csv` (Phase 2), `windows/fitness_trajectory.csv` (Phase 3) — each shipped no-bump; per maintainer decision the single coordinated bump for the whole EPIC landed as **v7** with Phase 4 (#119).)

### How to Bump Version

1. Update the version number in the `VERSION` file (e.g., from `v1` to `v2`)
2. Document all changes in `CHANGELOG.md` under a new version heading
3. Update `OUTPUT_SCHEMA.md` if the schema has changed
4. Update any affected documentation (README.md, this file)
5. Commit the changes
6. Create a git tag matching the version: `git tag v2`
7. Push the tag: `git push origin v2`

> **Also audit output-path defaults.** Any script with its own default output directory (currently `workflow/scripts/validate_test.py`) must track `VERSION` the same way the Snakefile does, or it will silently read a stale `results/<label>_v{N-1}` directory after a bump (see #88).

### Development with Versioning

When making changes that will affect outputs:
1. Work on a feature branch as usual
2. Update VERSION, CHANGELOG.md, and OUTPUT_SCHEMA.md as part of your changes
3. When the PR is ready to merge, the new version tag should be created

## Development Commands

### Running the Pipeline

Every run requires exactly one explicit `--configfile` (there is no default
`configfile:` directive; see the Snakefile header for why).

`config/config.yaml` lists **all four production trees**, one of which
(`gisaid-and-public_2023-01-01`) needs restricted GISAID-derived inputs. A
bare `rule all` run builds every tree; to build a single tree, target its
per-tree output path (see "Run specific target" below).

```bash
# Run with pixi (recommended)
pixi run snakemake -j 4 --use-conda --configfile config/config.yaml

# Dry run
pixi run snakemake -n --configfile config/config.yaml

# Run specific target (per-tree output path)
pixi run snakemake -j 4 --use-conda --configfile config/config.yaml \
  results/public_2022-06-01_v5/aamut_fitness/all_aamut_fitness.csv

# Clean all outputs
pixi run snakemake --delete-all-output --configfile config/config.yaml
```

### Testing
```bash
# Unit tests over frozen/synthetic fixtures (no pipeline run needed; CI runs
# all three before the pipeline steps)
pixi run python tests/test_glm_reference_equivalence.py
pixi run python tests/test_clade_founder_guards.py
pixi run python tests/test_founder_codon_background.py

# Run on test data (fast, uses max_clades=3)
pixi run snakemake -j 2 --use-conda \
  --configfile config/config_test_data.yaml

# Validate results
pixi run python workflow/scripts/validate_test.py

# Multi-tree smoke test: two trees in one config that resolve their
# get_clades checkpoints independently (the earlier max_leaf_date yields
# fewer clades). CI builds only each tree's clade_counts.csv and asserts
# the counts differ.
pixi run snakemake -j 2 --use-conda \
  --configfile config/config_test_data_multitree.yaml
```

## Architecture

### Directory Structure
```
mat-pcp/
├── config/
│   ├── config.yaml                  # Default config (all public/GISAID trees)
│   ├── config_test_data.yaml        # Test data configuration
│   └── config_test_data_multitree.yaml # Multi-tree smoke-test configuration
├── workflow/
│   ├── Snakefile            # Main workflow definition
│   ├── envs/                # Conda environments
│   │   ├── default.yml      # Python/pandas/biopython
│   │   ├── usher.yml        # UShER tools (matUtils)
│   │   └── bte.yml          # BTE library for tree traversal
│   ├── modules/             # Python modules
│   │   ├── glm.py           # General Linear Model for rate estimation
│   │   ├── rates.py         # Rate calculations
│   │   ├── probfit.py       # Probabilistic fitness estimation
│   │   ├── aamutfit.py      # Amino acid fitness aggregation
│   │   └── load.py          # Data loading utilities
│   └── scripts/             # Snakemake scripts
├── data/                    # Static data files
│   ├── usher_masked_sites.yaml
│   ├── lan_2022/           # RNA structure data
│   └── test_data/          # Test datasets and expected results
├── results/                 # Active pipeline outputs (gitignored)
│   ├── test_data/          # Test runs
│   └── {dataset}/          # Production runs
└── archives/                # Archived results (gitignored)
    └── README.md           # Archive documentation (tracked)
```

### Key Technologies
- **Snakemake**: Workflow management with checkpoints for dynamic clade processing
- **UShER/matUtils**: MAT tree manipulation and filtering
- **BTE**: Efficient tree traversal for PCP extraction
- **Pandas/NumPy**: Data processing and analysis
- **BioPython**: Sequence handling
- **SciPy/Statsmodels**: GLM rate estimation (in conda envs)

### Important Patterns

1. **Checkpoint pattern**: The `get_clades` rule is a checkpoint that determines the list of clades to process at runtime. All downstream per-clade rules depend on this checkpoint.

2. **Three conda environments**:
   - `default.yml` - Python, pandas, numpy, biopython
   - `usher.yml` - matUtils for tree manipulation
   - `bte.yml` - BTE library for tree traversal

3. **max_clades parameter**: Set `max_clades` in config to limit the number of clades processed. Useful for testing (CI uses `max_clades=3` for ~4 minute runs).

4. **Validation stage**: Cross-validates two independent mutation counting approaches (PCP-derived vs matUtils-derived) to ensure pipeline consistency.

5. **The GLM's synonymous-only filter lives in the loader, not in curation.**
   `glm/curated_counts.csv` deliberately holds *all* mutation classes (~22%
   synonymous on the v7/v8 test tree; the share is tree-specific) because
   downstream rules need the full grid.
   The neutral-rate restriction is applied when that file is *read*:
   `master_table.py` and `expected_counts_by_clade.py` both go through
   `modules/load.py::load_curated_muts`, an alias for `load_synonymous_muts`
   (`synonymous == True`), mirroring `SARS2-mut-fitness-v2`'s
   `master_tables.py.ipynb`. Inspecting `curated_counts.csv` therefore makes it
   *look* as though the GLM is fit on everything; check
   `log/master_table.log`'s "Loaded N curated mutations for training" instead.
   Replacing either loader call with a bare `pd.read_csv` would silently train
   the neutral model on mutations under purifying selection. Guarded by
   `tests/test_glm_reference_equivalence.py` (unit, frozen fixtures) and
   `validate_test.py::validate_glm_synonymous_training` (end-to-end). See
   [#126](https://github.com/matsengrp/mat-pcp/issues/126).

6. **An unresolvable clade founder is a hard failure, not a fallback.**
   `clade_founder_fasta.py` builds each clade's founder by applying that clade's
   defining mutations to the Wuhan reference. If the clade resolves in none of its
   three sources (Neher JSON → Roemer via `roemer_nextstrain_to_pango` → Roemer
   keyed directly), it used to apply *zero* mutations and write the unmodified
   reference — exit 0, logged as ordinary output. Nothing downstream caught it
   because the wrong founder propagated *consistently*: every consumer compared
   against the same wrong sequence, so no internal cross-check disagreed. It
   surfaced only when an outside consumer noticed implausible cross-clade
   correlations (24H/24I/25A/25B/25C in v8).
   Two guards now cover the two routes to that output — a missing *source*, and a
   resolved source that yields no *effective* mutations (empty substitution list,
   or every mutation dropped by the reference-base mismatch check). Both raise.
   The second check is deliberately **genome-wide, not per-gene**: at gene level
   reference-identity is legitimate (19B's founder equals 19A's for S, since
   C8782T and T28144C both fall outside spike), so a per-gene assertion would need
   an allowlist growing with every clade/gene pair — the same shape as the bug it
   guards. Guarded by `tests/test_clade_founder_guards.py`. See
   [#129](https://github.com/matsengrp/mat-pcp/issues/129).

7. **A clade's counts include only mutations off its founder's *codon*, not just
   its founder's nucleotide.** matUtils reports each branch's mutations relative
   to its parent, so the same nucleotide substitution also appears on codons that
   had already mutated — a different amino-acid change sharing one nucleotide
   event. `build_count_table.py` joined observed counts onto the founder grid on
   `(clade, nt_mutation)` alone, dropping `codon_change` (the only column that
   distinguishes the two) and using an unvalidated `how="left"` merge. Off-founder
   observations were therefore relabelled onto the founder codon and kept as
   *extra rows for the same mutation*: 6,415 duplicated keys / 13,146 rows in v8,
   every clade affected, concentrated in the highest-count mutations. Both
   estimators were corrupted by *different* mechanisms — `delta_fitness_bayesian`
   inverse-variance-averages the phantom row in, `delta_fitness_ratio` sums its
   expected count — so they never diverged enough to trip a cross-check, and the
   reported amino-acid uncertainty got *smaller*.
   The filter and the one-row-per-mutation invariant are ports of upstream's
   `merge_expected_and_actual_counts.py.ipynb`, lost when three rules were fused.
   Codon strings are `;`-joined across genes at overlap sites and compared
   whole, which requires the same gene order on both sides — true, and verified
   against the out-of-frame N/ORF9b overlap. A *structural* disagreement (one
   side coding, or a different component count) raises rather than reading as
   off-founder, because that would silently zero a legitimate count and the
   row-count invariant cannot see it. `aamutfit.aa_fitness` also rejects a
   repeated `nt_mutation` within one amino-acid group as independent defence.
   Guarded by `tests/test_founder_codon_background.py` (unit; the test dataset
   has no off-founder backgrounds, so it is the *only* coverage) and
   `validate_test.py::validate_founder_codon_background` (dataset invariant).
   See [#130](https://github.com/matsengrp/mat-pcp/issues/130).

### Configuration Sections

The configuration (`config/config.yaml`) is organized into:

1. **Input/Output**: `output_root` (optional, default `results`) and `trees` (a
   mapping of tree **label** → per-tree config). Each tree requires `mat` and
   `metadata`; `max_leaf_date` is the only per-tree-overridable key. Output dirs
   are derived as `{output_root}/{label}_{version}` — there is no `output_dir`
   key. A single run builds every tree in `trees`.
2. **Reference Genome**: `ref_fasta`, `ref_gtf`, `add_to_ref_gtf`
3. **Metadata QC**: `qc_filter_options` (date filtering via `min_leaf_date` /
   `max_leaf_date`, host validation, genome length). Global except that a tree
   may override `max_leaf_date`. `max_leaf_date` is optional; when unset it
   defaults to `None`, meaning no upper date cutoff.
4. **Tree Processing**: `optimize_tree`, `mat_subset_kwargs`
4b. **Node Dating (chronumental)**: `chronumental_kwargs` (passed verbatim to
   `chronumental`, e.g. `--steps`), `chronumental_max_median_abs_tip_residual_days`
   (the `chronumental_qc` gate — the run fails if the median absolute tip-date
   residual exceeds it). The QC rule runs on every tree.
4c. **Window Fitness (sliding windows)** (EPIC #115 Phases 2-3): `window_size`,
   `window_step`, `window_start_time`, `window_end_time`, `window_min_subtree_size`.
   The `window_subtree_counts` rule (Phase 2) buckets branch mutations into
   overlapping time windows (anchored at the end date, tiled backwards) and emits
   the `windows/` actual-count outputs. Its branch QC reuses the Mutation-Counting
   keys (`max_nt_mutations` / `max_reversions_to_ref` /
   `max_reversions_to_clade_founder`), so the clade-based and window-based paths
   filter branches identically. The `window_expected_counts` rule (Phase 3, #118,
   `modules/expected_counts.py` — a pure-pandas port of s2trajectory's
   `PossibleMutations`, no bte) then computes per-window expected counts from the
   GLM neutral model (`glm/master_table.csv`) and emits
   `windows/fitness_trajectory.csv` with
   `delta_fitness_ratio = log((actual + p) / (expected + p))`. Expected counts use
   the GLM only (no synonymous-average fallback) and **mat-pcp's own clade-path
   `evol_time` normalization applied per window** (`T = observed_syn / rate_syn`;
   `expected = T × rate`), NOT s2trajectory's `rate × subtree_mut_count` — so total
   expected synonymous == total observed synonymous per window by construction, and
   windowed fitness is directly comparable to the clade-based `aamut_fitness`.
5. **Clade Configuration**: `clade_column`, `min_clade_size`, `clades_to_exclude`, `max_clades`
6. **PCP Extraction**: `pcp_genes` (PCP branch QC reuses the shared Mutation-Counting keys `max_nt_mutations` / `max_reversions_to_ref` / `max_reversions_to_clade_founder`)
7. **Clade Founders**: `clade_founders` (Neher and Roemer JSON sources),
   `roemer_nextstrain_to_pango` (Nextstrain clade → Pango lineage, for clades not
   in the Neher JSON), `allow_reference_identical_founders` (optional; clades
   whose founder may legitimately equal the reference genome — the
   `reference_clade` is exempt automatically, so this is normally empty).
   **Every config carries its own copy of `roemer_nextstrain_to_pango`**, so
   adding a clade to one config does not add it to the others; that is how #129
   happened. `clade_founder_fasta.py` now raises rather than silently writing the
   reference as an unresolvable clade's founder — see the note under "Important
   Patterns" below.
8. **Mutation Counting**: `max_nt_mutations`, `max_reversions_*`, `site_include_range`
9. **GLM**: `rna_structure_data`, `glm_lightswitch_boundaries`
11. **Fitness**: `fitness_pseudocount`, `fitness_N_f`, `fitness_N_n_exp`
12. **Gene Configuration**: `orf1ab_to_nsps`, `gene_overlaps`, `genes`
13. **Site Masking**: `sites_to_exclude`, `site_include_range`
14. **Threads**: Per-rule thread allocation

### Output Organization

Results are stored per tree in `{output_root}/{label}_{version}/`:
- `mat/`: Filtered trees, clade subtrees
- `windows/`: Sliding-window subtree counts + founders (Phase 2) + per-window expected counts / `delta_fitness_ratio` trajectory (Phase 3)
- `pcps/`: Curated parent-child pairs
- `ref/`: Reference genome, coding sites
- `mutation_counts/`: Per-clade and aggregated counts
- `clade_founders/`: Founder sequences and mutations
- `glm/`: Annotated counts, curated counts, master table
- `fitness/`: Predicted counts, nucleotide fitness
- `aamut_fitness/`: Amino acid fitness
- `validation/`: PCP vs counts validation report
- `log/`: Processing logs

### Data Files
- `data/usher_masked_sites.yaml`: UShER masking information for problematic sites by clade
- `data/lan_2022/41467_2022_28603_MOESM11_ESM.txt`: RNA secondary structure data

## Pre-PR Quality Checklist

Before any pull request:

1. **Issue Compliance**: Verify completion of all requirements from the linked issue
2. **Versioning**: If changes affect outputs or downstream consumers:
   - Bump the version number in the `VERSION` file
   - Document changes in `CHANGELOG.md` under the new version
   - Update `OUTPUT_SCHEMA.md` if schema changed
   - Plan to tag the version in git when PR is ready to merge
3. **Documentation Sync**: If changes affect outputs or behavior:
   - Update `OUTPUT_SCHEMA.md` for any changes to output file formats
   - Update `README.md` for new outputs, changed behavior, or configuration options
   - Update this file (`CLAUDE.md`) for architectural changes or new patterns
4. **Test on Test Data**: Run the pipeline on test data to verify changes work:
   ```bash
   pixi run snakemake -j 2 --use-conda --configfile config/config_test_data.yaml
   ```
5. **Validate Results**: Run validation to check for regressions:
   ```bash
   pixi run python workflow/scripts/validate_test.py
   ```
6. **PR Description**: Ensure PR description accurately reflects the scope of changes
