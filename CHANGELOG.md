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

### v8

**Always run metadata QC; remove `skip_qc_filtering` (issue
[#121](https://github.com/matsengrp/mat-pcp/issues/121)).** The `skip_qc_filtering`
option — a development convenience that made `filter_metadata` a pass-through to match
legacy SARS2-mut-fitness behavior — is removed. Metadata QC now always runs, so every
tree is filtered by date cutoff (`min_leaf_date` / `max_leaf_date`), date resolution
(MISSING/MONTH/YEAR dropped), host, and genome length before tree filtering. This is a
breaking change: production configs previously ran with QC skipped, so the set of
sequences entering the tree — and therefore all downstream counts, fitness, and windowed
outputs — changes.

**Changes**:

- Removed the `skip_qc_filtering` config key and the corresponding branch in
  `workflow/scripts/filter_metadata.py` / the `filter_metadata` rule params. Leftover
  `skip_qc_filtering` keys in any config are now silently ignored (unknown global keys are
  not validated).
- `max_leaf_date` (per-tree) now actually takes effect; it was a no-op whenever
  `skip_qc_filtering: true`. `subset_frac` likewise now applies on every run.
- Fixes the `chronumental_qc` crash on production configs: the QC path emits the
  `date_resolution_str` column that `chronumental_qc` requires, which the skip path never
  created (`ValueError: Usecols do not match columns` on the missing column).
- Bumped `chronumental_kwargs` from `--steps 1000` to `--steps 5000` in the production
  configs (`config/config.yaml`, `config/config-2026-04-22-gisaid-and-public-S-ORF1a.yaml`)
  as convergence insurance for the dating stage; the test config stays at `--steps 100`.
- Removed the redundant `pcp_max_nt_mutations` / `pcp_max_reversions_to_ref` /
  `pcp_max_reversions_to_clade_founder` config keys. PCP extraction now reuses the shared
  Mutation-Counting keys (`max_nt_mutations` / `max_reversions_to_ref` /
  `max_reversions_to_clade_founder`), so PCP extraction and clade/window mutation counting
  filter branches with the same thresholds. (Config values were already identical, so this
  is a no-op on the shipped configs.)

### v7

**Time-resolved fitness (EPIC [#115](https://github.com/matsengrp/mat-pcp/issues/115)).**
The coordinated version bump for the three phases below. Each phase was purely
additive (new output files, reading no existing schema and read by no current
downstream consumer) and shipped no individual bump per maintainer decision; this
is the single bump that lands the whole EPIC. The new outputs are `mat/node_dates*`
+ `mat/chronumental_timetree.nwk` (Phase 1), `windows/subtree_counts.csv` +
`windows/founder_muts.csv` (Phase 2), and `windows/fitness_trajectory.csv`
(Phase 3). Everything that existed at v6 is byte-identical apart from the version
dir name.

**chronumental node dating** — Phase 1 of EPIC
[#115](https://github.com/matsengrp/mat-pcp/issues/115) (issue
[#116](https://github.com/matsengrp/mat-pcp/issues/116)). Adds time estimates for
every node of the filtered tree; nothing existing changes.

**Changes**:

- New per-tree outputs under `mat/`: `node_dates.tsv` (chronumental-inferred
  date for every tip and internal node), `chronumental_timetree.nwk` (the dated
  time-scaled tree), plus the always-on QC artifacts `node_dates_qc.tsv` and
  `node_dates_residuals.csv`.
- New rules `export_newick`, `chronumental`, `chronumental_qc` and a new conda
  env `workflow/envs/chronumental.yml`.
- New config keys `chronumental_kwargs` and
  `chronumental_max_median_abs_tip_residual_days`.

**sliding-window subtree counts** — Phase 2 of EPIC
[#115](https://github.com/matsengrp/mat-pcp/issues/115) (issue
[#117](https://github.com/matsengrp/mat-pcp/issues/117)). Buckets observed branch
mutations into overlapping time windows over the full dated tree and emits, per
window, each subtree founder plus per-mutation actual counts. No expected counts
or `delta_fitness_ratio` yet (Phase 3).

**Changes**:

- New per-tree outputs under `windows/`: `subtree_counts.csv` (per-window,
  per-mutation actual counts — coding only, with the same multi-gene `;`-collapse
  as `aggregated.csv`) and `founder_muts.csv` (per-window subtree founders +
  reconstructed haplotypes).
- New module `workflow/modules/windowing.py` (`WindowedMAT` + pure helpers) and
  thin wrapper `workflow/scripts/window_subtree_counts.py`; new rule
  `window_subtree_counts` (bte env). Windows are anchored at the end date and
  tiled backwards, so the most recent window is always the full `window_size`.
- Branch QC is identical to the clade path (`count_mutations.py`): it reuses the
  existing `max_nt_mutations` / `max_reversions_to_ref` /
  `max_reversions_to_clade_founder` keys and the unconditional
  drop-branches-with->1-nt-mutation-per-codon filter.
- New config keys `window_size`, `window_step`, `window_start_time`,
  `window_end_time`.
- New correctness harness: `workflow/scripts/validate_windowing.py` (pure +
  output-CSV checks + test-data count regression sentinels, base env) and the
  CI-only `validate_windowing` rule
  (`workflow/scripts/validate_windowing_tree.py`, bte env). The tree-level rule
  is **not** in `results_files`.

**windowed expected counts + `delta_fitness_ratio`** — Phase 3 of EPIC
[#115](https://github.com/matsengrp/mat-pcp/issues/115) (issue
[#118](https://github.com/matsengrp/mat-pcp/issues/118)). Turns the Phase-2
per-window founders + actual counts into a mutation × time-window fitness
trajectory: for each window's subtree founder it enumerates all possible coding
mutations from the founder sequence and looks up mat-pcp's GLM neutral rate for
each mutation's founder-sequence context. Expected counts use mat-pcp's **own
clade-path `evol_time` normalization** (`rates.py`) applied **per window** — the
same neutral-model application as the clade-based `aamut_fitness`, so the two are
directly comparable — rather than s2trajectory's `rate × subtree_mut_count`. Then
`delta_fitness_ratio = log((actual + p) / (expected + p))` per mutation per
window (matching `aamutfit.ratio_fitness`).

**Changes**:

- New per-tree output `windows/fitness_trajectory.csv` (per `(window_start,
  amino-acid ID)`: actual count, expected count, `delta_fitness_ratio`). Coding
  only; ORF1a dropped in favor of ORF1ab (aamutfit convention).
- New module `workflow/modules/expected_counts.py` (a pure-`pandas`/`numpy`
  founder-context enumeration ported from s2trajectory's `PossibleMutations` — no
  `bte`) and thin wrapper `workflow/scripts/window_expected_counts.py`; new rule
  `window_expected_counts` (default env, after the GLM stage).
- Expected counts come **only** from the GLM (via `glm/master_table.csv`); there
  is no synonymous-average fallback. Per window,
  `T = (observed synonymous count) / (rate summed over all possible synonymous
  mutations)` and `expected(ID) = T × Σ_founders rate(context)`, so total expected
  synonymous equals total observed synonymous **by construction** (the calibration
  that makes synonymous `delta_fitness_ratio` center at 0). Because `T` divides out
  the rate's absolute scale, mat-pcp's `rate` is used as-is (no `n_sites`/fraction
  renormalization).
- New shared helper `masking.global_masked_sites` — the single definition of the
  clade-independent global site mask (config `sites_to_exclude` ∪ ProblematicSites
  VCF ∪ outside `site_include_range`), reused as the windowed path's "bad sites".
- **Founder→reference reversions are masked** on the windowed path, mirroring the
  clade path's `to_nt == ref_nt` exclusion (`build_count_table`): a possible
  mutation whose destination is the reference base (arising only at
  founder-mutated sites) is dropped from both the expected enumeration and the
  actual counts, so these artifact-prone reversions are treated identically to
  the clade-based path and windowed fitness stays comparable to `aamut_fitness`.
- New config key `window_min_subtree_size` (default 0 — keep all subtrees).
- `validate_windowing.py` gains a per-window synonymous-calibration proof
  (`sum(actual_syn) == sum(expected_syn)` per window, exact by construction), a
  synonymous-centering sanity check, a formula check, and test-data regression
  sentinels for the trajectory.

### v6

**Multi-protein pipeline** — generalizes PCP extraction, founder slicing, DASM
codon masking, DASM-ready filtering, and PCP-vs-counts validation from a single
hardcoded gene (Spike) to a config-driven list of proteins. Tracks EPIC
[#84](https://github.com/matsengrp/mat-pcp/issues/84).

**Changes**:

- **`pcp_genes` list replaces the scalar `pcp_gene_region` + `masking_gene`
  config keys.** A single Snakemake invocation fans PCP extraction, gene-sliced
  founders, codon masks, DASM-ready filtering, and validation over every protein
  in `pcp_genes` (must be non-empty). Production `config/config.yaml`
  ships `pcp_genes: [S]` to keep v5->v6 Spike outputs byte-identical; the
  test configs exercise the multi-protein path with `[S, ORF1ab, E]`.
  Per-protein coordinates are derived from the reference GTF
  (`cds_blocks_for_gene`) rather than a hardcoded region.

- **Per-gene output layout.** Outputs that were previously Spike-only now live
  under a per-gene subdirectory or filename (gene symbol verbatim):
  - `pcps/{gene}/{clade}_pcps.csv`, `pcps/{gene}/all_pcps.csv`, and the
    `pcps/{gene}/*_pcps_dasm_ready.csv` / `*_dasm_filter_report.txt` /
    `*_dasm_filter_stats.json` / `*_filtered_branches.csv` /
    `filtered_branch_ids.csv` / `pcp_curation_report.txt` siblings.
  - `clade_founders/{gene}/{clade}.fa` (was `clade_founders/spike/{clade}.fa`;
    the `GENE_FOUNDER_DIRNAMES` `S`->`spike` mapping is retired).
  - `{gene}_codon_sites_to_mask.txt` (was `spike_codon_sites_to_mask.txt`).
  - `validation/{gene}/pcp_vs_counts_report.txt` (+ `pcp_aa_counts.csv`,
    `pcp_vs_counts_discrepancies.csv`), one set per protein.

- **ORF1ab is handled as one frameshift-reconstructed CDS**, not split into
  nsps: its two GTF CDS blocks share nt 13468, giving a single in-frame
  21288 nt = 7096 codon founder slice. Gene-sliced founders are now emitted
  natively by the pipeline, superseding the post-hoc `slice_spike_founders.sh`
  step ([#97](https://github.com/matsengrp/mat-pcp/issues/97)).

- **Delta-fitness columns renamed in the amino-acid fitness outputs**
  ([#89](https://github.com/matsengrp/mat-pcp/issues/89)): `naive_delta_fitness`
  -> **`delta_fitness_ratio`** (the primary per-mutation signal) and
  `delta_fitness` -> **`delta_fitness_bayesian`** (inverse-variance-weighted
  estimate). Pure rename (values unchanged), applied consistently across the
  per-clade `{clade}_aamut_fitness.csv` files and the aggregate
  `all_aamut_fitness.csv`.

- **Consolidated production config**
  ([#113](https://github.com/matsengrp/mat-pcp/issues/113)): the four per-tree
  production configs are merged into a single four-tree `config/config.yaml`.
  `min_clade_size` now defaults to `0` for every tree; the `public_2022-06-01`
  tree was previously `10000`, which biased its GLM rates
  ([#81](https://github.com/matsengrp/mat-pcp/issues/81)) — its v6 rate outputs
  therefore differ from v5. Removed the `masking_gene` and `clade_synonyms`
  config keys and the soft `min_clade_size` Snakefile guardrail.

- **`get_problematic_sites` uses exact gene-membership matching**
  ([#114](https://github.com/matsengrp/mat-pcp/issues/114)): compound `gene`
  cells (e.g. `S;ORF3a`) are split on `;` and matched exactly, so a gene symbol
  can no longer bleed into a longer one (e.g. `ORF1a` into `ORF1ab`) once fanned
  over `pcp_genes`. Latent-correctness fix; no shipped config had colliding
  symbols, so the S mask is byte-identical (17 sites).

- **Internal cleanups** (no output-schema impact): consolidated the by-path
  `load_workflow_module` helper into a single `workflow/modules/workflow_loader.py`
  ([#111](https://github.com/matsengrp/mat-pcp/issues/111)); removed the orphan
  `validate_pcp_vs_counts.py` script and a dead `ref_fasta` rule input
  ([#94](https://github.com/matsengrp/mat-pcp/issues/94)); `validate_test.py`'s
  default output dir is VERSION-aware
  ([#88](https://github.com/matsengrp/mat-pcp/issues/88)).

- **Spike outputs are byte-identical to v5** for every Spike-specific artifact
  (PCPs, gene-sliced founders, codon mask, DASM-ready PCPs) apart from the delta
  column renames above and the version dir name. Regression anchors held: 19A
  Spike founder md5 `b948f2a1f20edbce64cfff1dc540c119`, S DASM-ready aggregate
  md5 `5e538bd088e676b485a11a04745221e9`, 17-site S codon mask.

### v5

**Run multiple trees from one config** — tracks issue
[#104](https://github.com/matsengrp/mat-pcp/issues/104).

**Changes**:

- **New `trees:` config schema**: the pipeline now processes a mapping of named
  trees in a single Snakemake invocation. Shared parameters (reference genome,
  GLM, genes, masking, thresholds, fitness, clade settings, and all
  `qc_filter_options`) stay at the top level; a new top-level `trees:` mapping is
  keyed by a short tree **label**, and each entry supplies that tree's `mat` and
  `metadata` (both required). `max_leaf_date` is the only per-tree-overridable
  key; it is optional and defaults to `None` (no upper date cutoff). When a
  global `qc_filter_options.max_leaf_date` is set, a tree without its own value
  inherits it.

- **`output_dir` removed; derived automatically**: output paths are now
  `{output_root}/{label}_{version}` (e.g. `results/public_2022-06-01_v5`), where
  `output_root` is an optional global key defaulting to `results`. The scalar
  `output_dir`, `mat`, and `metadata` top-level keys no longer exist.

- **Parameter renames** (parallel naming for the two leaf-date bounds):
  `qc_filter_options.datetime_cutoff` → `max_leaf_date` and
  `qc_filter_options.min_date` → `min_leaf_date`.

- **Site masking now applied to the fitness outputs**
  ([#106](https://github.com/matsengrp/mat-pcp/issues/106)): masked mutations
  (`exclude == True` in the count table — the union of UShER/ProblematicSites-VCF
  masks, config `sites_to_exclude`, sites outside `site_include_range`, and
  founder→reference reversions) are now dropped in `ntmut_counts_by_clade.py`
  before nucleotide and amino-acid fitness are computed. Previously the fitness
  branch flagged but never removed these rows, so `ntmut_fitness_aggregated.csv`
  and `all_aamut_fitness.csv` included artifactual counts from masked sites; the
  masking had only ever been applied to the GLM training data (`curate_counts`).
  Rows for masked mutations no longer appear in the fitness outputs. The
  filtering is per `(clade, mutation)`, so the `"all"` cross-clade aggregate
  still counts a mutation from clades where it is not masked.

- **`"all"` cross-clade nucleotide aggregation now pools by mutation identity.**
  `ntmut_counts_by_clade` previously grouped the `"all"` cluster by *every*
  non-count column, including clade-founder-derived annotations (`motif`,
  `clade_founder_codon`/`aa`, `mutant_codon`/`aa`, `aa_mutation`, `synonymous`,
  `four_fold_degenerate`, `mut_class`) that can differ between clades for the
  same nucleotide mutation. That silently split a mutation into multiple `"all"`
  rows (one per context) and fragmented its pooled counts. The `"all"` cluster
  now keys strictly on `(nt_site, nt_mutation)`, summing counts across clades
  into one row per mutation; the differing annotation columns keep a
  representative (first contributing clade's) value. Per-clade rows are
  unchanged. In the test data this collapses 36 fragmented mutations (all
  zero-count) in `ntmut_fitness_aggregated.csv`.

- **`exclude` / `masked_in_usher` columns removed from the fitness outputs.**
  After masking is applied, every surviving row has `exclude == False`, so these
  now-constant flags are dropped in `ntmut_counts_by_clade` rather than carried
  as uninformative columns into `fitness/ntmut_counts/*`, `fitness/ntmut_fitness/*`,
  and `fitness/ntmut_fitness_aggregated.csv`. (They remain in `counts_all.csv`.)

- **`tau_squared` removed from the amino-acid fitness outputs**
  (`aamut_fitness/{clade}_aamut_fitness.csv`, `aamut_fitness/all_aamut_fitness.csv`).
  `tau_squared` is a nucleotide-level log-scale variance of the neutral expected
  count, consumed only by `probfit` when estimating nucleotide fitness. It never
  fed the amino-acid `delta_fitness` (an inverse-variance weighted average of the
  nucleotide estimates) or `naive_delta_fitness` (from summed counts); it was
  merely summed over the constituent nucleotide mutations and carried as an
  unused column. It remains in the nucleotide outputs
  (`fitness/ntmut_fitness_aggregated.csv`), where it is used.

- **Per-tree output layout is byte-identical to v4 for everything except the
  fitness outputs**, apart from the version dir name: a single-tree run writes
  the same files under `{output_root}/{label}_v5/...` that a v4 run wrote under
  `{output_dir}_v4/...`, except that the fitness files
  (`fitness/ntmut_counts/*`, `fitness/ntmut_fitness/*`,
  `fitness/ntmut_fitness_aggregated.csv`, `aamut_fitness/*`) now exclude masked
  mutations, drop the `exclude`/`masked_in_usher` columns, and pool the `"all"`
  cluster by nucleotide identity (see above).

- **Multi-tree CI coverage** (test-infra, no schema impact):
  `config/config_test_data_multitree.yaml` and a step in the Test Pipeline job
  build two trees from one config and assert their `get_clades` checkpoints
  resolve independently — guarding against a regression where a per-tree
  `max_leaf_date` override is ignored.

### v4

**Overhaul UShER masking, mask all founder→reference reversions, and make the codon mask
clade-agnostic** — tracks issues [#99](https://github.com/matsengrp/mat-pcp/issues/99) /
[#100](https://github.com/matsengrp/mat-pcp/pull/100) and
[#103](https://github.com/matsengrp/mat-pcp/issues/103).

**Changes**:

- **New UShER mask source**: `data/usher_masked_sites.yaml` is now the generated,
  Nextstrain-clade-keyed file produced by the standalone
  [`s2-clade-mask`](https://github.com/matsengrp/s2-clade-mask) repository (from the UCSC
  branch-specific mask), replacing the hand-maintained file. Each block carries whole-site
  `sites` masks and an optional per-clade `mutations` list of directional reversions;
  `masking.compute_exclusions` consumes both (site masks and mutation masks set
  `exclude` + `masked_in_usher`).

- **Mask all founder→reference reversions**: any count-grid mutation whose destination base
  equals the reference base at that site (`to_nt == reference_nt`) — i.e. a clade-defining
  base reverting to the reference — is now excluded (`exclude = True`, but **not**
  `masked_in_usher`). `build_count_table` gains the reference FASTA as an input to supply the
  per-site reference base. This is unconditional (no config toggle). The GLM is unaffected
  (these sites are clade-defining, hence non-conserved, and were never in GLM training).

- **Clade-agnostic codon mask**: `get_problematic_sites` now emits a codon site in
  `spike_codon_sites_to_mask.txt` only if it is excluded **across all clades** (every row at
  the `nt_site` has `exclude == True`). This keeps clade-independent masks
  (`site_include_range`, `sites_to_exclude`, ProblematicSites VCF) and drops clade-specific
  ones (UShER per-clade masks and founder→reference reversions), which should not
  clade-specifically mask the target gene.

- **Config consolidation**: `additional_sites_to_ignore` is removed; its sites are merged into
  `sites_to_exclude` (de-duplicated), so all globally-excluded sites follow one downstream
  path (`exclude = True`). Original upstream comments from `jbloomlab/SARS2-mut-fitness` are
  restored on `sites_to_exclude`. The merged sites now show `exclude = True` in
  `counts_all.csv` (previously absent from that flag). GLM training set and the codon mask are
  otherwise unchanged.

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
  - The now-unused `clade_founders/{clade}_ref_to_founder_muts.csv` output is removed;
    `clade_founder_fasta` still applies those mutations in-memory to build the founder FASTA.

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
