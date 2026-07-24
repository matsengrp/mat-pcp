# Output Schema Reference

This document describes the column schemas for the primary consumer-facing output files produced by the mat-pcp pipeline. These schemas are fixed and do not vary based on configuration.

For information on running the pipeline, see [README.md](README.md).

## Parent-Child Pairs (PCPs)

### `pcps/{gene}/all_pcps.csv` and `pcps/{gene}/{clade}_pcps.csv`

Curated parent-child sequence pairs extracted from the mutation-annotated tree. These represent evolutionary transitions along tree branches and are suitable for training sequence models (e.g., DASM).

| Column | Type | Description |
|--------|------|-------------|
| `clade` | string | Nextstrain clade identifier (e.g., "19A", "20B", "21K"). Only present in `all_pcps.csv`. |
| `parent_name` | string | Internal tree node identifier for the parent node. May be a sample accession (for leaf nodes) or an internal node ID. |
| `child_name` | string | Internal tree node identifier for the child node. |
| `parent` | string | Parent nucleotide sequence, sliced to this protein's coding sequence. PCP extraction is fanned over the `pcp_genes` config list; each protein's CDS blocks are derived from the reference GTF (`cds_blocks_for_gene`), so outputs live under a per-gene subdirectory (`pcps/{gene}/`). ORF1ab is reconstructed as one frameshifted CDS (its two blocks share nt 13468). |
| `child` | string | Child nucleotide sequence after applying branch mutations. Same length as `parent`. |
| `branch_length` | integer | Number of nucleotide mutations on this branch (within the gene region). |

**Notes:**
- Sequences are already sliced to the gene region of interest, so they can be used directly for training.
- PCPs have been curated to remove branches with identical parent/child sequences (in the gene region). As of v2, this filter is applied during extraction in the same pass that builds the table; there is no separate `mat/clades/{clade}_pcps_raw.csv` file.
- The filtering criteria (max mutations, max reversions) are controlled by the shared Mutation-Counting config parameters `max_nt_mutations`, `max_reversions_to_ref`, and `max_reversions_to_clade_founder` (as of v8; PCP extraction and clade/window mutation counting filter branches identically).

### `pcps/{gene}/all_pcps_dasm_ready.csv` and `pcps/{gene}/{clade}_pcps_dasm_ready.csv` (v2+)

DASM-training-ready sibling of `*_pcps.csv`. Same row schema as the parent file plus one column (`in_train`). The filter set + train/val split are applied in the pipeline so downstream consumers (e.g. `data-central` manifests) can read these files directly without any post-processing step.

| Column | Type | Description |
|--------|------|-------------|
| (all columns from `*_pcps.csv`) | — | Preserved verbatim. |
| `in_train` | boolean | `True` for the train split, `False` for the val split. The 90/10 split is realized by `np.random.default_rng(seed).random(N) < val_fraction` (with `val_fraction=0.10`), then `in_train = ~val_mask`. |

**Filters applied (in order)**:

1. Empty parent OR empty child sequence
2. Length mismatch parent vs child
3. Length not divisible by 3 (non-codon-aligned)
4. Ambiguous nucleotides (`N`) in either sequence
5. Stop codons (`TAA`, `TAG`, `TGA`) in parent or child
6. Every codon-level mutation falls in a codon position listed in `<gene>_codon_sites_to_mask.txt` (no-signal PCP)

**Seed convention**:

- **Per-clade**: `seed=1729` (matches `data-central/scripts/filter_viral_pcps.sh`).
- **Aggregate (`all_pcps_dasm_ready.csv`)**: `seed=42` (matches the explicit-seed guardrail in `filter_viral_pcps.sh`; the different seed prevents accidentally producing the aggregate with the per-clade default).

**Byte-equivalence to `dasm2 preprocess`**: this rule is a NumPy-only reimplementation of `dasm2.preprocess_command.preprocess` (which `filter_viral_pcps.sh` wraps). On the test pipeline, the aggregate file is byte-identical to the dasm2 reference output; md5 `5e538bd088e676b485a11a04745221e9` is asserted in `validate_test.py`.

**Companion files**:

- `pcps/{gene}/{clade}_dasm_filter_report.txt` — human-readable per-filter drop counts.
- `pcps/{gene}/{clade}_dasm_filter_stats.json` — machine-readable equivalent for downstream tooling.

### `pcps/{gene}/{clade}_filtered_branches.csv`

Branch identifiers that were excluded during PCP curation. Used internally to harmonize mutation counting with PCP extraction.

| Column | Type | Description |
|--------|------|-------------|
| `child_name` | string | Tree node identifier for branches that failed curation filters. |

### `pcps/{gene}/filtered_branch_ids.csv`

Aggregated list of all filtered branch IDs across all clades.

| Column | Type | Description |
|--------|------|-------------|
| `child_name` | string | Tree node identifier for filtered branches. |
| `clade` | string | Clade from which the branch was filtered. |

## Fitness Estimates

### `aamut_fitness/all_aamut_fitness.csv`

Amino acid-level fitness estimates aggregated across all clades. Fitness is computed as log(actual/expected) with Bayesian uncertainty quantification.

| Column | Type | Description |
|--------|------|-------------|
| `gene` | string | Gene name (e.g., "S" for Spike, "N" for Nucleocapsid). |
| `clade_founder_aa` | string | Single-letter amino acid in the clade founder sequence. |
| `mutant_aa` | string | Single-letter mutant amino acid. |
| `aa_site` | integer | 1-indexed amino acid position within the gene. |
| `aa_mutation` | string | Mutation string in standard format (e.g., "D614G"). |
| `subset_of_ORF1ab` | boolean | True if this entry is for an NSP (ORF1ab subset); False otherwise. |
| `expected_count` | float | GLM-predicted neutral count (sequence context + RNA structure + position); the baseline for both `delta_fitness_ratio` and `delta_fitness_bayesian`. |
| `actual_count` | float | Observed mutation count in the data. |
| `delta_fitness_ratio` | float | Pseudocount log-ratio fitness: `log((actual + 0.5) / (expected + 0.5))`. The primary per-mutation signal rendered by the dasm2 dashboard (renamed from `naive_delta_fitness` in v6). |
| `delta_fitness_bayesian` | float | Inverse-variance-weighted Bayesian fitness estimate accounting for uncertainty (renamed from `delta_fitness` in v6). |
| `uncertainty` | float | Standard error of the `delta_fitness_bayesian` estimate. |
| `source_clade` | string | Semicolon-separated list of clades that contributed to this aggregated estimate. |

**Notes:**
- Positive `delta_fitness_ratio` / `delta_fitness_bayesian` indicates the mutation is observed more often than expected (potentially beneficial).
- Negative values indicate the mutation is observed less often than expected (potentially deleterious).
- The `uncertainty` column can be used to filter for high-confidence estimates.
- **Masking:** masked mutations (`exclude == True` in the count table) are dropped before fitness is computed (in `ntmut_counts_by_clade`), so no row here reflects counts from masked sites. Because masking is per-clade, a mutation masked only in some clades still contributes to the aggregate from clades where it is not masked.

### `aamut_fitness/{clade}_aamut_fitness.csv`

Per-clade amino acid fitness estimates. Same schema as above but with `cluster` instead of aggregated fields:

| Column | Type | Description |
|--------|------|-------------|
| `cluster` | string | Clade name (e.g., "19A"). Note: named "cluster" for historical compatibility. |

(Other columns same as `all_aamut_fitness.csv`, excluding `source_clade`.)

> **Column names:** these two delta columns were renamed in v6 (`naive_delta_fitness` → `delta_fitness_ratio`, `delta_fitness` → `delta_fitness_bayesian`) and the new names are used consistently across the per-clade `{clade}_aamut_fitness.csv` files and the aggregate `all_aamut_fitness.csv`.

### `fitness/ntmut_fitness_aggregated.csv`

Nucleotide-level fitness estimates. Similar structure to amino acid fitness but at the nucleotide level.

| Column | Type | Description |
|--------|------|-------------|
| `clade` | string | Clade name for this entry. |
| `nt_mutation` | string | Nucleotide mutation (e.g., "A23403G"). |
| `nt_site` | integer | 1-indexed genomic position. |
| `clade_founder_nt` | string | Nucleotide in the clade founder. |
| `gene` | string | Gene containing this site (semicolon-separated if multiple). |
| `aa_mutation` | string | Amino acid mutation if nonsynonymous (semicolon-separated if multiple). |
| `synonymous` | boolean | True if mutation is synonymous. |
| `actual_count` | float | Observed mutation count. |
| `expected_count` | float | GLM-predicted neutral count (sequence context + RNA structure + position). |
| `tau_squared` | float | Variance component from GLM. |
| `f_max` | float | Maximum a posteriori (MAP) fitness estimate. |
| `left_conf_int` | float | Left 68% credible interval distance from `f_max`. |
| `right_conf_int` | float | Right 68% credible interval distance from `f_max`. |
| `f_mean` | float | Posterior mean fitness estimate. |
| `f_st_dev` | float | Posterior standard deviation of fitness estimate. |
| `source_clade` | string | Semicolon-separated list of clades that contributed to this entry. |

**Notes:**
- Additional technical columns (e.g., `motif`, `unpaired`, `codon_site`, `mut_class`) are present for internal pipeline use but are not documented here.
- **Masking:** as with the amino-acid fitness, masked mutations (`exclude == True`) are dropped in `ntmut_counts_by_clade` before nucleotide fitness is computed, so masked sites do not appear here. The `exclude` / `masked_in_usher` flag columns from `counts_all.csv` are not carried into this file.
- **`cluster == "all"` rows are pooled per nucleotide mutation.** The cross-clade aggregate keys on `(nt_site, nt_mutation)` and sums counts across clades, so there is exactly one `"all"` row per mutation. Annotation columns that derive from the clade founder and can differ between clades (`motif`, `clade_founder_codon`/`aa`, `mutant_codon`/`aa`, `aa_mutation`, `synonymous`, …) hold a representative (first contributing clade's) value on the `"all"` row; the per-clade rows retain the exact per-clade values.

## Mutation Counts

### `mutation_counts/aggregated.csv`

Raw mutation counts aggregated across all clades, with annotations. Masking is no longer
applied here; the `exclude` / `masked_in_usher` flags are computed once, on the full mutation
grid, in `counts_all.csv` (see below).

| Column | Type | Description |
|--------|------|-------------|
| `clade` | string | Clade name. |
| `nt_mutation` | string | Nucleotide mutation (e.g., "A23403G"). |
| `nt_site` | integer | 1-indexed genomic position. |
| `count` | integer | Total mutation count (terminal + non-terminal). |
| `count_terminal` | integer | Count on terminal (leaf) branches. |
| `count_non_terminal` | integer | Count on internal branches. |
| `synonymous` | boolean | True if the mutation is synonymous. |
| `noncoding` | boolean | True if the site is in a non-coding region. |
| `protein` | string | Protein name if in a coding region. |
| `aa_mutation` | string | Amino acid mutation if nonsynonymous. |
| `codon_change` | string | Codon change (e.g., "GAT>GGT"). |
| `reference_nt` | string | Nucleotide in Wuhan-Hu-1 reference. |
| `clade_founder_nt` | string | Nucleotide in the clade founder sequence. |
| `mean_log_size` | float | Mean log of subtree sizes for branches with this mutation. |

### `counts_all.csv`

The unified per-clade count table: every possible mutation in every processed clade (the
zero-`actual_count` rows are what the GLM trains on), with actual counts joined in and the
masking flags computed once via `workflow/modules/masking.py`.

| Column | Type | Description |
|--------|------|-------------|
| `clade` | string | Clade name. |
| `nt_site` | integer | 1-indexed genomic position. |
| `nt_mutation` | string | Nucleotide mutation (e.g., "A23403G"). |
| `mut_type` | string | Two-letter substitution type (e.g., "CT"). |
| `exclude` | boolean | True if the mutation is masked for any reason: config `sites_to_exclude`, the ProblematicSites VCF, the UShER yaml (per-clade site or mutation masks), a founder→reference reversion (destination base equals the reference base), or out of `site_include_range`. |
| `masked_in_usher` | boolean | True if excluded specifically by an artifactual-site mask (UShER yaml or the ProblematicSites VCF). Founder→reference reversions set `exclude` but not this flag. |
| `actual_count` | integer | Observed mutation count (0 if unobserved). |
| `synonymous` | boolean | True if the mutation is synonymous. |
| `noncoding` | boolean | True if the site is non-coding. |

**Notes:**
- Additional annotation columns (`gene`, `clade_founder_nt`, `clade_founder_codon`,
  `mutant_codon`, `mutant_aa`, `aa_mutation`, `codon_position`, `codon_site`,
  `four_fold_degenerate`, `count_terminal`, `count_non_terminal`, `subset`, `mean_log_size`)
  are present for internal pipeline use.
- Replaces the former `expected_vs_actual_mut_counts.csv`; there is no `expected_count` here —
  the neutral expectation is estimated downstream by the GLM.

## Clade Founder Sequences

### `clade_founders/{clade}.fa`

Whole-genome founder sequence (FASTA) for each clade, computed by applying the clade's founder mutations (sourced from Neher and/or Roemer) to the Wuhan-Hu-1 reference (NC_045512.2). One sequence record per file; the header is the clade name (e.g. `>19A`). Length: 29903 nt.

### `clade_founders/{gene}/{clade}.fa`

Gene-sliced founder sequence for each clade and each protein in the `pcp_genes` config list. Contains just that protein's coding sequence, taken from the reference GTF's CDS blocks (via `cds_blocks_for_gene`) and excluding the trailing stop codon; the subdirectory is the gene symbol verbatim (e.g. `clade_founders/S/`, `clade_founders/ORF1ab/`, `clade_founders/E/`). ORF1ab is reconstructed as one in-frame CDS across the ribosomal frameshift — its two GTF blocks share nt 13468 — giving 21288 nt = 7096 codons; Spike is positions 21563-25381 (1-based, inclusive), 3819 nt = 1273 codons. Each slice byte-aligns with the `parent`/`child` sequences in the matching `pcps/{gene}/{clade}_pcps.csv`, so it can be used directly as the wildtype for training ParentIndependent DASM models.

**Cross-version integrity sentinel**: the 19A Spike slice (`clade_founders/S/19A.fa`, Wuhan-Hu-1 reference founder) has md5 `b948f2a1f20edbce64cfff1dc540c119` and is byte-stable across all SARS-CoV-2 cutoff dates.

## Neutral Mutation Rates

### `rates_table.csv`

Estimated neutral mutation rates by sequence context and RNA structure, derived from synonymous mutations via GLM.

| Column | Type | Description |
|--------|------|-------------|
| `mut_type` | string | Mutation type (e.g., "C>T", "G>A"). |
| `motif` | string | Trinucleotide sequence context (e.g., "TCG" for CpG context). |
| `unpaired` | boolean | True if the site is in an unpaired (single-stranded) RNA region. |
| `rate` | float | Estimated neutral mutation rate for this context. |

**Notes:**
- These rates are used to compute expected mutation counts for fitness estimation.
- The `unpaired` column reflects RNA secondary structure from Lan et al. 2022.

## Node Dates (chronumental)

Time estimates for every node of the filtered tree, produced by the
`chronumental` rule and sanity-checked by the always-on `chronumental_qc` rule.
Added in Phase 1 of EPIC [#115](https://github.com/matsengrp/mat-pcp/issues/115).

### `mat/node_dates.tsv`

chronumental's `--dates_out`: an inferred date for every tip **and** internal
node of `mat/filtered_tree.nwk`.

| Column | Type | Description |
|--------|------|-------------|
| `strain` | string | Node identifier — a strain name for tips, an internal-node id (e.g. `node_123`) for internal nodes. Matches the `bte` `node.id` used downstream. |
| `predicted_date` | date | chronumental-inferred date (`YYYY-MM-DD`). |

### `mat/chronumental_timetree.nwk`

chronumental's `--tree_out`: the input tree rescaled so branch lengths are in
time, with the same topology and node labels as `mat/filtered_tree.nwk`.

### `mat/node_dates_qc.tsv`

Per-check pass/fail report from `chronumental_qc` (one row per check).

| Column | Type | Description |
|--------|------|-------------|
| `check` | string | `node_id_correspondence`, `time_monotonicity`, `date_coverage`, or `tip_residual`. |
| `status` | string | `PASS` or `FAIL`. |
| `detail` | string | Measured numbers for the check. |

The rule raises (fails the run) on any hard-invariant violation, and on a
`tip_residual` median above `chronumental_max_median_abs_tip_residual_days`.

### `mat/node_dates_residuals.csv`

Per-tip comparison of predicted vs. observed collection date (day-resolution tips
only).

| Column | Type | Description |
|--------|------|-------------|
| `strain` | string | Tip strain name. |
| `predicted` | date | chronumental-inferred date. |
| `observed` | date | Metadata collection date. |
| `residual_days` | int | `predicted − observed` in days. |

## Sliding-Window Subtree Counts (windows/)

EPIC [#115](https://github.com/matsengrp/mat-pcp/issues/115) Phase 2
([#117](https://github.com/matsengrp/mat-pcp/issues/117)). The
`window_subtree_counts` rule buckets observed branch mutations into overlapping
time windows over the full dated tree. Windows are anchored at the end date
(`window_end_time`, default = latest node date) and tiled backwards in
`window_step`-day strides, so the most recent window is always the full
`window_size`; the oldest window may extend before the start date. Branch QC is
identical to the clade path (`count_mutations.py`). Counts are **coding only**;
per-window expected counts and `delta_fitness_ratio` are added by
`window_expected_counts` in `windows/fitness_trajectory.csv` (Phase 3, below).

### `windows/subtree_counts.csv`

One row per `(window_start, founder_id, nt_mutation)`; genes overlapping a site
are `;`-collapsed exactly as in `aggregated.csv`.

| Column | Type | Description |
|--------|------|-------------|
| `window_start` | date | Window start (`YYYY-MM-DD`). |
| `founder_id` | string | `bte` node id of the subtree founder. |
| `protein` | string | `;`-joined gene name(s) for the site (e.g. `ORF1a;ORF1ab`). |
| `nt_site` | int | 1-based nucleotide position. |
| `nt_mutation` | string | Nucleotide mutation (e.g. `C8782T`). |
| `aa_mutation` | string | `;`-joined per-protein amino-acid change(s). |
| `codon_change` | string | `;`-joined `orig>alt` codon(s). |
| `synonymous` | bool | `all()` over the collapsed genes' synonymy. |
| `count` | int | # subtree branches carrying exactly this mutation (shared across collapsed genes, not summed). |

### `windows/founder_muts.csv`

One row per `(founder_id, window_start)` subtree founder.

| Column | Type | Description |
|--------|------|-------------|
| `founder_id` | string | `bte` node id of the subtree founder. |
| `window_start` | date | Window start (`YYYY-MM-DD`). |
| `subtree_size` | int | # nodes in the subtree, including the founder. |
| `branches_with_mutations` | int | # subtree branches with ≥1 mutation. |
| `founder_nt_mutations` | string | Space-joined founder haplotype vs. reference (e.g. `C241T C3037T`); empty when the founder equals the reference. |

### `windows/fitness_trajectory.csv`

EPIC [#115](https://github.com/matsengrp/mat-pcp/issues/115) Phase 3
([#118](https://github.com/matsengrp/mat-pcp/issues/118)). The
`window_expected_counts` rule turns the per-window founders + actual counts into
a mutation × time-window fitness trajectory. For each window's subtree founder it
enumerates all possible coding mutations from the founder sequence and looks up
mat-pcp's GLM neutral rate for each mutation's founder-sequence context. Expected
counts use mat-pcp's own clade-path `evol_time` normalization (`rates.py`) applied
per window: `T = (observed synonymous count) / (rate summed over all possible
synonymous mutations)` and `expected(ID) = T × Σ_founders rate(context)`, so total
expected synonymous equals total observed synonymous by construction. Expected
counts come only from the GLM (no synonymous-average fallback). Coding only; ORF1a
is dropped in favor of ORF1ab (matching `aamut_fitness`). **Founder→reference
reversions** (a possible mutation whose destination is the reference base, arising
only at founder-mutated sites) are masked from both the expected enumeration and
the actual counts, mirroring the clade path's `to_nt == ref_nt` exclusion. A window
with **no observed synonymous mutations** cannot be calibrated (its `T` is
undefined) and is omitted, so this file may cover fewer windows than
`subtree_counts.csv`.

One row per `(window_start, ID)`, where `ID` is the amino-acid mutation identity
`<gene>-<wt_aa><aa_site><mut_aa>`.

| Column | Type | Description |
|--------|------|-------------|
| `window_start` | date | Window start (`YYYY-MM-DD`). |
| `ID` | string | Amino-acid mutation identity, `<gene>-<wt_aa><aa_site><mut_aa>` (e.g. `S-D614G`). |
| `gene` | string | Gene name (ORF1a folded into ORF1ab). |
| `aa_site` | int | Codon (amino-acid) site within the gene. |
| `aa_mutation` | string | Amino-acid change `<wt_aa><aa_site><mut_aa>`. |
| `wt_aa` | string | Founder (wild-type) amino acid. |
| `mut_aa` | string | Mutant amino acid (`*` = stop). |
| `synonymous` | bool | `wt_aa == mut_aa`. |
| `actual_count` | float | Observed count summed over the window's subtrees. |
| `expected_count` | float | GLM-predicted count, `rate × subtree_mut_count` summed over subtrees (can be slightly negative for sparse contexts). |
| `delta_fitness_ratio` | float | `log((actual_count + p) / (expected_count + p))`, natural log, pseudocount `p = fitness_pseudocount`. Synonymous mutations center near 0. `NaN` in the rare case `expected_count + p ≤ 0` (a very negative sparse-context rate), matching the clade-path `aamut_fitness` `delta_fitness_ratio`. |
