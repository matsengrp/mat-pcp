# Output Schema Reference

This document describes the column schemas for the primary consumer-facing output files produced by the mat-pcp pipeline. These schemas are fixed and do not vary based on configuration.

For information on running the pipeline, see [README.md](README.md).

## Parent-Child Pairs (PCPs)

### `pcps/{gene}/all_pcps.csv.gz` and `pcps/{gene}/{clade}_pcps.csv.gz`

Curated parent-child sequence pairs extracted from the mutation-annotated tree. These represent evolutionary transitions along tree branches and are suitable for training sequence models (e.g., DASM).

> **Renamed in v9 (#134):** these were `.csv` through v8. Contents are unchanged — the column schema below is identical — but every path gained a `.csv.gz` suffix. `pandas.read_csv` and `dasm2.data.registry.load_dataframe` both handle it transparently; the latter matches on the extension, so the name must end in exactly `.csv.gz`. The motivation was size: on `gisaid-and-public_2026-04-22_v8`, `pcps/` was 504 GB apparent (470 GiB), 93% of the entire 542 GB output tree. Measured gzip level 9 gives ~125x on ORF1a and ~110x on S; on the v9 run `pcps/` is 7.7 GB apparent, across five genes where v8 had two.

| Column | Type | Description |
|--------|------|-------------|
| `clade` | string | Nextstrain clade identifier (e.g., "19A", "20B", "21K"). Only present in `all_pcps.csv.gz`. |
| `parent_name` | string | Internal tree node identifier for the parent node. May be a sample accession (for leaf nodes) or an internal node ID. |
| `child_name` | string | Internal tree node identifier for the child node. |
| `parent` | string | Parent nucleotide sequence, sliced to this protein's coding sequence. PCP extraction is fanned over the `pcp_genes` config list; each protein's CDS blocks come from the `pcp_gene_cds_blocks` config map if it names the gene, otherwise from the reference GTF (`gtf.resolve_cds_blocks`), so outputs live under a per-gene subdirectory (`pcps/{gene}/`). ORF1ab is reconstructed as one frameshifted CDS (its two blocks share nt 13468); `ORF1b` (v9+, nsp12–nsp16) is a config-supplied gene using the same two-block form, 8112 nt = 2704 codons. |
| `child` | string | Child nucleotide sequence after applying branch mutations. Same length as `parent`. |
| `branch_length` | integer | Number of nucleotide mutations on this branch (within the gene region). |

**Notes:**
- Sequences are already sliced to the gene region of interest, so they can be used directly for training.
- PCPs have been curated to remove branches with identical parent/child sequences (in the gene region). As of v2, this filter is applied during extraction in the same pass that builds the table; there is no separate `mat/clades/{clade}_pcps_raw.csv` file.
- The filtering criteria (max mutations, max reversions) are controlled by the shared Mutation-Counting config parameters `max_nt_mutations`, `max_reversions_to_ref`, and `max_reversions_to_clade_founder` (as of v8; PCP extraction and clade/window mutation counting filter branches identically).

### `pcps/{gene}/all_pcps_dasm_ready.csv.gz` and `pcps/{gene}/{clade}_pcps_dasm_ready.csv.gz` (v2+; gzipped in v9+)

DASM-training-ready sibling of `*_pcps.csv.gz`. Same row schema as the parent file plus one column (`in_train`). The filter set + train/val split are applied in the pipeline so downstream consumers (e.g. `data-central` manifests) can read these files directly without any post-processing step.

| Column | Type | Description |
|--------|------|-------------|
| (all columns from `*_pcps.csv.gz`) | — | Preserved verbatim. |
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
- **Aggregate (`all_pcps_dasm_ready.csv.gz`)**: `seed=42` (matches the explicit-seed guardrail in `filter_viral_pcps.sh`; the different seed prevents accidentally producing the aggregate with the per-clade default).

**Byte-equivalence to `dasm2 preprocess`**: this rule is a NumPy-only reimplementation of `dasm2.preprocess_command.preprocess` (which `filter_viral_pcps.sh` wraps). On the test pipeline, the aggregate file is byte-identical to the dasm2 reference output; md5 `5e538bd088e676b485a11a04745221e9` is asserted in `validate_test.py`. Since v9 our side is gzipped while the reference writes plain CSV, so that md5 is taken over the **decompressed** stream — the sentinel value is unchanged by compression. Writes pass `mtime=0`, so a re-run of the same file is byte-reproducible; note that gzip headers also embed the filename, so two identically-contented files under different names are not byte-comparable.

**Companion files**:

- `pcps/{gene}/{clade}_dasm_filter_report.txt` — human-readable per-filter drop counts.
- `pcps/{gene}/{clade}_dasm_filter_stats.json` — machine-readable equivalent for downstream tooling.

### `pcps/{gene}/{clade}_filtered_branches.csv`

Branch identifiers that were excluded during PCP curation. Deliberately **not** gzipped in v9 (#134): a single ID column of at most a few hundred MB, where staying readable with plain shell tools is worth more than the space.

| Column | Type | Description |
|--------|------|-------------|
| `child_name` | string | Tree node identifier for branches that failed curation filters. |

### `pcps/{gene}/filtered_branch_ids.csv`

Aggregated list of all filtered branch IDs across all clades. Not gzipped, as above.

> **Correction (v9):** earlier revisions of this document described the per-clade `{clade}_filtered_branches.csv` as "used internally to harmonize mutation counting with PCP extraction." That is stale. Nothing in the pipeline reads either file; harmonization comes from PCP extraction and mutation counting sharing the `max_nt_mutations` / `max_reversions_to_ref` / `max_reversions_to_clade_founder` config keys. Both files are retained for inspection only.

| Column | Type | Description |
|--------|------|-------------|
| `child_name` | string | Tree node identifier for filtered branches. |
| `clade` | string | Clade from which the branch was filtered. |

### `pcps/{gene}/censored_mutations.csv` and `pcps/{gene}/censored_mutations.clades.txt` (v11+)

Per-clade, per-nucleotide censor table for DASM training (#144), plus the clade roster it was built over. Additive: no existing output changed, and nothing in this repository reads it — the consumer is [matsengrp/dasm2-experiments#1079](https://github.com/matsengrp/dasm2-experiments/issues/1079).

UShER masking does not mark a site ambiguous; it **deletes** the mutations from the MAT, so the node inherits its parent's state. A masked site is therefore *censored, not ambiguous* — the state is known, the event was destroyed — which is why this is a table of coordinates for the loss rather than an `N` written into the PCP sequence.

| Column | Type | Description |
|--------|------|-------------|
| `clade` | string | Nextstrain clade. Matches the `clade` column of `pcps/{gene}/all_pcps_dasm_ready.csv.gz`; note that column exists **only** in the aggregate, so a per-clade consumer must key on the filename. |
| `nt_site` | int | 1-based genomic coordinate, as it appears in `data/usher_masked_sites.yaml`. |
| `nt_index` | int | 0-based index into this gene's PCP nucleotide sequence. |
| `from_nt` | string | Empty for a whole-site mask; the required parent base for a directional mask. |
| `to_nt` | string | Empty for a whole-site mask; the censored destination base for a directional mask. |

Sorted by `(clade, nt_index, to_nt)`. The sort is required, not cosmetic: `masking.load_usher_masks` returns unordered *sets*, and `data-central` manifests pin `md5`/`xxh3` per file.

**CSV contract.** A whole-site mask writes the empty string in both `from_nt` and `to_nt`; a row has both populated or neither. Read with `keep_default_na=False`, or pandas turns the empty cells into `NaN` and a `from_nt == ""` test silently matches nothing.

**`nt_index` is one-to-many.** Coordinates are gene-local because the consumer has no CDS knowledge — the PCP sequence *is* the gene. Where two CDS blocks share a nucleotide at a −1 ribosomal frameshift, `slice_cds` reads it twice, so one `nt_site` yields **two rows**. ORF1b's blocks `[[13442, 13468], [13468, 21552]]` put genomic 13468 at both index 26 and index 27, in different codons. A masked site outside the gene is dropped silently and counted in the per-gene log.

**This table does not replace `{gene}_codon_sites_to_mask.txt`.** It carries two of `masking.compute_exclusions`' five reasons — `usher_clade_sites` and `usher_clade_mutations`, the clade-specific ones. The ProblematicSites VCF, the config `sites_to_exclude` and the `site_include_range` bound remain in the codon mask at codon granularity, and founder→reference reversions are in neither (they are a selection filter, not an artifact mask). A training run must apply **both** files. The two are *not* disjoint by site: a site can be both clade-masked in the yaml and in the VCF (32 such `(clade, nt_site)` pairs on the v11 gisaid-and-public tree).

**The roster sidecar is part of the deliverable.** A clade with nothing censored contributes zero rows, so in the CSV alone it is indistinguishable from a mistyped clade name. `censored_mutations.clades.txt` lists every clade the run processed, one per line, sorted; a consumer materializes `{clade: [] for clade in roster}` before filling from the CSV so that absence-from-roster (an error) and absence-from-rows (normal) stay distinguishable.

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
| `actual_count` | float | Observed mutation count, summed over the nucleotide mutations giving this amino-acid change. Counts observations on the clade founder's codon background only (v9, #130). |
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
| `actual_count` | float | Observed mutation count, on the clade founder's codon background only (v9, #130). |
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
| `actual_count` | integer | Observed mutation count (0 if unobserved), counting **only** observations made on the clade founder's codon background. The same nucleotide substitution observed on an already-mutated codon is a different amino-acid change and is excluded (v9; see the note below). |
| `synonymous` | boolean | True if the mutation is synonymous. |
| `noncoding` | boolean | True if the site is non-coding. |

**Notes:**
- Additional annotation columns (`gene`, `clade_founder_nt`, `clade_founder_codon`,
  `mutant_codon`, `mutant_aa`, `aa_mutation`, `codon_position`, `codon_site`,
  `four_fold_degenerate`, `count_terminal`, `count_non_terminal`, `subset`, `mean_log_size`)
  are present for internal pipeline use.
- Replaces the former `expected_vs_actual_mut_counts.csv`; there is no `expected_count` here —
  the neutral expectation is estimated downstream by the GLM.
- **Exactly one row per `(clade, nt_mutation)`** — the table is the full enumeration of possible
  mutations per clade, so this is a structural invariant, checked at write time and by
  `validate_test.py`. Through v8 it did not hold: observations on non-founder codon backgrounds
  were relabelled onto the founder codon and kept as extra rows, so a mutation could appear
  2-5 times, with the copies byte-identical in every column except the counts. Consumers that
  grouped or joined on `(clade, nt_mutation)` silently double-counted those mutations, and both
  fitness estimators were affected (#130).

## Clade Founder Sequences

### `clade_founders/{clade}.fa`

Whole-genome founder sequence (FASTA) for each clade, computed by applying the clade's founder mutations (sourced from Neher and/or Roemer) to the Wuhan-Hu-1 reference (NC_045512.2). One sequence record per file; the header is the clade name (e.g. `>19A`). Length: 29903 nt.

### `clade_founders/{gene}/{clade}.fa`

Gene-sliced founder sequence for each clade and each protein in the `pcp_genes` config list. Contains just that protein's coding sequence, taken from that gene's CDS blocks (via `gtf.resolve_cds_blocks` — the `pcp_gene_cds_blocks` config map where it names the gene, otherwise the reference GTF) and excluding the trailing stop codon; the subdirectory is the gene symbol verbatim (e.g. `clade_founders/S/`, `clade_founders/ORF1ab/`, `clade_founders/E/`). ORF1ab is reconstructed as one in-frame CDS across the ribosomal frameshift — its two GTF blocks share nt 13468 — giving 21288 nt = 7096 codons; Spike is positions 21563-25381 (1-based, inclusive), 3819 nt = 1273 codons. Each slice byte-aligns with the `parent`/`child` sequences in the matching `pcps/{gene}/{clade}_pcps.csv.gz`, so it can be used directly as the wildtype for training ParentIndependent DASM models.

**`ORF1b` (v9+, #133).** A PCP-only gene, added through `pcp_gene_cds_blocks` rather than the reference GTF: blocks `[[13442, 13468], [13468, 21552]]`, 8112 nt = 2704 codons, covering nsp12–nsp16. It is exactly ORF1ab codons 4393–7096, so `clade_founders/ORF1b/{clade}.fa` is byte-identical to `clade_founders/ORF1ab/{clade}.fa[13176:21288]` for every clade — an invariant `validate_test.py` enforces. Note it spans the frameshift rather than starting after it: nsp12's first 9 residues are translated in ORF1a's frame, so the post-slippage block alone would clip them. ORF1b appears in **no** fitness output by design (`counts_all.csv`, the GLM, and `aamut_fitness` are unchanged); adding it to the GTF would have made an 8 kb stretch multi-gene and altered every fitness output.

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
- **384 rows, not 512.** `get_rates_table.py` drops the `nt_site_before_boundary
  == True` half of the grid, so consumers re-derive the boundary side per gene.
  The full grid is in `glm/master_table.csv` below.

### `rates_table_{era}.csv` (v10+)

One per era present on the tree — `rates_table_pre_omicron.csv` and/or
`rates_table_omicron.csv`. Identical schema and row count (384) to
`rates_table.csv`; only the training population differs. See
`glm/era_fit_summary.csv` for which eras a given run produced.

### `glm/master_table.csv`

The full GLM condition grid: **512 rows**, one per
`(mut_type, motif, unpaired, nt_site_before_boundary)` combination. Documented
here for the first time in v10 — it shipped from v1 with no schema entry.

| Column | Type | Description |
|--------|------|-------------|
| `mut_type` | string | Two-letter mutation type (e.g., `CT`, `GA`) — note this differs from `rates_table.csv`'s `C>T` rendering. |
| `motif` | string | Trinucleotide sequence context. |
| `unpaired` | int | 1 if the site is in an unpaired (single-stranded) RNA region. |
| `nt_site_boundary` | int | The configured light-switch boundary for this `mut_type`; 0 for the eight types that have none. |
| `nt_site_before_boundary` | boolean | Which side of the boundary this row serves. The boundary site itself belongs to the **before** bucket (`site <= boundary`, `modules/lightswitch.py`). |
| `rate` | float | Estimated neutral mutation rate for this condition. |
| `expected_count` | float | Predicted count for the training population. |
| `residual` | float | Empirical-Bayes-shrunk residual variance for this condition. |

**Row arithmetic:** 4 boundary types (`AT`, `CG`, `GC`, `CT`) × 16 motifs × 2
unpaired × 2 boundary sides = 256, plus 8 non-boundary types × 16 motifs × 2
unpaired × 1 = 256.

### `glm/master_table_{era}.csv` (v10+)

Same 512-row schema and grid as `glm/master_table.csv`, fitted on one era's
clades only. Both eras' tables share the grid exactly, so they are comparable
cell by cell.

**Three tables, not one with an `era` column.** DASM2's
`TabularRatesNeutralModel` assigns `wt_to_mut_rates[wt][mut]` in a row loop, so
an unrecognized extra column would be read *last-row-wins, silently*. The pooled
`master_table.csv` is retained with its v9 meaning as the control arm of the
per-era comparison, not as backwards compatibility.

### `glm/era_fit_summary.csv` (v10+)

Which eras this tree's run fitted, and which it skipped.

| Column | Type | Description |
|--------|------|-------------|
| `era` | string | `pre_omicron` or `omicron`. **Every** era appears, fitted or not. |
| `fitted` | boolean | Whether a rate table was produced for this era. |
| `n_clades` | int | Clades in this era on this tree. |
| `clades` | string | `;`-joined clade labels, in Nextstrain order. |
| `master_table` | string | Relative path to the era's master table, or empty if skipped. |
| `rates_table` | string | Relative path to the era's rates table, or empty if skipped. |

**Why this file exists.** An era with no clades on a tree is skipped — the GLM
cannot fit an empty partition, and single-era trees are normal (the test config's
three clades are all pre-Omicron; the BA.2/JN.1/XFG trees are single-era by
construction). But a skipped era leaves *no file*, which is indistinguishable
from an era that went missing through a wiring mistake. This records the absence
explicitly so it is auditable from the output tree alone.

On a single-era tree the era's tables are **byte-identical** to the pooled ones,
which `validate_test.py::validate_era_rate_tables` asserts.

### `windows/` expected counts use the pooled table

`window_expected_counts` reads `glm/master_table.csv`, **not** a per-era table.
The windowed path has no clade column to dispatch on (`subtree_counts.csv` and
`founder_muts.csv` key on `founder_id`, a bare node id), a per-window date label
is wrong exactly across the transition, and a per-founder label breaks the
shared-base-enumeration optimization. Since `T_window` divides out any constant
rate factor, only the split's *relative* structure leaks through — dominated by
G→T. Tracked as a follow-up.

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
| `check` | string | See the two lists below. |
| `status` | string | `PASS`, `FAIL`, or `REPORT`. |
| `detail` | string | Measured numbers for the check. |

**Gating checks** (`PASS`/`FAIL`): `node_id_correspondence`, `time_monotonicity`,
`date_coverage`, `tip_residual`. The rule raises on any hard-invariant violation, and
on a `tip_residual` median above `chronumental_max_median_abs_tip_residual_days`.

**Report-only checks** (`REPORT`, added by
[#139](https://github.com/matsengrp/mat-pcp/issues/139)) — these validate the
*internal* node dates, which `modules/windowing.py` uses to bucket every branch into
a time window. They **gate nothing**: `status` is always `REPORT` and they cannot
fail a run. Thresholds are deferred to a follow-up, because gating on a distribution
nobody has characterised would fail production runs uninterpretably.

| `check` | measures |
|---|---|
| `branch_dt_vs_mutations` | Δt binned by branch mutation count, plus days/mutation against a nominal ~15.3 (order-of-magnitude only: the MAT is masked and branch-filtered, so the realized slope runs higher). Also reports the **aggregate** rate — Σ mutations / Σ branch-years, exact over the per-branch arrays rather than reconstructed from the ranged bins below — as a ratio to that nominal and to chronumental's own pre-fit root-to-tip estimate parsed from `log/chronumental.stdout`, plus the zero-mutation share of total branch-*time*. The days/mutation figure is the **incremental** slope between adjacent single-count bins, which cancels the zero-mutation intercept by construction and is therefore blind to the duration inflation confirmed in [#145](https://github.com/matsengrp/mat-pcp/issues/145); the aggregate is the statistic that sees it ([#146](https://github.com/matsengrp/mat-pcp/issues/146)). When the stdout log is absent (it is a Snakemake `log:` entry, not a declared output) the detail says so and names the path it looked at, rather than silently omitting the comparison. |
| `zero_mut_leaf_predicted_dt` | Predicted Δt on zero-mutation leaf branches. **One-sided**: `time_monotonicity` is a hard gate already forcing Δt ≥ 0, so `frac<0` here is an *invariant rather than a measurement* and must not be quoted as a confirmed prediction. Read `min` instead: `EPS_DAYS = 1` admits a −1 day edge past that gate, and a 4-decimal fraction over 12 M observations hides up to ~600 negatives. |
| `zero_mut_leaf_vs_earliest_observed` | *Earliest observed* date among a parent's zero-mutation leaf children **minus** the parent's predicted date. Note this is the **opposite orientation** to `internal_date_vs_min_descendant_observed`: here a **negative** value is the anomaly. Also reports how many dated strains had no node in the tree and were excluded. |
| `root_date_plausible` | Root date against a deliberately wide window (2019-10-01 … 2020-02-29). A tripwire for a mis-scaled time axis. |
| `internal_date_vs_min_descendant_observed` | Predicted internal date minus the earliest observed descendant date. **Positive is the anomaly** (the ancestor postdates its earliest sample). Compares *inferred* against *observed*, so a violation is a three-way inconsistency between the inferred date, the topology and a collection date — and it is not independent of `tip_residual`: monotonicity forces `pred(ancestor) ≤ pred(tip)`, so a positive gap requires that tip's own residual to be at least as large. Quantifies downstream impact (how many internal nodes, hence branches, land in the wrong window) rather than testing the dating algorithm. |
| `internal_dt_to_closest_descendant` | Predicted Δt from an internal node to its **nearest descendant, both dates inferred** — so unlike the row above this *is* a test of the algorithm: given the topology it was handed, how much time does the fit put below each internal node? Equals min-over-*children* of the branch Δt, because monotonicity makes predicted dates nondecreasing along every path, so the nearest descendant is always an immediate child. Also reports how many zero-mutation branches are internal rather than terminal (0 on the v9 tree, so `zero_mut_leaf_predicted_dt` covers all of them). |

Because checks 6 and 8 use opposite sign conventions, `format_stats` prints **both**
`frac<0` and `frac>0`, and the two rows above additionally carry an explicit
`VIOLATING (<0)` / `VIOLATING (>0)` fraction and count. `frac<0 + frac>0` does not
sum to 1 — exact zeros are in neither. `tip_residual`'s detail also now reports its
**positive tail** (count, count beyond 30 d, and max), because the gate is a median
over ~17 M tips and cannot see a tail of any size; that tail is what check 8
re-expresses.

Why report-only rather than gating: `tip_residual` is the only pre-#139 check that
compares an inferred date to evidence, and it is aimed where there is least to
infer — tip collection dates are the priors chronumental is fit to. On
`gisaid-and-public_2026-04-22_v9` it reports `p95|res|=1.0d` over 16,979,190 tips
against a 30-day gate.

### `mat/node_dates_branch_dt.csv`

Long-form Δt distribution per mutation-count bin (#139). Binned rather than
per-branch: `node_dates_residuals.csv` is already 1.6 GB on the production tree, and
Δt in whole days is an integer, so a histogram gives **exact** medians and
percentiles at a few tens of kB.

| Column | Type | Description |
|--------|------|-------------|
| `mut_count_bin` | string | Mutation-count bin (`0`, `1`, …, `6-10`, `11-20`, `21-50`, `51+`). |
| `dt_days` | int | `date(child) − date(parent)` in whole days. |
| `count` | int | Number of branches in that bin with that Δt. |

### `mat/node_dates_branch_dt_summary.csv`

Wide-form per-bin summary of the same distributions — median, interquartile range and
p95 for branch Δt binned by branch length in mutations (#139). One row per
mutation-count bin plus a final `all` row.

Built from the *same* `DayHistogram` objects as `node_dates_branch_dt.csv`, so the two
files cannot disagree, and every statistic is **exact** rather than interpolated: Δt is
whole days and the histogram retains every observation. `iqr` is `p75 − p25` on those
same order statistics, not a separate estimate. The `all` row is the marginal of the
pooled data (via `DayHistogram.merge`), **not** an average of per-bin medians; the
script raises if its `n_branches` does not equal the number of branches with a defined
Δt, so a bin set that stopped partitioning them would fail rather than under-report.

| Column | Type | Description |
|--------|------|-------------|
| `mut_count_bin` | string | Mutation-count bin (`0`, `1`, …, `6-10`, `11-20`, `21-50`, `51+`), or `all` for the pooled row. |
| `n_branches` | int | Branches in the bin with a defined Δt. |
| `min`, `p05`, `p25`, `median`, `p75`, `p95`, `max` | int | Exact order statistics of Δt in whole days. |
| `iqr` | int | `p75 − p25`. |
| `mean` | float | Mean Δt in days, to 2 dp. |
| `frac_dt_zero` | float | Fraction of the bin with Δt exactly 0 — parent and child dated the same day. |

### `mat/node_dates_internal_checks.csv`

Long-form day distributions for the three internal-node checks that produce one
(#139), so a follow-up can choose thresholds from a distribution rather than a
summary line.

| Column | Type | Description |
|--------|------|-------------|
| `check` | string | `zero_mut_leaf_predicted_dt`, `zero_mut_leaf_vs_earliest_observed`, `internal_date_vs_min_descendant_observed`, or `internal_dt_to_closest_descendant`. |
| `days` | int | Signed day difference. Negative values are retained — for the latter two they are the informative direction. |
| `count` | int | Number of nodes with that value. |

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
