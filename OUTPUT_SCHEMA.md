# Output Schema Reference

This document describes the column schemas for the primary consumer-facing output files produced by the mat-pcp pipeline. These schemas are fixed and do not vary based on configuration.

For information on running the pipeline, see [README.md](README.md).

## Parent-Child Pairs (PCPs)

### `pcps/all_pcps.csv` and `pcps/{clade}_pcps.csv`

Curated parent-child sequence pairs extracted from the mutation-annotated tree. These represent evolutionary transitions along tree branches and are suitable for training sequence models (e.g., DASM).

| Column | Type | Description |
|--------|------|-------------|
| `clade` | string | Nextstrain clade identifier (e.g., "19A", "20B", "21K"). Only present in `all_pcps.csv`. |
| `parent_name` | string | Internal tree node identifier for the parent node. May be a sample accession (for leaf nodes) or an internal node ID. |
| `child_name` | string | Internal tree node identifier for the child node. |
| `parent` | string | Parent nucleotide sequence. Sliced to the genomic region specified by `pcp_gene_region` in config (typically Spike: positions 21563-25384). |
| `child` | string | Child nucleotide sequence after applying branch mutations. Same length as `parent`. |
| `branch_length` | integer | Number of nucleotide mutations on this branch (within the gene region). |

**Notes:**
- Sequences are already sliced to the gene region of interest, so they can be used directly for training.
- PCPs have been curated to remove branches with identical parent/child sequences (in the gene region). As of v2, this filter is applied during extraction in the same pass that builds the table; there is no separate `mat/clades/{clade}_pcps_raw.csv` file.
- The filtering criteria (max mutations, max reversions) are controlled by config parameters `pcp_max_nt_mutations`, `pcp_max_reversions_to_ref`, and `pcp_max_reversions_to_clade_founder`.

### `pcps/all_pcps_dasm_ready.csv` and `pcps/{clade}_pcps_dasm_ready.csv` (v2+)

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
6. Every codon-level mutation falls in a codon position listed in `spike_codon_sites_to_mask.txt` (no-signal PCP)

**Seed convention**:

- **Per-clade**: `seed=1729` (matches `data-central/scripts/filter_viral_pcps.sh`).
- **Aggregate (`all_pcps_dasm_ready.csv`)**: `seed=42` (matches the explicit-seed guardrail in `filter_viral_pcps.sh`; the different seed prevents accidentally producing the aggregate with the per-clade default).

**Byte-equivalence to `dasm2 preprocess`**: this rule is a NumPy-only reimplementation of `dasm2.preprocess_command.preprocess` (which `filter_viral_pcps.sh` wraps). On the test pipeline, the aggregate file is byte-identical to the dasm2 reference output; md5 `5e538bd088e676b485a11a04745221e9` is asserted in `validate_test.py`.

**Companion files**:

- `pcps/{clade}_dasm_filter_report.txt` — human-readable per-filter drop counts.
- `pcps/{clade}_dasm_filter_stats.json` — machine-readable equivalent for downstream tooling.

### `pcps/{clade}_filtered_branches.csv`

Branch identifiers that were excluded during PCP curation. Used internally to harmonize mutation counting with PCP extraction.

| Column | Type | Description |
|--------|------|-------------|
| `child_name` | string | Tree node identifier for branches that failed curation filters. |

### `pcps/filtered_branch_ids.csv`

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
| `expected_count` | float | GLM-predicted neutral count (sequence context + RNA structure + position); the baseline for both `naive_delta_fitness` and `delta_fitness`. |
| `actual_count` | float | Observed mutation count in the data. |
| `tau_squared` | float | Variance component from the hierarchical model. |
| `naive_delta_fitness` | float | Simple log-ratio fitness: `log((actual + 0.5) / (expected + 0.5))`. |
| `delta_fitness` | float | Bayesian fitness estimate accounting for uncertainty. |
| `uncertainty` | float | Standard error of the `delta_fitness` estimate. |
| `source_clade` | string | Semicolon-separated list of clades that contributed to this aggregated estimate. |

**Notes:**
- Positive `delta_fitness` indicates the mutation is observed more often than expected (potentially beneficial).
- Negative `delta_fitness` indicates the mutation is observed less often than expected (potentially deleterious).
- The `uncertainty` column can be used to filter for high-confidence estimates.

### `aamut_fitness/{clade}_aamut_fitness.csv`

Per-clade amino acid fitness estimates. Same schema as above but with `cluster` instead of aggregated fields:

| Column | Type | Description |
|--------|------|-------------|
| `cluster` | string | Clade name (e.g., "19A"). Note: named "cluster" for historical compatibility. |

(Other columns same as `all_aamut_fitness.csv`, excluding `source_clade`.)

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
| `exclude` | boolean | True if the mutation is masked for any reason (config sites, ProblematicSites VCF, UShER yaml, out of `site_include_range`, or a ref→founder clade-defining mutation). |
| `masked_in_usher` | boolean | True if excluded specifically by an artifactual-site mask (UShER yaml or the ProblematicSites VCF). |
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

### `clade_founders/{clade}_ref_to_founder_muts.csv`

Mutation list used to produce the founder FASTA, one row per substitution.

| Column | Type | Description |
|--------|------|-------------|
| `mutation` | string | Mutation string in standard format (e.g., `A23403G`). |
| `site` | integer | 1-indexed genomic position. |
| `from_nt` | string | Reference nucleotide. |
| `to_nt` | string | Founder nucleotide. |

### `clade_founders/spike/{clade}.fa` (v2+)

Gene-sliced founder sequence for each clade, produced when `pcp_gene_region` is set in the config. For Spike (the v2 default), this contains positions 21563-25381 (1-based, inclusive) of the founder genome — 3819 nt = 1273 codons, excluding the trailing stop codon. The slice is byte-aligned with the `parent`/`child` sequences in `pcps/{clade}_pcps.csv`, so it can be used directly as the wildtype for training ParentIndependent DASM models.

**Cross-version integrity sentinel**: the 19A spike slice (Wuhan-Hu-1 reference founder) has md5 `b948f2a1f20edbce64cfff1dc540c119` and is byte-stable across all SARS-CoV-2 cutoff dates.

**Notes:**
- The subdirectory name (`spike/`) is derived from the gene symbol in `pcp_gene_region.gene` via the `GENE_FOUNDER_DIRNAMES` mapping in the Snakefile (`S` → `spike`); other gene symbols fall back to a lowercased name.
- Multi-gene slicing is not yet supported (tracked for v3).

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
