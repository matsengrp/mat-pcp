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
| `expected_count` | float | Expected mutation count based on synonymous mutation rates. |
| `predicted_count` | float | GLM-predicted count incorporating sequence context and RNA structure. |
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
| `expected_count` | float | Expected count from synonymous rates (pre-GLM). |
| `predicted_count` | float | GLM-predicted count incorporating sequence context. |
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

Raw mutation counts aggregated across all clades, with annotations.

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
| `exclude` | boolean | True if this mutation should be excluded from fitness analysis. |
| `mean_log_size` | float | Mean log of subtree sizes for branches with this mutation. |

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
