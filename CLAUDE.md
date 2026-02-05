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
| 1-2 | Input Processing & Reference | `filter_metadata`, `filter_mat_tree`, `get_ref_fasta`, `ref_coding_sites` |
| 3 | Clade Processing & PCPs | `get_clades` (checkpoint), `mat_clade_subset`, `extract_pcps_by_clade` |
| 4 | Mutation Counting | `count_mutations_clade`, `aggregate_mutation_counts` |
| 5 | Expected Counts | `synonymous_mut_rates`, `expected_mut_counts`, `merge_expected_and_actual_counts` |
| 6 | GLM Rate Estimation | `annotate_counts`, `curate_counts`, `master_table`, `get_rates_table` |
| 7 | Fitness Estimation | `ntmut_fitness_by_clade`, `aamut_fitness_by_clade`, `aggregate_*_fitness` |
| 8 | Validation | `extract_mutations_from_pcps`, `validate_pcp_vs_counts` |

## Development Commands

### Running the Pipeline
```bash
# Run with pixi (recommended)
pixi run snakemake -j 4 --use-conda

# Dry run
pixi run snakemake -n

# Run specific target
pixi run snakemake -j 4 --use-conda results/aamut_fitness/all_aamut_fitness.csv

# Clean all outputs
pixi run snakemake --delete-all-output
```

### Testing
```bash
# Run on test data (fast, uses max_clades=3)
pixi run snakemake -j 2 --use-conda \
  --configfile config/config_test_data.yaml

# Validate results
pixi run python workflow/scripts/validate_test.py
```

## Architecture

### Directory Structure
```
mat-pcp/
├── config/
│   ├── config.yaml                  # Default pipeline configuration
│   ├── config_test_data.yaml        # Test data configuration
│   └── config-public-2022-06-01.yaml # Production configuration
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

### Configuration Sections

The configuration (`config/config.yaml`) is organized into:

1. **Input/Output**: `output_dir`, `mat`, `metadata`
2. **Reference Genome**: `ref_fasta`, `ref_gtf`, `add_to_ref_gtf`
3. **Metadata QC**: `qc_filter_options` (date filtering, host validation, genome length)
4. **Tree Processing**: `optimize_tree`, `mat_subset_kwargs`
5. **Clade Configuration**: `clade_column`, `min_clade_size`, `clades_to_exclude`, `max_clades`
6. **PCP Extraction**: `pcp_gene_region`, `pcp_max_nt_mutations`, `pcp_max_reversions_to_ref`
7. **Clade Founders**: `clade_founders` (Neher and Roemer JSON sources)
8. **Mutation Counting**: `max_nt_mutations`, `max_reversions_*`, `site_include_range`
9. **Synonymous Rates**: `synonymous_spectra_min_counts`
10. **GLM**: `rna_structure_data`, `glm_lightswitch_boundaries`
11. **Fitness**: `fitness_pseudocount`, `fitness_N_f`, `fitness_N_n_exp`
12. **Gene Configuration**: `orf1ab_to_nsps`, `gene_overlaps`, `genes`
13. **Site Masking**: `masking_gene`, `sites_to_exclude`, `additional_sites_to_ignore`
14. **Threads**: Per-rule thread allocation

### Output Organization

Results are stored in `{output_dir}/`:
- `mat/`: Filtered trees, clade subtrees
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
