# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`mat-pcp` is a unified Snakemake-based bioinformatics pipeline for analyzing SARS-CoV-2 genomic evolution. The pipeline extracts parent-child pairs (PCPs) from mutation-annotated trees (MATs), counts mutations, estimates mutation rates, and calculates fitness effects. It processes pre-built MAT trees, filters sequences based on quality criteria, optionally optimizes trees for maximum parsimony, and performs comprehensive evolutionary analysis.

## Development Commands

### Running the Pipeline
```bash
# Run with 4 cores (adjust -j based on available cores)
pixi run snakemake -j 4 --use-conda

# Dry run to see what would be executed
pixi run snakemake -n

# Run specific target
pixi run snakemake -j 4 --use-conda results/mat/all_pcps.csv

# Clean all outputs
pixi run snakemake --delete-all-output

# Run with custom config overrides
pixi run snakemake -j 4 --use-conda \
  --config output_dir="results_test" mat="test_data/tree.pb.gz" metadata="test_data/metadata.tsv.gz"
```

### Development and Testing
```bash
# Check workflow syntax
pixi run snakemake --lint

# Generate workflow visualization
pixi run snakemake --dag | dot -Tpdf > workflow.pdf

# Run with detailed output
pixi run snakemake -j 4 --use-conda --verbose

# Run test pipeline
pixi run snakemake -j 2 --use-conda \
  --config output_dir="results_test_data" mat="test_data/filtered_tree_raw.pb.gz" metadata="test_data/filtered_metadata.tsv.gz"
```

## Architecture and Key Components

### Pipeline Structure
- **Configuration**: `config/config.yaml` - comprehensive configuration with all pipeline parameters
- **Workflow**: `workflow/Snakefile` - defines processing rules and dependencies
- **Environments**: `workflow/envs/` - conda environment specifications for reproducibility
- **Scripts**: `workflow/scripts/` - Python scripts for specific processing tasks
- **Data**: `data/` - reference data files (masking sites, RNA structure data)

### Configuration Sections

The unified configuration (`config/config.yaml`) is organized into these sections:

1. **Input/Output**: `output_dir`, `mat`, `metadata`
2. **Reference Genome**: `ref_fasta`, `ref_gtf`, `add_to_ref_gtf`
3. **Metadata QC Filtering**: `qc_filter_options` (date filtering, host validation, genome length)
4. **Tree Processing**: `optimize_tree`, `mat_subset_kwargs`
5. **Clade Configuration**: `clade_column`, `min_clade_size`, `clades_to_exclude`
6. **PCP Extraction**: `pcp_gene_region`, `pcp_max_nt_mutations`, `pcp_max_reversions_to_ref`
7. **Clade Founders**: `clade_founders` (Neher and Roemer JSON sources)
8. **Mutation Counting**: `max_nt_mutations`, `max_reversions_*`, `site_include_range`, `sites_to_exclude`
9. **Synonymous Rate Estimation**: `synonymous_spectra_min_counts`
10. **GLM Rate Estimation**: `rna_structure_data`, `glm_lightswitch_boundaries`
11. **Fitness Estimation**: `fitness_pseudocount`, `fitness_N_f`, `fitness_N_n_exp`
12. **Gene Configuration**: `orf1ab_to_nsps`, `gene_overlaps`, `genes`
13. **Site Masking**: `masking_gene`, `additional_sites_to_ignore`
14. **Thread Configuration**: Per-rule thread allocation

### Key Processing Steps
1. **get_mat_tree**: Links input MAT trees and metadata
2. **filter_metadata**: QC-based filtering (invalid dates, hosts, genome length)
3. **filter_mat_tree**: Extracts valid sequences using matUtils
4. **optimize_mat_tree**: Optional optimization for maximum parsimony using matOptimize
5. **get_ref_fasta/gtf**: Downloads reference genome and annotations
6. **edit_ref_gtf**: Adds missing gene annotations (e.g., ORF9b)
7. **ref_coding_sites**: Identifies coding regions in the genome
8. **get_clades**: Checkpoint to identify clades with adequate sample counts
9. **samples_by_clade**: Extract sample names per clade
10. **mat_clade_subset**: Extract clade-specific MAT subtrees
11. **extract_pcps_by_clade**: Extract parent-child pairs per clade using BTE
12. **aggregate_clade_pcps**: Combine all clade PCPs into single output

### Data Files
- `data/usher_masked_sites.yaml`: UShER masking information for problematic sites by clade
- `data/lan_2022/41467_2022_28603_MOESM11_ESM.txt`: RNA secondary structure data from Lan et al. 2022

### Output Organization
Results are stored in `{output_dir}/` directories:
- `mat/`: Filtered trees and extracted PCPs
  - `clades/`: Clade-specific trees and PCP files
  - `all_pcps.csv`: Aggregated PCPs from all clades with clade annotations
  - `clade_counts.csv`: Summary of samples per clade
- `ref/`: Reference genome files and coding site information
- `log/`: Processing logs for each step

### Important Implementation Details
- Uses UCSC genome browser for SARS-CoV-2 reference data
- Handles compressed files (.gz) throughout
- Implements double extraction in `filter_mat_tree` to address matUtils behavior
- Tree optimization can be toggled via `optimize_tree` config parameter
- Thread allocation is configurable per rule with fallback to default
- Clade-level processing uses Snakemake checkpoints for dynamic workflow
