# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`mat-pcp` is a Snakemake-based bioinformatics pipeline for analyzing SARS-CoV-2 genomic evolution by extracting parent-child pairs (PCPs) from mutation-annotated trees (MATs). The pipeline processes pre-built MAT trees, filters sequences based on quality criteria, optionally optimizes trees for maximum parsimony, and extracts PCPs for evolutionary analysis.

## Development Commands

### Running the Pipeline
```bash
# Run with 4 cores (adjust -j based on available cores)
snakemake -j 4 --use-conda

# Dry run to see what would be executed
snakemake -n

# Run specific target
snakemake -j 4 --use-conda results_June-01-2022/mat/all_pcps.csv

# Clean all outputs
snakemake --delete-all-output
```

### Development and Testing
```bash
# Check workflow syntax
snakemake --lint

# Generate workflow visualization
snakemake --dag | dot -Tpdf > workflow.pdf

# Run with detailed output
snakemake -j 4 --use-conda --verbose
```

## Architecture and Key Components

### Pipeline Structure
- **Configuration**: `config/config.yaml` controls all pipeline parameters
- **Workflow**: `workflow/Snakefile` defines the processing rules and dependencies
- **Environments**: `workflow/envs/` contains conda environment specifications for reproducibility
- **Scripts**: `workflow/scripts/` contains Python scripts for specific processing tasks
- **Notebooks**: `workflow/notebooks/filter_metadata.ipynb` for interactive metadata filtering

### Key Processing Steps
1. **get_mat_tree**: Links input MAT trees and metadata
2. **filter_metadata**: QC-based filtering using Jupyter notebook (invalid dates, hosts, genome length)
3. **filter_mat_tree**: Extracts valid sequences using matUtils (double extraction for reliability)
4. **optimize_mat_tree**: Optional optimization for maximum parsimony using matOptimize
5. **get_ref_fasta/gtf**: Downloads reference genome and annotations
6. **edit_ref_gtf**: Adds missing gene annotations (e.g., ORF9b)
7. **ref_coding_sites**: Identifies coding regions in the genome
8. **extract_pcps**: Extracts parent-child pairs using BTE

### Configuration Options
Key parameters in `config/config.yaml`:
- `mat_trees`: Input trees and metadata paths
- `qc_filter_options`: Date validity, host species, genome length filters
- `optimize_tree`: Toggle tree optimization (default: false)
- `subset_frac`: Random sampling fraction for metadata (default: 1.0)
- `datetime_cutoff`: Date cutoff for sample inclusion
- `threads`: Per-rule thread allocation

### Output Organization
Results are stored in `results_{mat}/` directories:
- `mat/`: Filtered trees and extracted PCPs (`all_pcps.csv` is the main output)
- `ref/`: Reference genome files and coding site information
- `log/`: Processing logs for each step

### Important Implementation Details
- Uses UCSC genome browser for SARS-CoV-2 reference data
- Handles compressed files (.gz) throughout
- Implements double extraction in `filter_mat_tree` to address a specific issue
- Tree optimization can be toggled via `optimize_tree` config parameter
- Thread allocation is configurable per rule with fallback to default