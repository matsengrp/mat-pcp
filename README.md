# `mat-pcp`

A unified Snakemake pipeline for analyzing SARS-CoV-2 genomic evolution from mutation-annotated trees (MATs). This pipeline extracts parent-child pairs (PCPs), counts mutations, estimates mutation rates, and calculates fitness effects.

## Pipeline Overview

The pipeline performs the following key steps:

1. **Data Acquisition**:
   - Obtains pre-built mutation-annotated trees and corresponding metadata
   - Downloads reference genome (FASTA) and annotation (GTF) files

2. **Metadata Filtering**:
   - Filters sequences based on quality control parameters (valid dates, hosts, genome length)
   - Allows subsetting data by a specified datetime cutoff
   - Supports random sampling of sequences with adjustable sampling fraction

3. **Tree Processing**:
   - Extracts valid sequences from the MAT tree based on filtered metadata
   - Optionally optimizes the filtered tree for maximum parsimony using `matOptimize`

4. **Reference Genome Processing**:
   - Processes reference genome and annotation files
   - Identifies coding sites from the reference genome
   - Supports manual addition of gene annotations (e.g., ORF9b)

5. **Clade-level PCP Extraction**:
   - Identifies clades with adequate sample counts (configurable minimum)
   - Splits the tree by clade for independent analysis
   - Extracts parent-child pairs from each clade-specific tree
   - Aggregates all clade PCPs into a single output file with clade annotations

## Configuration

The pipeline is configured through `config/config.yaml`, which contains comprehensive parameters organized into sections:

### Input/Output
- `output_dir`: Directory for results
- `mat`: Path to mutation-annotated tree (.pb.gz)
- `metadata`: Path to metadata file (.tsv.gz)

### QC Filtering (`qc_filter_options`)
- Date resolution filtering (exclude MISSING, MONTH, YEAR resolutions)
- Date range filtering (`min_date`, `datetime_cutoff`)
- Host species validation
- Genome length constraints

### Tree Processing
- `optimize_tree`: Toggle tree optimization (default: false)
- `mat_subset_kwargs`: Additional matUtils extract arguments

### Clade Configuration
- `clade_column`: Metadata column for clade assignment
- `min_clade_size`: Minimum samples per clade
- `clades_to_exclude`: List of clades to skip

### Data Files
The pipeline uses reference data files in `data/`:
- `usher_masked_sites.yaml`: UShER masking information for problematic sites
- `lan_2022/`: RNA secondary structure data from Lan et al. 2022

## Running the Pipeline

You must have an environment with `snakemake` and `conda` available. We recommend using `pixi`:

```bash
# Install dependencies
pixi install

# Run pipeline with 4 cores
pixi run snakemake -j 4 --use-conda

# Run with custom input files
pixi run snakemake -j 4 --use-conda \
  --config output_dir="results" mat="path/to/tree.pb.gz" metadata="path/to/metadata.tsv.gz"
```

## Output Structure

Results are organized in the output directory:

- `mat/`: Filtered trees and extracted PCP data
  - `clades/`: Clade-specific trees and PCP files
  - `all_pcps.csv`: Aggregated PCPs from all clades with clade annotations
  - `clade_counts.csv`: Summary of samples per clade
- `ref/`: Reference genome files and coding site information
- `log/`: Log files for each step of the pipeline

The main output is `mat/all_pcps.csv`, containing the extracted parent-child pairs from all clades with clade annotations for evolutionary analysis.

## Testing and CI

The pipeline includes automated testing via GitHub Actions that runs on every push and pull request. The CI workflow:

1. Validates the configuration schema
2. Verifies required data files exist
3. Runs the pipeline on test data in `test_data/`
4. Validates outputs against expected results

### Running Tests Locally

```bash
# Install dependencies
pixi install

# Run test pipeline
pixi run snakemake -j 2 --use-conda \
  --config output_dir="results_test_data" mat="test_data/filtered_tree_raw.pb.gz" metadata="test_data/filtered_metadata.tsv.gz" subset_frac=0.1

# Validate results
pixi run python workflow/scripts/validate_test.py
```

The validation checks:
- Total number of PCPs extracted
- Number of clades processed
- Total sequences across all clades
- Presence of expected output files
