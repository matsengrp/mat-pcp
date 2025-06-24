# `mat-pcp`

A Snakemake pipeline for filtering, optimizing, and extracting parent-child pairs (PCPs) from a mutation-annotated tree (MAT). This pipeline processes SARS-CoV-2 genomic data to facilitate evolutionary analysis.

## Pipeline Overview

The pipeline performs the following key steps:

1. **Data Acquisition**: 
   - Obtains pre-built mutation-annotated trees and corresponding metadata
   - Downloads reference genome (FASTA) and annotation (GTF) files

2. **Metadata Filtering**:
   - Filters sequences based on quality control parameters (valid dates, hosts, genome length)
   - Allows subsetting data by a specified datetime cutoff (configurable in `config.yaml`)
   - Supports random sampling of sequences with adjustable sampling fraction

3. **Tree Processing**:
   - Extracts valid sequences from the MAT tree based on filtered metadata
   - Optimizes the filtered tree for maximum parsimony using `matOptimize`

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

The pipeline is configured through `config/config.yaml`, which allows you to specify:

- Input MAT trees and metadata
- QC filter options (date resolution, host filtering, genome length)
- Reference genome and annotation sources
- Datetime cutoff for sample inclusion
- Sampling fraction for metadata subsetting
- Clade extraction parameters (minimum clade size, clades to exclude)
- Thread allocation for parallel processing

## Running the pipeline 

You must have an environment with `snakemake` and `conda` available. The pipeline uses conda environments defined in the `workflow/envs/` directory for each step, ensuring reproducibility.

```bash
snakemake -j 4 --use-conda
```

Where `-j 4` specifies using 4 cores for parallel execution. Adjust based on your system capabilities.

## Output Structure

Results are organized by input tree in directories named `results_{mat}`, where `{mat}` is the tree identifier. Each results directory contains:

- `mat/`: Filtered trees and extracted PCP data
  - `clades/`: Clade-specific trees and PCP files
  - `all_pcps.csv`: Aggregated PCPs from all clades with clade annotations
  - `clade_counts.csv`: Summary of samples per clade
- `ref/`: Reference genome files and coding site information
- `log/`: Log files for each step of the pipeline

The main output is `mat/all_pcps.csv`, containing the extracted parent-child pairs from all clades with clade annotations for evolutionary analysis.

**NOTE**: Currently, this pipeline is focused on SARS-CoV-2 data but can be adapted for other organisms by modifying the configuration file.

## Testing and CI

The pipeline includes automated testing via GitHub Actions that runs on every push and pull request. The CI workflow:

1. Sets up the environment using pixi (providing snakemake and conda)
2. Runs the pipeline on test data in `test_data/`
3. Validates outputs against expected results
4. Reports any differences in pipeline behavior

### Running Tests Locally

To run the same tests that CI runs:

```bash
# Install dependencies
pixi install

# Run test pipeline
pixi run snakemake -j 2 --use-conda \
  --config mat_trees.test_data.mat=test_data/filtered_tree_raw.pb.gz \
           mat_trees.test_data.metadata=test_data/filtered_metadata.tsv.gz

# Validate results
python workflow/scripts/validate_test.py
```

The validation checks:
- Total number of PCPs extracted
- Number of clades processed
- Total sequences across all clades
- Presence of expected output files
- Checksums of key output files


