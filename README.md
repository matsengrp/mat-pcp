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

5. **PCP Extraction**:
   - Extracts parent-child pairs from the optimized tree
   - Outputs PCPs for downstream evolutionary analysis

## Configuration

The pipeline is configured through `config/config.yaml`, which allows you to specify:

- Input MAT trees and metadata
- QC filter options (date resolution, host filtering, genome length)
- Reference genome and annotation sources
- Datetime cutoff for sample inclusion
- Sampling fraction for metadata subsetting
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
- `ref/`: Reference genome files and coding site information
- `log/`: Log files for each step of the pipeline

The main output is `mat/all_pcps.csv`, containing the extracted parent-child pairs for evolutionary analysis.

**NOTE**: Currently, this pipeline is focused on SARS-CoV-2 data but can be adapted for other organisms by modifying the configuration file.


