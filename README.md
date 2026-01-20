# mat-pcp: Unified MAT Analysis Pipeline

A Snakemake pipeline for analyzing SARS-CoV-2 genomic evolution from mutation-annotated trees (MATs). This pipeline extracts parent-child pairs (PCPs), counts mutations, estimates neutral mutation rates, and calculates fitness effects.

## Overview

This pipeline combines functionality from three previously separate projects:
1. **mat-pcp** - Parent-child pair extraction from MATs
2. **SARS2-mut-fitness** - Expected vs actual mutation count computation
3. **SARS2-mut-fitness-v2** - Neutral mutation rate estimation via GLM

## Quick Start

### Prerequisites

- [pixi](https://pixi.sh/) for environment management
- Input MAT tree (`.pb` or `.pb.gz` format)
- Sample metadata (`.tsv` or `.tsv.gz` format)

### Installation

```bash
git clone <repository-url>
cd mat-pcp
pixi install
```

### Configuration

Edit `config/config.yaml` to specify your inputs:

```yaml
mat: "path/to/your/tree.pb.gz"
metadata: "path/to/your/metadata.tsv.gz"
output_dir: "results"
```

### Running the Pipeline

```bash
# Dry run to see what will be executed
pixi run snakemake -n

# Run with 4 cores
pixi run snakemake -j 4 --use-conda

# Run specific target
pixi run snakemake -j 4 --use-conda results/aamut_fitness/all_aamut_fitness.csv
```

## Pipeline Stages

| Stage | Description | Key Outputs |
|-------|-------------|-------------|
| 1. Input | Filter metadata and link input files | `filtered_metadata.tsv.gz` |
| 2. Reference | Fetch reference genome and annotations | `ref.fa`, `coding_sites.csv` |
| 3. Clades | Extract clade subtrees and PCPs | `clades/{clade}_pcps.csv`, `all_pcps.csv` |
| 4. Mutations | Count mutations per clade | `mutation_counts/aggregated.csv` |
| 5. Rates | Compute expected mutation counts | `expected_vs_actual_mut_counts.csv` |
| 6. GLM | Estimate neutral mutation rates | `rates_table.csv`, `master_table.csv` |
| 7. Fitness | Compute fitness estimates | `ntmut_fitness_aggregated.csv`, `all_aamut_fitness.csv` |
| 8. Validation | Validate pipeline consistency | `pcp_vs_counts_report.txt` |

## Primary Outputs

| Output | Description |
|--------|-------------|
| `mat/all_pcps.csv` | All parent-child pairs with clade annotations |
| `expected_vs_actual_mut_counts.csv` | Expected and actual mutation counts |
| `rates_table.csv` | Neutral mutation rates by context |
| `glm/master_table.csv` | Full GLM rate table with light-switch boundaries |
| `fitness/ntmut_fitness_aggregated.csv` | Nucleotide-level fitness estimates |
| `aamut_fitness/all_aamut_fitness.csv` | Amino acid-level fitness estimates |
| `validation/pcp_vs_counts_report.txt` | Validation report |

## Configuration Options

See `config/config.yaml` for all configuration options. Key parameters:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `output_dir` | Output directory | `results` |
| `qc_filter_options.min_date` | Earliest sample date | `2019-12-01` |
| `qc_filter_options.datetime_cutoff` | Latest sample date | `2022-06-01` |
| `min_clade_size` | Minimum samples per clade | `2` |
| `max_nt_mutations` | Maximum mutations per branch | `4` |
| `clade_column` | Metadata column for clades | `Nextstrain_clade` |
| `max_clades` | Limit number of clades (for testing) | `null` |

### QC Filtering Options

The `qc_filter_options` section controls metadata filtering:
- Date resolution filtering (exclude MISSING, MONTH, YEAR resolutions)
- Date range filtering
- Host species validation
- Genome length constraints

### Data Files

The pipeline uses reference data files in `data/`:
- `usher_masked_sites.yaml`: UShER masking information for problematic sites by clade
- `lan_2022/`: RNA secondary structure data from Lan et al. 2022

## Output Structure

Results are organized in the output directory:

```
results/
├── mat/
│   ├── clades/           # Per-clade trees and PCPs
│   ├── all_pcps.csv      # Aggregated PCPs
│   └── clade_counts.csv  # Clade sample counts
├── ref/
│   ├── ref.fa            # Reference genome
│   └── coding_sites.csv  # Coding site annotations
├── mutation_counts/      # Per-clade and aggregated counts
├── glm/                  # GLM rate estimation outputs
├── fitness/              # Nucleotide fitness outputs
├── aamut_fitness/        # Amino acid fitness outputs
├── validation/           # Validation reports
└── log/                  # Processing logs
```

## Testing

```bash
# Run pipeline on test data
pixi run snakemake -j 4 --use-conda \
  --config output_dir="results_test_data" \
           mat="test_data/filtered_tree_raw.pb.gz" \
           metadata="test_data/filtered_metadata.tsv.gz" \
           max_clades=3

# Validate results against expected
pixi run python workflow/scripts/validate_test.py
```

The validation checks:
- Total number of PCPs extracted
- Number of clades processed
- Total sequences across all clades
- Presence of expected output files
- PCP vs mutation count consistency

## CI/CD

The pipeline includes automated testing via GitHub Actions that runs on every push and pull request:
1. Validates configuration schema
2. Verifies required data files
3. Runs the full pipeline on test data
4. Validates outputs against expected results
5. Checks PCP vs mutation count validation

## License

See LICENSE file for details.
