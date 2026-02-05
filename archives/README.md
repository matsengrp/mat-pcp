# Archives Directory

This directory stores archived results from previous pipeline runs for comparison purposes. Archives are **not version controlled** but may be preserved locally for validation and reproducibility studies.

## Purpose

Use this directory to:
- Preserve results from important production runs
- Compare outputs between pipeline versions
- Validate changes to pipeline logic
- Retain historical results for publications

## Naming Convention

Archive directories should follow this format:

**Format:** `{dataset}_{date}`

### Examples

- `public_2022-06-01` - Full public dataset through June 2022
- `test_data_2026-02-05` - Test run from February 5, 2026
- `omicron_2023-12-15` - Omicron-specific analysis from Dec 2023

**For multi-clade subsets:** `{dataset}_{clade}_{date}`
- `public_19A_2022-06-01` - Only clade 19A from public dataset

### What Not to Use

Avoid these patterns:
- Generic names like `old/` or `backup/`
- Version numbers like `v1/` or `v2/` (doesn't indicate data source)
- Temporary suffixes like `_old/` or `_backup/` (delete instead of archiving)

## Creating an Archive

To archive a completed run:

```bash
# Move from results to archives with descriptive name
mv results/public_2022-06-01 archives/

# Or create a compressed archive for long-term storage
tar -czf archives/public_2022-06-01_$(date +%Y%m%d).tar.gz results/public_2022-06-01/
```

## Managing Archives

Archives can consume significant disk space. Consider:

1. **Compress old archives:**
   ```bash
   cd archives
   tar -czf public_2022-06-01.tar.gz public_2022-06-01/
   rm -rf public_2022-06-01/
   ```

2. **Keep only key outputs:**
   ```bash
   # Preserve final results but remove intermediate files
   cd archives/public_2022-06-01
   rm -rf mat/clades/ mutations/ mutation_counts/[0-9]*.csv log/ benchmarks/
   ```

3. **Document what's archived:**
   Add a brief note here about important archives:
   ```
   - public_2022-06-01: Results for Bloom et al. (2024) paper
   - test_data_2026-01-15: Validation run before major refactor
   ```

## Comparing Results

To compare outputs from different runs:

```bash
# Compare fitness estimates
pixi run python workflow/scripts/compare_results.py \
  archives/public_2022-06-01/aamut_fitness/all_aamut_fitness.csv \
  results/public_2022-06-01/aamut_fitness/all_aamut_fitness.csv
```

## Git Status

This directory is git-ignored (except this README) to prevent accidentally committing large result files. Only this README.md is tracked.
