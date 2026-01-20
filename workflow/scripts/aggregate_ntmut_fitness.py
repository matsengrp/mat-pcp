#!/usr/bin/env python3
"""Aggregate nucleotide fitness across all clades."""

import os
import pandas as pd

# Get snakemake variables
fitness_files = snakemake.input.fitness_files
output_csv = snakemake.output.csv
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log(f"Aggregating nucleotide fitness from {len(fitness_files)} files")

dfs = []
for fitness_file in fitness_files:
    if not os.path.exists(fitness_file):
        log(f"WARNING: File not found: {fitness_file}")
        continue

    df = pd.read_csv(fitness_file)
    # Extract clade from filename
    basename = os.path.basename(fitness_file)
    clade = basename.replace("_ntmut_fitness.csv", "")
    df["source_clade"] = clade
    dfs.append(df)
    log(f"  {clade}: {len(df)} mutations")

if len(dfs) == 0:
    log("WARNING: No fitness files found")
    aggregated = pd.DataFrame()
else:
    aggregated = pd.concat(dfs, ignore_index=True)

aggregated.to_csv(output_csv, index=False)

log(f"Aggregated {len(aggregated)} fitness entries")
log(f"Wrote output to {output_csv}")
log("Done!")
