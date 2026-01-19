#!/usr/bin/env python3
"""Aggregate mutation counts across all clades."""

import os

import pandas as pd

# Get snakemake variables
count_files = snakemake.input.counts
output_csv = snakemake.output.csv
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log(f"Aggregating {len(count_files)} count files")

# Aggregate all count files
dfs = []
for count_file in count_files:
    df = pd.read_csv(count_file)
    # Extract clade from filename (format: {clade}.csv)
    basename = os.path.basename(count_file).replace(".csv", "")
    clade = basename

    df["clade"] = clade
    dfs.append(df)
    log(f"  {clade}: {len(df)} mutations")

if len(dfs) == 0:
    log("No count files, creating empty output")
    # Create empty output with expected columns
    aggregated = pd.DataFrame(
        columns=[
            "nt_mutation",
            "nt_site",
            "count",
            "count_terminal",
            "count_non_terminal",
            "clade",
            "synonymous",
            "noncoding",
            "protein",
            "aa_mutation",
            "codon_change",
            "reference_nt",
            "clade_founder_nt",
            "exclude",
            "mean_log_size",
        ]
    )
else:
    aggregated = pd.concat(dfs, ignore_index=True)
    log(f"Total aggregated: {len(aggregated)} mutations")

    # Reorder columns to put clade first
    cols = aggregated.columns.tolist()
    if "clade" in cols:
        cols.remove("clade")
        cols = ["clade"] + cols
        aggregated = aggregated[cols]

aggregated.to_csv(output_csv, index=False)
log(f"Wrote aggregated counts to {output_csv}")
log("Done!")
