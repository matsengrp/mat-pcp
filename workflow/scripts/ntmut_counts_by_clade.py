#!/usr/bin/env python3
"""Extract nucleotide mutation counts for a specific clade or cluster."""

import pandas as pd

# Get snakemake variables
counts_csv = snakemake.input.counts
output_csv = snakemake.output.csv
clade = snakemake.params.clade
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log(f"Extracting nucleotide mutation counts for clade: {clade}")

# Load predicted counts
counts_df = pd.read_csv(counts_csv)
log(f"Loaded {len(counts_df)} count entries")

if clade == "all":
    # Aggregate counts across all clades
    log("Aggregating counts across all clades")

    # Columns to exclude from groupby (they get aggregated)
    agg_cols = [
        "clade",
        "actual_count",
        "count_terminal",
        "count_non_terminal",
        "expected_count",
        "tau_squared",
    ]
    groupby_cols = [col for col in counts_df.columns if col not in agg_cols]

    agg_dict = {}
    for col in [
        "actual_count",
        "count_terminal",
        "count_non_terminal",
        "expected_count",
    ]:
        if col in counts_df.columns:
            agg_dict[col] = "sum"

    # For tau_squared, take mean (residual variance is similar across clades
    # for the same mutation context)
    for col in ["tau_squared"]:
        if col in counts_df.columns:
            agg_dict[col] = "mean"

    clade_counts = counts_df.groupby(groupby_cols, as_index=False).agg(agg_dict)
    clade_counts["cluster"] = "all"
else:
    # Filter to specific clade
    log(f"Filtering to clade: {clade}")
    clade_counts = counts_df[counts_df["clade"] == clade].copy()
    clade_counts["cluster"] = clade

log(f"Extracted {len(clade_counts)} mutations for cluster {clade}")

clade_counts.to_csv(output_csv, index=False)
log(f"Wrote output to {output_csv}")
log("Done!")
