#!/usr/bin/env python3
"""Extract rates for sites after boundary (final output format)."""

import pandas as pd

# Get snakemake variables
master_table_csv = snakemake.input.master_table
output_csv = snakemake.output.rates_table
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log(f"Extracting final rates table from {master_table_csv}")

# Load master table
rates = pd.read_csv(master_table_csv)
log(f"Loaded {len(rates)} rate entries")

# Filter to rates after boundary
if "nt_site_before_boundary" in rates.columns:
    rates_after = rates[~rates["nt_site_before_boundary"]]
else:
    rates_after = rates

log(f"Rates after boundary: {len(rates_after)}")

# Select output columns
output_cols = ["mut_type", "motif", "unpaired", "rate"]
available_cols = [c for c in output_cols if c in rates_after.columns]

rates_after[available_cols].to_csv(output_csv, index=False)

log(f"Exported {len(rates_after)} rate entries to {output_csv}")
log("Done!")
