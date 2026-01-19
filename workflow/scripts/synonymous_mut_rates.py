#!/usr/bin/env python3
"""Compute synonymous mutation rates and spectra by clade."""

import numpy as np
import pandas as pd

# Get snakemake variables
mutation_counts_csv = snakemake.input.mutation_counts_csv
clade_founder_nts_csv = snakemake.input.clade_founder_nts_csv
output_csv = snakemake.output.rates_by_clade
min_counts = snakemake.params.synonymous_spectra_min_counts
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log("Computing synonymous mutation rates by clade")

# Load inputs
mutation_counts = pd.read_csv(mutation_counts_csv)
clade_founder_nts = pd.read_csv(clade_founder_nts_csv)

log(f"Loaded {len(mutation_counts)} mutation counts")
log(f"Loaded {len(clade_founder_nts)} clade founder nucleotide entries")

# Filter to synonymous mutations
synonymous = mutation_counts.query("synonymous == True").copy()
log(f"Synonymous mutations: {len(synonymous)}")

if len(synonymous) == 0:
    log("No synonymous mutations found, creating empty output")
    rates_df = pd.DataFrame(
        columns=["clade", "mut_type", "actual_count", "n_sites", "rate"]
    )
    rates_df.to_csv(output_csv, index=False)
    log("Done!")
    exit()

# Calculate mutation rates by clade and mutation type
rates_by_clade = []

for clade in synonymous["clade"].unique():
    clade_data = synonymous[synonymous["clade"] == clade]

    # Get founder info for this clade
    founder = clade_founder_nts[clade_founder_nts["clade"] == clade]

    # Count 4-fold degenerate sites by nucleotide
    fourfold = founder[founder["four_fold_degenerate"] == True]
    nt_counts = fourfold["nt"].value_counts().to_dict()

    log(f"  {clade}: {len(fourfold)} 4-fold degenerate sites")

    # Calculate rates for each mutation type
    for mut_type in [
        "AC",
        "AG",
        "AT",
        "CA",
        "CG",
        "CT",
        "GA",
        "GC",
        "GT",
        "TA",
        "TC",
        "TG",
    ]:
        from_nt = mut_type[0]
        to_nt = mut_type[1]

        # Count mutations of this type
        # Handle column name variations
        count_col = "count" if "count" in clade_data.columns else "actual_count"
        type_muts = clade_data[
            (clade_data["nt_mutation"].str[0] == from_nt)
            & (clade_data["nt_mutation"].str[-1] == to_nt)
        ]
        mut_count = type_muts[count_col].sum() if len(type_muts) > 0 else 0

        # Get number of possible sites
        n_sites = nt_counts.get(from_nt, 0)

        # Calculate rate (with minimum count threshold)
        if n_sites > 0 and mut_count >= min_counts:
            rate = mut_count / n_sites
        else:
            rate = np.nan

        rates_by_clade.append(
            {
                "clade": clade,
                "mut_type": mut_type,
                "actual_count": mut_count,
                "n_sites": n_sites,
                "rate": rate,
            }
        )

rates_df = pd.DataFrame(rates_by_clade)
rates_df.to_csv(output_csv, index=False)

log(f"Computed rates for {len(rates_df['clade'].unique())} clades")
log(f"Valid rates: {rates_df['rate'].notna().sum()} / {len(rates_df)}")
log("Done!")
