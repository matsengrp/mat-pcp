#!/usr/bin/env python3
"""Compute expected mutation counts from synonymous rates."""

import numpy as np
import pandas as pd

# Get snakemake variables
rates_by_clade_csv = snakemake.input.rates_by_clade
clade_founder_nts_csv = snakemake.input.clade_founder_nts_csv
output_csv = snakemake.output.expected_counts
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log("Computing expected mutation counts from synonymous rates")

# Load inputs
rates_by_clade = pd.read_csv(rates_by_clade_csv)
clade_founder_nts = pd.read_csv(clade_founder_nts_csv)

log(f"Loaded {len(rates_by_clade)} rate entries")
log(f"Loaded {len(clade_founder_nts)} clade founder entries")

# Create all possible mutations for each site using vectorized operations
# Expand each site to 3 mutations (to each other nucleotide)
nts = ["A", "C", "G", "T"]
mutation_dfs = []

for to_nt in nts:
    # Filter to sites where the founder nt differs from to_nt
    subset = clade_founder_nts[clade_founder_nts["nt"] != to_nt].copy()
    subset["to_nt"] = to_nt
    subset["mut_type"] = subset["nt"] + to_nt
    subset["nt_mutation"] = subset["nt"] + subset["site"].astype(str) + to_nt
    mutation_dfs.append(subset)

# Combine all mutations
all_mutations = pd.concat(mutation_dfs, ignore_index=True)
log(f"Generated {len(all_mutations)} mutation possibilities")

# Merge with rates to get expected counts
# Create rate lookup as a dataframe for efficient merge
rate_lookup = rates_by_clade[rates_by_clade["rate"].notna()][
    ["clade", "mut_type", "rate"]
].copy()

all_mutations = all_mutations.merge(
    rate_lookup, on=["clade", "mut_type"], how="left"
)

# Fill missing rates with 0
all_mutations["expected_count"] = all_mutations["rate"].fillna(0)

# Select and rename columns for output
expected_df = all_mutations[
    [
        "clade",
        "site",
        "nt_mutation",
        "mut_type",
        "expected_count",
        "gene",
        "codon_site",
        "codon_position",
    ]
].rename(columns={"site": "nt_site"})

# Log summary per clade
for clade in expected_df["clade"].unique():
    clade_data = expected_df[expected_df["clade"] == clade]
    n_with_rate = (clade_data["expected_count"] > 0).sum()
    log(f"  {clade}: {len(clade_data)} mutations, {n_with_rate} with rates")

expected_df.to_csv(output_csv, index=False)

log(f"Computed expected counts for {len(expected_df)} mutation possibilities")
log("Done!")
