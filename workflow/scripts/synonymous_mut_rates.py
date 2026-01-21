#!/usr/bin/env python3
"""Compute synonymous mutation rates and spectra by clade.

Uses the same normalization approach as SARS2-mut-fitness:
- rate = fraction / parent_nt_frac
- where fraction = count / total_count (fraction of all synonymous mutations)
- and parent_nt_frac = n_sites / total_sites (fraction of 4-fold sites with parent nt)
"""

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


log("Computing synonymous mutation rates by clade (SARS2-mut-fitness compatible)")

# Load inputs
mutation_counts = pd.read_csv(mutation_counts_csv)
clade_founder_nts = pd.read_csv(clade_founder_nts_csv)

log(f"Loaded {len(mutation_counts)} mutation counts")
log(f"Loaded {len(clade_founder_nts)} clade founder nucleotide entries")

# Get 4-fold degenerate sites from clade_founder_nts
fourfold_sites = clade_founder_nts.query("four_fold_degenerate == True")[
    ["clade", "site"]
].drop_duplicates()
fourfold_sites["is_4fold"] = True
log(f"Total 4-fold degenerate sites across clades: {len(fourfold_sites)}")

# Filter to synonymous mutations and merge with 4-fold info
synonymous = mutation_counts.query("synonymous == True and exclude == False").copy()
synonymous = synonymous.merge(
    fourfold_sites,
    left_on=["clade", "nt_site"],
    right_on=["clade", "site"],
    how="inner",
)
log(f"Synonymous 4-fold degenerate non-excluded mutations: {len(synonymous)}")

if len(synonymous) == 0:
    log("No synonymous mutations found, creating empty output")
    rates_df = pd.DataFrame(
        columns=[
            "clade",
            "mut_type",
            "count",
            "total_count",
            "fraction",
            "parent_nt",
            "parent_nt_frac",
            "rate",
        ]
    )
    rates_df.to_csv(output_csv, index=False)
    log("Done!")
    exit()

# Add mutation type column (e.g., "AtoC" format to match SARS2-mut-fitness)
synonymous["mut_type"] = synonymous["nt_mutation"].apply(
    lambda x: f"{x[0]}to{x[-1]}"
)

# Get count column name
count_col = "count" if "count" in synonymous.columns else "actual_count"

# Tally mutation type counts for each clade
mut_type_counts = (
    synonymous.groupby(["clade", "mut_type"], as_index=False)
    .agg({count_col: "sum"})
    .rename(columns={count_col: "count"})
)

# Calculate parent nucleotide composition at 4-fold degenerate sites
parent_composition = (
    clade_founder_nts.query("four_fold_degenerate == True")
    .groupby(["clade", "nt"], as_index=False)
    .agg(parent_nt_count=pd.NamedAgg("site", "nunique"))
    .assign(
        parent_nt_frac=lambda x: (
            x["parent_nt_count"]
            / x.groupby("clade")["parent_nt_count"].transform("sum")
        ),
    )
    .rename(columns={"nt": "parent_nt"})
)

# Add parent_nt column to mut_type_counts
mut_type_counts["parent_nt"] = mut_type_counts["mut_type"].str[0]

# Calculate total counts per clade and filter by minimum count threshold
mut_type_counts["total_count"] = mut_type_counts.groupby("clade")["count"].transform(
    "sum"
)

# Filter to clades with adequate counts
mut_type_counts = mut_type_counts.query("total_count >= @min_counts").copy()

# Calculate fraction of all mutations that are this type
mut_type_counts["fraction"] = (
    mut_type_counts["count"] / mut_type_counts["total_count"]
)

# Merge with parent nucleotide composition
rates_df = mut_type_counts.merge(
    parent_composition[["clade", "parent_nt", "parent_nt_frac"]],
    on=["clade", "parent_nt"],
    how="left",
)

# Calculate rate = fraction / parent_nt_frac (same as SARS2-mut-fitness)
rates_df["rate"] = rates_df["fraction"] / rates_df["parent_nt_frac"]

# Select and order columns to match SARS2-mut-fitness format
rates_df = rates_df[
    [
        "clade",
        "mut_type",
        "count",
        "total_count",
        "fraction",
        "parent_nt",
        "parent_nt_frac",
        "rate",
    ]
].sort_values(["clade", "mut_type"])

rates_df.to_csv(output_csv, index=False)

log(f"Computed rates for {len(rates_df['clade'].unique())} clades")
log(f"Valid rates: {rates_df['rate'].notna().sum()} / {len(rates_df)}")
log("Done!")
