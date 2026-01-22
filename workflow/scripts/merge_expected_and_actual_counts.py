#!/usr/bin/env python3
"""Merge expected and actual mutation counts with exclusion annotations."""

import pandas as pd

# Get snakemake variables
expected_csv = snakemake.input.expected
actual_csv = snakemake.input.actual
muts_to_exclude_csv = snakemake.input.muts_to_exclude
output_csv = snakemake.output.csv
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log("Merging expected and actual mutation counts")

# Load inputs
expected = pd.read_csv(expected_csv)
actual = pd.read_csv(actual_csv)
muts_to_exclude = pd.read_csv(muts_to_exclude_csv)

log(f"Expected counts: {len(expected)} entries")
log(f"Actual counts: {len(actual)} entries")
log(f"Mutations to exclude: {len(muts_to_exclude)} entries")

# Identify actual count column
count_col = "count" if "count" in actual.columns else "actual_count"
actual_cols = ["clade", "nt_mutation", count_col]
if "count_terminal" in actual.columns:
    actual_cols.append("count_terminal")
if "count_non_terminal" in actual.columns:
    actual_cols.append("count_non_terminal")

# Merge expected and actual counts
merged = expected.merge(
    actual[actual_cols].rename(
        columns={count_col: "actual_count"} if count_col != "actual_count" else {}
    ),
    on=["clade", "nt_mutation"],
    how="left",
)

# Fill missing actual counts with 0
merged["actual_count"] = merged["actual_count"].fillna(0).astype(int)
if "count_terminal" in merged.columns:
    merged["count_terminal"] = merged["count_terminal"].fillna(0).astype(int)
if "count_non_terminal" in merged.columns:
    merged["count_non_terminal"] = merged["count_non_terminal"].fillna(0).astype(int)

log(f"Merged data: {len(merged)} entries")

# Add exclusion annotations using vectorized merge instead of apply
# Prepare exclusion data with renamed column for merge
exclusion_df = muts_to_exclude[["clade", "mutation", "masked_in_usher"]].rename(
    columns={"mutation": "nt_mutation"}
)

# Merge to get exclusion info
merged = merged.merge(
    exclusion_df, on=["clade", "nt_mutation"], how="left", suffixes=("", "_excl")
)

# Set exclude flag based on whether the mutation was found in exclusion list
merged["exclude"] = merged["masked_in_usher"].notna()
merged["masked_in_usher"] = merged["masked_in_usher"].fillna(False)
merged["masked_in_usher"] = merged["masked_in_usher"].infer_objects(copy=False).astype(bool)

# Add v2-compatible columns
# subset: v2 uses geographic subsets (e.g., "USA"), mat-pcp uses "all" for all data
merged["subset"] = "all"
# mean_log_size: average log of descendant clade size, not tracked by mat-pcp (set to 0)
merged["mean_log_size"] = 0.0

# Reorder columns to match v2 expected_vs_actual_mut_counts.csv format
# Note: mut_type is needed by annotate_counts.py for downstream processing
output_columns = [
    "clade",
    "subset",
    "nt_site",
    "nt_mutation",
    "mut_type",
    "exclude",
    "masked_in_usher",
    "expected_count",
    "actual_count",
    "clade_founder_nt",
    "gene",
    "clade_founder_codon",
    "clade_founder_aa",
    "mutant_codon",
    "mutant_aa",
    "aa_mutation",
    "synonymous",
    "noncoding",
    "codon_position",
    "codon_site",
    "four_fold_degenerate",
    "count_terminal",
    "count_non_terminal",
    "mean_log_size",
]
# Only include columns that exist in the merged data
available_cols = [c for c in output_columns if c in merged.columns]
merged = merged[available_cols]

# Save merged data
merged.to_csv(output_csv, index=False)

log(f"Excluded mutations: {merged['exclude'].sum()}")
log(f"UShER masked: {merged['masked_in_usher'].sum()}")
log(f"Wrote {len(merged)} entries to {output_csv}")
log("Done!")
