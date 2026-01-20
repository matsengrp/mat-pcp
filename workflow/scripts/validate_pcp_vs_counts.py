#!/usr/bin/env python3
"""Validate that PCP-derived counts match aggregated mutation counts."""

import pandas as pd

# Get snakemake variables
pcp_counts_csv = snakemake.input.pcp_counts
mutation_counts_csv = snakemake.input.mutation_counts
output_report = snakemake.output.report
output_discrepancies = snakemake.output.discrepancies
pcp_gene_region = snakemake.params.pcp_gene_region
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log("Validating PCP-derived counts against aggregated mutation counts")

# Load inputs
pcp_counts = pd.read_csv(pcp_counts_csv)
mutation_counts = pd.read_csv(mutation_counts_csv)

log(f"PCP counts: {len(pcp_counts)} entries")
log(f"Aggregated mutation counts: {len(mutation_counts)} entries")

# Filter mutation_counts to PCP region if applicable
if pcp_gene_region:
    start = pcp_gene_region.get("start", 1)
    end = pcp_gene_region.get("end", len(mutation_counts))
    mutation_counts = mutation_counts[
        (mutation_counts["nt_site"] >= start) & (mutation_counts["nt_site"] <= end)
    ]
    log(f"Filtered to PCP region [{start}, {end}]: {len(mutation_counts)} entries")

# Identify the count column in mutation_counts
count_col = "count" if "count" in mutation_counts.columns else "actual_count"
mutation_counts = mutation_counts.rename(columns={count_col: "aggregated_count"})

# Merge PCP counts with aggregated counts
merged = pcp_counts.merge(
    mutation_counts[["clade", "nt_mutation", "aggregated_count"]],
    on=["clade", "nt_mutation"],
    how="outer",
)

# Fill NaN with 0
merged["pcp_count"] = merged["pcp_count"].fillna(0).astype(int)
merged["aggregated_count"] = merged["aggregated_count"].fillna(0).astype(int)

# Find discrepancies
merged["difference"] = merged["pcp_count"] - merged["aggregated_count"]
discrepancies = merged[merged["difference"] != 0].copy()

# Summary statistics
pcp_per_clade = pcp_counts[pcp_counts["clade"] != "all"]
total_pcp = (
    len(pcp_per_clade["nt_mutation"].unique()) if len(pcp_per_clade) > 0 else 0
)
total_agg = len(mutation_counts["nt_mutation"].unique())
n_discrepancies = len(discrepancies)
n_exact_matches = len(merged) - n_discrepancies

# Categorize discrepancies
only_in_pcp = discrepancies[discrepancies["aggregated_count"] == 0]
only_in_agg = discrepancies[discrepancies["pcp_count"] == 0]
count_mismatch = discrepancies[
    (discrepancies["pcp_count"] > 0) & (discrepancies["aggregated_count"] > 0)
]

log(f"Exact matches: {n_exact_matches}")
log(f"Discrepancies: {n_discrepancies}")

# Write report
with open(output_report, "w") as f:
    f.write("=" * 60 + "\n")
    f.write("PCP vs Mutation Counts Validation Report\n")
    f.write("=" * 60 + "\n\n")

    f.write("SUMMARY\n")
    f.write("-" * 40 + "\n")
    f.write(f"Unique mutations in PCPs: {total_pcp}\n")
    f.write(f"Unique mutations in aggregated counts: {total_agg}\n")
    f.write(f"Exact matches: {n_exact_matches}\n")
    f.write(f"Total discrepancies: {n_discrepancies}\n\n")

    if n_discrepancies == 0:
        f.write("SUCCESS: All mutation counts match exactly!\n\n")
        f.write("The PCP extraction and mutation counting pipelines are consistent.\n")
    else:
        f.write("DISCREPANCY BREAKDOWN\n")
        f.write("-" * 40 + "\n")
        f.write(f"Only in PCPs (not in aggregated): {len(only_in_pcp)}\n")
        f.write(f"Only in aggregated (not in PCPs): {len(only_in_agg)}\n")
        f.write(f"Count mismatches (both have counts): {len(count_mismatch)}\n\n")

        if len(count_mismatch) > 0:
            f.write("WARNING: Count mismatches indicate filtering differences.\n")
            f.write("Review filtering criteria in both pipelines.\n\n")

            # Show worst mismatches
            worst = count_mismatch.nlargest(10, "difference", keep="first")
            f.write("Top 10 count mismatches:\n")
            f.write(
                worst[
                    [
                        "clade",
                        "nt_mutation",
                        "pcp_count",
                        "aggregated_count",
                        "difference",
                    ]
                ].to_string()
            )
            f.write("\n")

        if len(only_in_pcp) > 0:
            f.write("\nMutations only in PCPs (first 10):\n")
            f.write(
                only_in_pcp.head(10)[["clade", "nt_mutation", "pcp_count"]].to_string()
            )
            f.write("\n")

        if len(only_in_agg) > 0:
            f.write("\nMutations only in aggregated counts (first 10):\n")
            f.write(
                only_in_agg.head(10)[
                    ["clade", "nt_mutation", "aggregated_count"]
                ].to_string()
            )
            f.write("\n")

    f.write("\n" + "=" * 60 + "\n")
    f.write("End of Report\n")
    f.write("=" * 60 + "\n")

# Write discrepancies to CSV
discrepancies.to_csv(output_discrepancies, index=False)

log(f"Wrote report to {output_report}")
log(f"Wrote {len(discrepancies)} discrepancies to {output_discrepancies}")
log("Validation complete!")
