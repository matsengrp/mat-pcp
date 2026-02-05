#!/usr/bin/env python3
"""Validate that PCP-derived AA mutation counts match aamut_fitness actual_count.

Compares amino acid mutations (e.g., D614G) between:
1. PCP-derived counts: Mutations from translating parent-child pairs
2. aamut_fitness actual_count: Aggregated from mutation counting pipeline
"""

import pandas as pd

# Get snakemake variables
pcp_counts_csv = snakemake.input.pcp_counts
aamut_fitness_files = snakemake.input.aamut_fitness_files
all_aamut_fitness_csv = snakemake.input.all_aamut_fitness
output_report = snakemake.output.report
output_discrepancies = snakemake.output.discrepancies
pcp_gene_name = snakemake.params.pcp_gene_name
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log("Validating PCP-derived AA mutation counts against aamut_fitness actual_count")


# Load PCP-derived counts
pcp_counts = pd.read_csv(pcp_counts_csv)
log(f"PCP-derived AA mutation counts: {len(pcp_counts)} entries")

# Load per-clade aamut_fitness files
aamut_per_clade = []
for f in aamut_fitness_files:
    df = pd.read_csv(f)
    aamut_per_clade.append(df)

if aamut_per_clade:
    aamut_clade_df = pd.concat(aamut_per_clade, ignore_index=True)
    # Rename cluster to clade for consistency
    if "cluster" in aamut_clade_df.columns:
        aamut_clade_df = aamut_clade_df.rename(columns={"cluster": "clade"})
else:
    aamut_clade_df = pd.DataFrame()

log(f"Per-clade aamut_fitness: {len(aamut_clade_df)} entries")

# Load aggregate aamut_fitness
aamut_all = pd.read_csv(all_aamut_fitness_csv)
log(f"Aggregate aamut_fitness: {len(aamut_all)} entries")

# Filter to gene of interest
if pcp_gene_name and len(aamut_clade_df) > 0:
    aamut_clade_df = aamut_clade_df[aamut_clade_df["gene"] == pcp_gene_name]
    log(f"Filtered per-clade to {pcp_gene_name}: {len(aamut_clade_df)} entries")

if pcp_gene_name and len(aamut_all) > 0:
    aamut_all = aamut_all[aamut_all["gene"] == pcp_gene_name]
    log(f"Filtered aggregate to {pcp_gene_name}: {len(aamut_all)} entries")

# Also filter PCP counts to the gene
if pcp_gene_name and len(pcp_counts) > 0 and "gene" in pcp_counts.columns:
    pcp_counts = pcp_counts[pcp_counts["gene"] == pcp_gene_name]
    log(f"Filtered PCP counts to {pcp_gene_name}: {len(pcp_counts)} entries")


def compare_counts(pcp_df, aamut_df, clade_col="clade", is_aggregate=False):
    """
    Compare PCP-derived counts with aamut_fitness actual_count by aa_mutation.

    Returns: (n_matches, n_discrepancies, discrepancies_df)
    """
    if len(pcp_df) == 0 and len(aamut_df) == 0:
        return 0, 0, pd.DataFrame()

    # Get mutations with actual_count > 0 from aamut_fitness
    if "actual_count" in aamut_df.columns:
        aamut_with_counts = aamut_df[aamut_df["actual_count"] > 0].copy()
    else:
        aamut_with_counts = pd.DataFrame()

    # Prepare merge keys - now using aa_mutation instead of just aa_site
    if is_aggregate:
        # For aggregate, we don't have clade column
        pcp_agg = pcp_df[pcp_df["clade"] == "all"].copy()
        merge_keys = ["gene", "aa_mutation"]

        if len(pcp_agg) > 0 and len(aamut_with_counts) > 0:
            merged = pcp_agg.merge(
                aamut_with_counts[["gene", "aa_mutation", "actual_count"]],
                on=merge_keys,
                how="outer",
            )
        elif len(pcp_agg) > 0:
            merged = pcp_agg.copy()
            merged["actual_count"] = 0
        elif len(aamut_with_counts) > 0:
            merged = aamut_with_counts[["gene", "aa_mutation", "actual_count"]].copy()
            merged["pcp_count"] = 0
        else:
            return 0, 0, pd.DataFrame()
    else:
        # For per-clade comparison
        pcp_per_clade = pcp_df[pcp_df["clade"] != "all"].copy()
        merge_keys = [clade_col, "gene", "aa_mutation"]

        if len(pcp_per_clade) > 0 and len(aamut_with_counts) > 0:
            merged = pcp_per_clade.merge(
                aamut_with_counts[[clade_col, "gene", "aa_mutation", "actual_count"]],
                on=merge_keys,
                how="outer",
            )
        elif len(pcp_per_clade) > 0:
            merged = pcp_per_clade.copy()
            merged["actual_count"] = 0
        elif len(aamut_with_counts) > 0:
            merged = aamut_with_counts[
                [clade_col, "gene", "aa_mutation", "actual_count"]
            ].copy()
            merged["pcp_count"] = 0
        else:
            return 0, 0, pd.DataFrame()

    # Fill NaN with 0
    merged["pcp_count"] = merged["pcp_count"].fillna(0).astype(int)
    merged["actual_count"] = merged["actual_count"].fillna(0).astype(int)

    # Find discrepancies
    merged["difference"] = merged["pcp_count"] - merged["actual_count"]
    discrepancies = merged[merged["difference"] != 0].copy()

    n_matches = len(merged) - len(discrepancies)
    n_discrepancies = len(discrepancies)

    return n_matches, n_discrepancies, discrepancies


# Compare per-clade
log("\n=== Per-clade comparison ===")
clade_matches, clade_disc_count, clade_disc = compare_counts(
    pcp_counts, aamut_clade_df, clade_col="clade", is_aggregate=False
)
log(f"Per-clade matches: {clade_matches}")
log(f"Per-clade discrepancies: {clade_disc_count}")

# Compare aggregate
log("\n=== Aggregate comparison ===")
agg_matches, agg_disc_count, agg_disc = compare_counts(
    pcp_counts, aamut_all, is_aggregate=True
)
log(f"Aggregate matches: {agg_matches}")
log(f"Aggregate discrepancies: {agg_disc_count}")

# Combine discrepancies
all_discrepancies = pd.concat([clade_disc, agg_disc], ignore_index=True)
total_discrepancies = len(all_discrepancies)

# Categorize discrepancies
if total_discrepancies > 0:
    only_in_pcp = all_discrepancies[all_discrepancies["actual_count"] == 0]
    only_in_aamut = all_discrepancies[all_discrepancies["pcp_count"] == 0]
    count_mismatch = all_discrepancies[
        (all_discrepancies["pcp_count"] > 0) & (all_discrepancies["actual_count"] > 0)
    ]
else:
    only_in_pcp = pd.DataFrame()
    only_in_aamut = pd.DataFrame()
    count_mismatch = pd.DataFrame()

# Write report
with open(output_report, "w") as f:
    f.write("=" * 60 + "\n")
    f.write("PCP vs aamut_fitness Validation Report\n")
    f.write("=" * 60 + "\n\n")

    f.write("SUMMARY\n")
    f.write("-" * 40 + "\n")
    f.write(f"Gene: {pcp_gene_name}\n")
    f.write(f"Per-clade exact matches: {clade_matches}\n")
    f.write(f"Per-clade discrepancies: {clade_disc_count}\n")
    f.write(f"Aggregate exact matches: {agg_matches}\n")
    f.write(f"Aggregate discrepancies: {agg_disc_count}\n")
    f.write(f"Total discrepancies: {total_discrepancies}\n\n")

    if total_discrepancies == 0:
        f.write("SUCCESS: All amino acid mutation counts match exactly!\n\n")
        f.write(
            "The PCP extraction and aamut_fitness pipelines are consistent.\n"
        )
    else:
        f.write("DISCREPANCY BREAKDOWN\n")
        f.write("-" * 40 + "\n")
        f.write(f"Only in PCPs (not in aamut_fitness): {len(only_in_pcp)}\n")
        f.write(f"Only in aamut_fitness (not in PCPs): {len(only_in_aamut)}\n")
        f.write(f"Count mismatches (both have counts): {len(count_mismatch)}\n\n")

        if len(only_in_pcp) > 0:
            f.write("\nAA mutations only in PCPs (first 10):\n")
            cols = ["clade", "gene", "aa_mutation", "pcp_count"] if "clade" in only_in_pcp.columns else ["gene", "aa_mutation", "pcp_count"]
            cols = [c for c in cols if c in only_in_pcp.columns]
            f.write(only_in_pcp.head(10)[cols].to_string())
            f.write("\n")

        if len(only_in_aamut) > 0:
            f.write("\nAA mutations only in aamut_fitness (first 10):\n")
            cols = ["clade", "gene", "aa_mutation", "actual_count"] if "clade" in only_in_aamut.columns else ["gene", "aa_mutation", "actual_count"]
            cols = [c for c in cols if c in only_in_aamut.columns]
            f.write(only_in_aamut.head(10)[cols].to_string())
            f.write("\n")

        if len(count_mismatch) > 0:
            f.write("\nCount mismatches (first 10):\n")
            cols = ["clade", "gene", "aa_mutation", "pcp_count", "actual_count", "difference"] if "clade" in count_mismatch.columns else ["gene", "aa_mutation", "pcp_count", "actual_count", "difference"]
            cols = [c for c in cols if c in count_mismatch.columns]
            f.write(count_mismatch.head(10)[cols].to_string())
            f.write("\n")

    f.write("\n" + "=" * 60 + "\n")
    f.write("End of Report\n")
    f.write("=" * 60 + "\n")

# Write discrepancies to CSV
all_discrepancies.to_csv(output_discrepancies, index=False)

log(f"\nWrote report to {output_report}")
log(f"Wrote {len(all_discrepancies)} discrepancies to {output_discrepancies}")

# Report validation status (don't exit with error - just report)
if total_discrepancies > 0:
    log(f"\nValidation found {total_discrepancies} discrepancies")
else:
    log("\nValidation passed - all counts match!")
