#!/usr/bin/env python3
"""Curate counts for GLM training."""

import pandas as pd

# Get snakemake variables
counts_csv = snakemake.input.counts
founder_csv = snakemake.input.clade_founder
output_csv = snakemake.output.curated
additional_sites_to_ignore = snakemake.params.additional_sites_to_ignore or []
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log("Curating counts for GLM training")

# Load inputs
counts_df = pd.read_csv(counts_csv)
founder_df = pd.read_csv(founder_csv)

log(f"Loaded {len(counts_df)} count entries")
log(f"Loaded {len(founder_df)} founder entries")
log(f"Additional sites to ignore: {len(additional_sites_to_ignore)}")

# Identify conserved sites (same codon and motif across all clades)
conserved_data = founder_df[
    (founder_df["codon"] == founder_df["ref_codon"])
    & (founder_df["motif"] == founder_df["ref_motif"])
]
site_counts = conserved_data["nt_site"].value_counts()
n_clades = len(founder_df["clade"].unique())
conserved_sites = site_counts[site_counts == n_clades].index.tolist()

log(f"Conserved sites (present in all {n_clades} clades): {len(conserved_sites)}")

# Sites to ignore: masked, excluded, or in additional list
sites_to_ignore = []

if "masked_in_usher" in counts_df.columns:
    masked_sites = (
        counts_df[counts_df["masked_in_usher"] == True]["nt_site"].unique().tolist()
    )
    sites_to_ignore.extend(masked_sites)
    log(f"UShER masked sites: {len(masked_sites)}")

if "exclude" in counts_df.columns:
    excluded_sites = (
        counts_df[counts_df["exclude"] == True]["nt_site"].unique().tolist()
    )
    sites_to_ignore.extend(excluded_sites)
    log(f"Excluded sites: {len(excluded_sites)}")

sites_to_ignore.extend(additional_sites_to_ignore)
sites_to_ignore = list(set(sites_to_ignore))
log(f"Total sites to ignore: {len(sites_to_ignore)}")

# Filter to conserved, non-ignored sites
curated = counts_df[
    counts_df["nt_site"].isin(conserved_sites)
    & ~counts_df["nt_site"].isin(sites_to_ignore)
].copy()

log(f"After filtering: {len(curated)} entries")

# Define columns to aggregate and columns to group by
sum_cols = ["expected_count", "actual_count", "count_terminal", "count_non_terminal"]
sum_cols = [col for col in sum_cols if col in curated.columns]

# Columns to not include in groupby (will be aggregated or dropped)
ignore_cols = sum_cols + ["clade", "exclude", "masked_in_usher"]
groupby_cols = [
    col for col in curated.columns if col not in ignore_cols and col in curated.columns
]

# Only aggregate if we have data
if len(curated) > 0:
    curated_agg = curated.groupby(groupby_cols, as_index=False).agg(
        {col: "sum" for col in sum_cols}
    )

    # Verify no duplicate mutations
    if "nt_mutation" in curated_agg.columns:
        n_dups = curated_agg["nt_mutation"].duplicated().sum()
        if n_dups > 0:
            log(f"WARNING: {n_dups} duplicate mutations found")
        else:
            log("No duplicate mutations found")

    # Drop exclusion columns (no longer needed after filtering)
    curated_agg.drop(columns=["exclude", "masked_in_usher"], errors="ignore", inplace=True)
else:
    curated_agg = pd.DataFrame()
    log("WARNING: No data remaining after curation")

# Save curated data
curated_agg.to_csv(output_csv, index=False)

log(f"Curated {len(curated_agg)} mutations from {len(counts_df)} original entries")
log(f"Wrote curated data to {output_csv}")
log("Done!")
