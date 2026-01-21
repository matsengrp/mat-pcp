#!/usr/bin/env python3
"""Validate PCP-derived counts match aamut_fitness actual_count (S gene).

This validation compares two independent mutation counting approaches:
1. PCP-derived counts: Mutations extracted from parent-child pairs
2. aamut_fitness actual_count: Aggregated from mutation counting pipeline

For amino acid mutations in the S gene (or configured region), both should
produce identical counts.
"""

import pandas as pd

# Get snakemake variables
pcp_derived_counts_csv = snakemake.input.pcp_derived_counts
aamut_fitness_files = snakemake.input.aamut_fitness_files
coding_sites_csv = snakemake.input.coding_sites
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


log("Validating PCP-derived counts against aamut_fitness actual_count")

# Load PCP-derived counts (nucleotide mutations)
pcp_counts = pd.read_csv(pcp_derived_counts_csv)
log(f"Loaded {len(pcp_counts)} PCP-derived count entries")

# Load coding sites for nt->aa mapping
coding_sites = pd.read_csv(coding_sites_csv)
log(f"Loaded {len(coding_sites)} coding site entries")

# Load and combine aamut_fitness files
aamut_dfs = []
for f in aamut_fitness_files:
    df = pd.read_csv(f)
    aamut_dfs.append(df)
aamut_fitness = pd.concat(aamut_dfs, ignore_index=True)
log(f"Loaded {len(aamut_fitness)} aamut_fitness entries from {len(aamut_fitness_files)} files")

# Filter to PCP region if specified
target_gene = None
if pcp_gene_region:
    start = pcp_gene_region.get("start", 1)
    end = pcp_gene_region.get("end", float("inf"))
    target_gene = pcp_gene_region.get("gene", None)
    log(f"Filtering to PCP region: {pcp_gene_region}")

    # Filter PCP counts to region
    pcp_counts = pcp_counts[
        (pcp_counts["nt_site"] >= start) & (pcp_counts["nt_site"] <= end)
    ]
    log(f"  PCP counts in region: {len(pcp_counts)}")

    # Filter aamut_fitness to gene if specified
    if target_gene:
        aamut_fitness = aamut_fitness[aamut_fitness["gene"] == target_gene]
        log(f"  aamut_fitness entries for {target_gene}: {len(aamut_fitness)}")

# Create mapping from site to codon info
# Note: coding_sites uses "site" column
coding_site_map = coding_sites.set_index("site")[
    ["gene", "codon_site", "codon_position"]
].to_dict("index")


def get_codon_site_and_gene(nt_site, coding_site_map, target_gene=None):
    """Get the codon site and gene for a nucleotide site."""
    if nt_site not in coding_site_map:
        return None, None

    info = coding_site_map[nt_site]
    gene = info["gene"]
    codon_site = info["codon_site"]

    if pd.isna(gene) or gene == "noncoding":
        return None, None

    # Handle overlapping genes (semicolon-separated)
    if ";" in str(gene):
        genes = str(gene).split(";")
        codon_sites = str(codon_site).split(";")

        # If target gene specified, try to find it
        if target_gene and target_gene in genes:
            idx = genes.index(target_gene)
            return int(codon_sites[idx]), target_gene
        else:
            # Take first gene
            return int(codon_sites[0]), genes[0]

    return int(codon_site), gene


# Aggregate PCP counts by clade and codon site (amino acid position)
# This doesn't require knowing the actual AA - just the position
pcp_aa_counts_list = []
for _, row in pcp_counts.iterrows():
    if row["clade"] == "all":
        continue  # Skip aggregate counts

    nt_site = row["nt_site"]
    pcp_count = row["pcp_count"]
    clade = row["clade"]

    codon_site, gene = get_codon_site_and_gene(nt_site, coding_site_map, target_gene)
    if codon_site is not None:
        pcp_aa_counts_list.append({
            "clade": clade,
            "gene": gene,
            "aa_site": codon_site,
            "pcp_count": pcp_count,
        })

pcp_aa_df = pd.DataFrame(pcp_aa_counts_list)
log(f"Converted to {len(pcp_aa_df)} amino acid site entries")

# Aggregate PCP counts by clade, gene, and aa_site
if len(pcp_aa_df) > 0:
    pcp_aa_counts = (
        pcp_aa_df.groupby(["clade", "gene", "aa_site"], as_index=False)
        .agg(pcp_count=("pcp_count", "sum"))
    )
else:
    pcp_aa_counts = pd.DataFrame(columns=["clade", "gene", "aa_site", "pcp_count"])

log(f"Aggregated to {len(pcp_aa_counts)} unique (clade, gene, aa_site) entries")

# Aggregate aamut_fitness actual_count by clade, gene, and aa_site
# (summing all AA changes at each position)
aamut_by_site = (
    aamut_fitness.groupby(["cluster", "gene", "aa_site"], as_index=False)
    .agg(actual_count=("actual_count", "sum"))
    .rename(columns={"cluster": "clade"})
)
log(f"Aggregated aamut_fitness to {len(aamut_by_site)} unique site entries")

# Merge PCP counts with aamut_fitness actual_count
merged = pcp_aa_counts.merge(
    aamut_by_site,
    on=["clade", "gene", "aa_site"],
    how="outer",
)

# Fill NaN with 0
merged["pcp_count"] = merged["pcp_count"].fillna(0).astype(int)
merged["actual_count"] = merged["actual_count"].fillna(0).astype(int)

# Find discrepancies
merged["difference"] = merged["pcp_count"] - merged["actual_count"]
discrepancies = merged[merged["difference"] != 0].copy()

# Summary statistics
n_total = len(merged)
n_exact_matches = (merged["difference"] == 0).sum()
n_discrepancies = len(discrepancies)

# Categorize discrepancies
only_in_pcp = discrepancies[discrepancies["actual_count"] == 0]
only_in_aamut = discrepancies[discrepancies["pcp_count"] == 0]
count_mismatch = discrepancies[
    (discrepancies["pcp_count"] > 0) & (discrepancies["actual_count"] > 0)
]

log(f"Total comparisons: {n_total}")
log(f"Exact matches: {n_exact_matches}")
log(f"Discrepancies: {n_discrepancies}")

# Write report
with open(output_report, "w") as f:
    f.write("=" * 60 + "\n")
    f.write("PCP vs aamut_fitness Validation Report\n")
    f.write("=" * 60 + "\n\n")

    f.write("SUMMARY\n")
    f.write("-" * 40 + "\n")
    f.write(f"PCP-derived site entries: {len(pcp_aa_counts)}\n")
    f.write(f"aamut_fitness site entries with actual_count > 0: {(aamut_by_site['actual_count'] > 0).sum()}\n")
    f.write(f"Total site comparisons: {n_total}\n")
    f.write(f"Exact matches: {n_exact_matches}\n")
    f.write(f"Total discrepancies: {n_discrepancies}\n\n")

    if n_discrepancies == 0:
        f.write("SUCCESS: All mutation counts match exactly!\n\n")
        f.write("The PCP extraction and mutation counting pipelines are consistent.\n")
    else:
        f.write("DISCREPANCY BREAKDOWN\n")
        f.write("-" * 40 + "\n")
        f.write(f"Only in PCPs (not in aamut_fitness): {len(only_in_pcp)}\n")
        f.write(f"Only in aamut_fitness (not in PCPs): {len(only_in_aamut)}\n")
        f.write(f"Count mismatches (both have counts): {len(count_mismatch)}\n\n")

        if len(count_mismatch) > 0:
            f.write("WARNING: Count mismatches indicate filtering differences.\n")
            f.write("Review filtering criteria in both pipelines.\n\n")

            # Show worst mismatches
            worst = count_mismatch.nlargest(10, "difference", keep="first")
            f.write("Top 10 count mismatches:\n")
            f.write(
                worst[
                    ["clade", "gene", "aa_site", "pcp_count", "actual_count", "difference"]
                ].to_string()
            )
            f.write("\n")

        if len(only_in_pcp) > 0:
            f.write("\nSites only in PCPs (first 10):\n")
            f.write(
                only_in_pcp.head(10)[["clade", "gene", "aa_site", "pcp_count"]].to_string()
            )
            f.write("\n")

        if len(only_in_aamut) > 0 and len(only_in_aamut) <= 20:
            f.write("\nSites only in aamut_fitness (first 10):\n")
            f.write(
                only_in_aamut.head(10)[["clade", "gene", "aa_site", "actual_count"]].to_string()
            )
            f.write("\n")

    # Per-clade summary
    f.write("\nPER-CLADE SUMMARY\n")
    f.write("-" * 40 + "\n")
    for clade in sorted(merged["clade"].unique()):
        clade_data = merged[merged["clade"] == clade]
        clade_disc = discrepancies[discrepancies["clade"] == clade]
        clade_pcp_total = clade_data["pcp_count"].sum()
        clade_actual_total = clade_data["actual_count"].sum()
        f.write(f"{clade}: {len(clade_data)} sites, {len(clade_disc)} discrepancies, "
                f"pcp_total={clade_pcp_total}, actual_total={clade_actual_total}\n")

    f.write("\n" + "=" * 60 + "\n")
    f.write("End of Report\n")
    f.write("=" * 60 + "\n")

# Write discrepancies to CSV
discrepancies.to_csv(output_discrepancies, index=False)

log(f"Wrote report to {output_report}")
log(f"Wrote {len(discrepancies)} discrepancies to {output_discrepancies}")
log("Validation complete!")
