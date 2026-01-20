#!/usr/bin/env python3
"""Extract mutation counts from PCPs for validation against aggregated counts."""

import pandas as pd
import Bio.SeqIO

# Get snakemake variables
pcps_csv = snakemake.input.pcps
coding_sites_csv = snakemake.input.coding_sites
ref_fasta = snakemake.input.ref_fasta
output_csv = snakemake.output.csv

pcp_gene_region = snakemake.params.pcp_gene_region
max_nt_mutations = snakemake.params.max_nt_mutations
max_reversions_to_ref = snakemake.params.max_reversions_to_ref
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log("Extracting mutations from PCPs for validation")

# Load inputs
pcps = pd.read_csv(pcps_csv)
coding_sites = pd.read_csv(coding_sites_csv)
ref_seq = str(Bio.SeqIO.read(ref_fasta, "fasta").seq)

log(f"Loaded {len(pcps)} PCPs")
log(f"Reference sequence length: {len(ref_seq)}")


def extract_mutations(parent, child, offset=1):
    """
    Extract mutations from parent-child sequence pair.

    Args:
        parent: Parent sequence string
        child: Child sequence string
        offset: Genomic position of first character (1-indexed)

    Returns:
        List of (from_nt, site, to_nt) tuples
    """
    mutations = []
    for i, (p, c) in enumerate(zip(parent, child)):
        if p != c and p in "ACGT" and c in "ACGT":
            site = offset + i
            mutations.append((p, site, c))
    return mutations


def count_reversions_to_ref(mutations, ref_seq):
    """Count mutations that revert to reference."""
    return sum(1 for (_, site, to_nt) in mutations if ref_seq[site - 1] == to_nt)


# Determine offset based on PCP region
if pcp_gene_region:
    offset = pcp_gene_region.get("start", 1)
    log(f"PCP gene region: {pcp_gene_region}")
else:
    offset = 1
    log("No PCP gene region specified, using full genome")

# Extract mutations from each PCP
all_mutations = []
filtered_count = 0
total_pairs = 0

for _, row in pcps.iterrows():
    total_pairs += 1
    parent = row.get("parent", "")
    child = row.get("child", "")
    clade = row.get("clade", "unknown")

    if pd.isna(parent) or pd.isna(child):
        continue

    mutations = extract_mutations(parent, child, offset)

    # Apply same filtering criteria as count_mutations
    n_mutations = len(mutations)
    n_reversions = count_reversions_to_ref(mutations, ref_seq)

    if n_mutations > max_nt_mutations:
        filtered_count += 1
        continue
    if n_reversions > max_reversions_to_ref:
        filtered_count += 1
        continue

    for from_nt, site, to_nt in mutations:
        all_mutations.append(
            {
                "clade": clade,
                "nt_mutation": f"{from_nt}{site}{to_nt}",
                "nt_site": site,
                "from_nt": from_nt,
                "to_nt": to_nt,
            }
        )

log(f"Extracted {len(all_mutations)} total mutations from {total_pairs} PCPs")
log(f"Filtered {filtered_count} PCPs due to filtering criteria")

# Count mutations by clade and mutation
if all_mutations:
    mut_df = pd.DataFrame(all_mutations)

    # Per-clade counts
    clade_counts = (
        mut_df.groupby(["clade", "nt_mutation", "nt_site"])
        .size()
        .reset_index(name="pcp_count")
    )

    # Also compute aggregate counts
    agg_counts = (
        mut_df.groupby(["nt_mutation", "nt_site"]).size().reset_index(name="pcp_count")
    )
    agg_counts["clade"] = "all"

    mutation_counts = pd.concat([clade_counts, agg_counts], ignore_index=True)
else:
    mutation_counts = pd.DataFrame(
        columns=["clade", "nt_mutation", "nt_site", "pcp_count"]
    )

mutation_counts.to_csv(output_csv, index=False)

log(f"Wrote {len(mutation_counts)} count entries to {output_csv}")
log("Done!")
