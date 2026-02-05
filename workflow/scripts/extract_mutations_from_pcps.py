#!/usr/bin/env python3
"""Extract amino acid mutation counts from PCPs by translating sequences.

Translates parent and child nucleotide sequences to amino acids, then
compares them to identify amino acid mutations (e.g., D614G).
"""

import pandas as pd
import Bio.SeqIO
from Bio.Seq import Seq

# Get snakemake variables
pcps_csv = snakemake.input.pcps
output_aa_csv = snakemake.output.aa_counts

pcp_gene_region = snakemake.params.pcp_gene_region
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log("Extracting AA mutations from PCPs by translating sequences")

# Load inputs
pcps = pd.read_csv(pcps_csv)

log(f"Loaded {len(pcps)} PCPs")

# Get gene info from pcp_gene_region
gene_name = None
if pcp_gene_region:
    gene_name = pcp_gene_region.get("gene", "S")
    offset = pcp_gene_region.get("start", 1)
    log(f"PCP gene region: {pcp_gene_region}")
    log(f"Gene name: {gene_name}")
else:
    offset = 1
    log("No PCP gene region specified")


def translate_codon(codon):
    """Translate a single codon to amino acid."""
    if len(codon) != 3:
        return "X"
    try:
        return str(Seq(codon).translate())
    except Exception:
        return "X"


def extract_aa_mutations_from_codons(parent_nt, child_nt):
    """
    Extract amino acid mutations by comparing codons and translating.

    This approach identifies mutations at the codon level first, then translates
    to determine amino acid changes. This captures ALL mutations including
    synonymous ones (e.g., D614D where the codon changed but AA stayed the same).

    Returns:
        List of dicts with aa_site (1-indexed), parent_aa, child_aa, aa_mutation
    """
    mutations = []

    # Ensure lengths are multiples of 3
    trim_len = min(len(parent_nt), len(child_nt))
    trim_len = (trim_len // 3) * 3
    parent_nt = parent_nt[:trim_len]
    child_nt = child_nt[:trim_len]

    # Compare codons
    for i in range(0, trim_len, 3):
        p_codon = parent_nt[i:i+3]
        c_codon = child_nt[i:i+3]

        # Only record if codons differ
        if p_codon != c_codon:
            # Translate both codons
            p_aa = translate_codon(p_codon)
            c_aa = translate_codon(c_codon)

            # Skip if translation failed or contains stop codons
            if p_aa == "*" or c_aa == "*" or p_aa == "X" or c_aa == "X":
                continue

            aa_site = (i // 3) + 1  # 1-indexed
            mutations.append({
                "aa_site": aa_site,
                "parent_aa": p_aa,
                "child_aa": c_aa,
                "aa_mutation": f"{p_aa}{aa_site}{c_aa}",
            })

    return mutations


# Extract AA mutations from each PCP
all_aa_mutations = []
total_pairs = 0

for _, row in pcps.iterrows():
    total_pairs += 1
    parent_nt = row.get("parent", "")
    child_nt = row.get("child", "")
    clade = row.get("clade", "unknown")

    if pd.isna(parent_nt) or pd.isna(child_nt):
        continue

    # Extract AA mutations at codon level (includes synonymous mutations)
    aa_mutations = extract_aa_mutations_from_codons(parent_nt, child_nt)

    for mut in aa_mutations:
        all_aa_mutations.append({
            "clade": clade,
            "gene": gene_name,
            "aa_site": mut["aa_site"],
            "aa_mutation": mut["aa_mutation"],
            "clade_founder_aa": mut["parent_aa"],
            "mutant_aa": mut["child_aa"],
        })

log(f"Extracted {len(all_aa_mutations)} AA mutations from {total_pairs} PCPs")

# Create output dataframe with counts
if all_aa_mutations:
    mut_df = pd.DataFrame(all_aa_mutations)

    # Per-clade counts by aa_mutation
    aa_clade_counts = (
        mut_df.groupby(["clade", "gene", "aa_site", "aa_mutation", "clade_founder_aa", "mutant_aa"])
        .size()
        .reset_index(name="pcp_count")
    )

    # Aggregate counts (across all clades)
    aa_agg_counts = (
        mut_df.groupby(["gene", "aa_site", "aa_mutation", "clade_founder_aa", "mutant_aa"])
        .size()
        .reset_index(name="pcp_count")
    )
    aa_agg_counts["clade"] = "all"

    aa_counts = pd.concat([aa_clade_counts, aa_agg_counts], ignore_index=True)
else:
    aa_counts = pd.DataFrame(columns=["clade", "gene", "aa_site", "aa_mutation", "clade_founder_aa", "mutant_aa", "pcp_count"])

# Write output
aa_counts.to_csv(output_aa_csv, index=False)

log(f"Wrote {len(aa_counts)} AA mutation count entries to {output_aa_csv}")
log("Done!")
