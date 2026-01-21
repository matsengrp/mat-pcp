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
max_nt_mutations = snakemake.params.max_nt_mutations
max_reversions_to_ref = snakemake.params.max_reversions_to_ref
ref_fasta = snakemake.input.ref_fasta
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
ref_seq = str(Bio.SeqIO.read(ref_fasta, "fasta").seq)

log(f"Loaded {len(pcps)} PCPs")
log(f"Reference sequence length: {len(ref_seq)}")

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


def translate_sequence(nt_seq):
    """Translate nucleotide sequence to amino acids."""
    # Ensure length is multiple of 3
    trim_len = len(nt_seq) - (len(nt_seq) % 3)
    nt_seq = nt_seq[:trim_len]
    try:
        return str(Seq(nt_seq).translate())
    except Exception:
        return ""


def extract_aa_mutations(parent_aa, child_aa):
    """
    Extract amino acid mutations by comparing parent and child AA sequences.

    Returns:
        List of dicts with aa_site (1-indexed), parent_aa, child_aa, aa_mutation
    """
    mutations = []
    for i, (p_aa, c_aa) in enumerate(zip(parent_aa, child_aa)):
        if p_aa != c_aa and p_aa != "*" and c_aa != "*" and p_aa != "X" and c_aa != "X":
            aa_site = i + 1  # 1-indexed
            mutations.append({
                "aa_site": aa_site,
                "parent_aa": p_aa,
                "child_aa": c_aa,
                "aa_mutation": f"{p_aa}{aa_site}{c_aa}",
            })
    return mutations


def count_nt_mutations(parent_nt, child_nt):
    """Count nucleotide differences between parent and child."""
    count = 0
    for p, c in zip(parent_nt, child_nt):
        if p != c and p in "ACGT" and c in "ACGT":
            count += 1
    return count


def count_reversions_to_ref(parent_nt, child_nt, ref_nt, offset):
    """Count mutations that revert to reference."""
    count = 0
    for i, (p, c) in enumerate(zip(parent_nt, child_nt)):
        if p != c and p in "ACGT" and c in "ACGT":
            ref_pos = offset - 1 + i  # 0-indexed position in ref
            if ref_pos < len(ref_nt) and ref_nt[ref_pos] == c:
                count += 1
    return count


# Extract AA mutations from each PCP
all_aa_mutations = []
filtered_count = 0
total_pairs = 0

for _, row in pcps.iterrows():
    total_pairs += 1
    parent_nt = row.get("parent", "")
    child_nt = row.get("child", "")
    clade = row.get("clade", "unknown")

    if pd.isna(parent_nt) or pd.isna(child_nt):
        continue

    # Apply filtering based on nucleotide mutations
    n_nt_mutations = count_nt_mutations(parent_nt, child_nt)
    n_reversions = count_reversions_to_ref(parent_nt, child_nt, ref_seq, offset)

    if n_nt_mutations > max_nt_mutations:
        filtered_count += 1
        continue
    if n_reversions > max_reversions_to_ref:
        filtered_count += 1
        continue

    # Translate sequences
    parent_aa = translate_sequence(parent_nt)
    child_aa = translate_sequence(child_nt)

    if not parent_aa or not child_aa:
        continue

    # Extract AA mutations
    aa_mutations = extract_aa_mutations(parent_aa, child_aa)

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
log(f"Filtered {filtered_count} PCPs due to filtering criteria")

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
