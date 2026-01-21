#!/usr/bin/env python3
"""Extract mutation counts from PCPs for validation against aggregated counts.

Outputs both nucleotide mutations and amino acid site counts for validation.
"""

import pandas as pd
import Bio.SeqIO

# Get snakemake variables
pcps_csv = snakemake.input.pcps
coding_sites_csv = snakemake.input.coding_sites
ref_fasta = snakemake.input.ref_fasta
output_nt_csv = snakemake.output.nt_counts
output_aa_csv = snakemake.output.aa_counts

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

# Create mapping from site to codon info
coding_site_map = coding_sites.set_index("site")[
    ["gene", "codon_site", "codon_position"]
].to_dict("index")


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


def get_aa_site_info(nt_site, coding_site_map, target_gene=None):
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


# Determine offset and target gene based on PCP region
target_gene = None
if pcp_gene_region:
    offset = pcp_gene_region.get("start", 1)
    target_gene = pcp_gene_region.get("gene", None)
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
        # Get amino acid site info
        aa_site, gene = get_aa_site_info(site, coding_site_map, target_gene)

        all_mutations.append(
            {
                "clade": clade,
                "nt_mutation": f"{from_nt}{site}{to_nt}",
                "nt_site": site,
                "from_nt": from_nt,
                "to_nt": to_nt,
                "gene": gene,
                "aa_site": aa_site,
            }
        )

log(f"Extracted {len(all_mutations)} total mutations from {total_pairs} PCPs")
log(f"Filtered {filtered_count} PCPs due to filtering criteria")

# Create output dataframes
if all_mutations:
    mut_df = pd.DataFrame(all_mutations)

    # === Nucleotide mutation counts ===
    # Per-clade counts
    nt_clade_counts = (
        mut_df.groupby(["clade", "nt_mutation", "nt_site"])
        .size()
        .reset_index(name="pcp_count")
    )

    # Aggregate counts
    nt_agg_counts = (
        mut_df.groupby(["nt_mutation", "nt_site"]).size().reset_index(name="pcp_count")
    )
    nt_agg_counts["clade"] = "all"

    nt_counts = pd.concat([nt_clade_counts, nt_agg_counts], ignore_index=True)

    # === Amino acid site counts ===
    # Filter to coding mutations only
    coding_mut_df = mut_df[mut_df["aa_site"].notna()].copy()

    if len(coding_mut_df) > 0:
        # Per-clade counts by aa_site
        aa_clade_counts = (
            coding_mut_df.groupby(["clade", "gene", "aa_site"])
            .size()
            .reset_index(name="pcp_count")
        )

        # Aggregate counts
        aa_agg_counts = (
            coding_mut_df.groupby(["gene", "aa_site"])
            .size()
            .reset_index(name="pcp_count")
        )
        aa_agg_counts["clade"] = "all"

        aa_counts = pd.concat([aa_clade_counts, aa_agg_counts], ignore_index=True)
    else:
        aa_counts = pd.DataFrame(columns=["clade", "gene", "aa_site", "pcp_count"])

else:
    nt_counts = pd.DataFrame(columns=["clade", "nt_mutation", "nt_site", "pcp_count"])
    aa_counts = pd.DataFrame(columns=["clade", "gene", "aa_site", "pcp_count"])

# Write outputs
nt_counts.to_csv(output_nt_csv, index=False)
aa_counts.to_csv(output_aa_csv, index=False)

log(f"Wrote {len(nt_counts)} nucleotide count entries to {output_nt_csv}")
log(f"Wrote {len(aa_counts)} amino acid site count entries to {output_aa_csv}")
log("Done!")
