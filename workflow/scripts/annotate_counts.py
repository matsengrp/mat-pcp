#!/usr/bin/env python3
"""Annotate counts with RNA structure, sequence context, and other features."""

import pandas as pd
import numpy as np

# Get snakemake variables
counts_csv = snakemake.input.counts
founder_csv = snakemake.input.clade_founder
rna_struct_file = snakemake.input.rna_struct
output_counts_csv = snakemake.output.counts_csv
output_founder_csv = snakemake.output.founder_csv
lightswitch_boundaries = snakemake.params.lightswitch_boundaries
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log("Annotating counts with RNA structure and sequence context")

# Load inputs
counts_df = pd.read_csv(counts_csv)
founder_df = pd.read_csv(founder_csv)

log(f"Loaded {len(counts_df)} count entries")
log(f"Loaded {len(founder_df)} founder entries")

# Build founder sequence dictionary
# Note: founder_df uses 'site' column, counts_df uses 'nt_site' column
founder_df = founder_df.sort_values(["clade", "site"])
founder_seq_dict = {}
for clade, data in founder_df.groupby("clade"):
    founder_seq_dict[clade] = "".join(data["nt"])
    log(f"  {clade}: {len(data)} sites")


def get_motif(site, clade, founder_seq_dict):
    """Get 3-mer motif centered on site."""
    if clade not in founder_seq_dict:
        return np.nan
    seq = founder_seq_dict[clade]
    if site < 2 or site > len(seq):
        return np.nan
    # Site is 1-indexed, so site-2 gives left neighbor, site-1 gives center
    return seq[site - 2 : site + 1]


# Add motif column to founder data
min_site = founder_df["site"].min()
max_site = founder_df["site"].max()
founder_df["motif"] = founder_df.apply(
    lambda row: (
        np.nan
        if row["site"] in [min_site, max_site]
        else get_motif(row["site"], row["clade"], founder_seq_dict)
    ),
    axis=1,
)

# Add reference codon and motif columns
# Get reference clade from config (default to 19A for backward compatibility)
reference_clade = snakemake.params.get("reference_clade", "19A")

# Validate that reference clade exists in the dataset
available_clades = sorted(founder_df["clade"].unique())
if reference_clade not in available_clades:
    raise ValueError(
        f"Reference clade '{reference_clade}' not found in dataset. "
        f"Available clades: {available_clades}. "
        f"Please update 'reference_clade' in config.yaml to a valid clade."
    )

# Extract reference founder and merge
ref_founder = founder_df[founder_df["clade"] == reference_clade]
log(f"Using clade '{reference_clade}' as reference for motif/codon annotation")

founder_df = founder_df.merge(
    ref_founder[["site", "codon", "motif"]].rename(
        columns={"codon": "ref_codon", "motif": "ref_motif"}
    ),
    on="site",
    how="left",
)

# Rename site to nt_site to match counts_df
founder_df.rename(columns={"site": "nt_site"}, inplace=True)

log("Added motif annotations to founder data")


def get_mut_class(row):
    """Determine mutation classification."""
    if row.get("synonymous", False):
        return "synonymous"
    elif row.get("noncoding", False):
        return "noncoding"
    elif "*" in str(row.get("mutant_aa", "")):
        return "nonsense"
    else:
        return "nonsynonymous"


# Add mutation classification if synonymous column exists
if "synonymous" in counts_df.columns or "noncoding" in counts_df.columns:
    counts_df["mut_class"] = counts_df.apply(get_mut_class, axis=1)
    log(f"Mutation classes: {counts_df['mut_class'].value_counts().to_dict()}")
else:
    # If no synonymous info, mark all as unknown
    counts_df["mut_class"] = "unknown"
    log("No synonymous column found, marking all as unknown")


def light_switch(mut_type, site, boundaries):
    """Determine if site is before light-switch boundary for mutation type."""
    if mut_type in ["AT", "CG", "GC"]:
        return site < boundaries.get("AT", 21562)
    elif mut_type == "CT":
        return site < boundaries.get("CT", 13467)
    return False


# Apply lightswitch boundaries
# counts_df already has mut_type column from earlier processing
counts_df["nt_site_before_boundary"] = counts_df.apply(
    lambda x: light_switch(x["mut_type"], x["nt_site"], lightswitch_boundaries), axis=1
)
log(f"Sites before boundary: {counts_df['nt_site_before_boundary'].sum()}")

# Add RNA secondary structure information
log(f"Loading RNA structure from {rna_struct_file}")
with open(rna_struct_file) as f:
    lines = [line.rstrip().split() for line in f]

# Parse structure file (format: site nucleotide i j k pairing_partner)
# Skip header line
paired = np.array([[int(x[0]), int(x[4])] for x in lines[1:]])
paired_dict = dict(zip(paired[:, 0], paired[:, 1]))
log(f"Loaded {len(paired_dict)} RNA structure entries")


def assign_ss_pred(site):
    """Assign secondary structure prediction."""
    if site not in paired_dict:
        return "nd"
    elif paired_dict[site] == 0:
        return "unpaired"
    else:
        return "paired"


counts_df["ss_prediction"] = counts_df["nt_site"].apply(assign_ss_pred)
counts_df["unpaired"] = counts_df["ss_prediction"].apply(
    lambda x: 1 if x == "unpaired" else 0
)

ss_counts = counts_df["ss_prediction"].value_counts()
log(f"Secondary structure: {ss_counts.to_dict()}")

# Merge with founder motif information
counts_df = counts_df.merge(
    founder_df[["nt_site", "clade", "motif", "ref_motif"]],
    on=["nt_site", "clade"],
    how="left",
)

# Handle edge cases for motif
counts_df.loc[counts_df.nt_site == 1, "motif"] = counts_df.loc[
    counts_df.nt_site == 1, "nt_mutation"
].apply(lambda m: f"A{m[0]}T" if pd.notna(m) and len(m) > 0 else "AAA")
max_site_in_counts = counts_df["nt_site"].max()
counts_df.loc[counts_df.nt_site == max_site_in_counts, "motif"] = "AAA"

log(f"Annotated {len(counts_df)} count entries")

# Save outputs
counts_df.to_csv(output_counts_csv, index=False)
founder_df.to_csv(output_founder_csv, index=False)

log(f"Wrote annotated counts to {output_counts_csv}")
log(f"Wrote annotated founder to {output_founder_csv}")
log("Done!")
