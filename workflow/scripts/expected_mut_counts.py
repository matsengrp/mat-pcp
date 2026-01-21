#!/usr/bin/env python3
"""Compute expected mutation counts from synonymous rates.

Uses the same normalization approach as SARS2-mut-fitness:

For each clade, we compute T that satisfies:
    N_s = T * sum_nt1(s_nt1 * sum_nt2(r_nt1->nt2))

where:
- N_s = total synonymous mutations at 4-fold sites (total_count)
- s_nt = number of 4-fold sites with nucleotide nt
- r_nt1->nt2 = rate of mutation from nt1 to nt2

The expected count for each site is then: T * rate
"""

import numpy as np
import pandas as pd

# Get snakemake variables
rates_by_clade_csv = snakemake.input.rates_by_clade
clade_founder_nts_csv = snakemake.input.clade_founder_nts_csv
output_csv = snakemake.output.expected_counts
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


# Standard genetic code
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate_codon(codon):
    """Translate a codon to amino acid."""
    if codon == "noncoding" or pd.isna(codon) or len(str(codon)) != 3:
        return "noncoding"
    return CODON_TABLE.get(codon.upper(), "X")


def compute_mutant_codon(founder_codon, codon_position, to_nt):
    """Compute mutant codon given a nucleotide change at a specific position.

    Args:
        founder_codon: The original codon (e.g., 'ATG')
        codon_position: 1-indexed position within codon (1, 2, or 3)
        to_nt: The nucleotide to mutate to (A, C, G, or T)
    """
    if founder_codon == "noncoding" or pd.isna(founder_codon):
        return "noncoding"
    # codon_position is 1-indexed in the data, convert to 0-indexed
    pos = codon_position - 1 if isinstance(codon_position, int) else -1
    if pos < 0 or pos > 2:
        return "noncoding"
    codon_list = list(founder_codon)
    codon_list[pos] = to_nt
    return "".join(codon_list)


log("Computing expected mutation counts from synonymous rates (SARS2-mut-fitness compatible)")

# Load inputs
rates_by_clade = pd.read_csv(rates_by_clade_csv)
clade_founder_nts = pd.read_csv(clade_founder_nts_csv)

log(f"Loaded {len(rates_by_clade)} rate entries")
log(f"Loaded {len(clade_founder_nts)} clade founder entries")

# Check if rates file has new format (with total_count column)
if "total_count" not in rates_by_clade.columns:
    log("WARNING: Using old rate format without total_count - expected counts may differ")
    # Old format: rate = count/n_sites, so expected_count = rate directly
    use_old_format = True
else:
    use_old_format = False

# Compute normalization factor T for each clade (SARS2-mut-fitness approach)
# T = total_count / sum_nt1(s_nt1 * sum_nt2(r_nt1->nt2))
if not use_old_format:
    # Get number of 4-fold sites per nucleotide per clade
    sites_by_nt = (
        clade_founder_nts.query("four_fold_degenerate == True")
        .groupby(["clade", "nt"], as_index=False)
        .agg(sites_4fold=pd.NamedAgg("site", "nunique"))
        .rename(columns={"nt": "parent_nt"})
    )

    # Compute normalization factor
    normalization_factor = (
        sites_by_nt.merge(rates_by_clade, on=["clade", "parent_nt"], validate="one_to_many")
        .groupby(["clade", "parent_nt", "total_count", "sites_4fold"], as_index=False)
        .agg(rate_sum=pd.NamedAgg("rate", "sum"))
        .assign(s_times_rate_sum=lambda x: x["sites_4fold"] * x["rate_sum"])
        .groupby(["clade", "total_count"], as_index=False)
        .agg(summand=pd.NamedAgg("s_times_rate_sum", "sum"))
        .assign(normalization_factor=lambda x: x["total_count"] / x["summand"])
        .drop(columns="summand")
    )

    log(f"Computed normalization factors for {len(normalization_factor)} clades")
    for _, row in normalization_factor.iterrows():
        log(f"  {row['clade']}: T = {row['normalization_factor']:.4f}")

# Create all possible mutations for each site using vectorized operations
# Expand each site to 3 mutations (to each other nucleotide)
nts = ["A", "C", "G", "T"]
mutation_dfs = []

for to_nt in nts:
    # Filter to sites where the founder nt differs from to_nt
    subset = clade_founder_nts[clade_founder_nts["nt"] != to_nt].copy()
    subset["to_nt"] = to_nt
    # Use "XtoY" format for rate lookup (SARS2-mut-fitness compatible)
    subset["mut_type_long"] = subset["nt"] + "to" + to_nt
    # Use "XY" format for output (mat-pcp GLM compatible)
    subset["mut_type"] = subset["nt"] + to_nt
    subset["nt_mutation"] = subset["nt"] + subset["site"].astype(str) + to_nt
    mutation_dfs.append(subset)

# Combine all mutations
all_mutations = pd.concat(mutation_dfs, ignore_index=True)
log(f"Generated {len(all_mutations)} mutation possibilities")

# Compute amino acid information for each mutation
log("Computing amino acid translations...")


def compute_aa_info(row):
    """Compute founder_aa, mutant_aa, aa_mutation, and synonymous for a row."""
    founder_codon = str(row["codon"])
    codon_position = str(row["codon_position"])
    codon_site = str(row["codon_site"])
    to_nt = row["to_nt"]

    # Handle noncoding
    if founder_codon == "noncoding" or pd.isna(row["codon"]):
        return pd.Series(
            {
                "clade_founder_aa": "noncoding",
                "mutant_aa": "noncoding",
                "aa_mutation": "noncoding",
                "synonymous": False,
                "noncoding": True,
            }
        )

    # Handle overlapping genes (semicolon-separated values)
    # Process each gene separately and combine results
    codons = founder_codon.split(";")
    positions = codon_position.split(";")
    sites = codon_site.split(";")

    founder_aas = []
    mutant_aas = []
    aa_mutations = []
    synonymous_flags = []

    for i, (codon, pos, site) in enumerate(zip(codons, positions, sites)):
        if codon == "noncoding" or pos == "noncoding":
            founder_aas.append("noncoding")
            mutant_aas.append("noncoding")
            aa_mutations.append("noncoding")
            synonymous_flags.append(False)
            continue

        # Compute mutant codon
        try:
            pos_int = int(pos)
        except ValueError:
            founder_aas.append("X")
            mutant_aas.append("X")
            aa_mutations.append("X")
            synonymous_flags.append(False)
            continue

        mutant_codon = compute_mutant_codon(codon, pos_int, to_nt)

        # Translate
        f_aa = translate_codon(codon)
        m_aa = translate_codon(mutant_codon)

        founder_aas.append(f_aa)
        mutant_aas.append(m_aa)
        aa_mutations.append(f"{f_aa}{site}{m_aa}")
        synonymous_flags.append(f_aa == m_aa)

    # Combine results with semicolons for overlapping genes
    return pd.Series(
        {
            "clade_founder_aa": ";".join(founder_aas),
            "mutant_aa": ";".join(mutant_aas),
            "aa_mutation": ";".join(aa_mutations),
            "synonymous": all(synonymous_flags),  # Only synonymous if all are
            "noncoding": False,
        }
    )


# Apply AA computation
aa_info = all_mutations.apply(compute_aa_info, axis=1)
all_mutations = pd.concat([all_mutations, aa_info], axis=1)
log(f"Computed AA info: {(~all_mutations['noncoding']).sum()} coding mutations")

# Merge with rates to get expected counts
# Create rate lookup as a dataframe for efficient merge
# Rates file uses "AtoC" format, so we match on mut_type_long
rate_lookup = rates_by_clade[rates_by_clade["rate"].notna()][
    ["clade", "mut_type", "rate"]
].copy()

all_mutations = all_mutations.merge(
    rate_lookup.rename(columns={"mut_type": "mut_type_long"}),
    on=["clade", "mut_type_long"],
    how="left",
)

# Compute expected count = T * rate (SARS2-mut-fitness approach)
if not use_old_format:
    all_mutations = all_mutations.merge(
        normalization_factor[["clade", "normalization_factor"]],
        on="clade",
        how="left",
    )
    all_mutations["expected_count"] = (
        all_mutations["normalization_factor"].fillna(0)
        * all_mutations["rate"].fillna(0)
    )
else:
    # Old format: expected_count = rate directly
    all_mutations["expected_count"] = all_mutations["rate"].fillna(0)

# Select and rename columns for output
expected_df = all_mutations[
    [
        "clade",
        "site",
        "nt_mutation",
        "mut_type",
        "expected_count",
        "gene",
        "codon_site",
        "codon_position",
        "clade_founder_aa",
        "mutant_aa",
        "aa_mutation",
        "synonymous",
        "noncoding",
    ]
].rename(columns={"site": "nt_site"})

# Log summary per clade
for clade in expected_df["clade"].unique():
    clade_data = expected_df[expected_df["clade"] == clade]
    n_with_rate = (clade_data["expected_count"] > 0).sum()
    log(f"  {clade}: {len(clade_data)} mutations, {n_with_rate} with rates")

expected_df.to_csv(output_csv, index=False)

log(f"Computed expected counts for {len(expected_df)} mutation possibilities")
log("Done!")
