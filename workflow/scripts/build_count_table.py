#!/usr/bin/env python3
"""Build the unified per-clade mutation count table.

This fuses what used to be three rules (``expected_mut_counts``,
``merge_expected_and_actual_counts``, ``aggregate_mutations_to_exclude``) into a
single step, and drops the synonymous-rate ``expected_count`` entirely:

1. Enumerate the full grid of all possible mutations per clade (from the clade
   founder nucleotides), with amino-acid / codon annotations -- this is what
   supplies the ``actual_count == 0`` rows the GLM trains on.
2. Join the aggregated actual counts onto the grid (fill 0 where unobserved).
3. Compute ``exclude`` / ``masked_in_usher`` once, via
   ``modules/masking.compute_exclusions``.

Output ``counts_all.csv`` replaces the old ``expected_vs_actual_mut_counts.csv``
(minus ``expected_count``); the GLM neutral expectation is estimated downstream.
"""

import os
from importlib import util as importlib_util

import Bio.SeqIO
import pandas as pd

# Snakemake inputs / params
actual_csv = snakemake.input.actual
clade_founder_nts_csv = snakemake.input.clade_founder_nts_csv
ref_fasta = snakemake.input.ref_fasta
site_mask_csv = snakemake.input.site_mask
usher_masked_sites_yaml = snakemake.input.usher_masked_sites
output_csv = snakemake.output.csv

sites_to_exclude = snakemake.params.sites_to_exclude
site_include_range = snakemake.params.site_include_range

log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


# The by-path module loader is defined once in
# workflow/modules/workflow_loader.py. Because that loader is itself the
# mechanism for importing shared workflow code, it cannot be imported via
# itself -- it is wired up here by absolute path (#111).
_loader_spec = importlib_util.spec_from_file_location(
    "workflow_loader",
    os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "modules",
        "workflow_loader.py",
    ),
)
workflow_loader = importlib_util.module_from_spec(_loader_spec)
_loader_spec.loader.exec_module(workflow_loader)
load_workflow_module = workflow_loader.load_workflow_module


masking = load_workflow_module("masking")


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
    """Compute mutant codon given a nucleotide change at a codon position (1-indexed)."""
    if founder_codon == "noncoding" or pd.isna(founder_codon):
        return "noncoding"
    pos = codon_position - 1 if isinstance(codon_position, int) else -1
    if pos < 0 or pos > 2:
        return "noncoding"
    codon_list = list(founder_codon)
    codon_list[pos] = to_nt
    return "".join(codon_list)


def compute_aa_info(row):
    """Compute founder_aa, mutant_aa, aa_mutation, mutant_codon, synonymous for a row."""
    founder_codon = str(row["codon"])
    codon_position = str(row["codon_position"])
    codon_site = str(row["codon_site"])
    to_nt = row["to_nt"]

    if founder_codon == "noncoding" or pd.isna(row["codon"]):
        return pd.Series(
            {
                "clade_founder_aa": "noncoding",
                "mutant_aa": "noncoding",
                "aa_mutation": "noncoding",
                "mutant_codon": "noncoding",
                "synonymous": False,
                "noncoding": True,
            }
        )

    # Handle overlapping genes (semicolon-separated values)
    codons = founder_codon.split(";")
    positions = codon_position.split(";")
    sites = codon_site.split(";")

    founder_aas, mutant_aas, aa_mutations, mutant_codons, synonymous_flags = (
        [], [], [], [], [])

    for codon, pos, site in zip(codons, positions, sites):
        if codon == "noncoding" or pos == "noncoding":
            founder_aas.append("noncoding")
            mutant_aas.append("noncoding")
            aa_mutations.append("noncoding")
            mutant_codons.append("noncoding")
            synonymous_flags.append(False)
            continue

        try:
            pos_int = int(pos)
        except ValueError:
            founder_aas.append("X")
            mutant_aas.append("X")
            aa_mutations.append("X")
            mutant_codons.append("X")
            synonymous_flags.append(False)
            continue

        m_codon = compute_mutant_codon(codon, pos_int, to_nt)
        f_aa = translate_codon(codon)
        m_aa = translate_codon(m_codon)

        founder_aas.append(f_aa)
        mutant_aas.append(m_aa)
        aa_mutations.append(f"{f_aa}{site}{m_aa}")
        mutant_codons.append(m_codon)
        synonymous_flags.append(f_aa == m_aa)

    return pd.Series(
        {
            "clade_founder_aa": ";".join(founder_aas),
            "mutant_aa": ";".join(mutant_aas),
            "aa_mutation": ";".join(aa_mutations),
            "mutant_codon": ";".join(mutant_codons),
            "synonymous": all(synonymous_flags),
            "noncoding": False,
        }
    )


log("Building unified count table")

# ---------------------------------------------------------------------------
# 1. Enumerate the full grid of all possible mutations, with AA annotations.
# ---------------------------------------------------------------------------
clade_founder_nts = pd.read_csv(clade_founder_nts_csv)
log(f"Loaded {len(clade_founder_nts)} clade founder nucleotide entries")

nts = ["A", "C", "G", "T"]
mutation_dfs = []
for to_nt in nts:
    subset = clade_founder_nts[clade_founder_nts["nt"] != to_nt].copy()
    subset["to_nt"] = to_nt
    subset["mut_type"] = subset["nt"] + to_nt
    subset["nt_mutation"] = subset["nt"] + subset["site"].astype(str) + to_nt
    mutation_dfs.append(subset)

grid = pd.concat(mutation_dfs, ignore_index=True)
log(f"Generated {len(grid)} mutation possibilities")

aa_info = grid.apply(compute_aa_info, axis=1)
grid = pd.concat([grid, aa_info], axis=1)
log(f"Computed AA info: {(~grid['noncoding']).sum()} coding mutations")

grid = grid.rename(
    columns={
        "site": "nt_site",
        "nt": "clade_founder_nt",
        "codon": "clade_founder_codon",
    }
)

# Reference (Wuhan-Hu-1) nucleotide per site, used to flag founder->reference
# reversions (grid rows whose destination base is the reference base). Mirrors
# the reference_nt computation in count_mutations.py.
ref_seq = str(Bio.SeqIO.read(ref_fasta, "fasta").seq)
grid["reference_nt"] = grid["nt_site"].map(lambda s: ref_seq[s - 1])
log(f"Loaded reference sequence: {len(ref_seq)} nt")

# ---------------------------------------------------------------------------
# 2. Join aggregated actual counts (fill 0 where unobserved).
# ---------------------------------------------------------------------------
actual = pd.read_csv(actual_csv)
count_col = "count" if "count" in actual.columns else "actual_count"
actual_cols = ["clade", "nt_mutation", count_col]
for extra in ("count_terminal", "count_non_terminal"):
    if extra in actual.columns:
        actual_cols.append(extra)

grid = grid.merge(
    actual[actual_cols].rename(
        columns={count_col: "actual_count"} if count_col != "actual_count" else {}
    ),
    on=["clade", "nt_mutation"],
    how="left",
)
grid["actual_count"] = grid["actual_count"].fillna(0).astype(int)
for extra in ("count_terminal", "count_non_terminal"):
    if extra in grid.columns:
        grid[extra] = grid[extra].fillna(0).astype(int)
log(f"Joined actual counts: {len(grid)} rows, {(grid['actual_count'] > 0).sum()} observed")

# ---------------------------------------------------------------------------
# 3. Compute exclude / masked_in_usher once, via the shared masking helper.
# ---------------------------------------------------------------------------
vcf_masked_sites = pd.read_csv(site_mask_csv)["site"].tolist()
usher_clade_sites, usher_clade_mutations = masking.load_usher_masks(usher_masked_sites_yaml)

exclude, masked_in_usher = masking.compute_exclusions(
    grid,
    sites_to_exclude=sites_to_exclude,
    vcf_masked_sites=vcf_masked_sites,
    usher_clade_sites=usher_clade_sites,
    usher_clade_mutations=usher_clade_mutations,
    site_include_range=site_include_range,
)
grid["exclude"] = exclude
grid["masked_in_usher"] = masked_in_usher
log(f"Excluded mutations: {grid['exclude'].sum()}  (masked_in_usher: {grid['masked_in_usher'].sum()})")

# ---------------------------------------------------------------------------
# 4. Add v2-compatible bookkeeping columns and write.
# ---------------------------------------------------------------------------
grid["subset"] = "all"
grid["mean_log_size"] = 0.0

output_columns = [
    "clade",
    "subset",
    "nt_site",
    "nt_mutation",
    "mut_type",
    "exclude",
    "masked_in_usher",
    "actual_count",
    "clade_founder_nt",
    "gene",
    "clade_founder_codon",
    "clade_founder_aa",
    "mutant_codon",
    "mutant_aa",
    "aa_mutation",
    "synonymous",
    "noncoding",
    "codon_position",
    "codon_site",
    "four_fold_degenerate",
    "count_terminal",
    "count_non_terminal",
    "mean_log_size",
]
available_cols = [c for c in output_columns if c in grid.columns]
grid[available_cols].to_csv(output_csv, index=False)

log(f"Wrote {len(grid)} entries to {output_csv}")
log("Done!")
