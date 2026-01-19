#!/usr/bin/env python3
"""Get nucleotide at each coding site for all clade founders."""

import os

import Bio.SeqIO
import pandas as pd

# Get snakemake variables
coding_sites_file = snakemake.input.coding_sites
fasta_files = snakemake.input.fastas
output_csv = snakemake.output.csv
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log(f"Processing {len(fasta_files)} clade founder files")

# Codons that start with these two nucleotides are 4-fold degenerate
codons_4fold = {
    "TC",  # serine
    "CT",  # leucine
    "CC",  # proline
    "CG",  # arginine
    "AC",  # threonine
    "GT",  # valine
    "GC",  # alanine
    "GG",  # glycine
}

# Load coding sites
sites = (
    pd.read_csv(coding_sites_file)
    .sort_values("site")
    .assign(
        codon_position=lambda x: x["codon_position"].str.split(";"),
        codon_site=lambda x: x["codon_site"].str.split(";"),
        gene=lambda x: x["gene"].str.split(";"),
    )
    .explode(["codon_position", "codon_site", "gene"])
)

log(f"Loaded {len(sites)} coding sites")


def get_codon(row, seq):
    """Get codon for a site given the sequence."""
    if row["codon_position"] == "noncoding":
        return "noncoding"
    r = row["site"]
    codon_pos = int(row["codon_position"])
    return f"{seq[r - codon_pos]}{seq[r - codon_pos + 1]}{seq[r - codon_pos + 2]}"


dfs = []
for f_name in fasta_files:
    seq = str(Bio.SeqIO.read(f_name, "fasta").seq)
    name = os.path.splitext(os.path.basename(f_name))[0]
    log(f"  Processing {name}: {len(seq)} bp")

    clade_df = (
        sites.assign(
            clade=name,
            nt=lambda x: x["site"].map(lambda r: seq[r - 1]),
            codon=lambda x: x.apply(lambda row: get_codon(row, seq), axis=1),
            four_fold_degenerate=lambda x: (
                (x["codon_position"] == "3") & x["codon"].str[:2].isin(codons_4fold)
            ),
        )
        .groupby(["clade", "site", "nt"], as_index=False)
        .aggregate(
            gene=pd.NamedAgg("gene", lambda s: ";".join(s)),
            codon=pd.NamedAgg("codon", lambda s: ";".join(s)),
            codon_position=pd.NamedAgg(
                "codon_position", lambda s: ";".join(map(str, s))
            ),
            codon_site=pd.NamedAgg("codon_site", lambda s: ";".join(s)),
            four_fold_degenerate=pd.NamedAgg("four_fold_degenerate", "all"),
        )
    )
    dfs.append(clade_df)

result = pd.concat(dfs)
result.to_csv(output_csv, index=False)

log(f"Wrote {len(result)} entries to {output_csv}")
log("Done!")
