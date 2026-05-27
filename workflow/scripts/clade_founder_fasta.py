#!/usr/bin/env python3
"""Generate FASTA and mutation list for each clade founder.

When ``pcp_gene_region`` is set in the pipeline config, a second FASTA file is
also written containing only the gene-region slice of the founder sequence
(e.g. Spike: 21563-25381, 1-based inclusive). The gene-sliced founder is
needed for training ParentIndependent DASM models.
"""

import json
import Bio.SeqIO
import pandas as pd

# Get inputs from snakemake
neher_json_path = snakemake.input.neher_json
roemer_json_path = snakemake.input.roemer_json
ref_fasta = snakemake.input.ref_fasta
roemer_nextstrain_to_pango = snakemake.params.roemer_nextstrain_to_pango
pcp_gene_region = getattr(snakemake.params, "pcp_gene_region", None)
clade = snakemake.wildcards.clade
output_fasta = snakemake.output.fasta
output_muts = snakemake.output.muts
# gene_fasta is only present when pcp_gene_region is configured. Guard with
# getattr so this script still runs if a future config omits the slice output.
output_gene_fasta = getattr(snakemake.output, "gene_fasta", None)
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log(f"Generating founder sequence for clade: {clade}")

# Load reference sequence
ref_record = Bio.SeqIO.read(ref_fasta, "fasta")
ref_seq = str(ref_record.seq)
log(f"Reference length: {len(ref_seq)} bp")

# Load founder data
with open(neher_json_path) as f:
    neher_data = json.load(f)
with open(roemer_json_path) as f:
    roemer_data = json.load(f)

log(f"Neher data: {len(neher_data)} clades")
log(f"Roemer data: {len(roemer_data)} lineages")


def apply_mutations(seq, muts):
    """
    Apply mutations to a sequence (1-indexed positions).

    Args:
        seq: Original sequence string
        muts: List of mutation strings like "A123T"

    Returns:
        Modified sequence string
    """
    seq_list = list(seq)
    applied = 0
    skipped = 0

    for mut in muts:
        try:
            from_nt = mut[0]
            to_nt = mut[-1]
            pos = int(mut[1:-1])

            if pos < 1 or pos > len(seq_list):
                log(f"  WARNING: Position {pos} out of range, skipping {mut}")
                skipped += 1
                continue

            if seq_list[pos - 1] == from_nt:
                seq_list[pos - 1] = to_nt
                applied += 1
            else:
                log(
                    f"  WARNING: Expected {from_nt} at position {pos}, "
                    f"found {seq_list[pos - 1]}, skipping {mut}"
                )
                skipped += 1
        except (ValueError, IndexError) as e:
            log(f"  WARNING: Could not parse mutation {mut}: {e}")
            skipped += 1

    log(f"  Applied {applied} mutations, skipped {skipped}")
    return "".join(seq_list)


def get_founder_mutations(clade_name):
    """
    Get mutations for clade founder from available sources.

    Tries Neher data first, then falls back to Roemer via Pango mapping.
    """
    mutations = []
    source = None

    # Try Neher data first (Nextstrain clade names)
    if clade_name in neher_data:
        founder_info = neher_data[clade_name]
        # Neher data uses "nuc" key for nucleotide mutations
        if "nuc" in founder_info:
            mutations = founder_info["nuc"]
            source = "neher"
        elif "mutations" in founder_info:
            mutations = founder_info["mutations"]
            source = "neher"
        elif "nuc_mutations" in founder_info:
            mutations = founder_info["nuc_mutations"]
            source = "neher"

    # Fall back to Roemer data via Pango mapping
    if not mutations and clade_name in roemer_nextstrain_to_pango:
        pango = roemer_nextstrain_to_pango[clade_name]
        if pango in roemer_data:
            roemer_info = roemer_data[pango]
            if "nucSubstitutions" in roemer_info:
                mutations = roemer_info["nucSubstitutions"]
                source = "roemer"
            elif "substitutions" in roemer_info:
                # Parse substitution format if different
                mutations = roemer_info["substitutions"]
                source = "roemer"

    # Also check if clade name is directly in Roemer (for Pango lineages)
    if not mutations and clade_name in roemer_data:
        roemer_info = roemer_data[clade_name]
        if "nucSubstitutions" in roemer_info:
            mutations = roemer_info["nucSubstitutions"]
            source = "roemer"

    return mutations, source


# Get mutations and generate founder sequence
mutations, source = get_founder_mutations(clade)
log(f"Found {len(mutations)} mutations from source: {source or 'none'}")

if mutations:
    log(f"Mutations: {mutations[:10]}{'...' if len(mutations) > 10 else ''}")

founder_seq = apply_mutations(ref_seq, mutations)

# Verify sequence length unchanged
assert len(founder_seq) == len(ref_seq), "Founder sequence length mismatch!"

# Count differences from reference
n_diffs = sum(1 for a, b in zip(ref_seq, founder_seq) if a != b)
log(f"Founder differs from reference at {n_diffs} positions")

# Write FASTA
with open(output_fasta, "w") as f:
    f.write(f">{clade}\n")
    # Write sequence in 80-character lines
    for i in range(0, len(founder_seq), 80):
        f.write(founder_seq[i : i + 80] + "\n")
log(f"Wrote founder FASTA to {output_fasta}")

# Write mutations CSV
if len(mutations) > 0:
    muts_data = []
    for mut in mutations:
        try:
            from_nt = mut[0]
            to_nt = mut[-1]
            site = int(mut[1:-1])
            muts_data.append(
                {"mutation": mut, "site": site, "from_nt": from_nt, "to_nt": to_nt}
            )
        except (ValueError, IndexError):
            continue
    muts_df = pd.DataFrame(muts_data)
else:
    muts_df = pd.DataFrame(columns=["mutation", "site", "from_nt", "to_nt"])

muts_df.to_csv(output_muts, index=False)
log(f"Wrote {len(muts_df)} mutations to {output_muts}")

# Optional gene-sliced founder FASTA
# ``pcp_gene_region`` boundaries are 1-based inclusive (start/end both included)
# relative to the Wuhan-Hu-1 (NC_045512.2) reference. For Spike, end=25381
# excludes the trailing stop codon by convention (positions 21563-25381 ->
# 3819 nt = 1273 codons). We use [start-1:end] for half-open Python slicing
# that yields the same inclusive nt range.
if output_gene_fasta and pcp_gene_region:
    gene_name = pcp_gene_region.get("gene") or "gene"
    start = pcp_gene_region["start"]
    end = pcp_gene_region["end"]
    gene_seq = founder_seq[start - 1 : end]
    expected_len = end - start + 1
    assert len(gene_seq) == expected_len, (
        f"Gene slice length {len(gene_seq)} != expected {expected_len} "
        f"for {gene_name} {start}-{end}"
    )
    with open(output_gene_fasta, "w") as f:
        f.write(f">{clade}\n")
        for i in range(0, len(gene_seq), 80):
            f.write(gene_seq[i : i + 80] + "\n")
    log(
        f"Wrote gene-sliced ({gene_name} {start}-{end}, {len(gene_seq)} nt) "
        f"founder FASTA to {output_gene_fasta}"
    )

log("Done!")
