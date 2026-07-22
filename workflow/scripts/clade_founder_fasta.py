#!/usr/bin/env python3
"""Generate the founder FASTA for each clade.

For every gene in the ``pcp_genes`` config list, a gene-sliced founder FASTA is
also written at ``clade_founders/<gene>/<clade>.fa``, containing just that
gene's coding sequence sliced from the reference GTF's CDS blocks (so the stop
codon is excluded and ORF1ab's ribosomal frameshift is reconstructed into a
single in-frame CDS). The gene-sliced founders are needed for training
ParentIndependent DASM models.
"""

import json
import os
import sys
from importlib import util as importlib_util

import Bio.SeqIO


def load_workflow_module(module_name):
    """Load a module from the workflow's modules directory."""
    workflow_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    module_path = os.path.join(workflow_dir, "modules", f"{module_name}.py")
    spec = importlib_util.spec_from_file_location(module_name, module_path)
    module = importlib_util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


sequences = load_workflow_module("sequences")
gtf_utils = load_workflow_module("gtf")


def write_fasta(filepath, header, seq, line_width=80):
    """Write a single-record FASTA with a wrapped sequence.

    The header is written verbatim (caller supplies any prefix); the sequence
    is wrapped at ``line_width`` characters per line. Used for both the
    whole-genome founder and the optional gene-sliced founder.
    """
    with open(filepath, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), line_width):
            f.write(seq[i : i + line_width] + "\n")


# Get inputs from snakemake
neher_json_path = snakemake.input.neher_json
roemer_json_path = snakemake.input.roemer_json
ref_fasta = snakemake.input.ref_fasta
gtf_path = snakemake.input.gtf
roemer_nextstrain_to_pango = snakemake.params.roemer_nextstrain_to_pango
# One gene-sliced founder FASTA is emitted per configured protein. ``pcp_genes``
# and the ``gene_fasta`` output list are built from the same config list in the
# same order (see ``_gene_sliced_founder_outputs`` in the Snakefile), so they
# zip position-for-position. Coordinates come from the GTF's CDS blocks, so
# ORF1ab's frameshift is reconstructed into a single in-frame CDS.
pcp_genes = list(snakemake.params.pcp_genes)
clade = snakemake.wildcards.clade
output_fasta = snakemake.output.fasta
gene_fastas = list(getattr(snakemake.output, "gene_fasta", []) or [])
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
write_fasta(output_fasta, clade, founder_seq)
log(f"Wrote founder FASTA to {output_fasta}")

# Gene-sliced founder FASTA, one per configured protein.
# CDS-block coordinates are 1-based inclusive relative to the Wuhan-Hu-1
# (NC_045512.2) reference and come from the GTF (so the stop codon is excluded
# and ORF1ab's frameshift is reconstructed). Slicing is delegated to
# ``sequences.slice_cds``, which byte-aligns with the PCP parent/child slices.
assert len(pcp_genes) == len(gene_fastas), (
    f"pcp_genes ({len(pcp_genes)}) and gene_fasta outputs ({len(gene_fastas)}) must match in length/order"
)
for gene, gene_fasta_path in zip(pcp_genes, gene_fastas):
    cds_blocks = gtf_utils.cds_blocks_for_gene(gtf_path, gene)
    gene_seq = sequences.slice_cds(founder_seq, cds_blocks)
    expected_len = sum(end - start + 1 for start, end in cds_blocks)
    assert len(gene_seq) == expected_len, (
        f"Gene slice length {len(gene_seq)} != expected {expected_len} "
        f"for {gene} {cds_blocks}"
    )
    write_fasta(gene_fasta_path, clade, gene_seq)
    log(
        f"Wrote gene-sliced ({gene} {cds_blocks}, {len(gene_seq)} nt) "
        f"founder FASTA to {gene_fasta_path}"
    )

log("Done!")
