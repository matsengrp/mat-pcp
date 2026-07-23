#!/usr/bin/env python3
"""Slice one gene's in-frame coding sequence out of a clade founder FASTA.

Reads a clade's whole-genome founder FASTA (from ``clade_founder_fasta``) and
writes the gene-sliced founder at ``clade_founders/<gene>/<clade>.fa``. The
slice contains just that gene's coding sequence, taken from the reference GTF's
CDS blocks, so the stop codon is excluded and ORF1ab's ribosomal frameshift is
reconstructed into a single in-frame CDS. Slicing is delegated to
``sequences.slice_cds`` on the GTF blocks from ``gtf.cds_blocks_for_gene``, so
the founder slice byte-aligns with the PCP parent/child slices produced by
``pcp_from_tree.py``. The gene-sliced founders are the wildtype sequences used
when training ParentIndependent DASM models.

One gene, one clade per invocation: the Snakemake ``slice_founder_gene`` rule
fans this over the ``{gene}`` x ``{clade}`` grid, so the per-protein
conditionality lives in the DAG rather than in a conditional output list.
"""

import os
from importlib import util as importlib_util

import Bio.SeqIO

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


sequences = load_workflow_module("sequences")
gtf_utils = load_workflow_module("gtf")


def write_fasta(filepath, header, seq, line_width=80):
    """Write a single-record FASTA with a wrapped sequence."""
    with open(filepath, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), line_width):
            f.write(seq[i : i + line_width] + "\n")


founder_fasta = snakemake.input.fasta
gtf_path = snakemake.input.gtf
gene = snakemake.wildcards.gene
clade = snakemake.wildcards.clade
output_fasta = snakemake.output.gene_fasta
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


founder_seq = str(Bio.SeqIO.read(founder_fasta, "fasta").seq)
cds_blocks = gtf_utils.cds_blocks_for_gene(gtf_path, gene)
gene_seq = sequences.slice_cds(founder_seq, cds_blocks)

expected_len = sum(end - start + 1 for start, end in cds_blocks)
assert len(gene_seq) == expected_len, (
    f"Gene slice length {len(gene_seq)} != expected {expected_len} "
    f"for {gene} {cds_blocks}"
)
# Independent structural invariant: a coding-sequence slice must be a whole
# number of codons. Unlike the length check above (derived from the same
# blocks), this catches a malformed GTF whose CDS blocks aren't in-frame.
assert len(gene_seq) % 3 == 0, (
    f"Gene slice for {gene} is {len(gene_seq)} nt, not a whole number of "
    f"codons (blocks {cds_blocks})"
)

write_fasta(output_fasta, clade, gene_seq)
log(
    f"Wrote gene-sliced ({gene} {cds_blocks}, {len(gene_seq)} nt) founder "
    f"FASTA for clade {clade} to {output_fasta}"
)
