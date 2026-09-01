#!/usr/bin/env python3
"""Extract parent-child pairs from a mutation-annotated tree.

This step also drops branches whose parent and child are identical within the
gene-region slice. These arise when a branch's only mutations fall outside the
sliced region - for non-spike branches in a spike-only extraction, this catches
the same case that the old `curate_pcps_by_clade` rule used to handle, but
fused into one pass so no `_pcps_raw.csv` intermediate is needed.

Outputs:
  pcps              CSV of curated PCPs (no identical parent/child within slice).
  filtered_branches CSV with a single `child_name` column listing branches that
                    were dropped because the slice was identical.
"""

import os
from collections import defaultdict
from importlib import util as importlib_util
from time import time

import bte
import pandas as pd
from pandas import DataFrame as DF


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
pcpio = load_workflow_module("pcpio")

# Get inputs from snakemake
tree_path = snakemake.input.tree
fasta_path = snakemake.input.fasta
gtf_path = snakemake.input.gtf
coding_sites_csv = snakemake.input.coding_sites
output_path = snakemake.output.pcps
filtered_branches_path = snakemake.output.filtered_branches
branch_qc_path = snakemake.output.branch_qc
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None

# The gene to slice to comes from the ``{gene}`` wildcard. Its CDS blocks come
# from ``pcp_gene_cds_blocks`` if that config names the gene, otherwise from the
# reference GTF (ORF1ab -> two frameshift blocks, everything else a single
# block). ORF1b is the only configured override: it is a PCP-only gene the
# shared GTF deliberately does not carry (#133).
gene = snakemake.wildcards.gene

# Get params
max_nt_mutations = snakemake.params.get("max_nt_mutations", 4)
max_reversions = snakemake.params.get("max_reversions_to_ref", 1)
max_reversions_to_founder = snakemake.params.get("max_reversions_to_clade_founder", 1)
# Attribute access, not .get: the Snakefile always declares this param (coerced
# to {}), so a default here would mask a dropped declaration as "no overrides
# configured" and silently fall back to the GTF for a PCP-only gene.
pcp_gene_cds_blocks = snakemake.params.pcp_gene_cds_blocks

# Clade founder sequence (optional, for filtering reversions to clade founder)
clade_founder_fasta = snakemake.input.get("clade_founder_fasta", None)


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


def load_reference_from_fasta(reference_file):
    """Load single sequence from FASTA file."""
    with open(reference_file, "r") as fh:
        next(fh)  # Skip header
        reference_seq = "".join(map(str.strip, fh))
    return reference_seq


def load_founder_sequence(founder_fasta):
    """Load clade founder sequence from FASTA file if provided."""
    if founder_fasta is None:
        return None
    return load_reference_from_fasta(founder_fasta)


def apply_muts(sequence, muts):
    """
    Return the sequence after applying mutations.

    Args:
        sequence: The sequence string.
        muts: Dict of {site: nucleotide} mutations (1-indexed sites).
    """
    if len(muts) == 0:
        return sequence
    sites = sorted(muts.keys())
    substrings = [sequence[: sites[0] - 1]]
    for start, stop in zip(sites[:-1], sites[1:]):
        tobase = muts[start]
        seq = sequence[start : stop - 1]
        substrings.append(tobase)
        substrings.append(seq)
    substrings.append(muts[sites[-1]])
    substrings.append(sequence[sites[-1] :])
    return "".join(substrings)


class PCPHelper:
    """Helper class to extract PCPs from a tree with filtering."""

    def __init__(
        self,
        tree_path,
        fasta_path,
        gtf_path,
        coding_sites_csv,
        max_nt_mutations=4,
        max_reversions=1,
        max_reversions_to_founder=1,
        founder_seq=None,
    ):
        log(f"Loading tree from {tree_path}")
        self.tree = bte.MATree(tree_path)
        log(f"Translating mutations using {gtf_path}")
        self.translations = defaultdict(list, self.tree.translate(gtf_path, fasta_path))
        self.nodes = self.tree.depth_first_expansion()
        self.int_nodes = [n for n in self.nodes if not n.is_leaf()]
        self.n_nodes = len(self.nodes)
        self.n_int_nodes = len(self.int_nodes)
        self.ref_seq = load_reference_from_fasta(fasta_path)
        self.founder_seq = founder_seq
        self.max_nt_mutations = max_nt_mutations
        self.max_reversions = max_reversions
        self.max_reversions_to_founder = max_reversions_to_founder

        # Load coding sites for filtering (exclude noncoding from quality criteria)
        log(f"Loading coding sites from {coding_sites_csv}")
        coding_sites_df = pd.read_csv(coding_sites_csv)
        self.coding_sites = set(
            coding_sites_df[coding_sites_df["gene"] != "noncoding"]["site"]
        )
        log(f"Loaded {len(self.coding_sites)} coding sites")

        log(f"Tree has {self.n_nodes} total nodes, {self.n_int_nodes} internal nodes")
        if founder_seq:
            log(f"Clade founder sequence loaded ({len(founder_seq)} bp)")

    def failure_reason(self, nt_muts, codon_muts):
        """Name the QC limit this branch fails, or ``None`` if it passes.

        Quality filters are applied to CODING mutations only, to align with
        SARS2-mut-fitness behavior and exclude error-prone noncoding regions
        from quality assessment.

        Returns a *reason* rather than a bool (#151) so
        ``subsample_summary.csv`` can report which limit cost an arm its
        branches; the limits are independent, so the early-return order below
        decides only the attribution of a branch failing several, never which
        branches survive. ``count_mutations.py`` attributes in this same order,
        or the counts and PCP halves of that summary would file one branch
        under two reasons.
        """
        # Extract coding mutations only (exclude noncoding mutations)
        # Uses pre-loaded coding_sites set from coding_sites.csv
        coding_nt_muts = [m for m in nt_muts if int(m[1:-1]) in self.coding_sites]

        # Too many CODING mutations
        if len(coding_nt_muts) > self.max_nt_mutations:
            return "max_nt_mutations"

        # No mutations (coding or noncoding)
        if len(nt_muts) == 0:
            return "no_mutations"

        # Too many reversions to reference (CODING only)
        n_reversions = sum(
            self.ref_seq[int(mut[1:-1]) - 1] == mut[-1]
            for mut in coding_nt_muts
        )
        if n_reversions > self.max_reversions:
            return "max_reversions_to_ref"

        # Too many reversions to clade founder (CODING only)
        if self.founder_seq is not None:
            n_reversions_to_founder = sum(
                self.founder_seq[int(mut[1:-1]) - 1] == mut[-1]
                for mut in coding_nt_muts
            )
            if n_reversions_to_founder > self.max_reversions_to_founder:
                return "max_reversions_to_clade_founder"

        # Multiple mutations at same codon position (ambiguous)
        codon_positions = {(mut.gene, mut.aa_index) for mut in codon_muts}
        if len(codon_positions) < len(codon_muts):
            return "duplicate_codon_position"

        return None

    def count_mutations_from_parent(self, parent):
        """Extract PCPs from a parent node."""
        parent_node_haplotype = self.tree.get_haplotype(parent.id)
        p_muts = {int(mut[1:-1]): mut[-1] for mut in parent_node_haplotype}
        parent_seq = apply_muts(self.ref_seq, p_muts)
        pcps = (parent.id, parent_seq, [])

        n_filtered = 0
        for node in parent.children:
            nt_muts = node.mutations
            codon_muts = self.translations[node.id]
            self.qc_hist[
                sum(1 for m in nt_muts if int(m[1:-1]) in self.coding_sites)
            ] += 1
            reason = self.failure_reason(nt_muts, codon_muts)
            if reason is not None:
                self.qc_drops[reason] += 1
                continue
            n_filtered += 1
            pcps[2].append((node.id, nt_muts))
        return n_filtered, pcps

    def count_mutations_on_tree(self):
        """Iterate over all nodes and extract PCPs (no DataFrame conversion yet)."""
        # Reset per call so a second invocation on one helper cannot double-count.
        self.qc_drops = defaultdict(int)
        self.qc_hist = defaultdict(int)
        n_passed = 0
        all_pcps = []
        t0 = time()
        for i, node in enumerate(self.int_nodes, 1):
            count, pcps = self.count_mutations_from_parent(node)
            n_passed += count
            if count != 0:
                all_pcps.append(pcps)
            if i % 1000 == 0:
                elapsed = time() - t0
                log(f"  {i}/{self.n_int_nodes} nodes processed in {elapsed:.1f}s")
        return n_passed, all_pcps

    def pcps_to_df(self, pcps, cds_blocks=None):
        """Convert PCP list to DataFrame with optional CDS slicing.

        ``cds_blocks`` is a list of ``(start, end)`` tuples (1-based, inclusive
        on both ends, relative to the Wuhan-Hu-1 reference NC_045512.2) for the
        gene being extracted, derived from the reference GTF. Single-block genes
        pass one tuple (e.g. Spike ``[(21563, 25381)]`` = 3819 nt / 1273 codons,
        stop codon excluded); ORF1ab passes its two frameshift blocks.

        Slicing is delegated to ``sequences.slice_cds`` (see
        ``workflow/modules/sequences.py``), which concatenates the blocks into a
        single in-frame coding sequence.
        """

        def mut_dict(muts):
            return {int(mut[1:-1]): mut[-1] for mut in muts}

        def child_seq(seq, muts):
            return apply_muts(seq, mut_dict(muts))

        def slice_seq(seq):
            if cds_blocks is None:
                return seq
            return sequences.slice_cds(seq, cds_blocks)

        cols = ["parent_name", "child_name", "parent", "child", "branch_length"]
        rows = (
            (
                p_id,
                c_id,
                slice_seq(p_seq),
                slice_seq(child_seq(p_seq, c_muts)),
                len(c_muts),
            )
            for (p_id, p_seq, child_entries) in pcps
            for c_id, c_muts in child_entries
        )
        return DF(rows, columns=cols)

    def make_and_write_pcps(
        self,
        output_path,
        filtered_branches_path,
        cds_blocks=None,
    ):
        """Extract, filter, and write PCPs + filtered-branch IDs."""
        log("Extracting PCPs from tree...")
        n_branches_passed, pcps = self.count_mutations_on_tree()
        log(f"Found {n_branches_passed} branches passing QC filters")

        if cds_blocks:
            log(f"Slicing sequences to CDS blocks {cds_blocks}")

        pcp_df = self.pcps_to_df(pcps, cds_blocks=cds_blocks)
        log(f"PCPs after slicing: {len(pcp_df)}")

        # Drop PCPs whose parent and child are identical within the (possibly
        # sliced) sequence. These were previously removed by the now-defunct
        # `curate_pcps_by_clade` rule; folding the check in here means callers
        # see a single curated output and no `_pcps_raw.csv` intermediate.
        identical_mask = pcp_df["parent"] == pcp_df["child"]
        n_identical = int(identical_mask.sum())
        filtered_df = pcp_df.loc[identical_mask, ["child_name"]].reset_index(drop=True)
        pcp_df = pcp_df.loc[~identical_mask].reset_index(drop=True)
        log(f"  identical: {n_identical}")
        log(f"Final PCPs: {len(pcp_df)} (dropped {n_identical} identical)")

        # gzipped: one of the four PCP families (#134). filtered_branches stays
        # plain -- it is a single ID column, kilobytes to a few hundred MB, and
        # staying readable with plain shell tools is worth more than the space.
        pcp_df.to_csv(
            output_path, index=False, compression=pcpio.PCP_CSV_COMPRESSION
        )
        log(f"Wrote PCPs to {output_path}")
        filtered_df.to_csv(filtered_branches_path, index=False)
        log(f"Wrote {len(filtered_df)} filtered branch IDs to {filtered_branches_path}")

        self.write_branch_qc(branch_qc_path, n_identical, len(pcp_df))

    def write_branch_qc(self, path, n_identical, n_kept):
        """Tally why each branch was kept or dropped, in long format (#151).

        Additive output: nothing existing reads it and no existing schema
        moves. ``subsample_summary.csv`` needs per-limit attribution because
        realized retention falls below the closed form for two separable
        reasons -- QC drops and same-site collapse -- and one aggregate
        "branches dropped" number cannot tell them apart.

        ``identical_slice`` is a *gene-level* outcome and has no counterpart on
        the counts arm: a branch whose only mutations fall outside this gene's
        CDS blocks is dropped here and kept there. That asymmetry is why the
        summary carries a ``gene`` column at all, empty on counts-arm rows.
        """
        rows = [("outcome", "kept", int(n_kept))]
        rows += [
            ("outcome", reason, int(self.qc_drops.get(reason, 0)))
            for reason in (
                "max_nt_mutations",
                "no_mutations",
                "max_reversions_to_ref",
                "max_reversions_to_clade_founder",
                "duplicate_codon_position",
            )
        ]
        rows.append(("outcome", "identical_slice", int(n_identical)))
        rows.append(
            ("total", "branches_scored", int(sum(self.qc_hist.values())))
        )
        rows += [
            ("n_nt_mutations_coding", str(int(k)), int(v))
            for k, v in sorted(self.qc_hist.items())
        ]
        DF(rows, columns=["category", "key", "n_branches"]).to_csv(path, index=False)
        log(f"Wrote branch QC tally to {path}")


def main():
    """Main entry point."""
    log("=" * 60)
    log("PCP Extraction")
    log("=" * 60)

    # Resolve this gene's CDS blocks: the pcp_gene_cds_blocks override if it
    # names the gene, otherwise the reference GTF.
    #
    # This re-reads the 30 kB reference that PCPHelper loads again below, on
    # purpose: resolving (and validating) the blocks here means a malformed
    # override fails in seconds, rather than after loading a multi-GB MAT. Do
    # not "optimize" it by moving the resolution after the helper.
    cds_blocks = gtf_utils.resolve_cds_blocks(
        gtf_path,
        gene,
        pcp_gene_cds_blocks,
        genome_length=len(load_reference_from_fasta(fasta_path)),
    )
    log(f"Gene {gene}: CDS blocks {cds_blocks}")

    log(f"Max mutations: {max_nt_mutations}")
    log(f"Max reversions to ref: {max_reversions}")
    log(f"Max reversions to founder: {max_reversions_to_founder}")

    # Load clade founder sequence if provided
    founder_seq = load_founder_sequence(clade_founder_fasta)
    if clade_founder_fasta:
        log(f"Loaded clade founder from {clade_founder_fasta}")

    helper = PCPHelper(
        tree_path,
        fasta_path,
        gtf_path,
        coding_sites_csv,
        max_nt_mutations=max_nt_mutations,
        max_reversions=max_reversions,
        max_reversions_to_founder=max_reversions_to_founder,
        founder_seq=founder_seq,
    )

    helper.make_and_write_pcps(
        output_path,
        filtered_branches_path,
        cds_blocks=cds_blocks,
    )

    log("=" * 60)
    log("Done!")


if __name__ == "__main__":
    main()
