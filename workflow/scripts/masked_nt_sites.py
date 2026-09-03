#!/usr/bin/env python3
"""Emit one gene's clade-independent nucleotide mask for DASM training.

The companion to #144's ``pcps/{gene}/censored_mutations.csv``, and the other
half of the same coordinate contract. That file carries the two *clade-specific*
reasons of ``masking.compute_exclusions``' five; this one carries the three
*clade-independent* ones, at nucleotide granularity:

  * ``sites_to_exclude``   -- the config's global site list
  * ``vcf_masked_sites``   -- the ProblematicSites VCF
  * ``site_include_range`` -- outside the configured genome bounds

Founder->reference reversions are deliberately in neither: they are a selection
filter, not an artifact mask. Together the two files give a consumer the
complete nucleotide mask with no codon widening anywhere.

**Why this file exists, and what replaced.** The nucleotide-level answer used
to be computed and thrown away by ``get_problematic_sites.py``:
``site_all_excluded = excluded.groupby(nt_site).transform("all")`` is exactly
this set, and that script widened it to codons and wrote only the codons, so one
suspect nucleotide took its two innocent neighbours with it -- a 2.85x over-mask
across the six v12 production genes (121 nt sites became 345 masked nt). #153
emitted this file as a separate rule rather than a second output there, because
``{gene}_codon_sites_to_mask.txt`` had seven consumers, one of them
``curate_pcps_for_dasm``, and rewriting it would have restated every PCP for a
file nothing yet read.

#158 then retired the codon file: ``curate_pcps_for_dasm`` reads this one, so
the drop filter and DASM's loss apply one mask at one granularity. Two
consequences. This is now the sole producer of the clade-independent mask, so
guard 1's reconciliation against ``counts_all.csv`` -- which reaches ``exclude``
through ``compute_exclusions`` rather than through
``global_masked_site_reasons`` -- is the only independent derivation left.
And the codon cross-check that made the retirement mechanical (guard 2 in #153:
widening these indices reproduced the codon file exactly, verified for every
configured gene on both v12 production trees and the test tree, with exact set
equality) is gone with the file it compared against; ``tests/test_masked_nt_sites.py``
is now the sole coverage of the CDS-block inversion.

Coordinates match ``censored_mutations.csv``: ``nt_site`` is the 1-based genomic
coordinate and ``nt_index`` is a 0-based index into this gene's PCP nucleotide
sequence, which the consumer sees as the whole sequence because it has no CDS
knowledge. The genomic->index map is **one-to-many** at a -1 ribosomal
frameshift (``sequences.genomic_to_pcp_indices``), so a mask on ORF1b's shared
nucleotide 13468 emits two rows, in different codons.
"""

import csv
import os
from collections import Counter
from importlib import util as importlib_util

import pandas as pd

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

gtf_utils = load_workflow_module("gtf")
masking = load_workflow_module("masking")
sequences = load_workflow_module("sequences")

counts_csv = snakemake.input.counts
site_mask_csv = snakemake.input.site_mask
gtf_path = snakemake.input.gtf
ref_fasta = snakemake.input.ref_fasta
output_table = snakemake.output.table
gene = snakemake.wildcards.gene
sites_to_exclude = snakemake.params.sites_to_exclude
site_include_range = snakemake.params.site_include_range
pcp_gene_cds_blocks = snakemake.params.pcp_gene_cds_blocks
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None

COLUMNS = ["nt_site", "nt_index", "reason"]


def log(msg):
    """Write to the log file if available, and always to stdout."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


def read_fasta(path):
    """Return a single-record FASTA's sequence, uppercased, newlines stripped."""
    with open(path) as f:
        return "".join(line.strip() for line in f if not line.startswith(">")).upper()


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

log(f"Building the global nucleotide mask for gene: {gene}")

# resolve_cds_blocks must be called WITH the overrides: ORF1b exists only in
# pcp_gene_cds_blocks, and the bare GTF lookup raises "No CDS blocks found".
cds_blocks = gtf_utils.resolve_cds_blocks(gtf_path, gene, pcp_gene_cds_blocks)
positions = sequences.cds_nt_positions(cds_blocks)
occurrences = Counter(positions)
gene_sites = sorted(occurrences)
log(
    f"{gene} CDS blocks {cds_blocks} -> {len(positions)} nt / "
    f"{len(positions) // 3} codons over {len(gene_sites)} distinct sites"
)

vcf_masked_sites = set(pd.read_csv(site_mask_csv)["site"].astype(int))
log(
    f"Global mask sources: {len(vcf_masked_sites)} ProblematicSites VCF sites, "
    f"{len(sites_to_exclude or [])} config sites_to_exclude, "
    f"site_include_range {site_include_range}"
)

# The reason vocabulary lives in modules/masking.py, beside the per-row
# authority that applies the same three tests inside compute_exclusions. One
# definition, so this file and counts_all.csv's `exclude` column cannot come to
# disagree about what "globally masked" means.
site_reasons = masking.global_masked_site_reasons(
    gene_sites,
    sites_to_exclude=sites_to_exclude,
    vcf_masked_sites=vcf_masked_sites,
    site_include_range=site_include_range,
)

# ---------------------------------------------------------------------------
# Rows
# ---------------------------------------------------------------------------
#
# One row per (nt_index, reason), NOT per nt_index: a site can carry two
# reasons, and 18 of the six v12 production genes' 121 masked sites do (both in
# the VCF and in the config list). Emitting the reasons as separate rows rather
# than a ";"-joined cell keeps the file long-form -- CLAUDE.md pattern 10
# records what positionally-paired ";"-joined cells cost when a consumer
# flattens them.
rows = []
for nt_site, reasons in sorted(site_reasons.items()):
    # One row per index: a site shared by two blocks at the frameshift is read
    # twice by slice_cds and so occupies two positions of the PCP sequence, in
    # different codons.
    for nt_index in sequences.genomic_to_pcp_indices(cds_blocks, nt_site):
        for reason in reasons:
            rows.append((nt_site, nt_index, reason))

# Sorted for byte-reproducibility: data-central manifests pin md5/xxh3 per file
# (CLAUDE.md pattern 9). (nt_index, reason) is a total order because the two
# together are unique -- one reason cannot be emitted twice at one index.
rows.sort(key=lambda row: (row[1], row[2]))

per_reason = Counter(reason for _, _, reason in rows)
log(
    f"{len(site_reasons)} globally masked nt sites in {gene} -> {len(rows)} rows "
    f"({dict(sorted(per_reason.items()))})"
)

# ---------------------------------------------------------------------------
# Guard 1: the emitted set reconciles against the count table, both directions
# ---------------------------------------------------------------------------
#
# A coordinate error that maps every site outside the gene writes a well-formed
# empty CSV at exit 0, and emptiness is *legitimate* for a short gene (E
# contributes exactly one site). So reconcile against a derivation that shares
# no code with the one above: `site_all_excluded` over counts_all.csv, which is
# the set get_problematic_sites.py computes and throws away, and which reaches
# `exclude` through compute_exclusions rather than through
# global_masked_site_reasons.
#
# Checked in BOTH directions with distinct messages, because the two failures
# mean different things -- see each branch.
# A callable ``usecols`` rather than a list, so a missing column produces the
# attributable message below instead of pandas' "Usecols do not match columns".
counts = pd.read_csv(counts_csv, usecols=lambda name: name in ("nt_site", "exclude"))
for column in ("nt_site", "exclude"):
    if column not in counts.columns:
        raise KeyError(
            f"{counts_csv} is missing the required {column!r} column, so the "
            f"independent reconciliation of {gene}'s mask cannot run. This is a "
            f"schema problem upstream in build_count_table, not an absence of "
            f"masked sites -- refusing to write an unreconciled file."
        )
all_excluded = (counts["exclude"] == True).groupby(  # noqa: E712
    counts["nt_site"].astype(int)
).all()
grid_sites = set(all_excluded.index.astype(int))
site_all_excluded = set(all_excluded.index[all_excluded].astype(int))
log(
    f"counts_all.csv: {len(grid_sites)} grid sites, "
    f"{len(site_all_excluded)} excluded in every clade genome-wide"
)

# A masked site the grid does not mention at all is a third failure mode, and it
# has to be separated out first: it makes the set comparison below meaningless
# in the "missing from the count table" direction.
off_grid = sorted(set(site_reasons) - grid_sites)
if off_grid:
    # min/max only when the grid is non-empty. A wholly empty count table is one
    # of the upstream failures this guard exists to name, and computing a range
    # over it would replace the diagnostic with a bare "min() arg is an empty
    # sequence" traceback pointing at this f-string.
    span = (
        f"the grid spans {min(grid_sites)}-{max(grid_sites)}"
        if grid_sites
        else "the grid is empty"
    )
    raise ValueError(
        f"{len(off_grid)} globally masked sites of {gene} are absent from "
        f"counts_all.csv's nt_site grid ({off_grid[:10]}...): {span}. "
        f"The count table and this gene's "
        f"CDS blocks {cds_blocks} disagree about which coordinates exist, so "
        f"neither the reconciliation below nor the codon cross-check means "
        f"anything."
    )

from_counts = site_all_excluded & set(gene_sites)
missing_here = sorted(from_counts - set(site_reasons))
if missing_here:
    # The count table says these sites are excluded in every clade, but none of
    # the three clade-independent reasons applies to them -- so they are
    # excluded clade-specifically in *every* clade (a UShER yaml block naming
    # them for all of them, or a reversion present in all of them). That is
    # possible in principle and true of nothing on any shipped tree (628
    # globally excluded sites genome-wide on both v12 production trees and the
    # test tree, and the two derivations agree exactly). It raises rather than
    # being logged because `reason` has no value that would describe such a
    # site honestly, and because the DASM-ready drop filter reads this file --
    # so a site the count table excludes everywhere and this file omits is a
    # position trained on as live.
    raise ValueError(
        f"{len(missing_here)} sites of {gene} are excluded in every clade of "
        f"counts_all.csv but match none of "
        f"{list(masking.GLOBAL_MASK_REASONS)}: {missing_here[:10]}. They are "
        f"masked clade-specifically in every clade, which this file's reason "
        f"vocabulary cannot express -- they belong in "
        f"pcps/{gene}/censored_mutations.csv, per clade."
    )

extra_here = sorted(set(site_reasons) - from_counts)
if extra_here:
    # The opposite direction, and it cannot happen by construction: every one of
    # the three reasons sets `exclude` on *every* grid row at the site, so a
    # site this script masks must be all-clade-excluded in the count table. A
    # failure here means build_count_table ran against a different config, VCF
    # or site_include_range than this rule did.
    raise ValueError(
        f"{len(extra_here)} sites of {gene} are globally masked here but are "
        f"NOT excluded in every clade of counts_all.csv: {extra_here[:10]}. "
        f"Each of {list(masking.GLOBAL_MASK_REASONS)} excludes every row at its "
        f"site, so this means the count table was built from a different "
        f"config, ProblematicSites VCF or site_include_range than this rule "
        f"was given."
    )

# The set equality above says nothing about row *multiplicity*, which is where
# the one-to-many frameshift index lives. Counted from the position list rather
# than from the emit loop, so the two derivations stay independent.
expected_rows = sum(
    occurrences[nt_site] * len(reasons) for nt_site, reasons in site_reasons.items()
)
if len(rows) != expected_rows:
    raise ValueError(
        f"the {gene} mask has {len(rows)} rows but its masked sites and CDS "
        f"blocks {cds_blocks} account for {expected_rows}. The emit loop and "
        f"the independent position count disagree, which means a row was "
        f"dropped or duplicated -- most likely the one-to-many frameshift index."
    )

if not rows:
    # Legitimate for a short gene, but it also means this gene's DASM training
    # applies no global mask at all -- the observable symptom of #133's ORF1b
    # defect, one file over. Say why rather than leaving a header-only file to
    # speak for itself.
    log(
        f"NOTE: {gene} has no globally masked nucleotide site, so this file is "
        f"header-only. That is legitimate -- {gene}'s CDS spans "
        f"{len(gene_sites)} sites, none of which is in the ProblematicSites "
        f"VCF ({len(vcf_masked_sites)} sites), in sites_to_exclude "
        f"({len(sites_to_exclude or [])} sites), or outside "
        f"{site_include_range} -- and it reconciles with counts_all.csv, which "
        f"excludes none of {gene}'s sites in every clade. Verify against the "
        f"gene's length before trusting it."
    )

# ---------------------------------------------------------------------------
# What codon granularity would have cost
# ---------------------------------------------------------------------------
#
# Reported, not checked: the counterparty is gone. Through v12 this was guard 2,
# an equality against {gene}_codon_sites_to_mask.txt reached by a different
# route, and it is what made #158's retirement of that file mechanical. The
# widening factor is still worth logging, because it is the size of the effect
# the DASM-ready drop filter no longer pays.
#
# It compares distinct masked *positions*, not rows: a two-reason site
# contributes two rows and one position, so dividing by the row count would
# understate the widening.
masked_indices = {nt_index for _, nt_index, _ in rows}
codon_widened = {nt_index // 3 + 1 for nt_index in masked_indices}
log(
    f"Granularity: {len(masked_indices)} masked nt positions would widen to "
    f"{len(codon_widened)} codons ({3 * len(codon_widened)} nt)"
    + (
        f" -- a {3 * len(codon_widened) / len(masked_indices):.2f}x over-mask "
        f"this file avoids"
        if masked_indices
        else ""
    )
)

# ---------------------------------------------------------------------------
# Guard 3: nt_index round-trips as a string identity
# ---------------------------------------------------------------------------
#
# Still numbered 3 with 2 retired above, deliberately: renumbering would make
# "guard 2" mean the codon cross-check in #153 and this check in #158, across
# ``tests/test_masked_nt_sites.py``, CLAUDE.md pattern 17 and both issues.
#
# Stated as a *sequence* identity rather than index arithmetic, so it also
# catches a blocks-vs-sequence mismatch (the wrong gene's blocks) that
# arithmetic alone would pass. Same guard as censor_table.py's, deliberately
# duplicated: it is three lines, and sharing it would mean the two files' only
# coordinate assertion had one implementation.
reference = read_fasta(ref_fasta)
gene_reference = sequences.slice_cds(reference, cds_blocks)
if len(gene_reference) != len(positions):
    raise ValueError(
        f"{gene}'s reference CDS slice is {len(gene_reference)} nt but its "
        f"position map has {len(positions)} entries; slice_cds and "
        f"cds_nt_positions disagree about blocks {cds_blocks}."
    )
for nt_site, nt_index, reason in rows:
    if gene_reference[nt_index] != reference[nt_site - 1]:
        raise ValueError(
            f"mask row (nt_site={nt_site}, nt_index={nt_index}, {reason}) does "
            f"not round-trip: {gene}'s reference CDS has "
            f"{gene_reference[nt_index]!r} at that index but the genome has "
            f"{reference[nt_site - 1]!r} at that site."
        )

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

os.makedirs(os.path.dirname(output_table) or ".", exist_ok=True)
with open(output_table, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(COLUMNS)
    writer.writerows(rows)

log(f"Wrote {len(rows)} rows to {output_table}")
log("Done!")
