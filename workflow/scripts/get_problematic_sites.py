#!/usr/bin/env python3
"""Identify problematic codon sites in a specific gene for masking.

Two ways a gene's codon sites are resolved:

1. **Genes the reference GTF carries** (S, ORF1ab, ORF1a, M, N, E, ...) are
   matched against ``counts_all.csv``'s ``gene`` column, whose ``codon_site`` is
   already in that gene's own numbering.
2. **PCP-only genes** given coordinates through ``pcp_gene_cds_blocks``
   (currently ORF1b) never appear in that column at all -- it is derived from the
   GTF via ``ref_coding_sites``. Matching on the name therefore selected **zero
   rows** and wrote an empty mask file at exit 0, so ORF1b's DASM-ready curation
   masked nothing while every other gene masked its problematic sites (35 codon
   sites on the v9 test tree). For these, the gene's codons are derived from its
   own CDS blocks: a fully-excluded ``nt_site`` marks whichever of this gene's
   codons covers it.

Which route applies is decided from the *config* (does the gene have a
``pcp_gene_cds_blocks`` entry), never from whether the count table happens to
mention the gene -- see the comment at the dispatch below.

Three routes to a silently empty mask all raise, rather than emitting one --
same reasoning as ``genemap.require_nonempty_fitness_slice`` (#133), since a
silent empty file is indistinguishable from "nothing needed masking": a gene
resolvable by neither route, blocks whose coordinates intersect nothing in the
count grid, and a count table missing its ``codon_site`` column.

An empty mask for a *resolvable* gene is left legal, on the theory that a short
gene might contain no fully-excluded site -- but note that **no shipped config
actually produces one**: every configured gene has a non-empty mask on both the
v9 test tree (S 17, ORF1ab 86, ORF1b 35, E 1) and the v9 production tree (S 17,
ORF1a 52, ORF1b 35, M 2, N 8, E 1). So the WARNING below fires in no CI run and
is untested. If a future gene legitimately masks nothing, confirm it against the
gene's length before trusting it: at E's 75 codons the observed per-codon
masking rate (~1.3%) predicts 0.95 masked codons, so an empty mask would be
entirely unremarkable -- and empty is exactly the signature of the #133 ORF1b
defect. Recompute rather than reason from the rate.
"""

import os
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

genemap = load_workflow_module("genemap")
gtf_utils = load_workflow_module("gtf")

# Get snakemake variables
expected_vs_actual_csv = snakemake.input.expected_vs_actual
gtf_path = snakemake.input.gtf
output_file = snakemake.output.sites_to_mask
masking_gene = snakemake.wildcards.gene
pcp_gene_cds_blocks = snakemake.params.pcp_gene_cds_blocks
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log(f"Identifying problematic codon sites in gene: {masking_gene}")

# Load expected vs actual counts
expected_v_actual = pd.read_csv(expected_vs_actual_csv)
log(f"Loaded {len(expected_v_actual)} count entries")

# Find problematic codon sites to mask in the target gene.
#
# A codon site is masked only if it is excluded across ALL clades -- i.e. every
# row at that ``nt_site`` has ``exclude == True``. This keeps the
# clade-independent masks (out of ``site_include_range``, config
# ``sites_to_exclude``, and ProblematicSites VCF sites -- each of which sets
# ``exclude`` on every row at the site) and drops the clade-specific ones (UShER
# per-clade ``sites``/``mutations`` masks and founder->reference reversions),
# which should not clade-specifically mask positions of the target gene.
if "exclude" not in expected_v_actual.columns:
    raise KeyError("counts table is missing the required 'exclude' column")

excluded = expected_v_actual["exclude"] == True  # noqa: E712
site_all_excluded = excluded.groupby(expected_v_actual["nt_site"]).transform("all")

def _numeric(text):
    """Parse one ``codon_site`` component, or None if it is not a number."""
    text = text.strip()
    if not text.replace(".", "").replace("-", "").isdigit():
        return None
    return int(float(text))


def codon_sites_from_counts_column(counts, keep_row, masking_gene):
    """Read ``masking_gene``'s own ``codon_site`` components (GTF-carried gene).

    ``ref_coding_sites`` joins *every* covering CDS feature at an overlap site,
    so ``gene`` and ``codon_site`` are both ``;``-joined and paired
    **positionally**: ``"N;ORF9b"`` with ``"4;1"``, and ``"ORF1a;ORF1ab;ORF1ab"``
    with ``"4401;4401;4402"`` at the frameshift nucleotide. They must therefore
    be zipped and filtered, not flattened: appending every component would file
    ORF9b's codon 1 as a codon of N. Both are valid N codon numbers (ORF9b runs
    1-98, well inside N's range), so no range check could see the contamination
    -- it would simply over-mask a real codon of the wrong gene.

    A row whose two cells disagree on component *count* raises rather than being
    zipped short, because that is a structural disagreement about the same site
    and silently dropping the tail is how the wrong codon gets masked. Same
    reasoning as ``build_count_table``'s structural check (#130).

    The gene name is matched as a whole component: a substring match would let
    ``ORF1a`` bleed into ``ORF1ab`` once fanned over ``pcp_genes`` (#114). NaN
    cells become the string ``nan``, which equals no real gene symbol.
    """
    rows = counts[keep_row]
    if "codon_site" not in rows.columns:
        # Raise, symmetric with the missing-'exclude' and missing-'gene' checks:
        # returning [] here wrote a 0-byte mask at exit 0 for *every* GTF-carried
        # gene, and the empty-mask WARNING below is indistinguishable from a
        # legitimately empty one. It compounds with build_count_table.py's
        # `available_cols = [c for c in output_columns if c in grid.columns]`,
        # which drops a missing column from counts_all.csv silently -- so one
        # upstream schema regression would unmask every gene with nothing to
        # attribute it to. This is the third such route found in this file; the
        # other two now raise as well.
        raise KeyError(
            f"counts table has no 'codon_site' column, so {masking_gene!r}'s "
            f"codon numbers cannot be read from it. This is a schema problem "
            f"upstream in build_count_table, not an absence of problematic "
            f"sites -- refusing to write an empty mask."
        )
    found = []
    for gene_cell, codon_cell in zip(rows["gene"], rows["codon_site"]):
        if pd.isna(codon_cell):
            continue
        genes = [g.strip() for g in str(gene_cell).split(";")]
        codons = [c.strip() for c in str(codon_cell).split(";")]
        if len(genes) != len(codons):
            raise ValueError(
                f"count table row has {len(genes)} 'gene' components "
                f"({gene_cell!r}) but {len(codons)} 'codon_site' components "
                f"({codon_cell!r}). They are positionally paired, so a mismatch "
                f"means the two columns disagree about this site and zipping "
                f"them would file another gene's codon number under "
                f"{masking_gene!r}."
            )
        for gene_component, codon_component in zip(genes, codons):
            if gene_component != masking_gene:
                continue
            codon_number = _numeric(codon_component)
            if codon_number is not None:
                found.append(codon_number)
    return found


def codon_sites_from_cds_blocks(counts, keep_site, cds_blocks, masking_gene):
    """Derive codon sites from ``nt_site`` alone (PCP-only gene).

    ``counts_all.csv``'s ``gene`` column comes from the reference GTF, so a gene
    added only through ``pcp_gene_cds_blocks`` (ORF1b) never appears there.
    Asking which of *this gene's* codons covers each fully-excluded ``nt_site``
    gives the right answer without depending on that column at all, and is
    immune to the cross-gene pairing above by construction.
    """
    site_to_codons = {}
    for codon_number, positions in enumerate(
        genemap.codon_genomic_positions(cds_blocks), start=1
    ):
        for position in positions:
            site_to_codons.setdefault(position, set()).add(codon_number)

    # A coordinate typo that puts the blocks outside the count grid reaches the
    # same silent empty mask the name-resolution guard exists to prevent, by a
    # different route -- every lookup below simply misses. Any real CDS has
    # sites in the grid, so requiring a non-empty intersection cannot false-alarm.
    grid_sites = set(counts["nt_site"].astype(int))
    covered = grid_sites & set(site_to_codons)
    if not covered:
        raise ValueError(
            f"none of {masking_gene!r}'s configured CDS blocks {cds_blocks} "
            f"intersect the count table's nt_site range "
            f"({min(grid_sites)}-{max(grid_sites)}). Every codon lookup would "
            f"miss and the mask would be silently empty, which is "
            f"indistinguishable from 'nothing needed masking'. Check the "
            f"pcp_gene_cds_blocks coordinates."
        )
    log(f"{masking_gene}: {len(covered)} of {len(site_to_codons)} CDS "
        f"nucleotides are present in the count grid")

    found = []
    for nt_site in counts.loc[keep_site, "nt_site"].astype(int).unique():
        found.extend(site_to_codons.get(int(nt_site), ()))
    return found


# Which derivation applies is a fact about the *config*, not about whichever
# genes the count table happens to mention. Branching on "did the table match
# any rows" would let a gene defined in both the GTF and pcp_gene_cds_blocks
# silently ignore its override here (bypassing resolve_cds_blocks' add-only
# guard, which then only fires because other rules call it too), and would
# misreport a GTF gene with an empty slice as "appears in neither".
gene_has_override = masking_gene in pcp_gene_cds_blocks
gene_in_gtf = masking_gene in gtf_utils.gtf_cds_genes(gtf_path)

if not gene_has_override and not gene_in_gtf:
    raise ValueError(
        f"gene {masking_gene!r} appears in neither the reference GTF "
        f"({gtf_path}, the source of the count table's 'gene' column) nor "
        f"'pcp_gene_cds_blocks'. Masking it would silently write an empty file, "
        f"which is indistinguishable from 'nothing needed masking' -- and would "
        f"leave this gene's DASM-ready PCPs unmasked. Add it to the GTF via "
        f"add_to_ref_gtf, or give it explicit coordinates in "
        f"pcp_gene_cds_blocks."
    )

if gene_has_override:
    # ---- PCP-only gene (pcp_gene_cds_blocks): derive from nt_site ----
    cds_blocks = gtf_utils.resolve_cds_blocks(
        gtf_path, masking_gene, pcp_gene_cds_blocks
    )
    log(f"{masking_gene} is not in the reference GTF; deriving codon sites from "
        f"its configured CDS blocks {cds_blocks}")
    codon_sites = codon_sites_from_cds_blocks(
        expected_v_actual, site_all_excluded, cds_blocks, masking_gene
    )
else:
    # ---- GTF-carried gene: codon_site is already in this gene's numbering ----
    if "gene" not in expected_v_actual.columns:
        raise KeyError(
            f"counts table has no 'gene' column, so {masking_gene!r} cannot be "
            f"matched against it"
        )
    in_gene = expected_v_actual["gene"].astype(str).apply(
        lambda cell: masking_gene in [g.strip() for g in cell.split(";")]
    )
    codon_sites = codon_sites_from_counts_column(
        expected_v_actual, site_all_excluded & in_gene, masking_gene
    )

codon_sites = sorted(set(codon_sites))

log(f"Found {len(codon_sites)} problematic codon sites")
if not codon_sites:
    # Legitimate for a small gene, but it means this gene's DASM-ready PCPs are
    # masked at zero positions -- the observable symptom of #132's ORF1b defect.
    # Say so at the top of the log rather than leaving a 0-byte file to speak.
    log(f"WARNING: {masking_gene} has no fully-excluded codon site, so its mask "
        f"is empty and its DASM-ready PCPs will be masked at no position. This "
        f"is expected only if {masking_gene} genuinely contains no excluded "
        f"site; verify against the gene's length before trusting it.")

# Write output
with open(output_file, "w") as f:
    f.write("\n".join([str(c) for c in codon_sites]))

log(f"Wrote {len(codon_sites)} sites to {output_file}")
log("Done!")
