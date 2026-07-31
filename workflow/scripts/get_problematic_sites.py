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

A gene resolvable by neither route raises, rather than emitting an empty mask --
same reasoning as ``genemap.require_nonempty_fitness_slice`` (#133): a silent
empty file here is indistinguishable from "nothing needed masking".
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

gene_members = None
gene_in_table = False
if "gene" in expected_v_actual.columns:
    # Match gene membership exactly. Compound ``gene`` cells from overlapping
    # ORFs (e.g. ``S;ORF3a``) are split on ``;`` so the masking gene must equal
    # one whole member -- a substring match would let a gene symbol bleed into a
    # longer one (e.g. ``ORF1a`` into ``ORF1ab``) once fanned over pcp_genes (#114).
    # NaN gene cells become the string ``nan``, which never equals a real gene
    # symbol, so they drop out of the mask -- matching the old ``na=False``.
    gene_members = expected_v_actual["gene"].astype(str).str.split(";")
    in_gene = gene_members.apply(
        lambda members: masking_gene in [g.strip() for g in members]
    )
    gene_in_table = bool(in_gene.any())

codon_sites = []

if gene_in_table:
    # ---- GTF-derived gene: codon_site is already in this gene's numbering ----
    problematic_sites = expected_v_actual[site_all_excluded & in_gene]
    if "codon_site" in problematic_sites.columns:
        raw_sites = problematic_sites["codon_site"].unique()
        # Handle potential semicolon-separated values (e.g., "123;456")
        for site in raw_sites:
            if pd.isna(site):
                continue
            site_str = str(site)
            if ";" in site_str:
                for s in site_str.split(";"):
                    s = s.strip()
                    if s.replace(".", "").replace("-", "").isdigit():
                        codon_sites.append(int(float(s)))
            elif site_str.replace(".", "").replace("-", "").isdigit():
                codon_sites.append(int(float(site_str)))
else:
    # ---- PCP-only gene (pcp_gene_cds_blocks): derive from nt_site ----
    # ``counts_all.csv``'s gene column comes from the reference GTF, so a gene
    # added only through pcp_gene_cds_blocks (ORF1b) never appears there and the
    # name match above selects zero rows. Resolving the gene's own CDS blocks and
    # asking which of *its* codons covers each fully-excluded nt_site gives the
    # right answer without depending on that column at all.
    if masking_gene not in (pcp_gene_cds_blocks or {}):
        raise ValueError(
            f"gene {masking_gene!r} appears in neither the count table's 'gene' "
            f"column (derived from the reference GTF) nor 'pcp_gene_cds_blocks'. "
            f"Masking it would silently write an empty file, which is "
            f"indistinguishable from 'nothing needed masking' -- and would leave "
            f"this gene's DASM-ready PCPs unmasked. Add it to the GTF via "
            f"add_to_ref_gtf, or give it explicit coordinates in "
            f"pcp_gene_cds_blocks."
        )
    cds_blocks = gtf_utils.resolve_cds_blocks(
        gtf_path, masking_gene, pcp_gene_cds_blocks
    )
    log(f"{masking_gene} is not in the reference GTF; deriving codon sites from "
        f"its configured CDS blocks {cds_blocks}")
    site_to_codons = {}
    for codon_number, positions in enumerate(
        genemap.codon_genomic_positions(cds_blocks), start=1
    ):
        for position in positions:
            site_to_codons.setdefault(position, set()).add(codon_number)
    excluded_sites = (
        expected_v_actual.loc[site_all_excluded, "nt_site"].astype(int).unique()
    )
    for nt_site in excluded_sites:
        codon_sites.extend(site_to_codons.get(int(nt_site), ()))

codon_sites = sorted(set(codon_sites))

log(f"Found {len(codon_sites)} problematic codon sites")

# Write output
with open(output_file, "w") as f:
    f.write("\n".join([str(c) for c in codon_sites]))

log(f"Wrote {len(codon_sites)} sites to {output_file}")
log("Done!")
