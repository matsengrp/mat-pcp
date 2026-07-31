#!/usr/bin/env python3
"""Validate that PCP-derived AA mutation counts match aamut_fitness actual_count.

Compares amino acid mutations (e.g., D614G) between:
1. PCP-derived counts: Mutations from translating parent-child pairs
2. aamut_fitness actual_count: Aggregated from mutation counting pipeline

A PCP gene need not appear in ``aamut_fitness`` under its own name or codon
numbering, so the two sides are aligned through ``pcp_gene_to_fitness_gene``
before merging (#133). Filtering the fitness table on ``gene == {gene}``
unconditionally, as this script used to, made the ORF1a comparison **vacuous**:
``aamutfit.get_coding`` drops ORF1a (``aamutfit.py:80``), so the filter matched
zero rows and the report announced no discrepancies because it had nothing to
compare. An empty fitness slice for a gene that has PCP counts is now a hard
failure -- that guard alone would have caught the gap.
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
pcp_counts_csv = snakemake.input.pcp_counts
aamut_fitness_files = snakemake.input.aamut_fitness_files
all_aamut_fitness_csv = snakemake.input.all_aamut_fitness
gtf_path = snakemake.input.gtf
output_report = snakemake.output.report
output_discrepancies = snakemake.output.discrepancies
pcp_gene_name = snakemake.wildcards.gene
pcp_gene_cds_blocks = snakemake.params.pcp_gene_cds_blocks
pcp_gene_to_fitness_gene = snakemake.params.pcp_gene_to_fitness_gene
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log("Validating PCP-derived AA mutation counts against aamut_fitness actual_count")


# Load PCP-derived counts
pcp_counts = pd.read_csv(pcp_counts_csv)
log(f"PCP-derived AA mutation counts: {len(pcp_counts)} entries")

# Load per-clade aamut_fitness files
aamut_per_clade = []
for f in aamut_fitness_files:
    df = pd.read_csv(f)
    aamut_per_clade.append(df)

if aamut_per_clade:
    aamut_clade_df = pd.concat(aamut_per_clade, ignore_index=True)
    # Rename cluster to clade for consistency
    if "cluster" in aamut_clade_df.columns:
        aamut_clade_df = aamut_clade_df.rename(columns={"cluster": "clade"})
else:
    aamut_clade_df = pd.DataFrame()

log(f"Per-clade aamut_fitness: {len(aamut_clade_df)} entries")

# Load aggregate aamut_fitness
aamut_all = pd.read_csv(all_aamut_fitness_csv)
log(f"Aggregate aamut_fitness: {len(aamut_all)} entries")

# -----------------------------------------------------------------------------
# Align the PCP gene with the fitness-table gene it is validated against.
#
# Genes present in aamut_fitness under their own name (S, M, N, E, ORF1ab) need
# no mapping and take the identity path. ORF1a and ORF1b do: see the module
# docstring and workflow/modules/genemap.py.
# -----------------------------------------------------------------------------
fitness_gene = pcp_gene_name
codon_map = None
unmapped_codons = []
mapping = (pcp_gene_to_fitness_gene or {}).get(pcp_gene_name)

if mapping:
    fitness_gene = mapping["gene"]
    declared_offset = mapping["codon_offset"]
    pcp_blocks = gtf_utils.resolve_cds_blocks(
        gtf_path, pcp_gene_name, pcp_gene_cds_blocks
    )
    fitness_blocks = gtf_utils.resolve_cds_blocks(
        gtf_path, fitness_gene, pcp_gene_cds_blocks
    )
    codon_map, unmapped_codons = genemap.align_to_fitness_gene(
        pcp_gene_name, fitness_gene, pcp_blocks, fitness_blocks, declared_offset
    )
    log(
        f"{pcp_gene_name} -> {fitness_gene}: {len(codon_map)} codons correspond "
        f"(offset {declared_offset}), {len(unmapped_codons)} have no counterpart"
    )


def to_pcp_codon_frame(df, label):
    """Restrict a fitness table to ``fitness_gene``, renumbered into PCP codons.

    Renumbering rebuilds ``aa_mutation`` as well as ``aa_site``, because the
    merge keys downstream are the ``aa_mutation`` strings and the site number is
    embedded in them.
    """
    if len(df) == 0:
        return df
    sub = df[df["gene"] == fitness_gene].copy()
    log(f"Filtered {label} to {fitness_gene}: {len(sub)} entries")
    if codon_map is None:
        return sub
    fitness_to_pcp = {fitness: pcp for pcp, fitness in codon_map.items()}
    sub = sub[sub["aa_site"].isin(fitness_to_pcp)].copy()
    sub["aa_site"] = sub["aa_site"].map(fitness_to_pcp)
    sub["gene"] = pcp_gene_name
    sub["aa_mutation"] = (
        sub["clade_founder_aa"] + sub["aa_site"].astype(str) + sub["mutant_aa"]
    )
    log(f"  renumbered into {pcp_gene_name} codons: {len(sub)} entries")
    return sub


aamut_clade_df = to_pcp_codon_frame(aamut_clade_df, "per-clade")
aamut_all = to_pcp_codon_frame(aamut_all, "aggregate")

# Filter PCP counts to the gene, then set aside codons with no counterpart in
# the fitness gene's frame. ORF1a's codons 4402-4405 are the real case: they run
# past the -1 ribosomal frameshift in ORF1a's own frame, so ORF1ab has no codon
# on those nucleotides. Left in the merge they would be reported as ordinary
# discrepancies, which is exactly the "explained, not generic" requirement.
pcp_counts_no_counterpart = pd.DataFrame()
if pcp_gene_name and len(pcp_counts) > 0 and "gene" in pcp_counts.columns:
    pcp_counts = pcp_counts[pcp_counts["gene"] == pcp_gene_name]
    log(f"Filtered PCP counts to {pcp_gene_name}: {len(pcp_counts)} entries")

if codon_map is not None and len(pcp_counts) > 0:
    in_map = pcp_counts["aa_site"].isin(codon_map)
    pcp_counts_no_counterpart = pcp_counts[~in_map].copy()
    pcp_counts = pcp_counts[in_map].copy()
    log(
        f"PCP counts on codons with a {fitness_gene} counterpart: "
        f"{len(pcp_counts)}; without: {len(pcp_counts_no_counterpart)}"
    )

genemap.require_nonempty_fitness_slice(
    pcp_gene_name,
    fitness_gene,
    len(pcp_counts),
    len(aamut_all),
    all_aamut_fitness_csv,
)


def compare_counts(pcp_df, aamut_df, clade_col="clade", is_aggregate=False):
    """
    Compare PCP-derived counts with aamut_fitness actual_count by aa_mutation.

    Returns: (n_matches, n_discrepancies, discrepancies_df)
    """
    if len(pcp_df) == 0 and len(aamut_df) == 0:
        return 0, 0, pd.DataFrame()

    # Get mutations with actual_count > 0 from aamut_fitness
    if "actual_count" in aamut_df.columns:
        aamut_with_counts = aamut_df[aamut_df["actual_count"] > 0].copy()
    else:
        aamut_with_counts = pd.DataFrame()

    # Prepare merge keys - now using aa_mutation instead of just aa_site
    if is_aggregate:
        # For aggregate, we don't have clade column
        pcp_agg = pcp_df[pcp_df["clade"] == "all"].copy()
        merge_keys = ["gene", "aa_mutation"]

        if len(pcp_agg) > 0 and len(aamut_with_counts) > 0:
            merged = pcp_agg.merge(
                aamut_with_counts[["gene", "aa_mutation", "actual_count"]],
                on=merge_keys,
                how="outer",
            )
        elif len(pcp_agg) > 0:
            merged = pcp_agg.copy()
            merged["actual_count"] = 0
        elif len(aamut_with_counts) > 0:
            merged = aamut_with_counts[["gene", "aa_mutation", "actual_count"]].copy()
            merged["pcp_count"] = 0
        else:
            return 0, 0, pd.DataFrame()
    else:
        # For per-clade comparison
        pcp_per_clade = pcp_df[pcp_df["clade"] != "all"].copy()
        merge_keys = [clade_col, "gene", "aa_mutation"]

        if len(pcp_per_clade) > 0 and len(aamut_with_counts) > 0:
            merged = pcp_per_clade.merge(
                aamut_with_counts[[clade_col, "gene", "aa_mutation", "actual_count"]],
                on=merge_keys,
                how="outer",
            )
        elif len(pcp_per_clade) > 0:
            merged = pcp_per_clade.copy()
            merged["actual_count"] = 0
        elif len(aamut_with_counts) > 0:
            merged = aamut_with_counts[
                [clade_col, "gene", "aa_mutation", "actual_count"]
            ].copy()
            merged["pcp_count"] = 0
        else:
            return 0, 0, pd.DataFrame()

    # Fill NaN with 0
    merged["pcp_count"] = merged["pcp_count"].fillna(0).astype(int)
    merged["actual_count"] = merged["actual_count"].fillna(0).astype(int)

    # Find discrepancies
    merged["difference"] = merged["pcp_count"] - merged["actual_count"]
    discrepancies = merged[merged["difference"] != 0].copy()

    n_matches = len(merged) - len(discrepancies)
    n_discrepancies = len(discrepancies)

    return n_matches, n_discrepancies, discrepancies


# Compare per-clade
log("\n=== Per-clade comparison ===")
clade_matches, clade_disc_count, clade_disc = compare_counts(
    pcp_counts, aamut_clade_df, clade_col="clade", is_aggregate=False
)
log(f"Per-clade matches: {clade_matches}")
log(f"Per-clade discrepancies: {clade_disc_count}")

# Compare aggregate
log("\n=== Aggregate comparison ===")
agg_matches, agg_disc_count, agg_disc = compare_counts(
    pcp_counts, aamut_all, is_aggregate=True
)
log(f"Aggregate matches: {agg_matches}")
log(f"Aggregate discrepancies: {agg_disc_count}")

# Combine discrepancies
all_discrepancies = pd.concat([clade_disc, agg_disc], ignore_index=True)
total_discrepancies = len(all_discrepancies)

# Categorize discrepancies
if total_discrepancies > 0:
    only_in_pcp = all_discrepancies[all_discrepancies["actual_count"] == 0]
    only_in_aamut = all_discrepancies[all_discrepancies["pcp_count"] == 0]
    count_mismatch = all_discrepancies[
        (all_discrepancies["pcp_count"] > 0) & (all_discrepancies["actual_count"] > 0)
    ]
else:
    only_in_pcp = pd.DataFrame()
    only_in_aamut = pd.DataFrame()
    count_mismatch = pd.DataFrame()

# Write report
with open(output_report, "w") as f:
    f.write("=" * 60 + "\n")
    f.write("PCP vs aamut_fitness Validation Report\n")
    f.write("=" * 60 + "\n\n")

    f.write("SUMMARY\n")
    f.write("-" * 40 + "\n")
    f.write(f"Gene: {pcp_gene_name}\n")
    if codon_map is not None:
        f.write(
            f"Compared against aamut_fitness gene: {fitness_gene} "
            f"(codon offset {declared_offset}; {pcp_gene_name} codon 1 == "
            f"{fitness_gene} codon {declared_offset + 1})\n"
        )
        f.write(f"Corresponding codons: {len(codon_map)}\n")
    else:
        f.write(f"Compared against aamut_fitness gene: {fitness_gene}\n")
    f.write(f"Per-clade exact matches: {clade_matches}\n")
    f.write(f"Per-clade discrepancies: {clade_disc_count}\n")
    f.write(f"Aggregate exact matches: {agg_matches}\n")
    f.write(f"Aggregate discrepancies: {agg_disc_count}\n")
    f.write(f"Total discrepancies: {total_discrepancies}\n\n")

    # Reported as its own category, never folded into the discrepancy counts:
    # these codons have no counterpart to disagree with. ORF1a's 4402-4405 run
    # past the -1 ribosomal frameshift in ORF1a's own frame.
    if unmapped_codons:
        f.write("CODONS WITH NO COUNTERPART IN THE FITNESS GENE\n")
        f.write("-" * 40 + "\n")
        shown = unmapped_codons if len(unmapped_codons) <= 20 else unmapped_codons[:20]
        f.write(
            f"{len(unmapped_codons)} of {len(codon_map) + len(unmapped_codons)} "
            f"{pcp_gene_name} codons occupy reference nucleotides that no "
            f"{fitness_gene} codon covers, so they are excluded from the\n"
            f"comparison rather than counted as discrepancies. This is expected "
            f"where a gene extends past the -1 ribosomal frameshift in its own "
            f"reading frame.\n"
        )
        f.write(f"Codons: {shown}")
        f.write(" (truncated)\n" if len(unmapped_codons) > len(shown) else "\n")
        f.write(
            f"PCP-derived count rows on those codons: "
            f"{len(pcp_counts_no_counterpart)}\n\n"
        )

    if total_discrepancies == 0:
        f.write("SUCCESS: All amino acid mutation counts match exactly!\n\n")
        f.write(
            "The PCP extraction and aamut_fitness pipelines are consistent.\n"
        )
    else:
        f.write("DISCREPANCY BREAKDOWN\n")
        f.write("-" * 40 + "\n")
        f.write(f"Only in PCPs (not in aamut_fitness): {len(only_in_pcp)}\n")
        f.write(f"Only in aamut_fitness (not in PCPs): {len(only_in_aamut)}\n")
        f.write(f"Count mismatches (both have counts): {len(count_mismatch)}\n\n")

        if len(only_in_pcp) > 0:
            f.write("\nAA mutations only in PCPs (first 10):\n")
            cols = ["clade", "gene", "aa_mutation", "pcp_count"] if "clade" in only_in_pcp.columns else ["gene", "aa_mutation", "pcp_count"]
            cols = [c for c in cols if c in only_in_pcp.columns]
            f.write(only_in_pcp.head(10)[cols].to_string())
            f.write("\n")

        if len(only_in_aamut) > 0:
            f.write("\nAA mutations only in aamut_fitness (first 10):\n")
            cols = ["clade", "gene", "aa_mutation", "actual_count"] if "clade" in only_in_aamut.columns else ["gene", "aa_mutation", "actual_count"]
            cols = [c for c in cols if c in only_in_aamut.columns]
            f.write(only_in_aamut.head(10)[cols].to_string())
            f.write("\n")

        if len(count_mismatch) > 0:
            f.write("\nCount mismatches (first 10):\n")
            cols = ["clade", "gene", "aa_mutation", "pcp_count", "actual_count", "difference"] if "clade" in count_mismatch.columns else ["gene", "aa_mutation", "pcp_count", "actual_count", "difference"]
            cols = [c for c in cols if c in count_mismatch.columns]
            f.write(count_mismatch.head(10)[cols].to_string())
            f.write("\n")

    # Add notes section about expected differences
    f.write("\n" + "=" * 60 + "\n")
    f.write("NOTES ON EXPECTED DIFFERENCES\n")
    f.write("=" * 60 + "\n")

    f.write("\nThe following types of discrepancies are EXPECTED:\n")
    f.write("1. Mutations only in PCP counts:\n")
    f.write("   - Likely reversions to clade founder (filtered from matUtils)\n")
    f.write("   - These are excluded from fitness estimation by design\n")

    f.write("\n2. Mutations only in matUtils counts:\n")
    f.write("   - Should be rare after stop codon filter fix\n")
    f.write("   - Investigate if any appear\n")

    f.write("\n3. Count mismatches where matUtils > PCP:\n")
    f.write("   - Likely due to terminal/non-terminal splits in matUtils\n")
    f.write("   - Differences typically ±1-2 for small number of mutations\n")

    f.write("\n4. Count mismatches where PCP > matUtils:\n")
    f.write("   - Likely synonymous mutations with bidirectional changes\n")
    f.write("   - One direction is reversion (filtered from matUtils)\n")

    f.write("\nSee README.md 'PCP vs matUtils Validation' section for details.\n")
    f.write("=" * 60 + "\n")

    f.write("\n" + "=" * 60 + "\n")
    f.write("End of Report\n")
    f.write("=" * 60 + "\n")

# Write discrepancies to CSV
all_discrepancies.to_csv(output_discrepancies, index=False)

log(f"\nWrote report to {output_report}")
log(f"Wrote {len(all_discrepancies)} discrepancies to {output_discrepancies}")

# Report validation status (don't exit with error - just report)
if total_discrepancies > 0:
    log(f"\nValidation found {total_discrepancies} discrepancies")
else:
    log("\nValidation passed - all counts match!")
