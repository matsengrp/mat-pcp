#!/usr/bin/env python3
"""Windowed expected counts + ``delta_fitness_ratio`` (EPIC #115 Phase 3, #118).

Thin Snakemake wrapper around :mod:`workflow.modules.expected_counts`. Consumes
the Phase-2 sliding-window outputs (``windows/founder_muts.csv`` +
``windows/subtree_counts.csv``) plus the GLM ``master_table.csv``, the reference
FASTA, ``ref/coding_sites.csv``, and the RNA-structure data, and emits a
mutation x time-window ``delta_fitness_ratio`` trajectory
(``windows/fitness_trajectory.csv``).

Expected counts come *only* from mat-pcp's GLM (via the master table), scaled per
window by the ``evol_time`` factor ``T`` (mat-pcp's clade-path normalization,
``rates.py``); there is no synonymous-average fallback. ``T`` divides out the
rate's absolute scale, so the raw ``rate`` is used and synonymous calibration
holds by construction.
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

def run(
    *,
    founder_muts_csv,
    subtree_counts_csv,
    master_table_csv,
    coding_sites_csv,
    ref_fasta,
    site_mask_csv,
    rna_structure_file,
    out_csv,
    sites_to_exclude,
    site_include_range,
    fitness_pseudocount,
    min_subtree_size,
    log=print,
):
    """Compute and write the windowed fitness trajectory."""
    workflow_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    expected_counts = load_workflow_module("expected_counts", workflow_dir)
    masking = load_workflow_module("masking", workflow_dir)

    ref_seq = expected_counts.load_reference_from_fasta(ref_fasta)
    genome_len = len(ref_seq)
    log(f"Reference genome length: {genome_len}")

    # mat-pcp's GLM per-context rates, used at native scale (the per-window
    # evol_time factor T divides out any constant rate factor).
    #
    # This is the POOLED master table, deliberately, while the clade path moved to
    # per-era tables in v10 (#137, EPIC #138 decision 7). The divergence is
    # documented rather than fixed because the windowed path has no clade column to
    # dispatch on: ``subtree_counts.csv`` / ``founder_muts.csv`` key on
    # ``founder_id``, a bare node id. A per-window date label is wrong exactly
    # across the pre-Omicron/Omicron transition -- the window that straddles it
    # contains both -- and a per-founder label would break the shared-base
    # enumeration that took this rule from ~20 h to 1.45 h in #123.
    #
    # What leaks through is bounded: T_window is fitted per window so it absorbs
    # any constant rate factor, leaving only the split's *relative* structure,
    # which is dominated by G->T (+0.62 mean log rate ratio between eras against
    # +0.27 for the next largest). Reading a per-era table here without a correct
    # era label would be worse than reading the pooled one, since it would apply
    # one era's relative structure to windows from the other. Era dispatch for the
    # windowed path is a follow-up, not a v10 item.
    rates_df = expected_counts.prepare_rates_df(master_table_csv)

    vcf_sites = pd.read_csv(site_mask_csv)["site"].tolist()
    bad_sites = masking.global_masked_sites(
        sites_to_exclude=sites_to_exclude,
        vcf_masked_sites=vcf_sites,
        site_include_range=site_include_range,
        genome_len=genome_len,
    )
    log(f"Globally masked (bad) sites: {len(bad_sites)}")

    basepair_df = expected_counts.build_site_unpaired_map(
        rna_structure_file, genome_len
    )
    ec = expected_counts.ExpectedCounts(
        ref_seq, coding_sites_csv, bad_sites, basepair_df, rates_df
    )

    founder_muts_df = pd.read_csv(founder_muts_csv).fillna(
        {"founder_nt_mutations": ""}
    )
    subtree_counts_df = pd.read_csv(subtree_counts_csv)
    log(
        f"Founder-windows: {len(founder_muts_df)}; "
        f"actual-count rows: {len(subtree_counts_df)}"
    )

    traj = expected_counts.compute_fitness_trajectory(
        founder_muts_df,
        subtree_counts_df,
        ec,
        bad_sites,
        pseudocount=fitness_pseudocount,
        min_subtree_size=min_subtree_size,
        log=log,
    )
    traj.to_csv(out_csv, index=False)

    if len(traj):
        syn = traj[traj["synonymous"]]
        syn_ratio = (
            syn["actual_count"].sum() / syn["expected_count"].sum()
            if syn["expected_count"].sum() != 0
            else float("nan")
        )
        log(
            f"Wrote {len(traj)} trajectory rows to {out_csv}. Synonymous "
            f"calibration sum(actual)/sum(expected)={syn_ratio:.4f} (should be ~1 by "
            f"construction); synonymous delta_fitness_ratio median="
            f"{syn['delta_fitness_ratio'].median():.4f} (n={len(syn)})."
        )
    else:
        log(f"Wrote empty trajectory to {out_csv}")


def main():
    """Entry point when run as a Snakemake ``script:``."""
    log_file = (
        snakemake.log[0]  # noqa: F821
        if hasattr(snakemake, "log") and snakemake.log  # noqa: F821
        else None
    )

    def log(msg):
        if log_file:
            with open(log_file, "a") as f:
                f.write(str(msg) + "\n")
        print(msg)

    run(
        founder_muts_csv=snakemake.input.founder_muts,  # noqa: F821
        subtree_counts_csv=snakemake.input.subtree_counts,  # noqa: F821
        master_table_csv=snakemake.input.master_table,  # noqa: F821
        coding_sites_csv=snakemake.input.coding_sites,  # noqa: F821
        ref_fasta=snakemake.input.ref_fasta,  # noqa: F821
        site_mask_csv=snakemake.input.site_mask,  # noqa: F821
        rna_structure_file=snakemake.input.rna_structure,  # noqa: F821
        out_csv=snakemake.output.trajectory,  # noqa: F821
        sites_to_exclude=snakemake.params.sites_to_exclude,  # noqa: F821
        site_include_range=snakemake.params.site_include_range,  # noqa: F821
        fitness_pseudocount=snakemake.params.fitness_pseudocount,  # noqa: F821
        min_subtree_size=snakemake.params.min_subtree_size,  # noqa: F821
        log=log,
    )


if __name__ == "__main__":
    main()
