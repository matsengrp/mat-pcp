#!/usr/bin/env python3
"""Compute amino acid mutation fitness estimates per clade."""

import os
import pandas as pd
from importlib import util as importlib_util

# Get snakemake variables
ntmut_fitness_csv = snakemake.input.ntmut_fitness
output_csv = snakemake.output.csv
orf1ab_to_nsps = snakemake.params.get("orf1ab_to_nsps", {})
gene_overlaps = snakemake.params.get("gene_overlaps", {"retain": [], "exclude": []})
genes = snakemake.params.get("genes", [])
fitness_pseudocount = snakemake.params.get("fitness_pseudocount", 0.5)
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


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


log("Computing amino acid mutation fitness")

# Determine workflow directory
script_dir = os.path.dirname(os.path.abspath(__file__))
workflow_dir = os.path.dirname(script_dir)

# Load workflow modules
aamutfit = load_workflow_module("aamutfit", workflow_dir)

# Load nucleotide fitness data
ntmut_fit = pd.read_csv(ntmut_fitness_csv)
log(f"Loaded {len(ntmut_fit)} nucleotide fitness entries")

# Columns to explode (semicolon-separated)
# Note: clade_founder_aa, mutant_aa are created in get_coding() from aa_mutation,
# so we only filter by what exists in input here, then add them after get_coding()
explode_cols = ["gene", "codon_site", "aa_mutation"]
explode_cols = [c for c in explode_cols if c in ntmut_fit.columns]

if len(ntmut_fit) == 0 or len(explode_cols) == 0:
    log("No data to process, creating empty output")
    aamut_fitness = pd.DataFrame(
        columns=[
            "cluster",
            "gene",
            "clade_founder_aa",
            "mutant_aa",
            "aa_site",
            "aa_mutation",
            "expected_count",
            "actual_count",
            "delta_fitness_ratio",
            "subset_of_ORF1ab",
            "delta_fitness_bayesian",
            "uncertainty",
        ]
    )
else:
    # Get coding mutations with gene overlaps handled
    ntmut_fit_coding = aamutfit.get_coding(ntmut_fit, gene_overlaps, explode_cols)
    log(f"Coding mutations after overlap handling: {len(ntmut_fit_coding)}")

    # Add clade_founder_aa and mutant_aa to explode_cols if they were created by get_coding()
    for col in ["clade_founder_aa", "mutant_aa"]:
        if col in ntmut_fit_coding.columns and col not in explode_cols:
            explode_cols.append(col)

    if len(ntmut_fit_coding) == 0:
        log("No coding mutations found")
        aamut_fitness = pd.DataFrame()
    else:
        # Aggregate counts for amino acid mutations
        aa_counts = aamutfit.aggregate_counts(ntmut_fit_coding, explode_cols)
        log(f"Aggregated to {len(aa_counts)} amino acid mutations")

        # Add the pseudocount log-ratio fitness estimate
        aamutfit.ratio_fitness(aa_counts, fitness_pseudocount=fitness_pseudocount)

        # Get refined fitness estimates (weighted average)
        aa_fit = aamutfit.aa_fitness(ntmut_fit_coding, explode_cols)
        log(f"Computed refined fitness for {len(aa_fit)} mutations")

        # Map ORF1ab to NSPs
        orf1ab_to_nsps_df = aamutfit.map_orf1ab_to_nsps(orf1ab_to_nsps)
        log(f"ORF1ab to NSP mapping: {len(orf1ab_to_nsps_df)} entries")

        # Add NSP entries
        aa_counts = aamutfit.add_nsps(aa_counts, orf1ab_to_nsps_df)
        aa_fit = aamutfit.add_nsps(aa_fit, orf1ab_to_nsps_df)

        # Merge counts and fitness
        aamut_fitness = aamutfit.merge_aa_df(aa_fit, aa_counts, explode_cols)

        # Order by gene and site if possible
        if genes and "gene" in aamut_fitness.columns:
            aamut_fitness["gene"] = pd.Categorical(
                aamut_fitness["gene"], categories=genes, ordered=True
            )
            if "aa_site" in aamut_fitness.columns:
                aamut_fitness = aamut_fitness.sort_values(
                    ["gene", "aa_site"]
                ).reset_index(drop=True)

aamut_fitness.to_csv(output_csv, index=False)

log(f"Computed fitness for {len(aamut_fitness)} amino acid mutations")
log(f"Wrote output to {output_csv}")
log("Done!")
