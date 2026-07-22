#!/usr/bin/env python3
"""Add GLM-expected counts to mutation counts for each clade.

Synchronized with SARS2-mut-fitness-v2 approach:
1. Train GLM on curated training data with L2 regularization
2. For each clade, compute clade-specific normalization factor T
3. Apply normalized expectations to all mutations in the clade
"""

import os
import pandas as pd
from importlib import util as importlib_util

# Get snakemake variables
counts_csv = snakemake.input.counts
curated_csv = snakemake.input.curated
output_csv = snakemake.output.csv
lightswitch_boundaries = snakemake.params.lightswitch_boundaries
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


log("Adding GLM-expected counts to mutation counts (v2-synchronized)")

# Determine workflow directory
script_dir = os.path.dirname(os.path.abspath(__file__))
workflow_dir = os.path.dirname(script_dir)

# Load workflow modules
rates_module = load_workflow_module("rates", workflow_dir)
load_module = load_workflow_module("load", workflow_dir)

# Load counts and curated training data
counts_df = pd.read_csv(counts_csv)
counts_curated = load_module.load_curated_muts(curated_csv)

log(f"Loaded {len(counts_df)} count entries")
log(f"Loaded {len(counts_curated)} curated mutations for training")

# Ensure counts_df has required columns for prediction
required_cols = ["mut_type", "motif", "unpaired", "nt_site_before_boundary", "clade"]
missing_cols = [col for col in required_cols if col not in counts_df.columns]
if missing_cols:
    log(f"WARNING: Missing columns for prediction: {missing_cols}")
    # Add default values for missing columns
    if "nt_site_before_boundary" not in counts_df.columns:
        counts_df["nt_site_before_boundary"] = False

# Get list of clades
clades = counts_df["clade"].unique()
log(f"Processing {len(clades)} clades: {list(clades)}")

# Use the v2-synchronized add_expected_count function
# This applies clade-specific normalization factor T
rates_module.add_expected_count(
    train_df=counts_curated,
    count_df=counts_df,
    clades=clades,
    lightswitch_boundaries=lightswitch_boundaries,
)

log(f"Added expected counts and tau_squared to {len(counts_df)} entries")
log(
    f"Expected count range: {counts_df['expected_count'].min():.2f} - "
    f"{counts_df['expected_count'].max():.2f}"
)
log(
    f"Tau squared range: {counts_df['tau_squared'].min():.4f} - "
    f"{counts_df['tau_squared'].max():.4f}"
)

counts_df.to_csv(output_csv, index=False)
log(f"Wrote output to {output_csv}")
log("Done!")
