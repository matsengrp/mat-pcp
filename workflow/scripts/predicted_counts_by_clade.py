#!/usr/bin/env python3
"""Add GLM-predicted counts to mutation counts for each clade."""

import sys
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


def load_workflow_module(module_name, workflow_dir=None):
    """Load a module from the workflow's modules directory."""
    if workflow_dir is None:
        workflow_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    module_path = os.path.join(workflow_dir, "modules", f"{module_name}.py")

    if not os.path.exists(module_path):
        raise ImportError(f"Cannot find module '{module_name}' at {module_path}")

    spec = importlib_util.spec_from_file_location(module_name, module_path)
    module = importlib_util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


log("Adding GLM-predicted counts to mutation counts")

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

# Initialize and train rates model
rate = rates_module.Rates(lightswitch_boundaries)
rate.populate_rates(counts_curated)

log("Trained GLM on curated mutations")

# Ensure counts_df has required columns for prediction
required_cols = ["nt_mutation", "motif", "unpaired", "nt_site_before_boundary"]
missing_cols = [col for col in required_cols if col not in counts_df.columns]
if missing_cols:
    log(f"WARNING: Missing columns for prediction: {missing_cols}")
    # Add default values for missing columns
    if "nt_site_before_boundary" not in counts_df.columns:
        counts_df["nt_site_before_boundary"] = False

# Calculate predicted counts for each row in counts_df
rate.predicted_counts_by_clade(counts_df)

# Add tau_squared from rate table lookup
rate_lookup = {}
for _, row in rate.rates.iterrows():
    key = (
        row["mut_type"],
        row["motif"],
        row["unpaired"],
        row["nt_site_before_boundary"],
    )
    rate_lookup[key] = row.get("residual", 0.1) if pd.notna(row.get("residual")) else 0.1


def get_tau_squared(row):
    key = (
        row.get("mut_type", ""),
        row.get("motif", ""),
        int(row.get("unpaired", 0)),
        bool(row.get("nt_site_before_boundary", False)),
    )
    return rate_lookup.get(key, 0.1)


counts_df["tau_squared"] = counts_df.apply(get_tau_squared, axis=1)

log(f"Added predicted counts and tau_squared to {len(counts_df)} entries")
log(
    f"Predicted count range: {counts_df['predicted_count'].min():.2f} - "
    f"{counts_df['predicted_count'].max():.2f}"
)
log(
    f"Tau squared range: {counts_df['tau_squared'].min():.4f} - "
    f"{counts_df['tau_squared'].max():.4f}"
)

counts_df.to_csv(output_csv, index=False)
log(f"Wrote output to {output_csv}")
log("Done!")
