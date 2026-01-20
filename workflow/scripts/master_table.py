#!/usr/bin/env python3
"""Create master table with predicted mutation rates via GLM."""

import sys
import os
import numpy as np
import pandas as pd
from importlib import util as importlib_util

# Get snakemake variables
curated_counts_csv = snakemake.input.curated_counts
output_csv = snakemake.output.master_table
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


log("Creating master rate table via GLM")

# Determine workflow directory
script_dir = os.path.dirname(os.path.abspath(__file__))
workflow_dir = os.path.dirname(script_dir)

# Load workflow modules
rates_module = load_workflow_module("rates", workflow_dir)
load_module = load_workflow_module("load", workflow_dir)

# Load curated counts (all curated mutations for training)
log(f"Loading curated counts from {curated_counts_csv}")
counts = load_module.load_curated_muts(curated_counts_csv)
log(f"Loaded {len(counts)} curated mutations for training")

# Check required columns
required_cols = ["nt_mutation", "actual_count", "motif", "unpaired"]
missing = [c for c in required_cols if c not in counts.columns]
if missing:
    log(f"ERROR: Missing required columns: {missing}")
    log(f"Available columns: {list(counts.columns)}")
    sys.exit(1)

# Ensure nt_site_before_boundary exists
if "nt_site_before_boundary" not in counts.columns:
    log("Adding nt_site_before_boundary column (defaulting to False)")
    counts["nt_site_before_boundary"] = False

# Initialize rates object with boundaries
log(f"Light-switch boundaries: {lightswitch_boundaries}")
rate = rates_module.Rates(lightswitch_boundaries)

# Fit GLM and populate rates
log("Training GLM...")
rate.populate_rates(counts)

# Calculate genome composition
log("Calculating genome composition...")
rate.genome_composition(counts)

# Calculate predicted counts for training data
rate.predicted_counts_by_clade(counts)

# Calculate residual variance (tau) per mutation type
log("Computing residual variance...")
tau = (
    counts.groupby("mut_type")
    .apply(
        lambda x: np.mean(
            (np.log(x.actual_count + 0.5) - np.log(x.predicted_count + 0.5)) ** 2
        ),
        include_groups=False,
    )
    .to_dict()
)

# Add residual variance to rates
rate.residual_variance(counts, tau)

# Add lightswitch boundary positions to output
rate.rates["nt_site_boundary"] = 0
for mut_type, boundary in lightswitch_boundaries.items():
    rate.rates.loc[rate.rates.mut_type == mut_type, "nt_site_boundary"] = int(boundary)

# Select output columns
output_cols = [
    "mut_type",
    "motif",
    "unpaired",
    "nt_site_boundary",
    "nt_site_before_boundary",
    "rate",
    "predicted_count",
    "residual",
]
available_cols = [c for c in output_cols if c in rate.rates.columns]

# Save master table
rate.rates[available_cols].to_csv(output_csv, index=False)

log(f"Generated master table with {len(rate.rates)} rate entries")
log(f"Mutation types: {sorted(rate.rates['mut_type'].unique())}")
log(f"Rate range: {rate.rates['rate'].min():.4f} - {rate.rates['rate'].max():.4f}")
log(f"Wrote master table to {output_csv}")
log("Done!")
