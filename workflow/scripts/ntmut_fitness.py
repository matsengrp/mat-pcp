#!/usr/bin/env python3
"""Compute probabilistic fitness estimates for nucleotide mutations."""

import sys
import os
import pandas as pd
from importlib import util as importlib_util

# Get snakemake variables
counts_csv = snakemake.input.counts
output_csv = snakemake.output.csv
N_f = snakemake.params.get("N_f", 300)
N_n_exp = snakemake.params.get("N_n_exp", 300)
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


log("Computing probabilistic fitness estimates")

# Determine workflow directory
script_dir = os.path.dirname(os.path.abspath(__file__))
workflow_dir = os.path.dirname(script_dir)

# Load workflow modules
probfit_module = load_workflow_module("probfit", workflow_dir)

# Load counts
counts_df = pd.read_csv(counts_csv)
log(f"Loaded {len(counts_df)} count entries")

# Check required columns
required_cols = ["actual_count", "expected_count", "tau_squared"]
missing_cols = [col for col in required_cols if col not in counts_df.columns]
if missing_cols:
    log(f"ERROR: Missing required columns: {missing_cols}")
    raise ValueError(f"Missing columns: {missing_cols}")

# Compute fitness estimates
log(f"Computing fitness with N_f={N_f}, N_n_exp={N_n_exp}")
probfit_module.add_probabilistic_estimates(counts_df, N_f=N_f, N_n_exp=N_n_exp)

# Log summary
log(
    f"Fitness range: f_mean [{counts_df['f_mean'].min():.2f}, "
    f"{counts_df['f_mean'].max():.2f}]"
)
log(
    f"Fitness std range: [{counts_df['f_st_dev'].min():.2f}, "
    f"{counts_df['f_st_dev'].max():.2f}]"
)

counts_df.to_csv(output_csv, index=False)
log(f"Computed fitness estimates for {len(counts_df)} mutations")
log(f"Wrote output to {output_csv}")
log("Done!")
