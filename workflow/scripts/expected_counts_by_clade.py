#!/usr/bin/env python3
"""Add GLM-expected counts to mutation counts for each clade.

Synchronized with SARS2-mut-fitness-v2 approach:
1. Train GLM on curated training data with L2 regularization
2. For each clade, compute clade-specific normalization factor T
3. Apply normalized expectations to all mutations in the clade

Since #137 there is **one fit per era**, not one global fit, and each clade is
served by its own era's model. The eras with training data arrive as named
inputs ``curated_<era>``; which of them exist is decided at DAG-construction
time (``Snakefile: eras_present``), so a single-era tree is fed exactly one
file and this script then reduces *exactly* to the v9 pooled behaviour --
same training frame, same fit, same output. That degeneracy is what makes the
test-data run (three pre-Omicron clades) a regression check on the refactor
rather than a new code path with no baseline.
"""

import os
import pandas as pd
from importlib import util as importlib_util

# Get snakemake variables
counts_csv = snakemake.input.counts
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


log("Adding GLM-expected counts to mutation counts (v2-synchronized, per era)")

# Determine workflow directory
script_dir = os.path.dirname(os.path.abspath(__file__))
workflow_dir = os.path.dirname(script_dir)

# Load workflow modules
rates_module = load_workflow_module("rates", workflow_dir)
load_module = load_workflow_module("load", workflow_dir)
era_mod = load_workflow_module("era", workflow_dir)

# Load counts
counts_df = pd.read_csv(counts_csv)
log(f"Loaded {len(counts_df)} count entries")

# Load one curated training frame per era with data. ERAS is iterated (rather
# than the input keys) so the order is the module's, not the DAG's.
# load_curated_muts filters to synonymous mutations -- the whole neutral-rate
# restriction. Do NOT replace it with pd.read_csv (#126).
input_keys = set(snakemake.input.keys())
train_by_era = {}
for era_name in era_mod.ERAS:
    key = f"curated_{era_name}"
    if key in input_keys:
        train_by_era[era_name] = load_module.load_curated_muts(snakemake.input[key])
        log(
            f"Loaded {len(train_by_era[era_name])} curated mutations for era "
            f"{era_name!r} from {snakemake.input[key]}"
        )
unexpected = sorted(
    k for k in input_keys
    if k.startswith("curated_") and k[len("curated_"):] not in era_mod.ERAS
)
if unexpected:
    raise ValueError(
        f"input(s) {unexpected} look like per-era curated tables but name no "
        f"known era ({era_mod.ERAS}). Silently ignoring one would drop an era "
        f"from the fit while still producing a full output file."
    )
if not train_by_era:
    raise ValueError(
        f"no per-era curated training inputs were supplied (inputs: "
        f"{sorted(input_keys)}). Expected at least one of "
        f"{[f'curated_{e}' for e in era_mod.ERAS]}."
    )
log(f"Eras with training data: {sorted(train_by_era)}")

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

# Report the dispatch before running it: a misrouted clade produces a plausible
# expected_count, so which model served which clade has to be inspectable.
for era_name, era_clades in sorted(era_mod.partition_clades(clades).items()):
    log(f"  era {era_name!r} serves {len(era_clades)} clades: {era_clades}")

# One fit per era, each clade served by its own era's model. era_for_clade
# raises on an unmapped clade -- never a positional slice (#137), never a
# default era.
rates_module.add_expected_count_by_era(
    train_by_era=train_by_era,
    count_df=counts_df,
    clades=clades,
    era_of_clade=era_mod.era_for_clade,
    lightswitch_boundaries=lightswitch_boundaries,
    log=log,
)

# Every row must have been served by some era's model. A clade whose rows were
# never touched keeps a null expected_count, which downstream reads as a
# legitimate zero-ish expectation rather than an omission (#129/#130 shape).
n_missing = int(counts_df["expected_count"].isna().sum())
if n_missing:
    unserved = sorted(
        counts_df.loc[counts_df["expected_count"].isna(), "clade"].unique().tolist()
    )
    raise ValueError(
        f"{n_missing} count rows have no expected_count after per-era dispatch "
        f"(clades: {unserved}). Every clade present must be served by exactly "
        f"one era's fitted model."
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
