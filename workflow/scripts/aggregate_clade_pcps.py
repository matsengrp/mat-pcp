#!/usr/bin/env python3
"""Aggregate clade-specific PCP files into a single output file."""

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


pcpio = load_workflow_module("pcpio")

# Get inputs from snakemake
pcp_files = snakemake.input.pcps
output_file = snakemake.output.all_pcps
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


def clade_from_pcp_path(pcp_file):
    """Derive a clade name from a per-clade PCP filename.

    Anchored on the suffix rather than the previous
    ``basename.replace("_pcps.csv", "")``. Since #134 these inputs are
    ``.csv.gz``, and the old form would have left the clade as ``"19A.gz"`` --
    a silently wrong label in every row of the aggregate, and in the ``clade``
    column that downstream per-clade grouping keys on. Nothing would have
    caught it: the row count is unchanged and every clade is wrong the same way.

    Raises rather than guessing if the name matches neither form.
    """
    basename = os.path.basename(pcp_file)
    for suffix in ("_pcps.csv.gz", "_pcps.csv"):
        if basename.endswith(suffix):
            return basename[: -len(suffix)]
    raise ValueError(
        f"PCP input {basename!r} ends in neither '_pcps.csv.gz' nor "
        f"'_pcps.csv'; cannot derive its clade"
    )


def main():
    """Aggregate all clade PCP files into a single file with clade annotations."""
    log("=" * 60)
    log("Aggregating PCP files")
    log("=" * 60)

    dfs = []
    total_pcps = 0

    for pcp_file in pcp_files:
        clade = clade_from_pcp_path(pcp_file)

        # Read PCP file
        df = pd.read_csv(pcp_file)
        df["clade"] = clade
        dfs.append(df)

        log(f"  {clade}: {len(df)} PCPs")
        total_pcps += len(df)

    if len(dfs) == 0:
        log("WARNING: No PCP files to aggregate!")
        all_pcps_df = pd.DataFrame(
            columns=[
                "parent_name",
                "child_name",
                "parent",
                "child",
                "branch_length",
                "clade",
            ]
        )
    else:
        all_pcps_df = pd.concat(dfs, ignore_index=True)

    # Reorder columns to put clade first
    cols = ["clade"] + [c for c in all_pcps_df.columns if c != "clade"]
    all_pcps_df = all_pcps_df[cols]

    # Save aggregated file (gzipped; see modules/pcpio.py, #134)
    all_pcps_df.to_csv(
        output_file, index=False, compression=pcpio.PCP_CSV_COMPRESSION
    )

    log("-" * 60)
    log(f"Aggregated {len(dfs)} clade files into {output_file}")
    log(f"Total PCPs: {total_pcps}")
    log(f"Clades processed: {sorted(all_pcps_df['clade'].unique())}")
    log("=" * 60)


if __name__ == "__main__":
    main()
