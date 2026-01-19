#!/usr/bin/env python3
"""Aggregate clade-specific PCP files into a single output file."""

import os
import pandas as pd

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


def main():
    """Aggregate all clade PCP files into a single file with clade annotations."""
    log("=" * 60)
    log("Aggregating PCP files")
    log("=" * 60)

    dfs = []
    total_pcps = 0

    for pcp_file in pcp_files:
        # Extract clade name from filename
        clade = os.path.basename(pcp_file).replace("_pcps.csv", "")

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

    # Save aggregated file
    all_pcps_df.to_csv(output_file, index=False)

    log("-" * 60)
    log(f"Aggregated {len(dfs)} clade files into {output_file}")
    log(f"Total PCPs: {total_pcps}")
    log(f"Clades processed: {sorted(all_pcps_df['clade'].unique())}")
    log("=" * 60)


if __name__ == "__main__":
    main()
