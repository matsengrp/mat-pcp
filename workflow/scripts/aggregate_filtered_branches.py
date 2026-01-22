#!/usr/bin/env python3
"""Aggregate per-clade filtered branch IDs and create curation report."""

import os
import pandas as pd

# Get inputs from snakemake
curated_pcps_file = snakemake.input.curated_pcps
raw_pcp_files = snakemake.input.raw_pcps
filtered_branch_files = snakemake.input.filtered_branches
output_report = snakemake.output.report
output_filtered = snakemake.output.filtered_branches
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


def main():
    """Aggregate filtered branches and create report."""
    log("=" * 60)
    log("Aggregating filtered branches")
    log("=" * 60)

    # Count raw and curated PCPs
    total_raw = 0
    for f in raw_pcp_files:
        df = pd.read_csv(f)
        total_raw += len(df)

    curated_df = pd.read_csv(curated_pcps_file)
    total_curated = len(curated_df)

    log(f"Total raw PCPs: {total_raw}")
    log(f"Total curated PCPs: {total_curated}")
    log(f"Total filtered: {total_raw - total_curated}")

    # Aggregate filtered branches
    filtered_dfs = []
    for f in filtered_branch_files:
        clade = os.path.basename(f).replace("_filtered_branches.csv", "")
        df = pd.read_csv(f)
        if len(df) > 0:
            df["clade"] = clade
            filtered_dfs.append(df)
            log(f"  {clade}: {len(df)} filtered branches")

    if filtered_dfs:
        all_filtered = pd.concat(filtered_dfs, ignore_index=True)
    else:
        all_filtered = pd.DataFrame(columns=["child_name", "clade"])

    # Save aggregated filtered branches
    all_filtered.to_csv(output_filtered, index=False)
    log(f"Saved {len(all_filtered)} filtered branch IDs to {output_filtered}")

    # Create report
    with open(output_report, "w") as f:
        f.write("PCP Curation Report\n")
        f.write("=" * 40 + "\n\n")
        f.write(f"Input PCPs: {total_raw}\n")
        f.write(f"Valid PCPs: {total_curated}\n")
        f.write(f"Removed: {total_raw - total_curated}\n\n")
        f.write("Filtered branches by clade:\n")
        for clade in sorted(all_filtered["clade"].unique()) if len(all_filtered) > 0 else []:
            count = len(all_filtered[all_filtered["clade"] == clade])
            f.write(f"  {clade}: {count}\n")

    log(f"Saved report to {output_report}")
    log("=" * 60)


if __name__ == "__main__":
    main()
