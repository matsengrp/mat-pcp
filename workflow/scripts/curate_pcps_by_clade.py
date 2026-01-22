#!/usr/bin/env python3
"""Filter PCPs for a single clade to remove invalid entries for DASM training."""

import pandas as pd
import sys
from pathlib import Path

# Get inputs from snakemake
pcps_file = snakemake.input.pcps
output_curated = snakemake.output.curated_pcps
output_filtered = snakemake.output.filtered_branches
clade = snakemake.wildcards.clade
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None

# Import pcp_validation module
sys.path.insert(0, str(Path(snakemake.scriptdir).parent / "modules"))
from pcp_validation import validate_pcp_df  # noqa: E402


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


def main():
    """Curate PCPs for a single clade."""
    log(f"Curating PCPs for clade: {clade}")

    # Load PCPs
    pcp_df = pd.read_csv(pcps_file)
    log(f"Loaded {len(pcp_df)} PCPs")

    if len(pcp_df) == 0:
        # Empty input - create empty outputs
        pcp_df.to_csv(output_curated, index=False)
        pd.DataFrame(columns=["child_name"]).to_csv(output_filtered, index=False)
        log("No PCPs to curate")
        return

    # Validate and filter
    valid_mask, invalid_reasons = validate_pcp_df(pcp_df)

    # Save valid PCPs
    valid_pcps = pcp_df[valid_mask].reset_index(drop=True)
    valid_pcps.to_csv(output_curated, index=False)

    # Log removal reasons
    log(f"Valid PCPs: {len(valid_pcps)}/{len(pcp_df)}")
    for reason, indices in invalid_reasons.items():
        if len(indices) > 0:
            log(f"  {reason}: {len(indices)}")

    # Save filtered branch IDs (child_name) for mutation counting
    filtered_df = pcp_df[~valid_mask]
    filtered_df[["child_name"]].to_csv(output_filtered, index=False)
    log(f"Filtered branches: {len(filtered_df)}")


if __name__ == "__main__":
    main()
