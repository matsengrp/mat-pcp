#!/usr/bin/env python3
"""Compute expected mutation counts from synonymous rates."""

import numpy as np
import pandas as pd

# Get snakemake variables
rates_by_clade_csv = snakemake.input.rates_by_clade
clade_founder_nts_csv = snakemake.input.clade_founder_nts_csv
output_csv = snakemake.output.expected_counts
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log("Computing expected mutation counts from synonymous rates")

# Load inputs
rates_by_clade = pd.read_csv(rates_by_clade_csv)
clade_founder_nts = pd.read_csv(clade_founder_nts_csv)

log(f"Loaded {len(rates_by_clade)} rate entries")
log(f"Loaded {len(clade_founder_nts)} clade founder entries")

expected_counts = []

for clade in clade_founder_nts["clade"].unique():
    founder = clade_founder_nts[clade_founder_nts["clade"] == clade]
    clade_rates = rates_by_clade[rates_by_clade["clade"] == clade]

    # Build rate lookup
    rate_lookup = {}
    for _, row in clade_rates.iterrows():
        if pd.notna(row["rate"]):
            rate_lookup[row["mut_type"]] = row["rate"]

    log(f"  {clade}: {len(rate_lookup)} valid rates, {len(founder)} sites")

    # Calculate expected counts for each possible mutation at each site
    for _, site_row in founder.iterrows():
        site = site_row["site"]
        from_nt = site_row["nt"]

        for to_nt in ["A", "C", "G", "T"]:
            if to_nt == from_nt:
                continue

            mut_type = from_nt + to_nt
            nt_mutation = f"{from_nt}{site}{to_nt}"

            # Get rate for this mutation type
            rate = rate_lookup.get(mut_type, np.nan)

            expected_counts.append(
                {
                    "clade": clade,
                    "nt_site": site,
                    "nt_mutation": nt_mutation,
                    "mut_type": mut_type,
                    "expected_count": rate if pd.notna(rate) else 0,
                    "gene": site_row["gene"],
                    "codon_site": site_row.get("codon_site", ""),
                    "codon_position": site_row.get("codon_position", ""),
                }
            )

expected_df = pd.DataFrame(expected_counts)
expected_df.to_csv(output_csv, index=False)

log(f"Computed expected counts for {len(expected_df)} mutation possibilities")
log("Done!")
