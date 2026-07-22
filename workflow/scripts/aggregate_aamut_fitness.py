#!/usr/bin/env python3
"""Aggregate amino acid fitness across all clades."""

import os
import pandas as pd
import numpy as np

# Get snakemake variables
fitness_files = snakemake.input.fitness_files
output_csv = snakemake.output.csv
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log(f"Aggregating amino acid fitness from {len(fitness_files)} files")

# Load all fitness files
dfs = []
for fitness_file in fitness_files:
    if not os.path.exists(fitness_file):
        log(f"WARNING: File not found: {fitness_file}")
        continue

    df = pd.read_csv(fitness_file)
    # Extract clade from filename
    basename = os.path.basename(fitness_file)
    clade = basename.replace("_aamut_fitness.csv", "")
    df["source_clade"] = clade
    dfs.append(df)
    log(f"  {clade}: {len(df)} amino acid mutations")

if len(dfs) == 0:
    log("WARNING: No fitness files found")
    # Create empty DataFrame with expected columns
    aggregated = pd.DataFrame(
        columns=[
            "cluster",
            "gene",
            "clade_founder_aa",
            "mutant_aa",
            "aa_site",
            "aa_mutation",
            "expected_count",
            "actual_count",
            "naive_delta_fitness",
            "subset_of_ORF1ab",
            "delta_fitness",
            "uncertainty",
            "source_clade",
        ]
    )
else:
    # Combine all data
    all_data = pd.concat(dfs, ignore_index=True)
    log(f"Total amino acid mutations across clades: {len(all_data)}")

    # Group by mutation identity for aggregation
    groupby_cols = ["gene", "clade_founder_aa", "mutant_aa", "aa_site", "aa_mutation"]
    if "subset_of_ORF1ab" in all_data.columns:
        groupby_cols.append("subset_of_ORF1ab")

    groupby_cols = [c for c in groupby_cols if c in all_data.columns]

    # Compute aggregated counts and weighted average fitness
    # tau_squared is intentionally not carried into the aa outputs (see
    # aamutfit.aggregate_counts / CHANGELOG v5): it is a nucleotide-level
    # quantity that never feeds the aa fitness estimates.
    agg_dict = {}
    for col in ["expected_count", "actual_count"]:
        if col in all_data.columns:
            agg_dict[col] = "sum"

    if "delta_fitness" in all_data.columns and "uncertainty" in all_data.columns:
        # Compute weighted average fitness
        all_data["weight"] = 1 / (all_data["uncertainty"] ** 2 + 0.001)
        all_data["weighted_fitness"] = all_data["delta_fitness"] * all_data["weight"]
        agg_dict["weight"] = "sum"
        agg_dict["weighted_fitness"] = "sum"

    agg_dict["source_clade"] = lambda x: ";".join(sorted(set(x)))

    if groupby_cols:
        aggregated = all_data.groupby(groupby_cols, as_index=False).agg(agg_dict)

        # Compute aggregated fitness
        if "weighted_fitness" in aggregated.columns and "weight" in aggregated.columns:
            aggregated["delta_fitness"] = (
                aggregated["weighted_fitness"] / aggregated["weight"]
            )
            aggregated["uncertainty"] = np.sqrt(1 / aggregated["weight"])
            aggregated.drop(columns=["weight", "weighted_fitness"], inplace=True)

        # Recompute naive fitness
        if "actual_count" in aggregated.columns and "expected_count" in aggregated.columns:
            aggregated["naive_delta_fitness"] = np.log(
                (aggregated["actual_count"] + 0.5)
                / (aggregated["expected_count"] + 0.5)
            )
    else:
        aggregated = all_data

# Rename the two delta-fitness columns at the final write boundary only (#89).
# All computation above uses the internal names (``naive_delta_fitness``,
# ``delta_fitness``); only the aggregate output CSV -- the downstream contract
# for the dasm2 dashboard -- carries the self-documenting names. Applied here so
# both the no-data and populated paths emit the same schema. Per-clade files and
# all internal code keep the original names (the aggregator re-reads them).
aggregated = aggregated.rename(
    columns={
        "naive_delta_fitness": "delta_fitness_pseudo",
        "delta_fitness": "delta_fitness_bayesian",
    }
)

aggregated.to_csv(output_csv, index=False)

log(f"Aggregated {len(aggregated)} unique amino acid mutations")
log(f"Wrote output to {output_csv}")
log("Done!")
