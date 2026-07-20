#!/usr/bin/env python3
"""Extract nucleotide mutation counts for a specific clade or cluster."""

import pandas as pd

# Get snakemake variables
counts_csv = snakemake.input.counts
output_csv = snakemake.output.csv
clade = snakemake.params.clade
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log(f"Extracting nucleotide mutation counts for clade: {clade}")

# Load predicted counts
counts_df = pd.read_csv(counts_csv)
log(f"Loaded {len(counts_df)} count entries")

# Apply masking before computing fitness. The `exclude` flag (built once in
# build_count_table via modules/masking.py) is a superset of `masked_in_usher`
# and also covers out-of-range sites, config `sites_to_exclude`, and
# founder->reference reversions. Dropping these per (clade, mutation) rows here
# keeps their artifactual actual counts out of both the per-clade fitness and
# the "all" cross-clade aggregation below (a site masked only in one clade still
# contributes from clades where it is not masked). The GLM/expected-count
# machinery is unaffected: it trains on separately curated counts.
#
# `exclude` is guaranteed present by build_count_table and must survive
# annotate_counts/expected_counts_by_clade untouched; its absence means masking
# would be silently skipped, so we fail hard rather than warn-and-continue.
if "exclude" not in counts_df.columns:
    raise KeyError(
        "ntmut_counts_by_clade: required 'exclude' column missing from the "
        "counts table. It is written by build_count_table (via "
        "modules/masking.py) and must survive annotate_counts and "
        "expected_counts_by_clade. Its absence would silently skip masking and "
        "reintroduce artifactual counts at masked sites into the fitness "
        "estimates."
    )

# `exclude` must be a clean boolean flag (only True/False, no NaN); otherwise the
# masking filter could silently keep rows it should drop. A stray NaN (e.g. from
# a bad upstream merge) or string values would slip past a naive comparison.
exclude_values = set(counts_df["exclude"].dropna().unique())
if counts_df["exclude"].isna().any() or not exclude_values.issubset({True, False}):
    raise ValueError(
        "ntmut_counts_by_clade: 'exclude' column must contain only booleans "
        f"(True/False); found values {sorted(exclude_values, key=str)} "
        f"(dtype={counts_df['exclude'].dtype}, "
        f"{int(counts_df['exclude'].isna().sum())} NaN). Cannot reliably mask."
    )

excluded_mask = counts_df["exclude"] == True
n_excluded = int(excluded_mask.sum())
counts_df = counts_df[~excluded_mask].copy()
# Drop the masking flags now that filtering is done: every surviving row has
# exclude == False (and therefore masked_in_usher == False), so carrying them
# downstream would only add two constant, uninformative columns to the fitness
# outputs.
counts_df = counts_df.drop(columns=["exclude", "masked_in_usher"], errors="ignore")
log(f"Dropped {n_excluded} excluded (masked) mutations; {len(counts_df)} remain")

if clade == "all":
    # Aggregate counts across all clades into one row per nucleotide mutation.
    #
    # We key strictly on mutation identity (nt_site, nt_mutation) and SUM the
    # counts / expected counts across clades. We deliberately do NOT group by the
    # annotation columns: several of them (motif, clade_founder_codon/aa,
    # mutant_codon/aa, aa_mutation, synonymous, four_fold_degenerate, mut_class)
    # derive from the clade *founder* sequence and so can differ between clades
    # for the very same nucleotide mutation. Grouping by them would split one
    # mutation into several "all" rows (one per context) and fragment its pooled
    # signal. For those columns we keep the first contributing clade's value as a
    # representative annotation; the pooled row is a nucleotide-level summary, and
    # the per-clade files retain the exact per-clade annotations.
    log("Aggregating counts across all clades")

    key_cols = ["nt_site", "nt_mutation"]

    agg_dict = {}
    for col in [
        "actual_count",
        "count_terminal",
        "count_non_terminal",
        "expected_count",
    ]:
        if col in counts_df.columns:
            agg_dict[col] = "sum"

    # For tau_squared, take mean (residual variance is similar across clades
    # for the same mutation context)
    for col in ["tau_squared"]:
        if col in counts_df.columns:
            agg_dict[col] = "mean"

    # Every other column (except the per-row `clade`, which is meaningless once
    # pooled) is a representative annotation carried through via "first".
    for col in counts_df.columns:
        if col in key_cols or col in agg_dict or col == "clade":
            continue
        agg_dict[col] = "first"

    clade_counts = counts_df.groupby(key_cols, as_index=False).agg(agg_dict)
    clade_counts["cluster"] = "all"
else:
    # Filter to specific clade
    log(f"Filtering to clade: {clade}")
    clade_counts = counts_df[counts_df["clade"] == clade].copy()
    clade_counts["cluster"] = clade

log(f"Extracted {len(clade_counts)} mutations for cluster {clade}")

clade_counts.to_csv(output_csv, index=False)
log(f"Wrote output to {output_csv}")
log("Done!")
