#!/usr/bin/env python3
"""Curate counts for GLM training, pooled over clades or restricted to one era.

Two rules run this script (#137):

``rule curate_counts``
    ``params.era`` is ``None`` -> the v9 behaviour, pooling every clade into
    ``glm/curated_counts.csv``. Byte-identical to pre-#137 output.
``rule curate_counts_era``
    ``params.era`` is ``"pre_omicron"`` / ``"omicron"`` -> the same curation
    restricted to that era's clades, written to
    ``glm/curated_counts_{era}.csv``.

**The conserved-site set is computed genome-globally over every clade, before
any era restriction**, matching upstream (``curate_counts.py:30-39`` there;
upstream computes ``conserved_sites`` from all clades and only then queries by
era). Computing conservation *within* an era would change the training site set,
so the two era fits would be trained on different sites from each other and from
the pooled fit -- which silently breaks the degeneracy guarantee that makes the
pooled table a usable control arm, and loses parity with upstream. The era filter
is therefore applied to the *rows*, never to the site-set derivation.

The aggregation drops ``clade`` (rows are pooled across clades by ``groupby``),
which is why the partition has to happen here rather than downstream: by the time
``curated_counts.csv`` exists there is no column left to partition on.
"""

import os
import pandas as pd
from importlib import util as importlib_util

# Get snakemake variables
counts_csv = snakemake.input.counts
founder_csv = snakemake.input.clade_founder
output_csv = snakemake.output.curated
# None for the pooled table; an era name for a per-era table.
era = snakemake.params.get("era", None)
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

script_dir = os.path.dirname(os.path.abspath(__file__))
workflow_dir = os.path.dirname(script_dir)
era_mod = load_workflow_module("era", workflow_dir)


if era is None:
    log("Curating counts for GLM training (pooled over all clades)")
else:
    log(f"Curating counts for GLM training (era={era!r})")

# Load inputs
counts_df = pd.read_csv(counts_csv)
founder_df = pd.read_csv(founder_csv)

log(f"Loaded {len(counts_df)} count entries")
log(f"Loaded {len(founder_df)} founder entries")

# Identify conserved sites (same codon and motif across all clades).
# Computed over EVERY clade, before any era restriction -- see the module
# docstring on why this must not be done within era.
conserved_data = founder_df[
    (founder_df["codon"] == founder_df["ref_codon"])
    & (founder_df["motif"] == founder_df["ref_motif"])
]
site_counts = conserved_data["nt_site"].value_counts()
n_clades = len(founder_df["clade"].unique())
conserved_sites = site_counts[site_counts == n_clades].index.tolist()

log(f"Conserved sites (present in all {n_clades} clades): {len(conserved_sites)}")

# Sites to ignore: masked in UShER or otherwise excluded (config sites_to_exclude,
# ProblematicSites VCF, out-of-range, and founder->reference reversions all set
# `exclude`, so they are covered by the excluded_sites branch below). Also
# derived pre-partition, for the same reason as the conserved sites.
sites_to_ignore = []

if "masked_in_usher" in counts_df.columns:
    masked_sites = (
        counts_df[counts_df["masked_in_usher"] == True]["nt_site"].unique().tolist()
    )
    sites_to_ignore.extend(masked_sites)
    log(f"UShER masked sites: {len(masked_sites)}")

if "exclude" in counts_df.columns:
    excluded_sites = (
        counts_df[counts_df["exclude"] == True]["nt_site"].unique().tolist()
    )
    sites_to_ignore.extend(excluded_sites)
    log(f"Excluded sites: {len(excluded_sites)}")

sites_to_ignore = list(set(sites_to_ignore))
log(f"Total sites to ignore: {len(sites_to_ignore)}")

# Filter to conserved, non-ignored sites
curated = counts_df[
    counts_df["nt_site"].isin(conserved_sites)
    & ~counts_df["nt_site"].isin(sites_to_ignore)
].copy()

log(f"After site filtering: {len(curated)} entries")

# Restrict to one era's clades, if asked. Every clade present is bucketed via
# era.era_for_clade, which RAISES on an unmapped clade rather than defaulting --
# a clade fitted into the wrong era yields a plausible rate table, not an error.
if era is not None:
    if era not in era_mod.ERAS:
        raise ValueError(
            f"params.era is {era!r}, which is not one of {era_mod.ERAS}"
        )
    clades_present = sorted(c for c in curated["clade"].unique() if pd.notna(c))
    partition = era_mod.partition_clades(clades_present)
    era_clades = partition.get(era, [])
    if not era_clades:
        # Not the legitimate empty-era skip: that is decided at DAG-construction
        # time by requesting only the non-empty eras (Snakefile: eras_present),
        # so being asked for an era with no clades means the two disagree.
        raise ValueError(
            f"era {era!r} has no clades among those present ({clades_present}); "
            f"eras present: {sorted(partition)}. An absent era is meant to be "
            f"skipped by not requesting its table at all, so this rule running "
            f"means the DAG requested a table that cannot be fitted."
        )
    log(f"Era {era!r} covers {len(era_clades)} of {len(clades_present)} clades: {era_clades}")
    curated = curated[curated["clade"].isin(era_clades)].copy()
    log(f"After era filtering: {len(curated)} entries")

# Define columns to aggregate and columns to group by
sum_cols = ["actual_count", "count_terminal", "count_non_terminal"]
sum_cols = [col for col in sum_cols if col in curated.columns]

# Columns to not include in groupby (will be aggregated or dropped)
ignore_cols = sum_cols + ["clade", "exclude", "masked_in_usher"]
groupby_cols = [
    col for col in curated.columns if col not in ignore_cols and col in curated.columns
]

# Only aggregate if we have data
if len(curated) > 0:
    curated_agg = curated.groupby(groupby_cols, as_index=False).agg(
        {col: "sum" for col in sum_cols}
    )

    # Verify no duplicate mutations
    if "nt_mutation" in curated_agg.columns:
        n_dups = curated_agg["nt_mutation"].duplicated().sum()
        if n_dups > 0:
            log(f"WARNING: {n_dups} duplicate mutations found")
        else:
            log("No duplicate mutations found")

    # Drop exclusion columns (no longer needed after filtering)
    curated_agg.drop(columns=["exclude", "masked_in_usher"], errors="ignore", inplace=True)
else:
    curated_agg = pd.DataFrame()
    log("WARNING: No data remaining after curation")

# Save curated data
curated_agg.to_csv(output_csv, index=False)

log(f"Curated {len(curated_agg)} mutations from {len(counts_df)} original entries")
log(f"Wrote curated data to {output_csv}")
log("Done!")
