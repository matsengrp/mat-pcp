"""File loading utilities for mutation count data."""

import numpy as np
import pandas as pd


def load_curated_muts(
    df_path,
    include_noncoding_muts=False,
    include_stop_tolerant_orfs=False,
    remove_orf9b=False,
):
    """
    Load curated mutations for GLM training.

    The curated counts file has already been filtered to conserved sites
    and aggregated across clades, so all entries are suitable for training.

    Args:
        df_path: Path to curated counts CSV
        include_noncoding_muts: Include non-coding mutations
        include_stop_tolerant_orfs: Include stop-tolerant ORF mutations
        remove_orf9b: Remove ORF9b mutations

    Returns:
        Filtered DataFrame of curated mutations
    """
    df = pd.read_csv(df_path)

    # Start with all curated mutations (already filtered by curation step)
    mask = np.ones(len(df), dtype=bool)

    # Optionally filter by gene
    if "gene" in df.columns and remove_orf9b:
        mask_orf9b = df["gene"].str.contains("ORF9b", na=False)
        mask = mask & ~mask_orf9b

    return df[mask]


def load_nonsynonymous_muts(df_path):
    """Load non-synonymous mutations."""
    df = pd.read_csv(df_path)
    if "mut_class" in df.columns:
        mask_nonsynonymous = (df["mut_class"] == "nonsynonymous").values
    elif "synonymous" in df.columns:
        mask_nonsynonymous = (
            (df["synonymous"] == False)  # noqa: E712
            & (df.get("noncoding", pd.Series([False] * len(df))) == False)  # noqa
        ).values
    else:
        mask_nonsynonymous = np.ones(len(df), dtype=bool)
    return df[mask_nonsynonymous]


def load_nonexcluded_muts(df_path):
    """Load all non-excluded mutations."""
    return pd.read_csv(df_path)
