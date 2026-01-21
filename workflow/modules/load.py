"""File loading utilities for mutation count data.

Synchronized with SARS2-mut-fitness-v2/modules/load.py to ensure consistent
data filtering between pipelines.
"""

import numpy as np
import pandas as pd


def load_synonymous_muts(
    df_path,
    include_noncoding_muts=False,
    include_stop_tolerant_orfs=False,
    remove_orf9b=False,
):
    """
    Load synonymous mutations for GLM training.

    Matches SARS2-mut-fitness-v2 load_synonymous_muts function exactly.
    The GLM should be trained on synonymous mutations only.

    Args:
        df_path: Path to curated counts CSV
        include_noncoding_muts: Include non-coding mutations
        include_stop_tolerant_orfs: Include stop-tolerant ORF mutations
        remove_orf9b: Remove ORF9b mutations

    Returns:
        Filtered DataFrame of synonymous mutations
    """
    df = pd.read_csv(df_path)

    # Create masks for synonymous mutations, non-coding mutations,
    # stop-codon tolerant ORFs, and ORF9b
    mask_synonymous = (df["synonymous"] == True).values  # noqa: E712
    mask_noncoding = (df["noncoding"] == True).values  # noqa: E712
    stop_tolerant_orfs = ["ORF6", "ORF7a", "ORF7b", "ORF8", "ORF10"]
    pattern = "|".join(stop_tolerant_orfs)
    mask_tolerant = df["gene"].str.contains(pattern, na=False).values

    # Handle ORF9b masking
    if "clade_founder_aa" in df.columns and "mutant_aa" in df.columns:
        mask_orf9b = (
            (df["gene"] == "N;ORF9b")
            & (
                df["clade_founder_aa"].apply(lambda x: str(x)[0] if pd.notna(x) else "")
                == df["mutant_aa"].apply(lambda x: str(x)[0] if pd.notna(x) else "")
            )
        ).values
    else:
        mask_orf9b = np.zeros(len(df), dtype=bool)

    # Apply masks
    if include_noncoding_muts:
        mask_synonymous = mask_synonymous | mask_noncoding
    if include_stop_tolerant_orfs:
        mask_synonymous = mask_synonymous | mask_tolerant
    if remove_orf9b:
        mask_synonymous = mask_synonymous | mask_orf9b

    return df[mask_synonymous]


def load_curated_muts(
    df_path,
    include_noncoding_muts=False,
    include_stop_tolerant_orfs=False,
    remove_orf9b=False,
):
    """
    Load curated mutations for GLM training.

    This is an alias for load_synonymous_muts to maintain backward compatibility
    while matching the v2 behavior of training on synonymous mutations only.

    Args:
        df_path: Path to curated counts CSV
        include_noncoding_muts: Include non-coding mutations
        include_stop_tolerant_orfs: Include stop-tolerant ORF mutations
        remove_orf9b: Remove ORF9b mutations

    Returns:
        Filtered DataFrame of synonymous mutations
    """
    return load_synonymous_muts(
        df_path, include_noncoding_muts, include_stop_tolerant_orfs, remove_orf9b
    )


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
