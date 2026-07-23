"""Module for converting nucleotide fitness to amino acid fitness."""

import pandas as pd
import numpy as np


def get_coding(df, gene_overlaps, explode_cols):
    """
    Get coding mutations, handling gene overlaps appropriately.

    Args:
        df: DataFrame with nucleotide fitness data
        gene_overlaps: Dict with "retain" and "exclude" lists of gene pair overlaps
        explode_cols: List of columns to explode (semicolon-separated values)

    Returns:
        DataFrame with coding mutations, exploded by gene
    """
    overlaps_to_retain = [set(overlap) for overlap in gene_overlaps.get("retain", [])]
    overlaps_to_exclude = [set(overlap) for overlap in gene_overlaps.get("exclude", [])]

    def exclude_overlap(gene):
        if ";" not in str(gene):
            return False
        genes = set(gene.split(";"))
        if genes in overlaps_to_exclude:
            return True
        elif genes in overlaps_to_retain:
            return False
        else:
            # Default: exclude unknown overlaps
            return True

    # Filter to coding mutations
    if "noncoding" in df.columns:
        df_coding = df[df["noncoding"] == False].copy()
    else:
        # If no noncoding column, assume all are coding
        df_coding = df.copy()

    if len(df_coding) == 0:
        return df_coding

    # Handle overlapping genes
    if "codon_site" in df_coding.columns:
        df_coding["is_overlapping"] = df_coding["codon_site"].apply(
            lambda x: len(set(str(x).split(";"))) > 1
        )
        df_coding["overlap_to_exclude"] = df_coding["is_overlapping"] & df_coding[
            "gene"
        ].apply(exclude_overlap)
        df_coding = df_coding[~df_coding["overlap_to_exclude"]]

    # Extract clade founder and mutant amino acids
    if "aa_mutation" in df_coding.columns:
        df_coding["clade_founder_aa"] = df_coding["aa_mutation"].apply(
            lambda x: x[0]
            if ";" not in str(x)
            else ";".join([y[0] for y in str(x).split(";")])
        )
        df_coding["mutant_aa"] = df_coding["aa_mutation"].apply(
            lambda x: x[-1]
            if ";" not in str(x)
            else ";".join([y[-1] for y in str(x).split(";")])
        )

    # Explode semicolon-separated columns
    # Include clade_founder_aa and mutant_aa if they were just created
    cols_to_explode = list(explode_cols)
    for col in ["clade_founder_aa", "mutant_aa"]:
        if col in df_coding.columns and col not in cols_to_explode:
            cols_to_explode.append(col)

    for col in cols_to_explode:
        if col in df_coding.columns:
            df_coding[col] = df_coding[col].astype(str).str.split(";")

    df_coding_exp = df_coding.explode(cols_to_explode)

    # Remove ORF1a (keep only ORF1ab)
    if "gene" in df_coding_exp.columns:
        df_coding_exp = df_coding_exp[df_coding_exp["gene"] != "ORF1a"]

    return df_coding_exp


def aggregate_counts(df, explode_cols):
    """
    Aggregate counts for amino acid mutations.

    Args:
        df: DataFrame with coding mutations
        explode_cols: List of grouping columns

    Returns:
        DataFrame with aggregated counts
    """
    groupby_cols = ["cluster"] + list(explode_cols)
    groupby_cols = [c for c in groupby_cols if c in df.columns]

    # tau_squared is deliberately not aggregated here: it is a nucleotide-level
    # log-scale variance used only by probfit (nt fitness). It never feeds the
    # amino-acid delta_fitness_bayesian (inverse-variance weighting of nt estimates) or
    # delta_fitness_pseudo (from summed counts), so summing it into aa rows would
    # only carry a meaningless column. See CHANGELOG v5.
    agg_dict = {}
    for col in ["expected_count", "actual_count"]:
        if col in df.columns:
            agg_dict[col] = "sum"

    if not agg_dict:
        return df

    df_agg = df.groupby(groupby_cols, as_index=False).aggregate(agg_dict)

    # Rename codon_site to aa_site
    if "codon_site" in df_agg.columns:
        df_agg = df_agg.rename(columns={"codon_site": "aa_site"})
        df_agg["aa_site"] = (
            pd.to_numeric(df_agg["aa_site"], errors="coerce").fillna(0).astype(int)
        )

    # Rebuild aa_mutation from components if not already present
    if (
        "aa_mutation" not in df_agg.columns
        and "clade_founder_aa" in df_agg.columns
        and "mutant_aa" in df_agg.columns
        and "aa_site" in df_agg.columns
    ):
        df_agg["aa_mutation"] = (
            df_agg["clade_founder_aa"].astype(str)
            + df_agg["aa_site"].astype(str)
            + df_agg["mutant_aa"].astype(str)
        )

    return df_agg


def map_orf1ab_to_nsps(orf1ab_to_nsps):
    """
    Create mapping from ORF1ab amino acid sites to NSP sites.

    Args:
        orf1ab_to_nsps: Dict mapping NSP names to (start, end) amino acid ranges

    Returns:
        DataFrame with ORF1ab_site, nsp_site, nsp columns
    """
    if not orf1ab_to_nsps:
        return pd.DataFrame(columns=["ORF1ab_site", "nsp_site", "nsp"])

    dfs = []
    for nsp, coords in orf1ab_to_nsps.items():
        if coords is None:
            # Skip entries with None values (used to override defaults)
            continue
        start, end = coords
        df = pd.DataFrame(
            {
                "ORF1ab_site": list(range(start, end + 1)),
                "nsp_site": list(range(1, end - start + 2)),
                "nsp": nsp,
            }
        )
        dfs.append(df)

    return pd.concat(dfs, ignore_index=True).drop_duplicates()


def add_nsps(df, orf1ab_to_nsps_df):
    """
    Add NSP entries for ORF1ab mutations.

    Args:
        df: DataFrame with amino acid data
        orf1ab_to_nsps_df: DataFrame with ORF1ab to NSP mapping

    Returns:
        DataFrame with original data plus NSP-mapped versions of ORF1ab mutations
    """
    if len(orf1ab_to_nsps_df) == 0:
        df["subset_of_ORF1ab"] = False
        return df

    # Original data
    df_orig = df.assign(subset_of_ORF1ab=False)

    # ORF1ab mutations mapped to NSPs
    if "gene" in df.columns and "aa_site" in df.columns:
        orf1ab_data = df[df["gene"] == "ORF1ab"].merge(
            orf1ab_to_nsps_df,
            left_on="aa_site",
            right_on="ORF1ab_site",
            validate="many_to_one",
            how="left",
        )

        if len(orf1ab_data) > 0:
            # Drop original ORF1ab columns and use NSP values
            orf1ab_data = orf1ab_data.drop(
                columns=["gene", "aa_site", "ORF1ab_site"], errors="ignore"
            )
            orf1ab_data = orf1ab_data.rename(columns={"nsp": "gene", "nsp_site": "aa_site"})
            orf1ab_data["subset_of_ORF1ab"] = True

            # Rebuild aa_mutation if possible
            if (
                "clade_founder_aa" in orf1ab_data.columns
                and "mutant_aa" in orf1ab_data.columns
            ):
                orf1ab_data["aa_mutation"] = (
                    orf1ab_data["clade_founder_aa"].astype(str)
                    + orf1ab_data["aa_site"].astype(str)
                    + orf1ab_data["mutant_aa"].astype(str)
                )

            return pd.concat([df_orig, orf1ab_data], ignore_index=True)

    return df_orig


def naive_fitness(aa_counts, fitness_pseudocount=0.5):
    """
    Add naive (log ratio) fitness estimate.

    Args:
        aa_counts: DataFrame with expected_count and actual_count columns
        fitness_pseudocount: Pseudocount for numerical stability
    """
    if "expected_count" in aa_counts.columns and "actual_count" in aa_counts.columns:
        aa_counts["delta_fitness_pseudo"] = np.log(
            (aa_counts["actual_count"] + fitness_pseudocount)
            / (aa_counts["expected_count"] + fitness_pseudocount)
        )


def aa_fitness(df, explode_cols):
    """
    Compute weighted average fitness for amino acid mutations.

    Weights by inverse variance (1/sigma^2) of nucleotide fitness estimates.

    Args:
        df: DataFrame with f_mean and f_st_dev columns
        explode_cols: Grouping columns

    Returns:
        DataFrame with delta_fitness_bayesian and uncertainty columns
    """
    if "f_st_dev" not in df.columns or "f_mean" not in df.columns:
        return pd.DataFrame()

    df = df.copy()
    df["weight"] = 1 / (df["f_st_dev"] ** 2 + 0.001)  # Add small constant for stability
    df["weighted_value"] = df["f_mean"] * df["weight"]

    groupby_cols = ["cluster"] + list(explode_cols)
    groupby_cols = [c for c in groupby_cols if c in df.columns]

    grouped_sums = df.groupby(groupby_cols).agg({"weight": "sum", "weighted_value": "sum"})

    grouped_sums["delta_fitness_bayesian"] = grouped_sums["weighted_value"] / grouped_sums["weight"]
    grouped_sums["uncertainty"] = np.sqrt(1 / grouped_sums["weight"])

    aa_fit = grouped_sums[["delta_fitness_bayesian", "uncertainty"]].reset_index()

    # Rename codon_site to aa_site
    if "codon_site" in aa_fit.columns:
        aa_fit = aa_fit.rename(columns={"codon_site": "aa_site"})
        aa_fit["aa_site"] = (
            pd.to_numeric(aa_fit["aa_site"], errors="coerce").fillna(0).astype(int)
        )

    # Rebuild aa_mutation from components if not already present
    if (
        "aa_mutation" not in aa_fit.columns
        and "clade_founder_aa" in aa_fit.columns
        and "mutant_aa" in aa_fit.columns
        and "aa_site" in aa_fit.columns
    ):
        aa_fit["aa_mutation"] = (
            aa_fit["clade_founder_aa"].astype(str)
            + aa_fit["aa_site"].astype(str)
            + aa_fit["mutant_aa"].astype(str)
        )

    return aa_fit


def merge_aa_df(aa_fit, aa_counts, explode_cols):
    """
    Merge amino acid fitness and counts DataFrames.

    Args:
        aa_fit: DataFrame with delta_fitness_bayesian and uncertainty
        aa_counts: DataFrame with counts
        explode_cols: Columns used for merging

    Returns:
        Merged DataFrame
    """
    merge_cols = ["cluster"] + list(explode_cols)

    # Handle aa_site rename
    if "codon_site" in merge_cols:
        merge_cols = [c if c != "codon_site" else "aa_site" for c in merge_cols]

    # Add subset_of_ORF1ab to merge columns if present
    if "subset_of_ORF1ab" in aa_fit.columns and "subset_of_ORF1ab" in aa_counts.columns:
        merge_cols.append("subset_of_ORF1ab")

    merge_cols = [c for c in merge_cols if c in aa_fit.columns and c in aa_counts.columns]

    aamut_fitness = pd.merge(
        aa_counts,
        aa_fit,
        on=merge_cols,
        how="inner",
    )
    return aamut_fitness
