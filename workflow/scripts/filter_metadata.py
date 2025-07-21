#!/usr/bin/env python3

import pandas as pd
from datetime import datetime
from enum import Enum, auto


class DateResolution(Enum):
    DAY = auto()
    MONTH = auto()
    YEAR = auto()
    MISSING = auto()


def determine_resolution(date_str) -> DateResolution:
    """
    Determine the resolution of the given date string.

    Args:
    - date_str: A string representation of a date

    Returns:
    - A DateResolution enum value indicating the resolution of the date
    """
    if pd.isna(date_str) or date_str == "?":
        return DateResolution.MISSING
    elif "-" in date_str:
        parts = date_str.split("-")
        if len(parts) == 3:
            return DateResolution.DAY
        elif len(parts) == 2:
            return DateResolution.MONTH
    elif date_str.isdigit():  # Only year is provided
        return DateResolution.YEAR

    return DateResolution.MISSING  # default


def robust_fromisoformat(date_string: str) -> pd.Timestamp:
    """
    Convert a string into a datetime object. This function is designed to handle common variations in date formatting.

    Args:
        date_string (str): A date in string format.

    Returns:
        pd.Timestamp: Corresponding datetime object.

    Notes:
    - If the string is "?", it returns a NaT (Not a Time).
    - If the string format represents a month only (e.g., "2023-03"), it returns the approximate midpoint of the month (e.g., "2023-03-15").
    - If the string format represents a year only (e.g., "2023"), it returns the approximate midpoint of the year (e.g., "2023-06-01").
    """

    if pd.isna(date_string) or date_string == "?":
        return pd.NaT
    elif "-" not in date_string and len(date_string) == 4:  # Only a year is provided
        return datetime.fromisoformat(f"{date_string}-06-01")
    elif date_string.count("-") == 1:  # Only year and month are provided
        mid_date = f"{date_string}-15"
        if len(mid_date) == 10:
            return datetime.fromisoformat(mid_date)
        else:
            padded_date_string = "-".join(i.zfill(2) for i in mid_date.split("-"))
            return datetime.fromisoformat(padded_date_string)
    elif len(date_string) != 10:
        padded_date_string = "-".join(i.zfill(2) for i in date_string.split("-"))
        return datetime.fromisoformat(padded_date_string)
    else:
        return datetime.fromisoformat(date_string)


def main():
    # Define variables provided by snakemake workflow
    META = snakemake.input["meta"]

    DATETIME_CUTOFF = snakemake.params["datetime_cutoff"]
    SUBSET_FRAC = snakemake.params["subset_frac"]
    INVALID_DATE_RES = snakemake.params.qc_filter_options["invalid_date_resolutions"]
    MIN_DATE = snakemake.params.qc_filter_options["min_date"]

    VALID_HOSTS = snakemake.params.qc_filter_options["valid_hosts"]
    MISSING_HOST_IS_VALID = snakemake.params.qc_filter_options["missing_host_is_valid"]

    MIN_GENOME_LENGTH = snakemake.params.qc_filter_options["min_genome_length"]
    MAX_GENOME_LENGTH = snakemake.params.qc_filter_options["max_genome_length"]

    # OUTPUT
    FILTERED_META_OUT = snakemake.output["filtered_metadata"]
    INVALID_STRAINS_OUT = snakemake.log["invalid_strains"]

    # Load metadata
    print("Loading metadata...")
    meta_df = pd.read_csv(META, compression="gzip", sep="\t", low_memory=False)

    print(f"Original dataset size: {len(meta_df)}")

    # Apply random subsetting based on subset_frac early
    if SUBSET_FRAC < 1.0:
        subset_size = int(len(meta_df) * SUBSET_FRAC)
        meta_df = meta_df.sample(n=subset_size, random_state=42)
        print(f"After random subsetting: {len(meta_df)}")

    # Process dates
    print("Processing dates...")
    meta_df["date_resolution"] = meta_df["date"].apply(determine_resolution)
    meta_df["date_resolution_str"] = meta_df["date_resolution"].apply(lambda x: x.name)
    meta_df.drop(columns=["date_resolution"], inplace=True)

    # Convert dates
    meta_df["date"] = meta_df["date"].apply(robust_fromisoformat)

    # Apply datetime cutoff filtering right after date conversion
    datetime_cutoff = datetime.fromisoformat(DATETIME_CUTOFF)
    valid_dates = ~pd.isna(meta_df["date"])
    date_filtered_df = meta_df[valid_dates & (meta_df["date"] <= datetime_cutoff)]

    # Keep track of dropped samples due to date cutoff
    dropped_by_cutoff = len(meta_df) - len(date_filtered_df)
    print(f"Dropped {dropped_by_cutoff} samples due to date cutoff or invalid dates")

    # Replace the original dataframe with the filtered one
    meta_df = date_filtered_df
    print(f"Dataset size after date cutoff filtering: {len(meta_df)}")

    # QC Filter #1 - Date Resolution
    print("Applying QC Filter #1 - Date Resolution...")
    invalid_resolution = meta_df["date_resolution_str"].isin(INVALID_DATE_RES)
    invalid_time = meta_df["date"] < robust_fromisoformat(MIN_DATE)
    invalid_date = meta_df[invalid_resolution | invalid_time]
    qc_1_strain_to_drop = invalid_date["strain"].to_list()
    qc_1_idx_to_drop = invalid_date.index.to_list()
    print(f"There are {len(qc_1_idx_to_drop)} samples to drop due to invalid date")

    meta_df.drop(qc_1_idx_to_drop, inplace=True)

    # QC Filter #2 - Host
    print("Applying QC Filter #2 - Host...")
    if "host" in meta_df.columns:
        missing_valid = lambda v: (
            not pd.isna(v) if MISSING_HOST_IS_VALID else pd.isna(v)
        )
        invalid_host_names = [
            v
            for v in meta_df["host"].unique()
            if v not in VALID_HOSTS and missing_valid(v)
        ]
        invalid_host = meta_df[meta_df["host"].isin(invalid_host_names)]
        qc_2_strain_to_drop = invalid_host["strain"].to_list()
        qc_2_idx_to_drop = invalid_host.index.to_list()
        print(f"There are {len(qc_2_idx_to_drop)} invalid hosts")

        meta_df.drop(qc_2_idx_to_drop, inplace=True)
    else:
        print(
            "Warning: 'host' column not found in metadata. Host filter will not be applied."
        )
        qc_2_strain_to_drop = []

    # QC Filter #3 - Length
    print("Applying QC Filter #3 - Length...")
    if "length" in meta_df.columns:
        print(
            f"There are {sum(meta_df['length'].isnull())} samples with missing length information"
        )

        invalid_length = meta_df[
            ~meta_df["length"].between(MIN_GENOME_LENGTH, MAX_GENOME_LENGTH)
        ]
        print(f"There are {len(invalid_length)} samples with invalid length attribute")
        qc_3_strain_to_drop = invalid_length["strain"].to_list()
        qc_3_idx_to_drop = invalid_length.index.to_list()
        meta_df.drop(qc_3_idx_to_drop, inplace=True)
    else:
        print(
            "Warning: 'length' column not found in metadata. Length filter will not be applied."
        )
        qc_3_strain_to_drop = []

    # Output results
    print("Saving outputs...")
    print(f"Final dataset size: {len(meta_df)}")
    print(f"Date range: {meta_df['date'].min()} to {meta_df['date'].max()}")

    if "host" in meta_df.columns:
        print(f"Number of unique hosts: {meta_df['host'].nunique()}")
    else:
        print("Host information not available in metadata")

    if "length" in meta_df.columns:
        print(
            f"Genome length range: {meta_df['length'].min()} to {meta_df['length'].max()}"
        )
    else:
        print("Length information not available in metadata")

    meta_df.to_csv(FILTERED_META_OUT, sep="\t", index=False)
    with open(INVALID_STRAINS_OUT, "w") as f:
        f.write(
            "\n".join(qc_1_strain_to_drop + qc_2_strain_to_drop + qc_3_strain_to_drop)
            + "\n"
        )

    print("Filtering complete!")


if __name__ == "__main__":
    main()
