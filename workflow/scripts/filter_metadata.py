#!/usr/bin/env python3
"""Filter metadata based on QC criteria (dates, hosts, genome length)."""

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
        date_str: A string representation of a date

    Returns:
        A DateResolution enum value indicating the resolution of the date
    """
    if pd.isna(date_str) or date_str == "?":
        return DateResolution.MISSING
    elif "-" in str(date_str):
        parts = str(date_str).split("-")
        if len(parts) == 3:
            return DateResolution.DAY
        elif len(parts) == 2:
            return DateResolution.MONTH
    elif str(date_str).isdigit():
        return DateResolution.YEAR
    return DateResolution.MISSING


def robust_fromisoformat(date_string: str) -> pd.Timestamp:
    """
    Convert a string into a datetime object with robust handling of variations.

    Args:
        date_string: A date in string format.

    Returns:
        Corresponding datetime object or NaT for invalid inputs.
    """
    if pd.isna(date_string) or date_string == "?":
        return pd.NaT

    date_string = str(date_string)

    if "-" not in date_string and len(date_string) == 4:
        return datetime.fromisoformat(f"{date_string}-06-01")
    elif date_string.count("-") == 1:
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
        try:
            return datetime.fromisoformat(date_string)
        except ValueError:
            return pd.NaT


def main():
    """Main filtering function."""
    # Get parameters from snakemake
    META = snakemake.input["metadata"]
    SUBSET_FRAC = snakemake.params["subset_frac"]
    SKIP_QC_FILTERING = snakemake.params.get("skip_qc_filtering", False)

    qc_opts = snakemake.params["qc_filter_options"]
    DATETIME_CUTOFF = qc_opts["datetime_cutoff"]
    INVALID_DATE_RES = qc_opts["invalid_date_resolutions"]
    MIN_DATE = qc_opts["min_date"]
    VALID_HOSTS = qc_opts["valid_hosts"]
    MISSING_HOST_IS_VALID = qc_opts["missing_host_is_valid"]
    MIN_GENOME_LENGTH = qc_opts["min_genome_length"]
    MAX_GENOME_LENGTH = qc_opts["max_genome_length"]

    FILTERED_META_OUT = snakemake.output["filtered_metadata"]
    INVALID_STRAINS_OUT = snakemake.log["invalid_strains"]

    # Load metadata
    print("Loading metadata...")
    meta_df = pd.read_csv(META, compression="gzip", sep="\t", low_memory=False)
    print(f"Original dataset size: {len(meta_df)}")

    # Skip QC filtering if requested (to match SARS2-mut-fitness behavior)
    if SKIP_QC_FILTERING:
        print("Skipping ALL QC filtering (skip_qc_filtering=true)")

        # Save outputs with no filtering at all
        print(f"Final dataset size: {len(meta_df)}")
        if FILTERED_META_OUT.endswith(".gz"):
            meta_df.to_csv(FILTERED_META_OUT, sep="\t", index=False, compression="gzip")
        else:
            meta_df.to_csv(FILTERED_META_OUT, sep="\t", index=False)

        with open(INVALID_STRAINS_OUT, "w") as f:
            f.write("# QC filtering skipped (skip_qc_filtering=true)\n")
            f.write("# No filtering applied\n")

        print("Filtering complete (all QC skipped)!")
        return

    # Apply random subsetting
    if SUBSET_FRAC < 1.0:
        subset_size = int(len(meta_df) * SUBSET_FRAC)
        meta_df = meta_df.sample(n=subset_size, random_state=42)
        print(f"After random subsetting ({SUBSET_FRAC}): {len(meta_df)}")

    qc_1_strain_to_drop = []
    qc_2_strain_to_drop = []
    qc_3_strain_to_drop = []

    # QC Filter #1 - Date
    if "date" in meta_df.columns:
        print("Applying QC Filter #1 - Date...")
        meta_df["date_resolution"] = meta_df["date"].apply(determine_resolution)
        meta_df["date_resolution_str"] = meta_df["date_resolution"].apply(
            lambda x: x.name
        )
        meta_df.drop(columns=["date_resolution"], inplace=True)
        meta_df["date"] = meta_df["date"].apply(robust_fromisoformat)

        datetime_cutoff = datetime.fromisoformat(DATETIME_CUTOFF)
        valid_dates = ~pd.isna(meta_df["date"])
        meta_df = meta_df[valid_dates & (meta_df["date"] <= datetime_cutoff)]

        invalid_resolution = meta_df["date_resolution_str"].isin(INVALID_DATE_RES)
        invalid_time = meta_df["date"] < robust_fromisoformat(MIN_DATE)
        invalid_date = meta_df[invalid_resolution | invalid_time]
        qc_1_strain_to_drop = invalid_date["strain"].to_list()
        meta_df.drop(invalid_date.index, inplace=True)
        print(f"  Dropped {len(qc_1_strain_to_drop)} samples due to invalid date")
    else:
        print("Warning: 'date' column not found. Date filters will not be applied.")

    # QC Filter #2 - Host
    if "host" in meta_df.columns:
        print("Applying QC Filter #2 - Host...")

        def is_invalid_host(v):
            if pd.isna(v):
                return not MISSING_HOST_IS_VALID
            return v not in VALID_HOSTS

        invalid_host_mask = meta_df["host"].apply(is_invalid_host)
        invalid_host = meta_df[invalid_host_mask]
        qc_2_strain_to_drop = invalid_host["strain"].to_list()
        meta_df.drop(invalid_host.index, inplace=True)
        print(f"  Dropped {len(qc_2_strain_to_drop)} samples due to invalid host")
    else:
        print("Warning: 'host' column not found. Host filter will not be applied.")

    # QC Filter #3 - Length
    if "length" in meta_df.columns:
        print("Applying QC Filter #3 - Length...")
        invalid_length = meta_df[
            ~meta_df["length"].between(MIN_GENOME_LENGTH, MAX_GENOME_LENGTH)
        ]
        qc_3_strain_to_drop = invalid_length["strain"].to_list()
        meta_df.drop(invalid_length.index, inplace=True)
        print(f"  Dropped {len(qc_3_strain_to_drop)} samples due to invalid length")
    else:
        print("Warning: 'length' column not found. Length filter will not be applied.")

    # Save outputs
    print(f"Final dataset size: {len(meta_df)}")

    # Ensure output is gzipped
    if FILTERED_META_OUT.endswith(".gz"):
        meta_df.to_csv(FILTERED_META_OUT, sep="\t", index=False, compression="gzip")
    else:
        meta_df.to_csv(FILTERED_META_OUT, sep="\t", index=False)

    # Write invalid strains log with summary statistics
    all_invalid = qc_1_strain_to_drop + qc_2_strain_to_drop + qc_3_strain_to_drop
    with open(INVALID_STRAINS_OUT, "w") as f:
        f.write("# Invalid strains removed during QC filtering\n")
        f.write(f"# Date filter: {len(qc_1_strain_to_drop)}\n")
        f.write(f"# Host filter: {len(qc_2_strain_to_drop)}\n")
        f.write(f"# Length filter: {len(qc_3_strain_to_drop)}\n")
        f.write(f"# Total: {len(all_invalid)}\n")
        for strain in all_invalid:
            f.write(f"{strain}\n")

    print("Filtering complete!")
    print(f"  Date-filtered: {len(qc_1_strain_to_drop)}")
    print(f"  Host-filtered: {len(qc_2_strain_to_drop)}")
    print(f"  Length-filtered: {len(qc_3_strain_to_drop)}")


if __name__ == "__main__":
    main()
