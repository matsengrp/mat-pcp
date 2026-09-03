#!/usr/bin/env python3
"""Filter metadata based on QC criteria (dates, hosts, genome length)."""

import os
from importlib import util as importlib_util

import pandas as pd
from datetime import datetime
from enum import Enum, auto

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
pcpio = workflow_loader.load_workflow_module("pcpio")


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

    Every branch is guarded, not just the exact-length-10 one. Through the
    2026-04-22 trees only that last branch caught ValueError, so any string the
    other three branches could not parse aborted the whole run -- which is what
    the 2026-08-06 metadata does: 231 rows of each tree carry the Solr-style
    range ``[2022 TO 2024]``, no ``-`` and length 14, so it reaches the
    ``len != 10`` branch and raises there. ``determine_resolution`` already
    classifies such a value MISSING, which every shipped config lists in
    ``invalid_date_resolutions``, so the intended disposition was always "drop
    it"; only the parser disagreed. Returning NaT lets the existing
    ``valid_dates`` filter do that, and ``main`` reports the count so the drop
    is visible rather than silent.
    """
    if pd.isna(date_string) or date_string == "?":
        return pd.NaT

    date_string = str(date_string)

    try:
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
            return datetime.fromisoformat(date_string)
    except ValueError:
        return pd.NaT


def config_date(date_string: str, key: str) -> pd.Timestamp:
    """Parse a *config-supplied* date, raising rather than returning NaT.

    ``robust_fromisoformat`` is lenient by design because it reads external
    per-row metadata, where an unrecognized value is data to be dropped. A date
    from the config is the opposite case: a typo there is a maintainer error, and
    a run must not proceed on it.

    Without this, widening ``robust_fromisoformat``'s guard to all four branches
    silently disabled a QC filter. ``min_leaf_date`` is used only as
    ``meta_df["date"] < robust_fromisoformat(MIN_LEAF_DATE)``, and every
    comparison against ``NaT`` is ``False`` -- so a malformed value made
    ``invalid_time`` all-False and dropped zero samples, at exit 0 with no
    warning. ``max_leaf_date`` is parsed with a bare ``datetime.fromisoformat``
    and so kept failing loudly; the two date knobs must not differ in strictness.
    The accepted formats are unchanged -- this only converts an unparseable value
    from a silent no-op back into an error.
    """
    parsed = robust_fromisoformat(date_string)
    if pd.isna(parsed):
        raise ValueError(
            f"qc_filter_options.{key} = {date_string!r} is not a parseable date. "
            f"Unlike a metadata row, a config date is not dropped when it fails "
            f"to parse: {key} would silently filter nothing. Fix the config."
        )
    return parsed


def main():
    """Main filtering function."""
    # Get parameters from snakemake
    META = snakemake.input["metadata"]
    SUBSET_FRAC = snakemake.params["subset_frac"]

    qc_opts = snakemake.params["qc_filter_options"]
    MAX_LEAF_DATE = qc_opts.get("max_leaf_date")  # None -> no upper date cutoff
    INVALID_DATE_RES = qc_opts["invalid_date_resolutions"]
    MIN_LEAF_DATE = qc_opts["min_leaf_date"]
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
        raw_dates = meta_df["date"]
        meta_df["date"] = raw_dates.apply(robust_fromisoformat)

        # Values that were present but could not be parsed at all. They are dropped
        # by `valid_dates` below like any other invalid date; reported here so the
        # drop is attributable, since an unparseable date means a metadata format
        # this parser has not seen rather than an ordinary missing value.
        unparseable = raw_dates.notna() & (raw_dates != "?") & meta_df["date"].isna()
        if unparseable.any():
            examples = raw_dates[unparseable].value_counts()
            print(
                f"Unparseable date values: {int(unparseable.sum())} rows across "
                f"{len(examples)} distinct value(s); dropped as invalid dates. "
                f"Most common: {examples.head(5).to_dict()}"
            )

        valid_dates = ~pd.isna(meta_df["date"])
        if MAX_LEAF_DATE is not None:
            max_leaf_date = datetime.fromisoformat(MAX_LEAF_DATE)
            meta_df = meta_df[valid_dates & (meta_df["date"] <= max_leaf_date)]
        else:
            meta_df = meta_df[valid_dates]

        invalid_resolution = meta_df["date_resolution_str"].isin(INVALID_DATE_RES)
        invalid_time = meta_df["date"] < config_date(MIN_LEAF_DATE, "min_leaf_date")
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

    # Ensure output is gzipped. The compression policy comes from pcpio rather
    # than being spelled here: pandas' default embeds mtime=now, so an unchanged
    # re-run produces a different md5, and `data-central` manifests pin md5/xxh3
    # per file. That is not hypothetical for this file -- v12 shipped it with
    # mtime=now baked in (#134's policy, applied to the one writer that missed
    # it).
    if FILTERED_META_OUT.endswith(".gz"):
        meta_df.to_csv(
            FILTERED_META_OUT,
            sep="\t",
            index=False,
            compression=pcpio.REPRODUCIBLE_GZIP,
        )
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
