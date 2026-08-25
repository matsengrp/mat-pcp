#!/usr/bin/env python3
"""Date-cutoff subtree selection for the truncated PCP sets (#149).

A *cutoff* names a date and asks for the PCPs that would have been extractable
from the tree using only samples collected **at or before** it -- a temporal
holdout, so DASM can be trained on what was knowable then.

Three decisions live here because each has a wrong answer that looks right.

**The cutoff selects from the clade's own sample list, never from the metadata's
clade column.** ``rule samples_by_clade`` matches a clade with the awk regex
``^clade( |$)`` against values like ``21L (Omicron)``, so re-deriving clade
membership in Python would be a second implementation of that match. It would
agree on ``19A`` and disagree on exactly the values the regex exists for --
and the disagreement would be invisible, because both lists are plausible clade
rosters and the resulting PCP counts merely differ. ``select_upto`` therefore
takes the *already-written* ``mat/clades/{clade}_samples.txt`` and filters it by
date alone. A cutoff at or after every sample's date then returns that file
unchanged, which is what makes the degenerate arm's byte-equality hold by
construction rather than by luck.

**Order is preserved.** ``matUtils extract -s`` is handed this list, and
returning a set (or a sorted copy of a list that was not sorted) would change
the file for no reason. ``data-central`` manifests pin ``md5``/``xxh3`` per file
(pattern 9), so a gratuitous reordering is a spurious diff on every rebuild --
and on the degenerate arm it would break byte-equality with the full-tree
outputs, the one check that says the cutoff machinery left the ordinary path
alone.

**The comparison is against the *coerced* date, not the raw metadata string.**
``filter_metadata.py`` rewrites the ``date`` column through
``robust_fromisoformat``, which maps a month-resolution value to the 15th and a
year-resolution value to June 1st. Every shipped config lists ``MONTH``/``YEAR``
in ``invalid_date_resolutions``, so in practice every surviving row is
day-resolution and the distinction is moot -- but a config that admitted them
would compare the coerced date, and that is the behaviour to reason about rather
than a coincidence to rely on.
"""

import datetime
import os

import pandas as pd

# The ``{cutoff}`` wildcard constraint, and the directory prefix under
# ``pcps/{gene}/``. Strict ISO: a cutoff is a *config* date, so an unparseable
# value is a maintainer error that must stop the run rather than something to
# coerce (the same argument ``filter_metadata.config_date`` makes for
# ``min_leaf_date``).
CUTOFF_PATTERN = r"[0-9]{4}-[0-9]{2}-[0-9]{2}"
CUTOFF_DIR_PREFIX = "upto_"


def cutoff_dir(cutoff):
    """The per-cutoff subdirectory name, e.g. ``upto_2023-01-01``."""
    return f"{CUTOFF_DIR_PREFIX}{cutoff}"


def parse_cutoff(value, where="pcp_date_cutoffs"):
    """Parse one configured cutoff, raising on anything but strict ISO.

    ``datetime.date.fromisoformat`` accepts ``2023-1-1`` on Python 3.11+, which
    would not match ``CUTOFF_PATTERN`` and so would fail later as a baffling
    ``MissingInputException`` on a path nobody asked for. Checking the length
    and the separators here means the error names the config key.
    """
    text = str(value)
    if len(text) != 10 or text[4] != "-" or text[7] != "-":
        raise ValueError(
            f"{where} entry {value!r} is not strict ISO YYYY-MM-DD. A cutoff "
            f"becomes a directory name and the '{{cutoff}}' wildcard, so a "
            f"loose form would fail later as a missing-input error on a path "
            f"nothing requested."
        )
    try:
        return datetime.date.fromisoformat(text)
    except ValueError as exc:
        raise ValueError(f"{where} entry {value!r} is not a real date: {exc}") from None


def validate_cutoffs(cutoffs, min_leaf_date, max_leaf_date, label):
    """Return the parsed, de-duplicated cutoffs for one tree, in config order.

    Each must postdate ``min_leaf_date`` (a cutoff at or before it selects
    nothing, and an empty sample list makes ``matUtils extract`` fail with an
    error naming neither the cutoff nor the config) and must not exceed
    ``max_leaf_date`` (a cutoff past the tree's own horizon reproduces the full
    tree while its directory name claims a truncation -- the most expensive kind
    of wrong, because every output looks right).

    A cutoff *equal* to ``max_leaf_date`` is allowed on purpose: it is the
    degenerate case whose outputs must match the full tree's byte-for-byte, and
    ``config_test_data.yaml`` carries it for exactly that reason.
    """
    parsed, seen = [], set()
    lo = parse_cutoff(min_leaf_date, "qc_filter_options.min_leaf_date")
    hi = None if max_leaf_date is None else parse_cutoff(
        max_leaf_date, "qc_filter_options.max_leaf_date"
    )
    for value in cutoffs:
        date = parse_cutoff(value, f"trees.{label}.pcp_date_cutoffs")
        if date <= lo:
            raise ValueError(
                f"Tree '{label}': pcp_date_cutoffs entry {value!r} is at or "
                f"before qc_filter_options.min_leaf_date ({min_leaf_date}), so "
                f"it selects no samples at all."
            )
        if hi is not None and date > hi:
            raise ValueError(
                f"Tree '{label}': pcp_date_cutoffs entry {value!r} is after the "
                f"tree's max_leaf_date ({max_leaf_date}), so it would select "
                f"every sample and reproduce the full tree under a directory "
                f"name claiming a cutoff. Drop it, or use {max_leaf_date!r} if "
                f"the degenerate arm is what you want."
            )
        if str(value) in seen:
            raise ValueError(
                f"Tree '{label}': pcp_date_cutoffs lists {value!r} twice."
            )
        seen.add(str(value))
        parsed.append(str(value))
    return parsed


def load_strain_dates(metadata_path):
    """``strain`` -> collection date, from a filtered-metadata TSV.

    Reads only the two columns it needs: the production metadata is ~8.4 M rows
    and this is called once per (cutoff, clade).
    """
    df = pd.read_csv(
        metadata_path,
        sep="\t",
        usecols=["strain", "date"],
        dtype={"strain": "string"},
        low_memory=False,
    )
    dates = pd.to_datetime(df["date"], errors="coerce")
    missing = int(dates.isna().sum())
    if missing:
        raise ValueError(
            f"{metadata_path}: {missing} of {len(df)} rows have an unparseable "
            f"`date` after filtering. filter_metadata.py drops unparseable "
            f"dates, so this means the wrong metadata file was wired in -- the "
            f"raw input rather than mat/filtered_metadata.tsv.gz."
        )
    return pd.Series(dates.values, index=df["strain"].values)


def read_sample_list(path):
    """The strain names in a ``{clade}_samples.txt``, in file order."""
    with open(path) as handle:
        return [line.strip() for line in handle if line.strip()]


def select_upto(strains, strain_dates, cutoff):
    """The subsequence of ``strains`` collected at or before ``cutoff``.

    Inclusive of the cutoff date itself -- "at or before", so a sample collected
    on 2023-01-01 is in the ``upto_2023-01-01`` set. An off-by-one here would be
    invisible in every summary statistic.

    A strain absent from ``strain_dates`` raises. Both files derive from the same
    ``mat/filtered_metadata.tsv.gz``, so a miss means they do not: most likely
    the raw metadata was wired in, in which case silently dropping the strain
    would shrink every cutoff set by an unreported amount.
    """
    bound = pd.Timestamp(parse_cutoff(cutoff, "cutoff"))
    unknown = [s for s in strains if s not in strain_dates.index]
    if unknown:
        raise ValueError(
            f"{len(unknown)} strain(s) in the clade sample list are absent from "
            f"the metadata date index, e.g. {unknown[:3]}. The sample list and "
            f"the metadata must come from the same filtered metadata."
        )
    return [s for s in strains if strain_dates[s] <= bound]
