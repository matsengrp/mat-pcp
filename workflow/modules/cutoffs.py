#!/usr/bin/env python3
"""Date-cutoff subtree selection for the truncated PCP sets (#149).

A *cutoff* names a date and asks for the PCPs that would have been extractable
from the tree using only samples collected **at or before** it -- a temporal
holdout, so DASM can be trained on what was knowable then.

Four decisions live here because each has a wrong answer that looks right.

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

**A sample the tree dates well past the cutoff does not belong in that
cutoff's set, whatever its metadata says.** Some metadata dates are wrong by
years, so a cutoff answered from the ``date`` column alone admits post-cutoff
genomes. ``pcp_max_predicted_days_past_cutoff`` drops them; see the comment
block above ``select_for_cutoff`` for the measurements behind the value and for
why the filter is confined to the cutoff path.

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
from typing import NamedTuple

import pandas as pd

# The ``{cutoff}`` wildcard constraint, and the directory prefix under
# ``pcps/{gene}/``. Strict ISO: a cutoff is a *config* date, so an unparseable
# value is a maintainer error that must stop the run rather than something to
# coerce (the same argument ``filter_metadata.config_date`` makes for
# ``min_leaf_date``).
CUTOFF_PATTERN = r"[0-9]{4}-[0-9]{2}-[0-9]{2}"
CUTOFF_DIR_PREFIX = "upto_"

# A clade with fewer than two samples at a cutoff is not a small clade -- it is a
# clade that cannot be made into a tree, so the configured ``min_clade_size`` is
# *floored* here rather than trusted directly. Two separate failures sit below
# this line, and only the second one fixes the threshold at 2 rather than 1.
#
# ``0 >= 0``. Every production config ships ``min_clade_size: 0``, which is
# harmless on the full-tree path -- a clade only reaches the roster because
# metadata rows carry it, so it always has samples, and ``adequate_sample_counts``
# is vacuously true there. At a cutoff it is not harmless: on the v12 2026-08-06
# trees it marked all 26 zero-sample (cutoff, clade) pairs ``adequate``, the
# checkpoint committed the fan-out, and ``samples_upto_cutoff``'s empty-list guard
# then correctly raised on each -- taking every downstream aggregate with it. So
# the whole skip-and-record mechanism this module exists to implement was inert on
# exactly the configs it was written for.
#
# ``bte segfaults on a single-leaf MAT``. This is the reason the floor is 2. A
# one-sample cutoff subtree writes a well-formed ``.pb`` that ``matUtils``
# produces without complaint, and ``bte.MATree`` then dies with SIGSEGV loading
# it -- measured on 25C at 2023-01-01 (1 sample), six times, once per gene, while
# the 5-leaf 23E subtree loads in 4 ms. Pinned by construction rather than
# inferred: subtrees cut to 2 and 3 leaves from that same 23E tree both load
# cleanly, so the boundary is exactly one leaf. A signal-11 in a C++ extension
# yields no traceback and no partial output, so nothing above it can attribute
# the failure to clade size.
#
# The floor is structural, not a threshold: ``min_clade_size`` says how big a
# clade the user wants, and this says what can be built at all. Note the same
# latent crash exists on the full-tree path -- ``get_clades`` applies
# ``min_clade_size`` unfloored -- but no shipped tree has a one-sample clade
# (the smallest on the released v11 public tree is 25I at 20), so it is recorded
# here rather than changed there, where it would move released semantics.
#
# Only the sub-2 cases are floored. A clade with a *handful* of samples is kept,
# even though it may contribute zero PCP rows for some genes (23E has 5 samples
# at 2023-01-01 on the public tree and 0 rows for S), because that is exactly
# what the full-tree path already ships under the same ``min_clade_size: 0``:
# 21E has 47 samples and 6 S rows in the released v11 tree. Dropping them here
# would make the cutoff roster disagree with the full roster for reasons the
# config does not state. A consumer that cares reads ``n_samples_upto`` from
# ``pcps/cutoff_summary.csv``.
MIN_SAMPLES_AT_CUTOFF = 2


def adequate_at_cutoff(n_kept, min_clade_size):
    """Whether a clade with ``n_kept`` samples at a cutoff gets a PCP set.

    The single definition of the adequacy verdict, so the checkpoint that writes
    ``adequate_sample_counts`` and anything that re-derives it cannot drift. See
    ``MIN_SAMPLES_AT_CUTOFF`` for why the configured floor is raised to 2.
    """
    return int(n_kept) >= effective_min_clade_size(min_clade_size)


def effective_min_clade_size(min_clade_size):
    """The floored ``min_clade_size`` actually applied at a cutoff.

    Exposed separately so log lines and error messages quote the threshold that
    was *enforced*, not the one that was configured -- a message reading
    ``>= 0`` next to a skipped one-sample clade is a contradiction the reader
    has to resolve from the source.
    """
    return max(int(min_clade_size), MIN_SAMPLES_AT_CUTOFF)


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


# A cutoff asks which samples had been *collected* by a date, and answers with
# the metadata's ``date`` column. Some of those dates are wrong by years, so the
# question is answered with post-cutoff genomes unless they are removed.
#
# On the v12 2026-08-06 public tree, 16 leaves in the 2023-01-01 cutoff set
# belong to Nextstrain clades 24A (JN.1) and 24E (KP.3.1.1) -- lineages that did
# not exist in 2023 -- with reported dates from 2020-01-10 to 2022-12-09 and
# chronumental predicted dates 337 to 570 days *past* the cutoff. Their pango
# calls (JN.1.4, KP.3.1.1, MC.24) agree with the tree and not with the date.
# They arrive in batches sharing one placeholder value: 20-odd
# ``Spain/AR-IACS_LaSARD-*`` records all report exactly 2021-07-16 or
# 2022-10-18, nine ``USA/ND-NDDH-260xx`` all report 2021-09-02. That is a
# submission date written into the collection field, not independent mistakes.
#
# ``pcp_max_predicted_days_past_cutoff`` therefore drops a strain from a cutoff
# set when its *predicted* date lands more than that many days after the cutoff.
# Four things about the shape are load-bearing.
#
# *The test is against the cutoff, not against the reported date.* Comparing the
# two dates to each other -- dropping when they disagree by more than N -- is the
# obvious alternative and it is a different, broader filter: it also removes
# records whose dates disagree while *both* sit comfortably before the cutoff,
# which are untrustworthy but not leaks. What a temporal holdout can actually be
# wrong about is admitting a genome from after the cutoff, so that is what is
# tested. It also means the predicate depends on the cutoff: the same strain can
# be kept at a later cutoff and dropped at an earlier one, which is correct.
#
# *The tolerance exists because the predicted date is fitted, not observed.* A
# sample genuinely collected days before the cutoff can predict past it on
# ordinary fit noise: over every leaf of both v12 trees the median
# |predicted - reported| is 10 d, p95 is 35/37 d and p99 is 56/62 d under the
# ``--variance_dates 30`` prior (pattern 13). So the bound is a config value
# rather than zero and rather than a constant here -- it is a property of the
# date prior, which pattern 13 has already moved once, and a run that retunes
# the prior must revisit it.
#
# *This filter is applied at a cutoff only, never on the full-tree path.* Not a
# cost compromise: the full-tree outputs do not read the reported date at all --
# clade rosters come from the clade column and the windowed path keys on
# predicted dates -- so filtering there would move every output while fixing
# nothing. The reported date is load-bearing in exactly two places, the
# ``min_leaf_date``/``max_leaf_date`` metadata filter and this cutoff selection.
#
# *A strain chronumental did not date is kept, and counted.* 154 leaves of the
# public tree's 2023-01-01 cutoff set have a metadata row and a place in the
# tree but no ``node_dates.tsv`` entry. Dropping them would remove data on the
# strength of a missing prediction rather than a contradicting one; keeping them
# silently would hide a chronumental output that had gone half-missing. So they
# are kept and reported as ``n_undated_kept`` in ``pcps/cutoff_summary.csv``.
DISCORDANCE_CONFIG_KEY = "pcp_max_predicted_days_past_cutoff"


class CutoffSelection(NamedTuple):
    """The outcome of selecting one clade's samples for one cutoff.

    Returned as a record rather than a bare list so that the checkpoint's
    counting and the sample-list writer cannot drift: ``n_samples_upto`` in
    ``pcps/cutoff_summary.csv`` is ``len(kept)`` from this same call, and
    ``validate_cutoff_pcps`` reconciles that number against the written file's
    line count.
    """

    kept: list
    n_after_cutoff: int
    n_predicted_late: int
    n_undated: int


def load_predicted_dates(node_dates_path):
    """``strain`` -> chronumental predicted date, from ``mat/node_dates.tsv``.

    Internal nodes are in this file too and are simply never looked up; the
    caller indexes it by the strain names in a clade sample list. Values carry a
    time component (``2023-12-10 16:57:02.752815``) which ``to_datetime`` keeps
    and the day-granular comparison below tolerates.
    """
    df = pd.read_csv(
        node_dates_path,
        sep="\t",
        usecols=["strain", "predicted_date"],
        dtype={"strain": "string"},
        low_memory=False,
    )
    dates = pd.to_datetime(df["predicted_date"], errors="coerce")
    return pd.Series(dates.values, index=df["strain"].values)


def select_for_cutoff(
    strains,
    strain_dates,
    cutoff,
    predicted_dates=None,
    max_days_past_cutoff=None,
):
    """One clade's samples for one cutoff: date-truncated, then tree-checked.

    The single entry point for both callers (the ``get_clades_by_cutoff``
    checkpoint and ``samples_upto_cutoff``), because they must agree exactly --
    a filter applied in one and not the other yields a roster whose
    ``n_samples_upto`` disagrees with the file it describes, and only the
    validator's line-count check would notice.

    ``predicted_dates=None`` or ``max_days_past_cutoff=None`` disables the
    predicted-date filter, which is what makes a config without
    ``pcp_max_predicted_days_past_cutoff`` behave exactly as before.
    """
    upto = select_upto(strains, strain_dates, cutoff)
    n_after_cutoff = len(strains) - len(upto)
    if predicted_dates is None or max_days_past_cutoff is None:
        return CutoffSelection(upto, n_after_cutoff, 0, 0)

    grace = int(max_days_past_cutoff)
    if grace < 0:
        raise ValueError(
            f"pcp_max_predicted_days_past_cutoff must be >= 0, got "
            f"{max_days_past_cutoff!r}. Use null to disable the filter; 0 would "
            f"demand that no predicted date exceed the cutoff at all, which "
            f"ordinary fit noise violates for samples collected just before it."
        )
    # Day granularity on both sides. ``node_dates.tsv`` carries a time
    # (``2023-12-10 16:57:02.752815``) while the reported dates and the cutoff
    # are dates, so comparing raw timestamps would make membership depend on the
    # time of day chronumental happened to emit: a prediction of exactly
    # cutoff+X at 12:00 would be excluded and the same prediction at 00:00 kept,
    # which is a boundary nobody could reason about from the config.
    bound = pd.Timestamp(parse_cutoff(cutoff, "cutoff")) + pd.Timedelta(days=grace)
    kept, n_predicted_late, n_undated = [], 0, 0
    for strain in upto:
        predicted = predicted_dates.get(strain)
        if predicted is None or pd.isna(predicted):
            # Kept deliberately -- see ``n_undated_kept`` in the block above.
            n_undated += 1
            kept.append(strain)
            continue
        if pd.Timestamp(predicted).normalize() > bound:
            n_predicted_late += 1
        else:
            kept.append(strain)
    return CutoffSelection(kept, n_after_cutoff, n_predicted_late, n_undated)
