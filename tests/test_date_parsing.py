#!/usr/bin/env python3
"""Pin the metadata date parser and the config/row strictness split (v11).

``filter_metadata.py`` classifies every leaf's date, and both of its failure
modes are silent. Neither is reachable from any shipped config today, which is
exactly why they need a unit test rather than a pipeline assertion.

**Why the parser is lenient at all.** Through v10 only the last of
``robust_fromisoformat``'s four branches carried a ``try/except``, so any value
the other three could not parse aborted the whole run. The 2026-08-06 metadata
does that: 231 rows per tree carry the Solr-style range ``[2022 TO 2024]``, which
has no ``-`` and length 14, so it reaches the ``len != 10`` branch and raises
there. ``determine_resolution`` already classifies such a value ``MISSING``,
which every shipped config lists in ``invalid_date_resolutions``, so the intended
disposition was always "drop it" -- only the parser disagreed. All four branches
are now guarded.

**Why widening that guard needed a second function.** ``robust_fromisoformat``
reads two unrelated kinds of input. Per-row metadata is external and messy, where
an unrecognized value is data to be dropped. A date from the *config* is a
maintainer's typo, and ``min_leaf_date`` is consumed only as
``meta_df["date"] < robust_fromisoformat(MIN_LEAF_DATE)`` -- every comparison
against ``NaT`` is ``False``, so returning ``NaT`` there made ``invalid_time``
all-False and dropped **zero** samples, at exit 0, with no warning: a QC filter
silently switched off by a one-character config error. ``max_leaf_date`` is parsed
with a bare ``datetime.fromisoformat`` and so kept failing loudly, so the two
date knobs also disagreed on strictness. ``config_date`` restores the raise
without changing which formats are accepted.

Both halves are asserted against the **real script**, exec'd from disk, because a
reimplementation of a four-branch string parser here would be the thing under
test rather than a check on it. The two directions are deliberately separate:
a test that only checked "the Solr range no longer raises" would pass on a parser
that returned ``NaT`` for *everything*, and a test that only checked round-trips
would pass on the v10 parser that aborted the run.
"""

import importlib.util as importlib_util
import sys
import traceback
from datetime import datetime
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "workflow" / "scripts" / "filter_metadata.py"

# (raw value, expected ISO date). One per branch of robust_fromisoformat, so a
# branch that stops parsing is caught by value and not just by "did not raise".
#   year-only          -> midpoint 06-01
#   year-month         -> midpoint day 15
#   unpadded           -> zero-filled components
#   exact length 10    -> passed through
PARSEABLE = [
    ("2021", "2021-06-01", "year-only, mapped to the year midpoint"),
    ("2021-03", "2021-03-15", "year-month, mapped to the month midpoint"),
    ("2021-3-7", "2021-03-07", "unpadded month and day"),
    ("2021-3", "2021-03-15", "unpadded year-month"),
    ("2021-03-07", "2021-03-07", "already exact ISO"),
    # No '-' and length 8, so this reaches the len != 10 branch, where the
    # zero-fill is a no-op and Python 3.11+ fromisoformat accepts ISO 8601 basic
    # format directly. Pinned because it parses by way of a branch written for a
    # different shape, so it would be easy to break while "tidying" that branch.
    ("20210307", "2021-03-07", "compact ISO basic format, no separators"),
]

# Values that must become NaT rather than raising. The Solr range is the one that
# killed the 2026-08-06 run; the rest are the neighbouring shapes that reach the
# same branches.
UNPARSEABLE = [
    ("[2022 TO 2024]", "the Solr-style range that aborted the 2026-08-06 run"),
    ("2021-99-99", "in-format but out-of-range components"),
    ("not-a-date", "free text with the right number of separators"),
    ("2021-", "trailing separator, empty component"),
]

# Missing-value sentinels, which were already NaT before the widening. Listed
# separately because they must NOT be counted as "unparseable" by the diagnostic:
# conflating an ordinary blank with a novel metadata format is what the report
# exists to distinguish.
MISSING = [None, float("nan"), "?", pd.NA]

# Every min_leaf_date / max_leaf_date in a shipped config must keep parsing.
SHIPPED_CONFIG_DATES = ["2019-12-01", "2026-04-22", "2026-08-06"]


def _load_script():
    """Exec the real filter_metadata.py.

    It touches the ``snakemake`` global only inside ``main()``, so importing it
    outside a pipeline run is safe.
    """
    spec = importlib_util.spec_from_file_location("filter_metadata", SCRIPT)
    module = importlib_util.module_from_spec(spec)
    sys.modules["filter_metadata"] = module
    spec.loader.exec_module(module)
    return module


def check_parseable_values_round_trip():
    """Every branch still parses, and to the same date as before the widening.

    The direction that catches a parser broken *open*: a bare
    ``return pd.NaT`` would satisfy "the Solr range does not raise" while
    silently dropping every leaf in the tree.
    """
    failures = []
    fm = _load_script()
    for raw, expected, why in PARSEABLE:
        got = fm.robust_fromisoformat(raw)
        if pd.isna(got):
            failures.append(f"{raw!r} ({why}) returned NaT; expected {expected}")
        elif got != datetime.fromisoformat(expected):
            failures.append(f"{raw!r} ({why}) parsed to {got}, expected {expected}")
    if not failures:
        print(f"  {len(PARSEABLE)} values parse, one per branch")
    return failures


def check_unparseable_values_are_nat_not_raises():
    """The widened guard: unrecognized values return NaT instead of aborting."""
    failures = []
    fm = _load_script()
    for raw, why in UNPARSEABLE:
        try:
            got = fm.robust_fromisoformat(raw)
        except Exception as exc:
            failures.append(
                f"{raw!r} ({why}) raised {type(exc).__name__}: {exc}. This aborts "
                f"the whole run in filter_metadata; it must return NaT and be "
                f"dropped by the valid_dates filter"
            )
            continue
        if not pd.isna(got):
            failures.append(f"{raw!r} ({why}) parsed to {got}; expected NaT")
    if not failures:
        print(f"  {len(UNPARSEABLE)} unparseable values return NaT without raising")
    return failures


def check_missing_sentinels_still_nat():
    """The pre-existing missing values keep their disposition."""
    failures = []
    fm = _load_script()
    for raw in MISSING:
        try:
            got = fm.robust_fromisoformat(raw)
        except Exception as exc:
            failures.append(f"{raw!r} raised {type(exc).__name__}: {exc}")
            continue
        if not pd.isna(got):
            failures.append(f"{raw!r} parsed to {got}; expected NaT")
    if not failures:
        print(f"  {len(MISSING)} missing-value sentinels still NaT")
    return failures


def check_unparseable_is_distinguishable_from_missing():
    """The reporting mask separates a novel format from an ordinary blank.

    ``main()`` reports unparseable values so the drop is attributable. That mask
    is ``raw.notna() & (raw != "?") & parsed.isna()``; this asserts it selects the
    Solr range and none of the missing sentinels, since a mask that caught both
    would bury the signal it exists to surface.
    """
    failures = []
    fm = _load_script()
    raw = pd.Series(["[2022 TO 2024]", "2021-03-07", "?", None, float("nan")])
    parsed = raw.apply(fm.robust_fromisoformat)
    mask = raw.notna() & (raw != "?") & parsed.isna()
    selected = raw[mask].tolist()
    if selected != ["[2022 TO 2024]"]:
        failures.append(
            f"the unparseable-reporting mask selected {selected}; expected only "
            f"['[2022 TO 2024]'] -- a missing value must not be reported as an "
            f"unfamiliar metadata format, and a real one must not be hidden"
        )
    if not failures:
        print("  reporting mask selects the novel format only, not blanks or '?'")
    return failures


def check_config_date_raises_on_unparseable():
    """A malformed config date must stop the run, not silently filter nothing.

    This is the direction that regressed when the guard was widened, and the one
    no shipped config can exercise.
    """
    failures = []
    fm = _load_script()
    if not hasattr(fm, "config_date"):
        return [
            "filter_metadata.py has no config_date; min_leaf_date is presumably "
            "parsed with robust_fromisoformat again, which returns NaT and "
            "silently disables the min-date filter (every comparison against NaT "
            "is False)"
        ]
    for raw, why in UNPARSEABLE:
        try:
            got = fm.config_date(raw, "min_leaf_date")
        except ValueError:
            continue
        except Exception as exc:
            failures.append(
                f"config_date({raw!r}) raised {type(exc).__name__}, expected "
                f"ValueError: {exc}"
            )
            continue
        failures.append(
            f"config_date({raw!r}) ({why}) returned {got!r} instead of raising. "
            f"If that is NaT, min_leaf_date compares False against every row and "
            f"the QC filter drops nothing, at exit 0"
        )
    if not failures:
        print(f"  config_date raises on all {len(UNPARSEABLE)} unparseable values")
    return failures


def check_config_date_accepts_shipped_values():
    """The strict wrapper must not narrow the accepted formats.

    Paired with the check above: a ``config_date`` that raised on everything
    would satisfy it while breaking every config in the repo.
    """
    failures = []
    fm = _load_script()
    if not hasattr(fm, "config_date"):
        return []  # already reported by the previous check
    for raw in SHIPPED_CONFIG_DATES + [p[0] for p in PARSEABLE]:
        try:
            got = fm.config_date(raw, "min_leaf_date")
        except Exception as exc:
            failures.append(
                f"config_date({raw!r}) raised {type(exc).__name__}: {exc}. This "
                f"format is accepted by robust_fromisoformat, so the strict "
                f"wrapper has narrowed what a config may say"
            )
            continue
        if pd.isna(got):
            failures.append(f"config_date({raw!r}) returned NaT instead of raising")
    if not failures:
        print("  config_date accepts every format robust_fromisoformat does")
    return failures


def check_min_leaf_date_uses_the_strict_parser():
    """``main()`` must route min_leaf_date through config_date, not the lenient one.

    Asserted against the source text because reaching the call needs a full
    ``snakemake`` object and a metadata TSV. Without this, both parsers can be
    perfectly correct while ``main()`` still calls the wrong one -- the shape of
    CLAUDE.md patterns 8 and 11, where a test named for a behaviour never runs the
    code that implements it.
    """
    failures = []
    source = SCRIPT.read_text()
    if 'config_date(MIN_LEAF_DATE, "min_leaf_date")' not in source:
        failures.append(
            "filter_metadata.py's main() does not parse MIN_LEAF_DATE with "
            "config_date(MIN_LEAF_DATE, 'min_leaf_date'). If it uses "
            "robust_fromisoformat, a malformed min_leaf_date yields NaT and the "
            "min-date QC filter silently drops nothing"
        )
    if not failures:
        print("  main() routes min_leaf_date through the strict parser")
    return failures


CHECKS = [
    ("all four branches still parse", check_parseable_values_round_trip),
    ("unparseable values return NaT", check_unparseable_values_are_nat_not_raises),
    ("missing sentinels unchanged", check_missing_sentinels_still_nat),
    ("reporting mask is specific", check_unparseable_is_distinguishable_from_missing),
    ("config dates fail loudly", check_config_date_raises_on_unparseable),
    ("config dates stay permissive", check_config_date_accepts_shipped_values),
    ("main() uses the strict parser", check_min_leaf_date_uses_the_strict_parser),
]


def run_all_checks():
    failures = []
    for label, check in CHECKS:
        print(f"Checking {label}...")
        try:
            failures.extend(check())
        except Exception:
            failures.append(f"{label} raised:\n{traceback.format_exc()}")
    return failures


def test_date_parsing():
    """Pytest entry point.

    This repo runs tests as standalone scripts (``pixi run python tests/<f>.py``)
    and does not depend on pytest. This wrapper exists so that anyone who does
    reach for ``pytest tests/`` gets a real result rather than a green run that
    collected nothing.
    """
    failures = run_all_checks()
    assert not failures, "\n".join(failures)


def main():
    failures = run_all_checks()
    if failures:
        print(f"\nFAILED ({len(failures)} problem(s)):")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    print(
        "\nPASSED: every branch parses, unrecognized row values become NaT "
        "instead of aborting the run, and a malformed config date still raises."
    )


if __name__ == "__main__":
    main()
