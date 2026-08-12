#!/usr/bin/env python3
"""Pin the pre-Omicron / Omicron era mapping and its year-prefix audit (#137).

``modules/era.py`` restores the ``pre_omicron_or_omicron`` covariate that mat-pcp
lost when it ported ``SARS2-mut-fitness-v2``'s GLM. The mapping is 58
hand-maintained entries, and its failure mode is a clade in the *wrong list*,
which produces a plausible rate table rather than an error. Nothing downstream can
see it: both estimators would be consistently wrong, exactly the #129/#130/#133
shape. So the mapping is audited against the Nextstrain naming convention, and
this file pins both halves.

What each check exists for:

1. **The mapping covers the tree the release actually runs.** Read from a frozen
   copy of ``gisaid-and-public_2026-04-22``'s ``mat/clade_counts.csv``, so adding a
   clade to the tree without mapping it fails here rather than at fit time. Also
   pins the 22/30 split, since a partition can be total and still be wrong.

2. **An unmapped clade raises.** Upstream's ``else 'omicron'`` fallback would
   label all ten pre-Omicron clades it never analyzed as Omicron. There is no
   default here, in either direction.

3. **The audit passes on the shipped mapping**, and

4. **the audit FAILS on each single-entry mutation** -- which is the check that
   stops (3) from being vacuous. Six mutations are tested, including all three
   2021-named Omicron clades.

5. **The two-sidedness is load-bearing, demonstrated rather than asserted.** This
   check builds the *one-sided* predicate someone would naturally write
   (``year <= 21 -> pre_omicron``), confirms it accepts the shipped mapping, and
   then confirms it **fails to catch** ``21K`` moved into ``pre_omicron`` -- the
   single most likely misassignment in the table, since 21K is BA.1 and sorts
   among the Delta-era clades. The real predicate catches it. Without this check,
   a future simplification of ``era_from_label`` that dropped the ``OMICRON_21``
   exclusion would leave checks 3 and 4 green for five of six mutations.

6. **The label format is asserted, not assumed.** ``checkpoint get_clades``
   normalizes ``"21K (BA.1)"`` to ``"21K"`` by accident of awk field splitting
   (``split($0, a, " ")`` then ``a[2]``, ``Snakefile:768-770``), which passes
   single-token non-clades through intact -- ``config_test_data.yaml``'s clade
   column carries ``recombinant`` and an empty value. Those must raise, not bucket.

7. **The two-digit year assumption is pinned.** Whether ``era_from_label`` parses
   the year as an integer or compares it as a string is *not observable* while the
   regex requires exactly two digits, because all two-digit strings order
   identically either way. What is observable, and is what makes the int parse
   safe, is that three- and one-digit years are rejected; that plus a sweep over
   every prefix ``00``-``99`` (pinning a single threshold between 21 and 22) is
   the honest coverage. A regex widened to ``\\d{2,}`` alongside a string
   comparison would fail the sweep.

8. **``partition_clades`` keys only the eras present.** This is what lets a
   single-era tree skip the absent era's rule at DAG-construction time (#138
   decision 9) rather than emitting a placeholder table. Exercised on the real
   test-config clade set (19A/19B/20A, all pre-Omicron).

Usage:
    python tests/test_era_mapping.py
"""

import csv
import sys
import traceback
from importlib import util as importlib_util
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
ERA_MODULE = REPO_ROOT / "workflow" / "modules" / "era.py"
PRODUCTION_CLADES = (
    REPO_ROOT / "tests" / "data" / "production_clade_counts_2026-04-22.csv"
)
# The 2026-08-06 trees add 25D-25I. Frozen separately rather than replacing the
# 2026-04-22 fixture, so the "CLADE_ERA maps nothing extra" direction of the check
# stays two-sided: a clade on neither shipped tree is still a drift failure.
# Frozen from public_2026-08-06's mat/clade_counts.csv. Both 2026-08-06 trees
# carry identical clade *labels* (verified; only the counts differ), so one
# fixture covers them -- but name the source tree, since a future regeneration
# from the other tree would still pass and the provenance would be lost.
PRODUCTION_CLADES_NEW = (
    REPO_ROOT / "tests" / "data" / "production_clade_counts_2026-08-06.csv"
)

# The 2026-04-22 tree's split, stated here so a mapping edit that keeps the
# partition total but moves a clade across the boundary still fails.
EXPECTED_N_PRE_OMICRON = 22
EXPECTED_N_OMICRON = 30
# The 2026-08-06 trees add six Omicron clades and no pre-Omicron ones.
EXPECTED_N_OMICRON_NEW = 36

# Single-entry mutations the audit must reject. All three 2021-named Omicron
# clades are included because the year prefix alone cannot separate them from
# Delta-era clades, so they are where a one-sided predicate goes blind.
MUTATIONS = [
    ("21K", "pre_omicron", "BA.1 relabelled Delta-era"),
    ("21L", "pre_omicron", "BA.2 relabelled Delta-era"),
    ("21M", "pre_omicron", "basal B.1.1.529 relabelled Delta-era"),
    ("21J", "omicron", "Delta relabelled Omicron"),
    ("20A", "omicron", "a 2020 clade relabelled Omicron"),
    ("22A", "pre_omicron", "BA.4 relabelled pre-Omicron"),
]

# The clade set config_test_data.yaml actually produces (max_clades: 3).
TEST_CONFIG_CLADES = ["19A", "19B", "20A"]


def _load_era():
    spec = importlib_util.spec_from_file_location("era", ERA_MODULE)
    module = importlib_util.module_from_spec(spec)
    sys.modules["era"] = module
    spec.loader.exec_module(module)
    return module


def _production_clades(path=PRODUCTION_CLADES):
    with open(path, newline="") as handle:
        return [row["Nextstrain_clade_usher"] for row in csv.DictReader(handle)]


def check_mapping_covers_production_tree():
    """Every clade on every shipped tree must be mapped, with the right split."""
    failures = []
    era = _load_era()
    clades = _production_clades()
    clades_new = _production_clades(PRODUCTION_CLADES_NEW)

    for label, got, expected in (
        ("2026-04-22", clades, EXPECTED_N_PRE_OMICRON + EXPECTED_N_OMICRON),
        ("2026-08-06", clades_new, EXPECTED_N_PRE_OMICRON + EXPECTED_N_OMICRON_NEW),
    ):
        if len(got) != expected:
            failures.append(
                f"the frozen {label} clade_counts.csv has {len(got)} clades but "
                f"this test expects {expected}; if the tree changed, update both "
                f"the fixture and the expected split"
            )

    for label, got in (("2026-04-22", clades), ("2026-08-06", clades_new)):
        missing = sorted(set(got) - set(era.CLADE_ERA))
        if missing:
            failures.append(f"clades on {label} with no era: {missing}")

    # Two-sided: a CLADE_ERA entry on neither shipped tree is drift, even though a
    # superset is harmless at runtime. Checked against the union so that adding a
    # newer tree does not require deleting the older tree's fixture.
    extra = sorted(set(era.CLADE_ERA) - set(clades) - set(clades_new))
    if extra:
        failures.append(
            f"CLADE_ERA maps clades absent from every shipped tree: {extra} "
            f"(harmless, but the fixture and mapping have drifted)"
        )

    for label, got, exp_pre, exp_omi in (
        ("2026-04-22", clades, EXPECTED_N_PRE_OMICRON, EXPECTED_N_OMICRON),
        ("2026-08-06", clades_new, EXPECTED_N_PRE_OMICRON, EXPECTED_N_OMICRON_NEW),
    ):
        n_pre = sum(1 for c in got if era.CLADE_ERA.get(c) == era.PRE_OMICRON)
        n_omi = sum(1 for c in got if era.CLADE_ERA.get(c) == era.OMICRON)
        if (n_pre, n_omi) != (exp_pre, exp_omi):
            failures.append(
                f"{label} era split is {n_pre} pre-Omicron / {n_omi} Omicron, "
                f"expected {exp_pre}/{exp_omi} -- a clade has moved across the "
                f"boundary"
            )

    bad = sorted(c for c, e in era.CLADE_ERA.items() if e not in era.ERAS)
    if bad:
        failures.append(f"clades mapped to an era outside ERAS: {bad}")

    if not failures:
        # Per tree, not once: reusing the loop's last n_pre/n_omi alongside the
        # first fixture's length printed "all 52 clades mapped (22 / 36)", which
        # does not add up and invites a maintainer to "fix" a real count.
        print(
            f"  all clades mapped on both shipped trees: "
            f"2026-04-22 {len(clades)} "
            f"({EXPECTED_N_PRE_OMICRON} pre-Omicron / {EXPECTED_N_OMICRON} Omicron), "
            f"2026-08-06 {len(clades_new)} "
            f"({EXPECTED_N_PRE_OMICRON} pre-Omicron / {EXPECTED_N_OMICRON_NEW} Omicron)"
        )
    return failures


def check_unmapped_clade_raises():
    """era_for_clade must raise rather than defaulting to either era."""
    failures = []
    era = _load_era()
    for clade in ("21N", "26A", "recombinant", "", "21K (BA.1)", None):
        try:
            result = era.era_for_clade(clade)
        except KeyError:
            continue
        except TypeError:
            continue  # unhashable/None -- still not a silent bucket
        failures.append(
            f"era_for_clade({clade!r}) returned {result!r} instead of raising; "
            f"an unmapped clade fitted into the wrong era produces a plausible "
            f"rate table, not an error"
        )
    if not failures:
        print("  era_for_clade raises on every unmapped label (no default bucket)")
    return failures


def check_audit_passes_on_shipped_mapping():
    failures = []
    era = _load_era()
    try:
        n = era.audit_clade_era_mapping()
    except ValueError as exc:
        failures.append(f"the shipped CLADE_ERA fails its own year-prefix audit:\n{exc}")
    else:
        if n != len(era.CLADE_ERA):
            failures.append(
                f"audit reported {n} entries but CLADE_ERA has "
                f"{len(era.CLADE_ERA)}; it is not checking every entry"
            )
        else:
            print(f"  year-prefix audit passes over all {n} entries")
    return failures


def check_audit_rejects_each_mutation():
    """The audit must fail on every single-entry mutation.

    This is what stops check 3 from being vacuous: an audit that returned
    unconditionally would pass there and fail here.
    """
    failures = []
    era = _load_era()
    for clade, wrong_era, description in MUTATIONS:
        mutated = dict(era.CLADE_ERA)
        if mutated[clade] == wrong_era:
            failures.append(
                f"{clade} is already mapped {wrong_era!r}, so this mutation is a "
                f"no-op and tests nothing"
            )
            continue
        mutated[clade] = wrong_era
        try:
            era.audit_clade_era_mapping(mutated)
        except ValueError as exc:
            if clade not in str(exc):
                failures.append(
                    f"audit rejected the {description} mutation but did not name "
                    f"{clade} in the message: {exc}"
                )
            continue
        failures.append(
            f"audit ACCEPTED {clade} -> {wrong_era!r} ({description}); the "
            f"year-prefix check is not two-sided"
        )

    # A non-era value must also be rejected, not silently compared.
    mutated = dict(era.CLADE_ERA)
    mutated["20A"] = "Omicron"  # right idea, wrong literal
    try:
        era.audit_clade_era_mapping(mutated)
    except ValueError:
        pass
    else:
        failures.append("audit accepted the era value 'Omicron', which is not in ERAS")

    # Every mutation above is a SINGLE wrong entry, so none of them can tell
    # whether the audit reports all disagreements or only the first. That gap is
    # real: truncating the message to ``disagreements[:1]`` leaves this whole file
    # green. It matters because the audit's purpose is to be read by whoever
    # mis-edited the table -- surfacing one of three wrong clades sends them off
    # to fix one and re-run into the next, and the natural reading of a
    # one-line failure is that there is one problem.
    multi = dict(era.CLADE_ERA)
    wrong = {"21K": era.PRE_OMICRON, "22A": era.PRE_OMICRON, "20A": era.OMICRON}
    for clade, wrong_era in wrong.items():
        if multi[clade] == wrong_era:
            failures.append(f"{clade} is already {wrong_era!r}; multi-entry case is a no-op")
        multi[clade] = wrong_era
    try:
        era.audit_clade_era_mapping(multi)
    except ValueError as exc:
        missing = sorted(c for c in wrong if c not in str(exc))
        if missing:
            failures.append(
                f"audit rejected a 3-entry mutation but named only "
                f"{sorted(set(wrong) - set(missing))}, omitting {missing}. It must "
                f"list every disagreement, not stop at the first."
            )
    else:
        failures.append(f"audit ACCEPTED three simultaneous wrong entries: {wrong}")

    if not failures:
        print(
            f"  audit rejects all {len(MUTATIONS)} single mutations, a bad era "
            f"value, and names all 3 clades in a multi-entry mutation"
        )
    return failures


def check_two_sidedness_is_load_bearing():
    """Show that a one-sided predicate would miss the 21K mutation.

    Documents *why* the OMICRON_21 exclusion exists by exhibiting the weaker
    predicate and the mutation it cannot see. If a future edit drops the
    exclusion, ``era_from_label`` becomes this predicate and the final assertion
    below fails.
    """
    failures = []
    era = _load_era()

    def one_sided_accepts(mapping):
        """The weaker check: an allowlist per era, with no exclusion.

        This is the natural first draft -- "pre-Omicron means a 19/20/21 prefix,
        Omicron means 22+ or one of 21K/21L/21M" -- and it is a *validity check*
        rather than a total function, so unlike a total predicate it accepts the
        shipped mapping exactly. What it cannot do is reject a 2021 clade on the
        pre-Omicron side, because the ``21`` prefix is allowed there.
        """
        for clade, mapped in mapping.items():
            year = int(clade[:2])
            if mapped == era.PRE_OMICRON:
                ok = year in (19, 20, 21)
            else:
                ok = year >= 22 or clade in era.OMICRON_21
            if not ok:
                return False
        return True

    # Premise: the weaker check accepts the shipped mapping, which is why it
    # looks adequate.
    if not one_sided_accepts(era.CLADE_ERA):
        failures.append(
            "the one-sided check rejects the shipped mapping, so this check's "
            "premise no longer holds"
        )

    # The mutation the one-sided check cannot see...
    mutated = dict(era.CLADE_ERA)
    mutated["21K"] = era.PRE_OMICRON
    if not one_sided_accepts(mutated):
        failures.append(
            "the one-sided check caught 21K -> pre_omicron, so this check no "
            "longer demonstrates anything"
        )

    # ...while it does still catch a mutation the year prefix alone can see, so
    # the demonstration is about the exclusion specifically and not about the
    # weaker check being useless.
    across_years = dict(era.CLADE_ERA)
    across_years["22A"] = era.PRE_OMICRON
    if one_sided_accepts(across_years):
        failures.append(
            "the one-sided check also missed 22A -> pre_omicron; it is weaker "
            "than this check assumes and the comparison is not isolating the "
            "OMICRON_21 exclusion"
        )

    # ...must be caught by the real one.
    try:
        era.audit_clade_era_mapping(mutated)
    except ValueError:
        pass
    else:
        failures.append(
            "the shipped audit missed 21K -> pre_omicron, the one mutation a "
            "one-sided predicate is blind to"
        )

    if not failures:
        print(
            "  the one-sided predicate accepts 21K -> pre_omicron; the shipped "
            "audit rejects it"
        )
    return failures


def check_label_format_is_asserted():
    failures = []
    era = _load_era()

    malformed = [
        "recombinant",  # a real config value
        "",  # a real config value
        "21K (BA.1)",  # pre-normalization Nextstrain_clade_usher
        "21",  # no letter
        "A21",  # transposed
        "21k",  # lowercase
        "2A",  # one-digit year
        "210A",  # three-digit year
        "21AB",  # two letters
        " 21K",  # leading space
    ]
    for label in malformed:
        try:
            result = era.era_from_label(label)
        except ValueError:
            continue
        failures.append(
            f"era_from_label({label!r}) returned {result!r} instead of raising"
        )

    for label in ("19A", "20J", "21K", "25C"):
        try:
            era.era_from_label(label)
        except ValueError as exc:
            failures.append(f"era_from_label rejected the valid label {label!r}: {exc}")

    if not failures:
        print(f"  {len(malformed)} malformed labels raise; valid labels do not")
    return failures


def check_year_threshold_sweep():
    """One threshold between 21 and 22, across every two-digit prefix.

    Whether the year is parsed as an int or compared as a string is unobservable
    while the regex pins two digits (see the module docstring, check 7). The
    observable part is that the ordering has exactly one transition, in the right
    place, and that non-two-digit years are rejected.
    """
    failures = []
    era = _load_era()

    for year in range(0, 100):
        label = f"{year:02d}A"
        try:
            result = era.era_from_label(label)
        except ValueError as exc:
            failures.append(f"era_from_label({label!r}) raised: {exc}")
            continue
        expected = era.PRE_OMICRON if year <= 21 else era.OMICRON
        if result != expected:
            failures.append(
                f"era_from_label({label!r}) == {result!r}, expected {expected!r}"
            )

    # The 21K/21L/21M exception applies only within year 21.
    for label in sorted(era.OMICRON_21):
        if era.era_from_label(label) != era.OMICRON:
            failures.append(f"{label} should be Omicron by the OMICRON_21 exception")
    if era.era_from_label("22K") != era.OMICRON:
        failures.append("22K should be Omicron by year, independent of the letter")
    if era.era_from_label("20K") != era.PRE_OMICRON:
        failures.append(
            "20K should be pre-Omicron -- the OMICRON_21 exception must not leak "
            "across years by matching on the letter"
        )

    if not failures:
        print("  years 00-99 sweep cleanly with one transition between 21 and 22")
    return failures


def check_partition_keys_only_present_eras():
    failures = []
    era = _load_era()

    single = era.partition_clades(TEST_CONFIG_CLADES)
    if set(single) != {era.PRE_OMICRON}:
        failures.append(
            f"partition_clades({TEST_CONFIG_CLADES}) keyed {sorted(single)}; the "
            f"test config is single-era, so the Omicron rule must never be built "
            f"(#138 decision 9). An empty list under an 'omicron' key would make "
            f"populate_rates raise on every CI run."
        )
    elif single[era.PRE_OMICRON] != TEST_CONFIG_CLADES:
        failures.append(
            f"partition_clades did not preserve input order: "
            f"{single[era.PRE_OMICRON]} != {TEST_CONFIG_CLADES}"
        )

    both = era.partition_clades(["21J", "21K", "19A", "24A"])
    if set(both) != {era.PRE_OMICRON, era.OMICRON}:
        failures.append(f"a two-era clade list keyed {sorted(both)}")
    elif both[era.PRE_OMICRON] != ["21J", "19A"] or both[era.OMICRON] != ["21K", "24A"]:
        failures.append(f"two-era partition mis-grouped or reordered: {both}")

    if era.partition_clades([]) != {}:
        failures.append("partition_clades([]) should key no eras at all")

    try:
        era.partition_clades(["19A", "21N"])
    except KeyError:
        pass
    else:
        failures.append(
            "partition_clades accepted the unmapped clade 21N; an empty era is a "
            "legitimate skip but an unknown clade never is"
        )

    if not failures:
        print("  partition_clades keys only the eras present and raises on unknowns")
    return failures


CHECKS = [
    ("the mapping covers the release tree", check_mapping_covers_production_tree),
    ("an unmapped clade raises", check_unmapped_clade_raises),
    ("the audit passes on the shipped mapping", check_audit_passes_on_shipped_mapping),
    ("the audit rejects each single-entry mutation", check_audit_rejects_each_mutation),
    ("two-sidedness is load-bearing", check_two_sidedness_is_load_bearing),
    ("the label format is asserted", check_label_format_is_asserted),
    ("the year threshold sweeps cleanly", check_year_threshold_sweep),
    ("partition_clades keys only present eras", check_partition_keys_only_present_eras),
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


def test_era_mapping():
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
        "\nPASSED: the era mapping covers the release tree, raises on unmapped "
        "clades, and is audited by a two-sided year-prefix predicate."
    )


if __name__ == "__main__":
    main()
