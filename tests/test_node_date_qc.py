#!/usr/bin/env python3
"""Pin the internal-node date check primitives (#139).

The checks themselves are **report-only**, so a pipeline run can never fail on
them -- which means a pipeline run also cannot tell you they are computing the
right thing. Everything load-bearing about them is therefore in
``modules/node_date_qc.py`` and asserted here, against the real functions.

What each check exists for:

1. **Exact order statistics from the histogram.** The whole memory argument for
   #139 -- 21 M nodes, 349 s, 38.2 GB, and a per-tip residual file already at
   1.6 GB -- rests on Delta-t being whole days, so a histogram loses *nothing*
   versus sorting every observation. This check asserts that against
   ``numpy.percentile(..., method="lower")`` on the expanded array, over random
   integer data including negatives. If it ever became an approximation, every
   number the follow-up picks a threshold from would be subtly wrong.

2. **Negative values survive.** For checks 6 and 8 a negative value is *the
   interesting case*: an ancestor postdating its earliest observed descendant, or a
   parent dated after a leaf identical to it. Clipping to non-negative would hide
   exactly the anomaly the checks are for, and would still produce a plausible
   distribution.

3. **Bin labels partition the mutation counts.** Every non-negative count lands in
   exactly one bin, bins are contiguous, and the last one is open -- a gap would
   silently drop branches from check 5's denominator.

4. **``propagate_min_descendant`` is a real post-order accumulation.** Tested on a
   tree where the minimum has to travel more than one level, where a subtree has no
   observed leaf at all (must stay ``nan``, not become 0 or +inf), and where the
   minimum arrives at a parent from its *second* child -- the case a naive
   "first non-nan wins" would get wrong. Also checks that it does not mutate its
   input, since the caller reuses ``own_observed``.

5. **The reverse-iteration precondition is asserted, not assumed.** The function
   relies on parents preceding children, which is what ``depth_first_expansion``
   gives. A shuffled ordering must produce a wrong answer -- if it did not, the
   function would not be doing the propagation it claims to.

Usage:
    python tests/test_node_date_qc.py
"""

import sys
import traceback
from importlib import util as importlib_util
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
MODULE = REPO_ROOT / "workflow" / "modules" / "node_date_qc.py"


def _load():
    spec = importlib_util.spec_from_file_location("node_date_qc", MODULE)
    module = importlib_util.module_from_spec(spec)
    sys.modules["node_date_qc"] = module
    spec.loader.exec_module(module)
    return module


def check_histogram_matches_exact_quantiles():
    """The histogram must give the same order statistics as the full array."""
    failures = []
    qc = _load()
    rng = np.random.default_rng(20260806)

    for trial, (size, low, high) in enumerate(
        [(1, 0, 1), (2, -5, 5), (37, 0, 3), (500, -30, 400), (5000, -2, 900)]
    ):
        values = rng.integers(low, high + 1, size=size)
        hist = qc.DayHistogram()
        hist.update(values)

        if hist.n != len(values):
            failures.append(f"trial {trial}: n={hist.n}, expected {len(values)}")
            continue
        stats = hist.stats()
        expected = {
            "min": int(values.min()),
            "max": int(values.max()),
            "median": int(np.percentile(values, 50, method="lower")),
            "p05": int(np.percentile(values, 5, method="lower")),
            "p95": int(np.percentile(values, 95, method="lower")),
        }
        for key, want in expected.items():
            if stats[key] != want:
                failures.append(
                    f"trial {trial} ({size} values in [{low},{high}]): "
                    f"{key}={stats[key]}, exact answer is {want}. The histogram is "
                    f"approximating, which invalidates the memory argument for "
                    f"binning at all."
                )
        if abs(stats["mean"] - values.mean()) > 1e-9:
            failures.append(
                f"trial {trial}: mean={stats['mean']}, expected {values.mean()}"
            )
        if sum(c for _, c in hist.as_rows()) != len(values):
            failures.append(f"trial {trial}: as_rows() does not sum to n")

    empty = qc.DayHistogram().stats()
    if empty["n"] != 0 or empty["median"] is not None:
        failures.append(f"an empty histogram must report n=0 and None: {empty}")

    if not failures:
        print("  histogram order statistics are exact over 5 random trials")
    return failures


def check_negatives_survive():
    """Negative day values are the anomaly signal and must not be clipped."""
    failures = []
    qc = _load()
    hist = qc.DayHistogram()
    hist.update([-40, -1, 0, 0, 3])
    stats = hist.stats()
    if stats["min"] != -40:
        failures.append(
            f"min={stats['min']}, expected -40; negative values are being clipped, "
            f"which would hide an ancestor postdating its earliest sample"
        )
    if abs(stats["frac_negative"] - 0.4) > 1e-12:
        failures.append(
            f"frac_negative={stats['frac_negative']}, expected 0.4 (2 of 5)"
        )
    if abs(stats["mean"] - (-38 / 5)) > 1e-12:
        failures.append(f"mean={stats['mean']}, expected {-38 / 5}")
    if not failures:
        print("  negative day values are retained and counted")
    return failures


def check_bins_partition_mutation_counts():
    failures = []
    qc = _load()
    seen = {}
    for n in range(0, 200):
        label = qc.mutation_count_bin(n)
        seen.setdefault(label, []).append(n)
    for label, members in seen.items():
        if members != list(range(members[0], members[-1] + 1)):
            failures.append(f"bin {label!r} is not contiguous: {members[:5]}...")
    covered = sorted(n for members in seen.values() for n in members)
    if covered != list(range(200)):
        failures.append("mutation counts 0-199 are not each in exactly one bin")
    if qc.mutation_count_bin(0) != "0":
        failures.append(
            f"count 0 is labelled {qc.mutation_count_bin(0)!r}, but check 6 reads "
            f"the '0' bin by name"
        )
    if qc.mutation_count_bin(10**6) != qc.mutation_count_bin(MAX_PROBE := 10**7):
        failures.append("the last bin is not open-ended")
    if not failures:
        print(f"  {len(seen)} bins partition mutation counts 0-199, last is open")
    return failures


def check_min_descendant_propagation():
    """A real post-order accumulation, including the cases a shortcut gets wrong."""
    failures = []
    qc = _load()

    #  0 root
    #  |-- 1 internal
    #  |    |-- 2 leaf (obs 50)
    #  |    +-- 3 leaf (obs 20)   <- minimum arrives from the SECOND child
    #  +-- 4 internal
    #       |-- 5 internal
    #       |    +-- 6 leaf (obs 90)  <- must travel two levels to node 4
    #       +-- 7 leaf (no observation)
    #  8 internal, 9 leaf: a subtree with NO observed leaf -> must stay nan
    parent = np.array([-1, 0, 1, 1, 0, 4, 5, 4, 0, 8], dtype="int64")
    own = np.array(
        [np.nan, np.nan, 50.0, 20.0, np.nan, np.nan, 90.0, np.nan, np.nan, np.nan]
    )
    own_before = own.copy()
    got = qc.propagate_min_descendant(parent, own, len(parent))

    expected = {0: 20.0, 1: 20.0, 2: 50.0, 3: 20.0, 4: 90.0, 5: 90.0, 6: 90.0}
    for node, want in expected.items():
        if got[node] != want:
            failures.append(f"node {node}: got {got[node]}, expected {want}")
    for node in (7, 8, 9):
        if not np.isnan(got[node]):
            failures.append(
                f"node {node} has no observed descendant, so it must stay nan; "
                f"got {got[node]}. A 0 or +inf here would be compared against a "
                f"real date and reported as an anomaly."
            )
    if not np.array_equal(own, own_before, equal_nan=True):
        failures.append("propagate_min_descendant mutated its leaf_values input")

    try:
        qc.propagate_min_descendant(parent, own[:-1], len(parent))
    except ValueError:
        pass
    else:
        failures.append("mismatched array lengths did not raise")

    if not failures:
        print("  min-descendant propagation handles multi-level, second-child, "
              "and unobserved-subtree cases")
    return failures


def check_ordering_precondition_matters():
    """A shuffled ordering must give a wrong answer.

    The function's correctness rests on parents preceding children, which
    ``depth_first_expansion`` guarantees. If a shuffled order produced the same
    result, the single reverse pass would not be doing the accumulation it claims
    and the check would be vacuous.
    """
    failures = []
    qc = _load()

    # A chain 0 -> 1 -> 2 -> 3, observed only at the deepest leaf. In DFS order the
    # value must reach the root; with children BEFORE parents it cannot.
    chain_parent = np.array([-1, 0, 1, 2], dtype="int64")
    chain_own = np.array([np.nan, np.nan, np.nan, 7.0])
    forward = qc.propagate_min_descendant(chain_parent, chain_own, 4)
    if forward[0] != 7.0:
        failures.append(
            f"in parent-before-child order the leaf value must reach the root; "
            f"got {forward[0]}"
        )

    # Same tree indexed in reverse (3 -> 2 -> 1 -> 0 becomes 0 <- 1 <- 2 <- 3 with
    # parent indices pointing FORWARD), which violates the precondition.
    bad_parent = np.array([1, 2, 3, -1], dtype="int64")
    bad_own = np.array([7.0, np.nan, np.nan, np.nan])
    reverse = qc.propagate_min_descendant(bad_parent, bad_own, 4)
    if reverse[3] == 7.0:
        failures.append(
            "a child-before-parent ordering still propagated to the root, so the "
            "single reverse pass is not order-dependent and this test proves "
            "nothing about the real call site"
        )
    if not failures:
        print("  propagation is order-dependent, as the single reverse pass requires")
    return failures


def check_constants_are_sane():
    failures = []
    qc = _load()
    if not 10 <= qc.DAYS_PER_MUTATION_NOMINAL <= 25:
        failures.append(
            f"DAYS_PER_MUTATION_NOMINAL is {qc.DAYS_PER_MUTATION_NOMINAL:.2f}, "
            f"expected ~15 from 8e-4 subs/site/yr x 29903 nt"
        )
    lo, hi = qc.ROOT_DATE_PLAUSIBLE
    if not (lo < "2019-12-09" < hi):
        failures.append(
            f"the plausible root window {lo}..{hi} excludes v9's actual root date "
            f"2019-12-09, so check 7 would report OUTSIDE on a good run"
        )
    if not failures:
        print(f"  nominal clock {qc.DAYS_PER_MUTATION_NOMINAL:.1f} d/mut; "
              f"root window brackets v9's 2019-12-09")
    return failures


CHECKS = [
    ("histogram order statistics are exact", check_histogram_matches_exact_quantiles),
    ("negative day values survive", check_negatives_survive),
    ("bins partition mutation counts", check_bins_partition_mutation_counts),
    ("min-descendant propagation", check_min_descendant_propagation),
    ("propagation depends on node ordering", check_ordering_precondition_matters),
    ("constants are sane", check_constants_are_sane),
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


def test_node_date_qc():
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
        "\nPASSED: the internal-node date primitives give exact statistics, keep "
        "negative values, and propagate descendant minima in true post-order."
    )


if __name__ == "__main__":
    main()
