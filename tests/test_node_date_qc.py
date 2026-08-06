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
   gives. A violated ordering must **raise**. It used to return silently, and
   what it returned was plausible -- ``nan`` at the root, indistinguishable from a
   genuinely unobserved subtree -- so a traversal-order change in ``bte`` would
   have degraded checks 6 and 8 with no symptom.

6. **The quantile boundary.** Item 1's random trials do not pin the ``seen >
   target`` comparison: mutating it to ``seen >= target`` leaves them all green,
   because the bug needs a cumulative count to land *exactly* on an integer
   target. That coincidence is constructed on purpose in
   ``check_quantile_boundary_is_exact`` rather than hoped for.

7. **``min_child_delta`` is the nearest-*descendant* Delta-t only because dates are
   monotone.** Check 9 is named for descendants but computes a minimum over
   immediate children; that is valid only because ``time_monotonicity`` forces
   predicted dates to be nondecreasing along every edge. Both halves are asserted:
   agreement with an independent brute-force min-over-all-descendants on a monotone
   tree, and **disagreement** on a tree with one decreasing edge. Without the second
   half the first is just two spellings of the same loop.

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
            "median": int(np.quantile(values, 0.5, method="lower")),
            "p05": int(np.quantile(values, 0.05, method="lower")),
            # The quartiles feed node_dates_branch_dt_summary.csv's median/IQR/p95
            # table, so they need the same exactness guarantee as the median.
            "p25": int(np.quantile(values, 0.25, method="lower")),
            "p75": int(np.quantile(values, 0.75, method="lower")),
            "p95": int(np.quantile(values, 0.95, method="lower")),
        }
        if stats["iqr"] != stats["p75"] - stats["p25"]:
            failures.append(
                f"trial {trial}: iqr={stats['iqr']} but p75-p25="
                f"{stats['p75'] - stats['p25']}; iqr must be derived from the same "
                f"exact order statistics, not estimated separately"
            )
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

        # update_from_array is the production path (the per-node Python loop cost
        # minutes at 21 M nodes), so it must agree with update() exactly rather
        # than merely approximately. Fed floats with nan padding, as the script
        # does, to pin that nan is dropped and not counted as a day value.
        bulk = qc.DayHistogram()
        padded = np.concatenate([values.astype("float64"), [np.nan, np.nan]])
        bulk.update_from_array(padded)
        if bulk.as_rows() != hist.as_rows():
            failures.append(
                f"trial {trial}: update_from_array disagrees with update(). "
                f"n={bulk.n} vs {hist.n}. The vectorized bulk loader is the path "
                f"the pipeline actually takes."
            )
        if bulk.stats() != stats:
            failures.append(
                f"trial {trial}: update_from_array stats differ from update stats"
            )

    empty = qc.DayHistogram().stats()
    if empty["n"] != 0 or empty["median"] is not None:
        failures.append(f"an empty histogram must report n=0 and None: {empty}")

    if not failures:
        print("  histogram order statistics are exact over 5 random trials")
    return failures


def check_quantile_boundary_is_exact():
    """Pin the ``seen > target`` comparison at targets that land on an integer.

    ``quantile`` walks the sorted day values accumulating counts and returns the
    first value whose cumulative count exceeds ``q * (n - 1)``. The comparison has
    to be strict: when the target is an exact integer *and* a cumulative count
    equals it, ``>=`` returns the value one rank too low.

    The random trials above cannot see this. Mutating the shipped comparison to
    ``seen >= target`` leaves every one of them green -- verified by running this
    file against a mutated copy of ``node_date_qc.py`` -- because hitting the bug
    needs a cumulative count to fall *exactly* on an integer target, which random
    integers in a wide range essentially never do at the three quantiles sampled.
    So this check constructs that coincidence on purpose instead of hoping for it.

    ``[0]*k + [1]*(n-k)`` with ``q = k/(n-1)`` gives ``target == k`` exactly, and
    the cumulative count after the zeros is exactly ``k`` -- the boundary, every
    time. Exactness is the whole justification for binning rather than keeping
    21 M per-branch rows (see the module docstring), so it deserves a test that
    fails when it breaks.
    """
    failures = []
    qc = _load()

    # merge() builds the summary table's "all" row from the per-bin histograms, so
    # its order statistics must equal those of the pooled data -- not an average of
    # per-bin medians. Checked here because a wrong merge produces a plausible table.
    a, b = qc.DayHistogram(), qc.DayHistogram()
    a.update([1, 1, 5]), b.update([2, 9, 9, 100])
    pooled = qc.DayHistogram()
    pooled.update([1, 1, 5, 2, 9, 9, 100])
    if qc.DayHistogram().merge(a).merge(b).stats() != pooled.stats():
        failures.append(
            "merge() does not reproduce the pooled distribution, so the summary "
            "table's 'all' row would not be the marginal of its bins"
        )

    for n in range(2, 60):
        for k in range(1, n):
            q = k / (n - 1)
            if q > 1.0:
                continue
            values = [0] * k + [1] * (n - k)
            hist = qc.DayHistogram()
            hist.update(values)
            got = hist.quantile(q)
            # np.quantile, not np.percentile: percentile takes q*100 and divides
            # by 100 internally, and that round trip is lossy. At q=1/3, n=4 it
            # yields a virtual index of 0.9999999999999999 which floors to 0,
            # while quantile() computes q*(n-1) == 1.0 exactly. Comparing against
            # the percentile form reports a failure against correct code -- the
            # reference has to do the same arithmetic as the implementation.
            want = int(np.quantile(values, q, method="lower"))
            if got != want:
                failures.append(
                    f"n={n}, k={k}, q={q!r}: quantile={got}, exact answer is {want}. "
                    f"target = q*(n-1) = {q * (n - 1)!r} and the cumulative count "
                    f"after the {k} zero(s) is exactly {k}, so the comparison "
                    f"against the target must be strict."
                )
                if len(failures) > 5:
                    failures.append("... further boundary failures suppressed")
                    return failures

    # The minimal hand-checked case, stated explicitly so the intent survives a
    # future refactor of the sweep above.
    hist = qc.DayHistogram()
    hist.update([1, 1, 1, 1, 2, 2, 2, 2, 2])
    if hist.quantile(0.5) != 2:
        failures.append(
            f"median of [1]*4+[2]*5 is {hist.quantile(0.5)}, must be 2 "
            f"(numpy lower); target is 0.5*(9-1) = 4.0 and the count of 1s is 4"
        )

    if not failures:
        print(
            "  quantile is exact at integer targets (1711 constructed boundary "
            "cases), where the random trials are blind"
        )
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
    """A violated ordering must RAISE, not return a wrong answer.

    The function's correctness rests on parents preceding children, which
    ``depth_first_expansion`` guarantees. Two things are asserted here, and the
    second is the one that matters:

    1. In a valid order the deepest value reaches the root -- otherwise the single
       reverse pass is not doing the accumulation it claims.
    2. In a violated order the function **raises**. It used to return silently,
       and what it returned was plausible: for the reversed chain below, ``nan``
       at the root instead of ``7.0``. Nothing downstream could distinguish that
       from a genuinely unobserved subtree, so a change in ``bte``'s traversal
       order would have degraded checks 6 and 8 with no symptom -- the
       #129/#130/#133 shape. Asserting the raise, rather than asserting the wrong
       answer, is what keeps that from coming back.
    """
    failures = []
    qc = _load()

    # A chain 0 -> 1 -> 2 -> 3, observed only at the deepest leaf. In DFS order the
    # value must reach the root.
    chain_parent = np.array([-1, 0, 1, 2], dtype="int64")
    chain_own = np.array([np.nan, np.nan, np.nan, 7.0])
    forward = qc.propagate_min_descendant(chain_parent, chain_own, 4)
    if forward[0] != 7.0:
        failures.append(
            f"in parent-before-child order the leaf value must reach the root; "
            f"got {forward[0]}"
        )

    # Same tree indexed in reverse (parent indices point FORWARD), violating the
    # precondition. This must raise rather than under-propagate.
    bad_parent = np.array([1, 2, 3, -1], dtype="int64")
    bad_own = np.array([7.0, np.nan, np.nan, np.nan])
    try:
        reverse = qc.propagate_min_descendant(bad_parent, bad_own, 4)
    except ValueError as exc:
        if "depth-first" not in str(exc):
            failures.append(
                f"a child-before-parent ordering raised, but the message does not "
                f"explain the ordering requirement: {exc}"
            )
    else:
        failures.append(
            f"a child-before-parent ordering returned {reverse.tolist()} instead "
            f"of raising. Whether or not the value is wrong, returning it means a "
            f"traversal-order change in bte degrades checks 6 and 8 with no "
            f"symptom at all."
        )

    # A single bad edge in an otherwise valid ordering must also be caught -- the
    # realistic failure is one misplaced node, not a wholesale reversal.
    one_bad = np.array([-1, 0, 3, 1], dtype="int64")
    try:
        qc.propagate_min_descendant(one_bad, np.array([np.nan, np.nan, 5.0, 6.0]), 4)
    except ValueError:
        pass
    else:
        failures.append(
            "a single out-of-order edge (node 2's parent at index 3) was accepted; "
            "the guard must be per-edge, not a whole-array heuristic"
        )

    if not failures:
        print("  a violated node ordering raises instead of under-propagating")
    return failures


def check_min_child_delta_equals_min_descendant():
    """``min_child_delta`` is the nearest-*descendant* Delta-t only if dates are monotone.

    Check 9 reports "Delta-t to the nearest descendant" but computes a minimum over
    immediate *children*. That shortcut is valid because ``time_monotonicity``
    forces predicted dates to be nondecreasing along every edge, so the nearest
    descendant is always a child. This asserts both halves of that:

    * on a monotone tree it agrees with the independent brute-force route
      (``propagate_min_descendant`` over each subtree, excluding the node itself);
    * on a tree with one *decreasing* edge it disagrees -- so the monotonicity
      precondition is load-bearing, not decoration. Without this second half the
      first is just two spellings of the same loop.
    """
    failures = []
    qc = _load()

    def min_descendant_delta(parent, values, n):
        """Brute force: min over ALL strict descendants, via child lists."""
        kids = {i: [] for i in range(n)}
        for i, p in enumerate(parent):
            if p >= 0:
                kids[p].append(i)
        out = np.full(n, np.nan)
        for i in range(n):
            stack, best = list(kids[i]), np.inf
            while stack:
                j = stack.pop()
                best = min(best, values[j] - values[i])
                stack.extend(kids[j])
            out[i] = best if np.isfinite(best) else np.nan
        return out

    # Monotone: a root with two children, one of which has a deep chain. The deep
    # chain is what makes "nearest descendant" differ from "any descendant".
    parent = np.array([-1, 0, 0, 1, 3, 2], dtype="int64")
    dates = np.array([0.0, 5.0, 30.0, 6.0, 40.0, 31.0])
    got = qc.min_child_delta(parent, dates, 6)
    want = min_descendant_delta(parent, dates, 6)
    if not np.allclose(got[~np.isnan(want)], want[~np.isnan(want)]):
        failures.append(
            f"on a monotone tree min_child_delta={got.tolist()} disagrees with the "
            f"min over all descendants {want.tolist()}; the check-9 shortcut is "
            f"not valid"
        )
    if not np.array_equal(np.isnan(got), np.isnan(want)):
        failures.append(
            f"nan placement differs: {np.isnan(got).tolist()} vs "
            f"{np.isnan(want).tolist()} (leaves must be nan, not 0 or inf)"
        )
    # Node 0's nearest descendant is child 1 at +5, NOT the deeper node 4 at +40.
    if got[0] != 5.0:
        failures.append(f"root's nearest-descendant delta should be 5.0, got {got[0]}")

    # Non-monotone: node 4 moves BEFORE its parent 3. min-over-children for node 3
    # is now negative, but node 0's min over children is unchanged at +5 while its
    # min over all descendants is -1. The two routes must part company.
    bad_dates = np.array([0.0, 5.0, 30.0, 6.0, -1.0, 31.0])
    got_bad = qc.min_child_delta(parent, bad_dates, 6)
    want_bad = min_descendant_delta(parent, bad_dates, 6)
    if got_bad[0] == want_bad[0]:
        failures.append(
            "with a decreasing edge the min over children still equals the min over "
            "descendants, so this test cannot show that monotonicity is what makes "
            "check 9's shortcut valid"
        )

    # Undated children must not be counted as a 0-day delta.
    with_nan = np.array([0.0, np.nan, 12.0, np.nan, np.nan, np.nan])
    got_nan = qc.min_child_delta(np.array([-1, 0, 0, 1, 3, 2], dtype="int64"), with_nan, 6)
    if got_nan[0] != 12.0:
        failures.append(
            f"a node with one undated and one dated child must report the dated "
            f"delta (12.0), got {got_nan[0]}"
        )

    if len(qc.min_child_delta(np.array([-1], dtype="int64"), np.array([1.0]), 1)) != 1:
        failures.append("a single-node tree must still return an array of length 1")
    try:
        qc.min_child_delta(np.array([-1, 0], dtype="int64"), np.array([1.0]), 2)
    except ValueError:
        pass
    else:
        failures.append("a length mismatch must raise")

    if not failures:
        print("  min_child_delta is the nearest-descendant delta iff dates are monotone")
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
    ("quantile is exact at integer targets", check_quantile_boundary_is_exact),
    ("negative day values survive", check_negatives_survive),
    ("bins partition mutation counts", check_bins_partition_mutation_counts),
    ("min-descendant propagation", check_min_descendant_propagation),
    ("a violated node ordering raises", check_ordering_precondition_matters),
    ("min_child_delta vs min descendant", check_min_child_delta_equals_min_descendant),
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
