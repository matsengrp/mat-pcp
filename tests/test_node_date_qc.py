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



def check_aggregate_sees_what_the_slope_cannot():
    """The whole point of #146's addition, on a tree built to fool the old report.

    Branch durations here are constructed so the **incremental** slope between
    adjacent single-count bins is 15 d/mutation -- indistinguishable from the
    nominal clock, exactly what check 5 reported on the v11 trees -- while a large
    zero-mutation intercept inflates total elapsed time. If the aggregate were
    computed the way the slope is (differences, intercept cancelled) or over only
    the mutation-carrying branches, it would come out at 10.4 subs/yr and agree that
    the tree is fine. Over all branches it is 6.64, a 1.57x deficit that is entirely
    the intercept.

    So this pins the discrimination, not just the arithmetic: the two statistics
    must **disagree** on this input. A version of ``implied_clock_summary`` that
    dropped the zero-mutation branches from either sum would leave every other
    assertion in this file green.
    """
    failures = []
    qc = _load()
    dt = np.array([30.0] * 1000 + [45.0] * 500 + [60.0] * 500)
    muts = np.array([0] * 1000 + [1] * 500 + [2] * 500)
    summary = qc.implied_clock_summary(dt, muts)

    # The old report's statistic, computed the way check 5 computes it.
    incremental_slope = 60.0 - 45.0
    if not 10 <= incremental_slope <= 25:
        failures.append(
            "the fixture is not adversarial: its incremental slope must look "
            "healthy against the nominal clock for this check to mean anything"
        )

    total_days = 1000 * 30 + 500 * 45 + 500 * 60
    expected_rate = 1500 / (total_days / qc.DAYS_PER_YEAR)
    if abs(summary["subs_per_year"] - expected_rate) > 1e-9:
        failures.append(
            f"aggregate rate is {summary['subs_per_year']:.4f}, expected "
            f"{expected_rate:.4f} = 1500 mutations / {total_days} branch-days; the "
            f"sums must run over ALL branches including zero-mutation ones"
        )
    rate_excluding_intercept = 1500 / ((500 * 45 + 500 * 60) / qc.DAYS_PER_YEAR)
    if abs(summary["subs_per_year"] - rate_excluding_intercept) < 1e-6:
        failures.append(
            "the aggregate matches the intercept-excluded rate, so it is blind to "
            "the same term the incremental slope cancels -- the addition is vacuous"
        )
    nominal = qc.DAYS_PER_YEAR / qc.DAYS_PER_MUTATION_NOMINAL
    if summary["subs_per_year"] / nominal > 0.5:
        failures.append(
            f"the fixture should read as badly depressed against nominal; got "
            f"{summary['subs_per_year'] / nominal:.2f}x"
        )
    expected_time_share = (1000 * 30) / total_days
    if abs(summary["frac_branch_time_zero_mut"] - expected_time_share) > 1e-9:
        failures.append(
            f"zero-mutation branch-time share is "
            f"{summary['frac_branch_time_zero_mut']:.4f}, expected "
            f"{expected_time_share:.4f}"
        )
    # Share of TIME, not of branches: the two differ here (36.4% vs 50%), and only
    # the former indicts the fit -- most branches carrying no mutations is normal
    # under dense sampling.
    if abs(summary["frac_branches_zero_mut"] - 0.5) > 1e-9:
        failures.append(
            f"zero-mutation branch COUNT share is "
            f"{summary['frac_branches_zero_mut']:.4f}, expected 0.5"
        )
    if not failures:
        print(
            f"  slope {incremental_slope:.0f} d/mut looks nominal while the "
            f"aggregate reads {summary['subs_per_year']:.2f} subs/yr "
            f"({summary['subs_per_year'] / nominal:.2f}x nominal), "
            f"{100 * summary['frac_branch_time_zero_mut']:.1f}% of time on "
            f"zero-mutation branches"
        )
    return failures


def check_aggregate_is_exact_not_binned():
    """Exact per-branch sums, not a reconstruction from the ranged bins.

    ``node_dates_branch_dt_summary.csv`` reports ``6-10``, ``11-20``, ``21-50`` and
    ``51+``, so recovering the aggregate from it requires inventing a representative
    count per bin -- which moved the answer by ~+/-7% on the production tree. The
    arrays are in hand at the call site, so the exact sum is free; this asserts the
    exact value and that it is distinguishable from the midpoint reconstruction.
    """
    failures = []
    qc = _load()
    muts = np.array([7, 10, 60])
    dt = np.array([100.0, 100.0, 100.0])
    summary = qc.implied_clock_summary(dt, muts)
    if summary["total_mutations"] != 77.0:
        failures.append(
            f"total_mutations is {summary['total_mutations']}, expected 77 = 7+10+60"
        )
    # Both 7 and 10 fall in the "6-10" bin, so a reconstruction has one number to
    # spend on them; 8 is the midpoint. The exact sum differs by 1.
    midpoint_total = 8 + 8 + 60
    if summary["total_mutations"] == midpoint_total:
        failures.append(
            "total_mutations equals the binned reconstruction, so the fixture "
            "cannot tell exact from binned; pick counts that differ from midpoints"
        )
    if not failures:
        print(f"  exact 77 mutations, distinguishable from binned {midpoint_total}")
    return failures


def check_undefined_branches_and_zero_span():
    """NaN Delta-t is excluded, and a zero-span tree yields no rate at all."""
    failures = []
    qc = _load()
    summary = qc.implied_clock_summary([10.0, np.nan, 20.0], [1, 99, 1])
    if summary["n_branches"] != 2:
        failures.append(
            f"n_branches is {summary['n_branches']}, expected 2: a NaN Delta-t "
            f"branch has no duration and must leave both sums"
        )
    if summary["total_mutations"] != 2.0:
        failures.append(
            f"total_mutations is {summary['total_mutations']}, expected 2 -- the "
            f"NaN branch's 99 mutations must not enter the numerator while its "
            f"(undefined) duration stays out of the denominator"
        )
    collapsed = qc.implied_clock_summary([0.0, 0.0], [1, 2])
    if collapsed["subs_per_year"] is not None:
        failures.append(
            f"a zero-span tree reported a rate of {collapsed['subs_per_year']}; it "
            f"must be None, since inf or NaN would render as a number in the report"
        )
    note = qc.implied_clock_note(collapsed, 39.36, "irrelevant.stdout")
    if "undefined" not in note:
        failures.append(f"the zero-span note must say so; got {note!r}")
    if not failures:
        print("  NaN branches excluded from both sums; zero-span gives no rate")
    return failures


def check_root_to_tip_reference_is_never_silently_dropped():
    """Absent reference must be *stated*. This is the non-vacuity guard.

    ``log/chronumental.stdout`` is a Snakemake ``log:`` entry, not a declared
    output, so it is legitimately absent on an archived or re-run tree. The failure
    mode is not the parse returning ``None`` -- it is a report line that omits the
    comparison, which reads identically to one where the comparison was fine.
    """
    failures = []
    qc = _load()
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        good = Path(tmp) / "chronumental.stdout"
        good.write_text(
            "Doing a rough initial fit\n"
            f"{qc.ROOT_TO_TIP_PREFIX} 39.3631189049384\n"
            "step:0\tloss:1.0\tlength_cor:0.9949\n"
        )
        if qc.parse_root_to_tip_rate(str(good)) != 39.3631189049384:
            failures.append(
                f"failed to parse the real stdout line format; got "
                f"{qc.parse_root_to_tip_rate(str(good))}"
            )
        empty = Path(tmp) / "no-such-line.stdout"
        empty.write_text("step:0\tloss:1.0\n")
        for label, path in (
            ("a missing file", str(Path(tmp) / "absent.stdout")),
            ("a log without the line", str(empty)),
        ):
            if qc.parse_root_to_tip_rate(path) is not None:
                failures.append(f"{label} must parse to None")
            summary = qc.implied_clock_summary([100.0], [1])
            note = qc.implied_clock_note(summary, None, path)
            if "unavailable" not in note or path not in note:
                failures.append(
                    f"with {label} the note must say the reference is unavailable "
                    f"and name the path it looked at; got {note!r}"
                )
            if "x chronumental's own" in note:
                failures.append(
                    f"with {label} the note claims a root-to-tip ratio: {note!r}"
                )
        present = qc.implied_clock_note(
            qc.implied_clock_summary([365.25], [1]), 4.0, str(good)
        )
        if "0.25x chronumental's own pre-fit root-to-tip 4.00" not in present:
            failures.append(
                f"with a reference available the note must carry the ratio; got "
                f"{present!r}"
            )
    if not failures:
        print("  reference present -> ratio; absent -> stated, with the path")
    return failures



def check_rtt_divergence_accumulates_and_excludes_the_root():
    """Root-to-tip divergence sums the path, and the root's own mutations are not on it.

    Three properties, and the middle one is the load-bearing surprise:

    1. Divergence accumulates down a path (a two-branch chain sums both branches),
       and a sibling's count does not leak across.
    2. The root's ``branch_mutations`` entry is **ignored**. On a MAT the root's
       mutation list is its difference from the reference, not evolution inside the
       tree, so counting it would offset every tip by the same constant -- which is
       invisible in the plot (the cloud shifts up as a block, the slope is
       unchanged) and wrong in the regression's intercept. A constant offset is
       exactly the kind of error a scatter cannot show you, so it is asserted here.
    3. A child-before-parent ordering **raises** rather than returning a truncated
       sum, via the guard shared with ``propagate_min_descendant``.
    """
    failures = []
    qc = _load()

    # 0 (root) -> 1 -> 3, and 0 -> 2. Root carries 7 "mutations" that must not count.
    parent = np.array([-1, 0, 0, 1], dtype="int64")
    muts = np.array([7, 2, 5, 3], dtype="int64")
    div = qc.root_to_tip_divergence(parent, muts, 4)
    expected = [0, 2, 5, 5]
    if div.tolist() != expected:
        failures.append(
            f"root-to-tip divergence must be {expected} (root excluded, path "
            f"summed, siblings independent); got {div.tolist()}"
        )

    # Same tree with a different root mutation count: every value must be unchanged.
    div_other_root = qc.root_to_tip_divergence(
        parent, np.array([999, 2, 5, 3], dtype="int64"), 4
    )
    if div_other_root.tolist() != div.tolist():
        failures.append(
            f"the root's own mutation count must not enter any divergence; "
            f"changing it moved {div.tolist()} to {div_other_root.tolist()}"
        )

    bad_parent = np.array([1, 2, 3, -1], dtype="int64")
    try:
        got = qc.root_to_tip_divergence(bad_parent, muts, 4)
    except ValueError as exc:
        if "depth-first" not in str(exc):
            failures.append(
                f"a child-before-parent ordering raised without explaining the "
                f"ordering requirement: {exc}"
            )
    else:
        failures.append(
            f"a child-before-parent ordering returned {got.tolist()} instead of "
            f"raising. The forward pass truncates silently without the "
            f"precondition, which would understate divergence for a whole subtree "
            f"and read as a flatter clock."
        )

    if not failures:
        print("  divergence sums the path, ignores the root's own mutations, "
              "and raises on a bad ordering")
    return failures


def check_rtt_histogram_is_lossless():
    """The binned scatter must give the SAME regression as the per-tip cloud.

    This is the property the whole design rests on: the figure ships an exact joint
    histogram (a few tens of kB) instead of 17 M points, and that is only legitimate
    if every statistic drawn from it is identical rather than approximate. Both axes
    are integers, so it should be exact -- but "should be" covers a key-encoding
    error that collides two cells, and a collision would show up as a *slightly*
    wrong slope, which is unfalsifiable by eye.

    So: build a raw point cloud with deliberate duplicates, regress it unweighted,
    regress the histogram weighted, and require agreement to floating-point
    tolerance. The fixture also puts two tips on the same date with different
    divergences and two tips with the same divergence on different dates, which is
    what a bad key encoding mixes up.
    """
    failures = []
    qc = _load()

    days = np.array([10, 10, 10, 12, 12, 40, 40, 55], dtype="int64")
    muts = np.array([0, 0, 3, 0, 4, 3, 9, 12], dtype="int64")

    h_days, h_muts, h_counts = qc.root_to_tip_histogram(days, muts)
    if int(h_counts.sum()) != len(days):
        failures.append(
            f"the histogram must account for every tip: {int(h_counts.sum())} "
            f"counted vs {len(days)} given"
        )
    # (10, 0) appears twice; every other pair once -> 7 cells.
    if len(h_days) != 7:
        failures.append(f"expected 7 non-empty cells, got {len(h_days)}")
    duplicated = [
        int(c) for d, m, c in zip(h_days, h_muts, h_counts) if d == 10 and m == 0
    ]
    if duplicated != [2]:
        failures.append(f"cell (10, 0) must hold 2 tips, got {duplicated}")

    binned = qc.weighted_linregress(h_days, h_muts, h_counts)
    raw = qc.weighted_linregress(days, muts, np.ones(len(days)))
    for key in ("slope", "intercept", "r"):
        if not np.isclose(binned[key], raw[key], rtol=1e-12, atol=1e-12):
            failures.append(
                f"binning changed the regression's {key}: {binned[key]} from the "
                f"histogram vs {raw[key]} from the raw points. The histogram is "
                f"only a legitimate stand-in for the cloud if these are equal."
            )
    if binned["n"] != float(len(days)):
        failures.append(f"weighted n must be {len(days)}, got {binned['n']}")

    # A NaN date (a tip with no parseable observed date) is dropped, not counted as
    # day zero -- which would drag the slope toward 1970.
    with_nan_days = np.array([10.0, np.nan, 12.0])
    dropped = qc.root_to_tip_histogram(with_nan_days, np.array([1, 2, 3]))
    if int(dropped[2].sum()) != 2:
        failures.append(
            f"a tip with no parseable date must be dropped, not binned at the "
            f"epoch; got {int(dropped[2].sum())} tips instead of 2"
        )

    if not failures:
        print("  the exact histogram reproduces the raw cloud's regression bit for bit")
    return failures


def check_rtt_slope_is_per_year_and_reference_is_never_dropped():
    """The reported slope is subs/**year**, and a missing reference is stated.

    A per-day slope reads as ~0.065 where the reference is ~40, so a units slip does
    not look like a units slip -- it looks like a catastrophically slow clock, which
    is precisely the finding this whole check exists to report. The fixture is
    therefore a cloud with an exactly known rate.

    The second half mirrors ``check_root_to_tip_reference_is_never_silently_dropped``
    for the new note: with chronumental's slope available the note carries the ratio,
    and without it the note says so **and names the path it looked at**.
    """
    failures = []
    qc = _load()

    # Exactly 24 mutations per 365.25 days, over 5 dates.
    days = np.arange(0, 5, dtype="int64") * 365
    muts = (24.0 / qc.DAYS_PER_YEAR) * days
    reg = qc.weighted_linregress(days, muts, np.ones(len(days)))
    per_year = reg["slope"] * qc.DAYS_PER_YEAR
    if not np.isclose(per_year, 24.0, rtol=1e-9):
        failures.append(
            f"a cloud built at exactly 24 subs/yr must regress to 24 subs/yr "
            f"after the DAYS_PER_YEAR conversion; got {per_year}"
        )

    note = qc.root_to_tip_note(reg, 3, 30.0, "/nonexistent/chronumental.stdout")
    if "24.0" not in note or "ratio" not in note:
        failures.append(
            f"with a reference available the note must carry the recomputed slope "
            f"and the ratio; got {note!r}"
        )
    if "3 excluded" not in note:
        failures.append(
            f"the note must say how many tips were excluded, so a regression over "
            f"a small subset is not read as the whole tree; got {note!r}"
        )
    # The tripwire is two-sided, and both sides are checked: a units error or an
    # off-by-a-tree-level divergence can land either high or low, and a guard that
    # only fires one way is half a guard.
    for logged, label in ((12.0, "2.00 (high)"), (40.0, "0.60 (low)")):
        flagged = qc.root_to_tip_note(reg, 0, logged, "/nonexistent")
        if "tripwire" not in flagged:
            failures.append(
                f"a ratio of {label} is outside RTT_SLOPE_REL_TOL "
                f"({qc.RTT_SLOPE_REL_TOL}) and must be flagged; got {flagged!r}"
            )
    # 24 / 30 = 0.80: a 20% gap is INSIDE the tolerance and must stay quiet, because
    # the two regressions run over different tip sets and are not expected to agree
    # exactly. This is the case that would make the check cry wolf on a correct run.
    if "tripwire" in note:
        failures.append(
            f"a 0.80 ratio is inside the {qc.RTT_SLOPE_REL_TOL} tolerance and must "
            f"NOT be flagged; got {note!r}"
        )

    absent = qc.root_to_tip_note(reg, 0, None, "/some/rotated/away.stdout")
    if "unavailable" not in absent or "/some/rotated/away.stdout" not in absent:
        failures.append(
            f"with no reference the note must say so and name the path it looked "
            f"at; got {absent!r}"
        )
    if "ratio" in absent:
        failures.append(
            f"with no reference the note must not imply a comparison was made; "
            f"got {absent!r}"
        )

    # A single-date cloud has no slope; that must be reported, not divided by zero.
    flat = qc.weighted_linregress(np.zeros(4), np.arange(4.0), np.ones(4))
    if flat["slope"] is not None:
        failures.append(f"a single-date cloud must have no slope; got {flat}")
    flat_note = qc.root_to_tip_note(flat, 0, 40.0, "/x")
    if "undefined" not in flat_note or "share one collection date" not in flat_note:
        failures.append(
            f"a single-date cloud must report as undefined AND name that as the "
            f"cause; got {flat_note!r}"
        )

    # An empty cloud is also slopeless, but for a different reason -- a tree whose
    # metadata carries no day-resolution dates at all, which is a normal state and
    # not a degenerate cloud. Reporting the wrong cause sends the reader after the
    # wrong problem.
    empty = qc.weighted_linregress(np.array([]), np.array([]), np.array([]))
    empty_note = qc.root_to_tip_note(empty, 0, 40.0, "/x")
    if "no day-resolution tips" not in empty_note:
        failures.append(
            f"an empty cloud must be reported as having no day-resolution tips, "
            f"not as a single-date cloud; got {empty_note!r}"
        )

    if not failures:
        print("  slope is subs/yr; the reference is compared, flagged, or its "
              "absence stated")
    return failures


def check_fitted_rate_comes_from_the_last_logged_step():
    """The fitted clock must be read from the LAST progress line, not the first.

    This is the whole reason the parse is tested rather than eyeballed. chronumental
    logs ``step:0`` with the optimizer's starting value, which is its own pre-fit
    root-to-tip estimate after one update -- so a first-match parse returns a number
    within a few percent of the pre-fit rate. That reads as a converged clock
    agreeing with its starting point, which is the *opposite* of the finding on every
    real tree (fitted 9.36 against pre-fit 40.36 on ``public_2026-08-06_v11``), and
    nothing downstream could tell the two apart: both are plausible clock rates in
    subs/yr and both plot as a straight line.
    """
    import os
    import tempfile

    qc = _load()
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        good = os.path.join(tmp, "chronumental.stdout")
        with open(good, "w") as f:
            f.write("Loading tree\n")
            f.write(f"{qc.ROOT_TO_TIP_PREFIX} 40.357335447508774\n")
            # Step 0 near the pre-fit rate, the last step far below it: exactly the
            # shape of a real run.
            f.write("step:0\tloss:4e8\tlength_cor:0.9942\tmutation_rate:39.1646\n")
            f.write("step:10\tloss:2e8\tlength_cor:0.9754\tmutation_rate:20.0000\n")
            f.write("step:20\tloss:1e8\tlength_cor:0.3964\tmutation_rate:9.3575\n")
            f.write("Fit completed. Extracting parameters.\n")
            f.write("Writing tree to file\n")

        rate, step = qc.parse_fitted_mutation_rate(good)
        if rate is None or abs(rate - 9.3575) > 1e-9 or step != 20:
            failures.append(
                f"expected the last logged step (9.3575, 20), got ({rate}, {step})"
            )
        # State the trap as an assertion, not only as a comment: a first-match parse
        # returns 39.1646, which is 0.97x the pre-fit rate and would pass any
        # sanity check written in terms of "is this a plausible clock".
        pre_fit = qc.parse_root_to_tip_rate(good)
        if rate is not None and abs(rate / pre_fit - 1) < 0.5:
            failures.append(
                f"parsed rate {rate} sits within 50% of the pre-fit {pre_fit}: this "
                f"fixture is built so the last step is 4.3x below it, so a value "
                f"this close means step 0 was read"
            )

        # A truncated final line -- a run killed mid-write -- must leave the last
        # good pair standing rather than blanking the reference or pairing the new
        # step with the old rate.
        truncated = os.path.join(tmp, "truncated.stdout")
        with open(truncated, "w") as f:
            f.write("step:0\tmutation_rate:39.1646\n")
            f.write("step:10\tmutation_rate:9.3575\n")
            f.write("step:20\tloss:1e8\tmutation_ra")
        rate, step = qc.parse_fitted_mutation_rate(truncated)
        if rate is None or abs(rate - 9.3575) > 1e-9 or step != 10:
            failures.append(
                f"a truncated final line should leave (9.3575, 10) standing, got "
                f"({rate}, {step})"
            )

        # No progress lines at all, and no file at all: both are "no rate", and both
        # must be reportable rather than crash.
        empty = os.path.join(tmp, "empty.stdout")
        with open(empty, "w") as f:
            f.write("Loading tree\nFit completed. Extracting parameters.\n")
        for label, path in [("a log with no step lines", empty),
                            ("a missing log", os.path.join(tmp, "nope.stdout"))]:
            rate, step = qc.parse_fitted_mutation_rate(path)
            if (rate, step) != (None, None):
                failures.append(f"{label} should give (None, None), got ({rate}, {step})")

    if not failures:
        print("  fitted rate is the last logged step; truncation and absence are safe")
    return failures


def check_model_implied_line_anchors_at_the_root_and_is_per_year():
    """The fitted-clock line passes through zero on the root's date, in subs/**year**.

    Two ways to get this wrong, and only one of them is visible on the figure. A
    per-*day* reading of the rate puts the line 365x too steep and shoots it off the
    top of the plot, which anyone would notice. But anchoring it anywhere other than
    the root -- at the cloud's centroid, say, as the pre-fit slope legitimately is --
    keeps it inside the axes and looking reasonable while destroying the only thing
    it is for: the vertical gap to the cloud, which is the duration inflation
    measured in mutations.
    """
    qc = _load()
    failures = []
    root_day = 18_000.0  # arbitrary epoch day for the root
    rate = 12.0

    at_root = qc.model_implied_divergence([root_day], root_day, rate)
    if at_root is None or abs(float(at_root[0])) > 1e-12:
        failures.append(f"divergence at the root must be exactly 0, got {at_root}")

    one_year = qc.model_implied_divergence(
        [root_day + qc.DAYS_PER_YEAR], root_day, rate
    )
    if one_year is None or abs(float(one_year[0]) - rate) > 1e-9:
        failures.append(
            f"one year past the root must be exactly one rate ({rate}), got {one_year}"
        )

    # Ten years out, the difference between per-year and per-day is not subtle;
    # asserted anyway because it is the cheap half of the pair above.
    ten_years = qc.model_implied_divergence(
        [root_day + 10 * qc.DAYS_PER_YEAR], root_day, rate
    )
    if ten_years is None or abs(float(ten_years[0]) - 10 * rate) > 1e-9:
        failures.append(f"ten years past the root should be {10 * rate}, got {ten_years}")

    # A date BEFORE the root gives negative predicted divergence rather than being
    # clamped: clamping would hide a root dated after its own samples, which is
    # exactly what check 7 and check 8 exist to catch.
    before = qc.model_implied_divergence([root_day - 365.25], root_day, rate)
    if before is None or float(before[0]) >= 0:
        failures.append(f"a pre-root date should predict negative divergence, got {before}")

    # Either input missing means "no line", which the caller states on the figure.
    for label, args in [
        ("no fitted rate", ([root_day], root_day, None)),
        ("no root anchor", ([root_day], None, rate)),
    ]:
        if qc.model_implied_divergence(*args) is not None:
            failures.append(f"{label} should give None, not a line")

    if not failures:
        print("  model line is anchored at the root, in subs/yr, and unclamped")
    return failures


def check_per_bin_rate_is_exact_and_reconciles():
    """Per-bin rates use real per-branch counts and sum back to the aggregate.

    The fixture is built so a **midpoint** reading of the ranged bins gives a visibly
    different answer: the two branches in the ``21-50`` bin carry 21 and 50
    mutations, whose midpoint (35) is wrong for both. That matters more than the
    small error it causes, because the ranged bins are exactly where the punctuated
    jumps live -- the branches this decomposition exists to characterize -- so a
    midpoint is least defensible precisely where it is most tempting.

    The reconciliation is the other half: the bins must partition the branches, so
    summing ``mutations`` and ``branch_years`` across them has to reproduce
    :func:`implied_clock_summary` exactly. A binning bug that dropped or
    double-counted a bin would otherwise show up only as a slightly odd rate in one
    row, which is unfalsifiable by eye.
    """
    qc = _load()
    failures = []
    year = qc.DAYS_PER_YEAR

    # 0-mutation branches holding real time, ordinary branches, and two jumps.
    dt = np.array([year, year, year, 2 * year, year, year, np.nan], dtype="float64")
    muts = np.array([0, 0, 1, 4, 21, 50, 7], dtype="float64")
    by_bin = qc.implied_rate_by_bin(dt, muts)

    if "0" not in by_bin:
        failures.append(
            "the zero-mutation bin must be reported, not dropped as undefined: "
            "0 subs/yr over a large share of the tree's time is the finding"
        )
    elif by_bin["0"]["subs_per_year"] != 0.0:
        failures.append(f"zero-mutation bin should be 0.0 subs/yr, got {by_bin['0']}")

    # 21 + 50 = 71 mutations over 2 branch-years. A midpoint of 35 would give 70/2.
    ranged = by_bin.get("21-50")
    if ranged is None or abs(ranged["mutations"] - 71) > 1e-9:
        failures.append(f"the 21-50 bin must sum real counts (71), got {ranged}")
    elif abs(ranged["subs_per_year"] - 35.5) > 1e-9:
        failures.append(f"21-50 rate should be 71/2 = 35.5, got {ranged}")

    # The NaN-duration branch carries 7 mutations. Crediting them without their time
    # would inflate the 6-10 bin (and the aggregate) out of nothing.
    if "6-10" in by_bin:
        failures.append(
            f"a branch with an undefined Delta-t must be excluded from both sums, "
            f"but its bin appeared: {by_bin['6-10']}"
        )

    aggregate = qc.implied_clock_summary(dt, muts)
    tot_m = sum(e["mutations"] for e in by_bin.values())
    tot_y = sum(e["branch_years"] for e in by_bin.values())
    if abs(tot_m - aggregate["total_mutations"]) > 1e-6:
        failures.append(
            f"per-bin mutations sum to {tot_m} but the aggregate says "
            f"{aggregate['total_mutations']}: the bins are not a partition"
        )
    if abs(tot_y - aggregate["branch_years"]) > 1e-6:
        failures.append(
            f"per-bin branch-years sum to {tot_y} but the aggregate says "
            f"{aggregate['branch_years']}: the bins are not a partition"
        )

    if qc.implied_rate_by_bin([], []) != {}:
        failures.append("no branches should give an empty mapping, not a raise")
    try:
        qc.implied_rate_by_bin([1.0, 2.0], [1.0])
        failures.append("mismatched array lengths should raise")
    except ValueError:
        pass

    if not failures:
        print("  per-bin rates are exact, keep the zero bin, and sum to the aggregate")
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
    ("aggregate rate sees the intercept", check_aggregate_sees_what_the_slope_cannot),
    ("aggregate is exact, not binned", check_aggregate_is_exact_not_binned),
    ("NaN branches and zero-span trees", check_undefined_branches_and_zero_span),
    ("root-to-tip absence is reported", check_root_to_tip_reference_is_never_silently_dropped),
    ("divergence accumulates, root excluded", check_rtt_divergence_accumulates_and_excludes_the_root),
    ("the binned scatter is lossless", check_rtt_histogram_is_lossless),
    ("slope units and reference handling", check_rtt_slope_is_per_year_and_reference_is_never_dropped),
    ("fitted rate is the last step", check_fitted_rate_comes_from_the_last_logged_step),
    ("model line anchors at the root", check_model_implied_line_anchors_at_the_root_and_is_per_year),
    ("per-bin rate is exact and reconciles", check_per_bin_rate_is_exact_and_reconciles),
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
