"""Internal-node date sanity checks for chronumental output (#139).

The existing ``chronumental_qc`` checks are node-id correspondence, monotonicity,
coverage, and **tip residuals**. That last one is the only check that compares an
inferred date against evidence, and it is aimed at the one place there is almost
nothing to infer: tip collection dates are the priors chronumental is *fit to*, so
reproducing them is close to circular. On ``gisaid-and-public_2026-04-22_v9`` it
reports ``median|res|=1.0d, p95|res|=1.0d`` across 16,979,190 tips against a 30-day
gate.

**The internal dates are what the pipeline consumes.** ``modules/windowing.py`` keys
every node id into ``self.times`` (``:350``), asserts all nodes are dated
(``:351-355``), and buckets each branch into a time window by its endpoints -- so
internal dates alone decide which window a branch's mutations land in, straight into
``windows/subtree_counts.csv`` and ``windows/fitness_trajectory.csv``. Nothing
checked them. It is a live risk rather than a theoretical one:
``chronumental_kwargs`` is ``--steps 1000`` and that config's own comment records
production step tuning as deferred, and under-convergence surfaces first in exactly
the parameters with the weakest priors.

Four checks, all **report-only**: they emit rows into ``node_dates_qc.tsv`` and
binned artifacts, and gate nothing. Gating on a distribution nobody has looked at
would fail production runs uninterpretably -- the same measure-then-decide shape as
#135's audit. Thresholds are chosen in a follow-up, once the distributions exist.

**A. Δt versus branch mutation count.** Binned; the clock should make these roughly
proportional. Order-of-magnitude expectation ~15 days/mutation
(8e-4 subs/site/yr × 29,903 nt ≈ 23.9 subs/yr), and no tighter, because the MAT is
masked and filtered so the realized rate on *counted* mutations is lower. Note this
check is itself partly circular -- chronumental fits a clock, so agreement shows the
fit **converged to something clock-like**, which is precisely the ``--steps 1000``
risk, not that the clock is right.

**B. Zero-mutation leaf/parent pairs.** A leaf identical to its parent (no mutations
on the connecting branch) should have a similar date. Two refinements over the naive
form: the predicted-date version is **one-sided** -- the existing monotonicity check
already forces ``date(child) >= date(parent)``, so it is a right tail rather than a
distribution centred on zero -- and identical leaves genuinely carry different
collection dates, so per-leaf equality is the wrong bar. The informative comparison
is the parent's predicted date against the **earliest observed date** among its
zero-mutation leaf children, which is the only one of the four that is not
circular.

**C. Root date plausibility.** v9 gives 2019-12-09, so this passes today. It is here
because it is the cheapest possible tripwire for a mis-scaled time axis.

**D. Internal date versus minimum descendant observed date.** An ancestor cannot
postdate the earliest sample beneath it. Computed by propagating the per-node
minimum observed descendant date in one reverse pass over the depth-first ordering.
Note this compares an *inferred* date against an *observed* one, so a violation is a
three-way inconsistency between the inferred date, the topology and a collection
date -- and it is **bounded by the tip residuals**, because monotonicity forces
``pred(ancestor) <= pred(tip)``. It quantifies downstream impact (how many internal
nodes, hence branches, land in the wrong window) but it is not a clean test of the
dating algorithm.

**E. Internal node to its nearest descendant, inferred versus inferred.** The
algorithm test D is not: both endpoints are chronumental outputs, so this asks how
much time the fit puts between an internal node and its closest descendant given the
topology it was handed. Because monotonicity makes predicted dates nondecreasing
along every path, the nearest descendant is always an immediate *child*, so this is
:func:`min_child_delta` -- one scatter-min over the edges already built for check A,
with no propagation pass and no dependence on node ordering.

**Everything is integer days in histograms, never per-branch rows.** The rule
already costs 349 s at 38.2 GB max RSS over 20,961,614 nodes, and
``node_dates_residuals.csv`` is already 1.6 GB. Δt in whole days is an integer, so a
histogram over it gives **exact** medians and percentiles at a few tens of kB --
there is no accuracy/memory tradeoff being made here.
"""

import os

import numpy as np

#: Mutation-count bins for check A. Upper bounds inclusive; the last is open.
MUT_COUNT_BINS = [0, 1, 2, 3, 4, 5, 10, 20, 50]

#: The leading bins that hold exactly one mutation count, so consecutive medians
#: differ by exactly one mutation and their difference estimates the clock slope.
#: The ranged bins above (``6-10``, ``11-20``, ``21-50``, ``51+``) have no
#: well-defined increment, and on the production tree they also saturate -- the ~650
#: branches with >20 mutations sit far below any linear extrapolation.
SINGLE_COUNT_BINS = [b for b in MUT_COUNT_BINS if b == 0 or b - 1 in MUT_COUNT_BINS]

#: Rough clock expectation, days per mutation. 8e-4 subs/site/yr x 29903 nt
#: ~= 23.9 subs/yr -> ~15.3 days/substitution. An order-of-magnitude reference
#: only: the MAT is masked and branch-filtered, so the realized rate on counted
#: mutations is lower and the observed slope is expected to exceed this.
DAYS_PER_MUTATION_NOMINAL = 365.25 / (8e-4 * 29903)

#: Plausible window for the SARS-CoV-2 MRCA. Deliberately wide -- this is a
#: tripwire for a mis-scaled time axis, not an estimate.
ROOT_DATE_PLAUSIBLE = ("2019-10-01", "2020-02-29")

#: Days per year for converting summed branch Delta-t into branch-years.
DAYS_PER_YEAR = 365.25

#: The line chronumental prints before fitting, carrying its own root-to-tip
#: regression estimate of the clock. That number is the ideal reference for
#: :func:`implied_clock_summary`: chronumental computes it from the same tree and the
#: same dates, before any duration is optimized, so comparing the two isolates what
#: the fit did to the durations from every other difference between this tree and a
#: nominal expectation.
ROOT_TO_TIP_PREFIX = "Root to tip regression: got rate of:"


def implied_clock_summary(branch_dt_days, branch_mutations):
    """The clock rate the timetree implies about itself: sum(mut) / sum(branch-years).

    Check A reports the **incremental** slope between adjacent single-count bins,
    which is deliberately intercept-free (see the note at its call site). That makes
    it blind to the failure mode #145 confirmed: chronumental's durations are free
    parameters under a flat prior while every day-resolved leaf date is asserted to
    sd 0.3 days, so the fit can buy leaf-date accuracy by inflating durations -- and
    the inflation lands almost entirely in the intercept. On
    ``public_2026-08-06_v11``, 56.7% of branches carry zero mutations yet hold 45.9%
    of the tree's total branch-time. Differencing adjacent bin medians cancels
    exactly that term: the incremental slope read 11 d/mutation against a ~15.3 d
    nominal and looked healthy while the aggregate rate sat 4.3x below chronumental's
    own root-to-tip estimate.

    So this is the aggregate the slope cannot see. Two properties matter:

    * **It is exact, not binned.** The sums run over the per-branch arrays, not over
      :data:`MUT_COUNT_BINS` midpoints. Reconstructing it from
      ``node_dates_branch_dt_summary.csv`` requires guessing representative counts
      for the ranged bins (``6-10``, ``11-20``, ``21-50``, ``51+``), which moves the
      answer by about +/-7% -- enough to argue about, and unnecessary here since the
      arrays are already in hand.
    * **The zero-mutation share is reported as a share of branch-*time*, not of
      branches.** A tree can legitimately have most branches carry no mutations
      (dense sampling); what indicts the fit is those branches holding a large
      fraction of the total elapsed time, which is duration with no mutational
      support.

    ``branch_dt_days`` may carry NaN where Delta-t is undefined; those branches leave
    both sums, and ``n_branches`` counts what was actually used.
    """
    dt = np.asarray(branch_dt_days, dtype=float)
    muts = np.asarray(branch_mutations, dtype=float)
    if dt.shape != muts.shape:
        raise ValueError(
            f"branch Delta-t and mutation counts must be parallel arrays; got "
            f"{dt.shape} and {muts.shape}"
        )
    defined = ~np.isnan(dt)
    dt = dt[defined]
    muts = muts[defined]
    total_days = float(dt.sum())
    total_muts = float(muts.sum())
    branch_years = total_days / DAYS_PER_YEAR
    zero = muts == 0
    return {
        "n_branches": int(dt.size),
        "total_mutations": total_muts,
        "total_branch_days": total_days,
        "branch_years": branch_years,
        # None rather than inf or NaN: a zero-span tree is a different failure (a
        # collapsed time axis), and inf would render as a number in the report.
        "subs_per_year": (total_muts / branch_years) if branch_years > 0 else None,
        "n_zero_mut_branches": int(zero.sum()),
        "frac_branches_zero_mut": (float(zero.mean()) if dt.size else None),
        "frac_branch_time_zero_mut": (
            float(dt[zero].sum() / total_days) if total_days > 0 else None
        ),
    }


def implied_clock_note(summary, root_to_tip_rate, stdout_path):
    """The report fragment for :func:`implied_clock_summary`.

    Lives here rather than at the call site so the one property that makes this
    diagnostic trustworthy is unit-testable: **the root-to-tip comparison is either
    made or its absence is stated**, never quietly dropped. A detail line that just
    omits the ratio when the log is missing reads exactly like one where the ratio
    was fine, which is the failure shape this whole check exists to end.
    """
    nominal = DAYS_PER_YEAR / DAYS_PER_MUTATION_NOMINAL
    rate = summary["subs_per_year"]
    if rate is None:
        return "; aggregate rate undefined: the dated branches span zero total time"
    note = (
        f"; aggregate {summary['total_mutations']:.0f} mutations over "
        f"{summary['branch_years']:.0f} branch-years = {rate:.2f} subs/yr "
        f"({rate / nominal:.2f}x the nominal {nominal:.1f})"
    )
    if root_to_tip_rate:
        note += (
            f", {rate / root_to_tip_rate:.2f}x chronumental's own pre-fit "
            f"root-to-tip {root_to_tip_rate:.2f} subs/yr"
        )
    else:
        note += (
            f", root-to-tip reference unavailable (no parsable "
            f"'{ROOT_TO_TIP_PREFIX}' line at {stdout_path})"
        )
    if summary["frac_branch_time_zero_mut"] is not None:
        note += (
            f"; zero-mutation branches hold "
            f"{100 * summary['frac_branch_time_zero_mut']:.1f}% of total branch-time "
            f"({100 * summary['frac_branches_zero_mut']:.1f}% of branches)"
        )
    return note


def parse_root_to_tip_rate(stdout_path):
    """chronumental's own pre-fit clock estimate from its stdout, or ``None``.

    ``None`` means "no reference available", which the caller must **report** rather
    than pass over: the log is a Snakemake ``log:`` entry rather than a declared
    output (it is not in ``OUTPUT_SCHEMA.md``), so it can legitimately be absent on
    a re-run or an archived tree, and a detail line that silently drops the
    comparison is indistinguishable from one where the comparison passed. Same shape
    as the era skip in ``Snakefile: eras_present``: absent is fine, absent-and-quiet
    is not.
    """
    if not stdout_path or not os.path.exists(stdout_path):
        return None
    with open(stdout_path) as f:
        for line in f:
            if line.startswith(ROOT_TO_TIP_PREFIX):
                return float(line[len(ROOT_TO_TIP_PREFIX):].strip())
    return None


def mutation_count_bin(n_mutations):
    """Label for ``n_mutations`` under :data:`MUT_COUNT_BINS`."""
    for i, upper in enumerate(MUT_COUNT_BINS):
        if n_mutations <= upper:
            lower = MUT_COUNT_BINS[i - 1] + 1 if i else 0
            return f"{lower}" if lower == upper else f"{lower}-{upper}"
    return f"{MUT_COUNT_BINS[-1] + 1}+"


class DayHistogram:
    """Counts of integer day values, for exact order statistics at fixed memory.

    Values are whole days, so nothing is approximated: the median and percentiles
    below are the same numbers a sorted array of every observation would give.
    Negative values are kept -- for checks B and D a negative value is the
    interesting case and must not be clipped away.
    """

    def __init__(self):
        self._counts = {}
        self.n = 0
        self.total = 0

    def add(self, days):
        self._counts[days] = self._counts.get(days, 0) + 1
        self.n += 1
        self.total += days

    def update(self, values):
        for v in values:
            self.add(int(v))

    def update_from_array(self, values):
        """Bulk-load a NumPy array of day values.

        ``np.unique`` collapses 21 M branches to a few thousand distinct day values
        in C, which is what keeps the whole-tree checks off the Python-level
        per-node loop they started as (that cost minutes on the production tree).
        NaNs are dropped; the caller is responsible for having masked out nodes with
        no date.
        """
        values = np.asarray(values)
        values = values[~np.isnan(values)]
        if not values.size:
            return
        uniques, counts = np.unique(values.astype("int64"), return_counts=True)
        for value, count in zip(uniques.tolist(), counts.tolist()):
            self._counts[value] = self._counts.get(value, 0) + count
            self.n += count
            self.total += value * count

    def merge(self, other):
        """Fold ``other``'s counts into this histogram.

        Used to build the marginal "all branches" row from the per-bin histograms.
        Doing it this way rather than re-scanning the tree means the marginal and the
        bins cannot disagree, and it keeps the marginal's order statistics exact --
        averaging per-bin medians would not be the median of the pooled data.
        """
        for value, count in other._counts.items():
            self._counts[value] = self._counts.get(value, 0) + count
            self.n += count
            self.total += value * count
        return self

    def quantile(self, q):
        """Lower-interpolated quantile, matching ``numpy.percentile(..., 'lower')``."""
        if not self.n:
            return None
        target = q * (self.n - 1)
        seen = 0
        for value in sorted(self._counts):
            seen += self._counts[value]
            if seen > target:
                return value
        return max(self._counts)

    def stats(self):
        if not self.n:
            return {
                "n": 0,
                "min": None,
                "p05": None,
                "p25": None,
                "median": None,
                "p75": None,
                "iqr": None,
                "p95": None,
                "max": None,
                "mean": None,
                "frac_negative": None,
                "frac_positive": None,
            }
        negatives = sum(c for v, c in self._counts.items() if v < 0)
        positives = sum(c for v, c in self._counts.items() if v > 0)
        p25 = self.quantile(0.25)
        p75 = self.quantile(0.75)
        return {
            "n": self.n,
            "min": min(self._counts),
            "p05": self.quantile(0.05),
            "p25": p25,
            "median": self.quantile(0.5),
            "p75": p75,
            "iqr": p75 - p25,
            "p95": self.quantile(0.95),
            "max": max(self._counts),
            "mean": self.total / self.n,
            "frac_negative": negatives / self.n,
            "frac_positive": positives / self.n,
        }

    def as_rows(self):
        """``(days, count)`` pairs, sorted -- the binned artifact."""
        return sorted(self._counts.items())


def propagate_min_descendant(parent_index, leaf_values, n_nodes):
    """Minimum ``leaf_values`` over each node's descendant leaves, including itself.

    ``parent_index[i]`` is the index of node ``i``'s parent, or ``-1`` for the root.
    Nodes must be ordered so that a parent always precedes its children (which is
    what ``depth_first_expansion`` gives), so iterating **backwards** visits every
    child before its parent -- a valid post-order accumulation without recursion or
    an explicit stack, which matters at 21 M nodes.

    ``leaf_values`` is a float array with ``nan`` where a node has no observed date
    of its own. Returns a float array with ``nan`` where no descendant is observed.

    The loop runs over Python **lists** with ``inf`` as the missing sentinel rather
    than over NumPy arrays with ``nan``. At 21 M nodes the NumPy version spent
    minutes building a scalar object per element access; lists plus a sentinel make
    the body two indexed reads and a comparison. It cannot be vectorized outright
    because each parent depends on its children, so the sequential pass stays --
    only its constant factor changes.
    """
    if len(parent_index) != n_nodes or len(leaf_values) != n_nodes:
        raise ValueError(
            f"parent_index ({len(parent_index)}) and leaf_values "
            f"({len(leaf_values)}) must both have length n_nodes ({n_nodes})"
        )
    parents = np.asarray(parent_index, dtype="int64")
    # The precondition is checked, not merely documented. Violating it does not
    # raise on its own: the reverse pass simply stops propagating past the
    # offending edge and returns a smaller-but-plausible minimum (a nan where a
    # real date belonged), at exit 0 -- the #129/#130/#133 silent-drop shape. The
    # loop below cannot be vectorized because each parent depends on its
    # children, but this guard can: one O(n) comparison, which is free next to
    # the pass it protects.
    positions = np.arange(n_nodes)
    has_parent = parents >= 0
    out_of_order = has_parent & (parents >= positions)
    if out_of_order.any():
        bad = np.flatnonzero(out_of_order)[:5].tolist()
        raise ValueError(
            f"parent_index must list every parent before its children (a "
            f"depth-first pre-order); {int(out_of_order.sum())} node(s) violate "
            f"this, first at indices {bad} with parents "
            f"{parents[bad].tolist()}. The single reverse pass below relies on "
            f"this ordering, and silently under-propagates without it rather "
            f"than failing -- so it is checked instead of assumed."
        )
    parents = parents.tolist()
    values = np.asarray(leaf_values, dtype="float64")
    best = np.where(np.isnan(values), np.inf, values).tolist()
    for i in range(n_nodes - 1, 0, -1):
        parent = parents[i]
        if parent < 0:
            continue
        value = best[i]
        if value < best[parent]:
            best[parent] = value
    out = np.array(best, dtype="float64")
    out[np.isinf(out)] = np.nan
    return out


def min_child_delta(parent_index, values, n_nodes):
    """Minimum ``values[child] - values[parent]`` over each node's own children.

    Returns a float array with ``nan`` for nodes that have no child with a defined
    delta (all leaves, plus any parent whose children are all undated).

    **Why this is the Delta-t to the nearest descendant, not merely to the nearest
    child.** When ``values`` is nondecreasing along every parent->child edge -- which
    is what the ``time_monotonicity`` gate enforces for predicted dates -- then for
    any deeper descendant ``d`` reached through child ``c``, ``values[d] >=
    values[c]``. So the minimum over the whole descendant set is always attained at
    an immediate child, and a one-pass ``np.minimum.at`` over the edges replaces a
    propagation over the tree. Unlike :func:`propagate_min_descendant` this needs no
    particular node ordering, because scatter-min does not care which edge it sees
    first.

    That equivalence is a *theorem about monotone values*, not about this function:
    fed values that decrease somewhere, it still returns the minimum over children,
    which is then **not** the minimum over descendants.
    ``tests/test_node_date_qc.py`` asserts both halves -- agreement with
    :func:`propagate_min_descendant` on a monotone tree, and disagreement on a
    non-monotone one -- so the precondition cannot quietly stop being load-bearing.
    """
    parents = np.asarray(parent_index, dtype="int64")
    vals = np.asarray(values, dtype="float64")
    if len(parents) != n_nodes or len(vals) != n_nodes:
        raise ValueError(
            f"parent_index ({len(parents)}) and values ({len(vals)}) must both "
            f"have length n_nodes ({n_nodes})"
        )
    out = np.full(n_nodes, np.inf, dtype="float64")
    child = np.flatnonzero(parents >= 0)
    parent_of = parents[child]
    delta = vals[child] - vals[parent_of]
    defined = ~np.isnan(delta)
    np.minimum.at(out, parent_of[defined], delta[defined])
    out[np.isinf(out)] = np.nan
    return out


def format_stats(stats, unit="d", violation=None):
    """Compact one-line rendering for the report's ``detail`` column.

    ``violation`` names the sign that indicates a problem -- ``"negative"`` or
    ``"positive"`` -- and makes the report state the violating fraction as such.
    It is required rather than optional-by-default because checks 6 and 8 use
    **opposite** sign conventions (``min_observed - predicted`` versus
    ``predicted - min_observed``), and printing only ``frac<0`` for both rendered
    one as 23.6% *bad* and the other as 96.1% *good* in the same column. A reader
    comparing those two numbers draws the opposite of the right conclusion.

    Both fractions are printed regardless, so the raw quantity stays available;
    ``violation`` only adds the labelled one. ``frac_negative`` +
    ``frac_positive`` does not sum to 1 -- exact zeros are in neither.
    """
    if not stats["n"]:
        return "no observations"
    text = (
        f"n={stats['n']}, median={stats['median']}{unit}, "
        f"mean={stats['mean']:.1f}{unit}, p05={stats['p05']}{unit}, "
        f"p95={stats['p95']}{unit}, min={stats['min']}{unit}, "
        f"max={stats['max']}{unit}, frac<0={stats['frac_negative']:.4f}, "
        f"frac>0={stats['frac_positive']:.4f}"
    )
    if violation is None:
        return text
    if violation not in ("negative", "positive"):
        raise ValueError(
            f"violation must be 'negative', 'positive' or None, not {violation!r}"
        )
    key = "frac_negative" if violation == "negative" else "frac_positive"
    sign = "<0" if violation == "negative" else ">0"
    count = round(stats[key] * stats["n"])
    return f"{text}; VIOLATING ({sign}) = {stats[key]:.4f} ({count} of {stats['n']})"
