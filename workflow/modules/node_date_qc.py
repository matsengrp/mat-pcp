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


def implied_rate_by_bin(branch_dt_days, branch_mutations):
    """Per-bin ``sum(mutations) / sum(branch-years)``, keyed by bin label.

    :func:`implied_clock_summary` gives one aggregate for the whole tree; this
    decomposes it by branch mutation count, which is what distinguishes the two
    stories a low aggregate could tell. If the fit were tracking gradual evolution
    and letting saltational branches sit as outliers, the ordinary low-count bins
    would come out near the clock and the high-count bins far below it. What the v11
    trees actually show is the reverse: bins 1-5 climb 13.5 -> 25.1 subs/yr, arriving
    at the ~23.9 nominal, the >20-mutation bins run *above* it (57.9 and 72.0, i.e. the fit
    compresses the jumps rather than stretching them), and the entire deficit is the
    zero-mutation bin -- 5.77 M branches averaging 22.5 days, contributing 46% of the
    tree's time and no mutations at all. Under ``--variance_dates 30`` that same bin
    averages 0.24 days and the aggregate moves 9.47 -> 27.25.

    Two things about the arithmetic matter.

    **Exact counts, not bin midpoints.** ``mutations`` sums the real per-branch
    counts, so the ranged bins (``6-10``, ``21-50``, ``51+``) are as exact as the
    single-count ones. Reconstructing this from
    ``node_dates_branch_dt_summary.csv``'s ``mean`` column requires inventing a
    representative count per ranged bin, which moved the aggregate 9.47 -> 9.59 in
    the analysis that prompted this function -- small, but it is the ranged bins
    that carry the punctuation, so that is exactly where a midpoint is least
    defensible.

    **The zero-mutation bin is reported, not dropped.** Its rate is 0 by
    construction, and a caller tempted to omit it as "undefined" would delete the
    finding: 0 subs/yr over 46% of the tree's elapsed time is the whole diagnosis.
    Branches with an undefined ``Delta-t`` are excluded from both sums, since a NaN
    duration cannot be credited with its mutations either.
    """
    dt = np.asarray(branch_dt_days, dtype="float64")
    muts = np.asarray(branch_mutations, dtype="float64")
    if len(dt) != len(muts):
        raise ValueError(
            f"branch arrays disagree: {len(dt)} durations vs {len(muts)} counts"
        )
    finite = np.isfinite(dt) & np.isfinite(muts)
    dt, muts = dt[finite], muts[finite]
    out = {}
    edges = np.asarray(MUT_COUNT_BINS)
    bin_pos = np.searchsorted(edges, muts, side="left")
    for pos in np.unique(bin_pos):
        mask = bin_pos == pos
        count = int(mask.sum())
        label = mutation_count_bin(
            int(edges[pos]) if pos < len(edges) else int(edges[-1]) + 1
        )
        years = float(dt[mask].sum()) / DAYS_PER_YEAR
        mutations = float(muts[mask].sum())
        out[label] = {
            "n_branches": count,
            "mutations": mutations,
            "branch_years": years,
            "subs_per_year": (mutations / years) if years > 0 else None,
        }
    return out


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


#: chronumental logs one progress line per ``--log_every`` steps, each carrying the
#: current value of every fitted parameter as whitespace-separated ``name:value``
#: fields. The *fitted* clock is the ``mutation_rate`` field of the **last** such
#: line: the run ends with "Fit completed. Extracting parameters." and no
#: final-parameter line, so the last logged step is the closest thing the log has to
#: the rate the run actually wrote out.
FITTED_RATE_STEP_PREFIX = "step:"
FITTED_RATE_FIELD = "mutation_rate"


def parse_fitted_mutation_rate(stdout_path):
    """chronumental's *fitted* clock, as ``(subs_per_year, step)``.

    ``(None, None)`` means "no fitted rate available", which the caller must
    **state** rather than pass over, for the same reason as
    :func:`parse_root_to_tip_rate`.

    **The last progress line, not the first.** Step 0 reports the optimizer's
    starting value, which is the pre-fit root-to-tip estimate after one update --
    39.16 against a pre-fit 40.36 on ``public_2026-08-06_v11``. A first-match parse
    therefore returns a number within 3% of the pre-fit rate: it looks like a
    converged clock agreeing with its own starting point, and it is the exact
    opposite of the finding, since the rate that run actually fitted is 9.36, 4.3x
    lower. Nothing downstream can tell the two apart, so the last-line semantics are
    unit-tested rather than trusted to the loop reading correctly.

    **The step is returned, not discarded.** chronumental logs every
    ``--log_every`` steps, so the last *logged* step need not be the last step
    taken, and the rate can still be moving there: over the final third of the #146
    ``--variance_dates 30`` run it went 21.9 -> 27.2. Printing the rate without the
    step implies a convergence the log does not show.
    """
    if not stdout_path or not os.path.exists(stdout_path):
        return None, None
    rate, step = None, None
    with open(stdout_path) as f:
        for line in f:
            if not line.startswith(FITTED_RATE_STEP_PREFIX):
                continue
            fields = {}
            for token in line.split():
                name, sep, value = token.partition(":")
                if sep:
                    fields[name] = value
            if FITTED_RATE_FIELD not in fields:
                continue
            try:
                # Assigned together, and only on success, so a truncated final line
                # (a killed run mid-write) leaves the last good pair standing rather
                # than blanking the reference or pairing a rate with a stale step.
                rate, step = float(fields[FITTED_RATE_FIELD]), int(fields["step"])
            except ValueError:
                continue
    return rate, step


def model_implied_divergence(day_ordinals, root_day_ordinal, fitted_rate):
    """Root-to-tip divergence the fitted clock predicts at each date.

    This is the model's own prediction, not a slope drawn near the cloud for
    comparison, and the reason is worth stating because it is what makes the line
    legitimate. chronumental parameterizes each branch's *duration*, and a node's
    date is the root date plus the cumulative sum along its path -- so the durations
    on any root-to-tip path sum to ``tip_date - root_date`` by construction. The
    per-branch Poisson expectation ``rate * t`` therefore telescopes: expected
    root-to-tip divergence is ``rate * (tip_date - root_date)``, exactly. Tip fitted
    dates equal their collection dates to a median 0.01 d on the shipped prior, so
    evaluating that at the *collection* date is the same statement.

    Which means the gap between this line and the cloud is not a plotting artifact
    to be tuned away -- it is the #145 finding drawn to scale. On
    ``public_2026-08-06_v11`` the fitted clock predicts 62 mutations at the last
    collection date, where the cloud's own regression reads 247: the fit reproduces
    every leaf date
    while under-predicting every leaf's divergence by ~4x, having parked the missing
    time on interior branches carrying no mutations.

    Anchoring at zero divergence on the root's date is exact for the same reason
    :func:`root_to_tip_divergence` excludes the root's own mutation list: on a MAT
    those are the root's difference from the *reference*, not evolution along a
    branch, so divergence at the root is 0 and the line passes through
    ``(root_date, 0)`` rather than through a fitted intercept.
    """
    if fitted_rate is None or root_day_ordinal is None:
        return None
    days = np.asarray(day_ordinals, dtype="float64")
    return fitted_rate * (days - float(root_day_ordinal)) / DAYS_PER_YEAR


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
    # The precondition is checked, not merely documented. Violating it does not
    # raise on its own: the reverse pass simply stops propagating past the
    # offending edge and returns a smaller-but-plausible minimum (a nan where a
    # real date belonged), at exit 0 -- the #129/#130/#133 silent-drop shape. The
    # loop below cannot be vectorized because each parent depends on its
    # children, but this guard can: one O(n) comparison, which is free next to
    # the pass it protects. Shared with root_to_tip_divergence's forward pass,
    # which needs the identical ordering for the identical reason.
    parents = require_parent_before_child(parent_index, n_nodes).tolist()
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


# =============================================================================
# Root-to-tip divergence (#146): the reference itself, as a picture.
#
# :func:`implied_clock_summary` compares the timetree's own implied rate against
# chronumental's pre-fit root-to-tip slope, and that slope is the whole reason the
# comparison means anything (see :data:`ROOT_TO_TIP_PREFIX`). But chronumental
# prints one number and discards the cloud it fitted, so nothing showed *whether the
# regression itself is any good*: a rate of 40.36 subs/yr summarizes a tight line
# and a fan of recombinants and metadata errors equally well, and only one of those
# is a reference worth dividing by. These primitives rebuild that cloud.
# =============================================================================

#: Relative gap between the root-to-tip slope recomputed here and the one
#: chronumental logs, past which the report calls the two inconsistent. Wide on
#: purpose: the two regressions are over **different tip sets** (see
#: :func:`root_to_tip_histogram`), so they are not expected to agree exactly, and a
#: tight tolerance here would fire on a correct pipeline. It is a tripwire for a
#: gross disagreement -- a units error, a divergence array off by a tree level, a
#: date axis in the wrong epoch -- not an equivalence test.
RTT_SLOPE_REL_TOL = 0.25


def require_parent_before_child(parent_index, n_nodes):
    """Validate the pre-order precondition shared by the single-pass accumulators.

    ``parent_index[i]`` is node ``i``'s parent, or ``-1`` for a root. Returns the
    parents as an ``int64`` array.

    Both the reverse pass in :func:`propagate_min_descendant` and the forward pass
    in :func:`root_to_tip_divergence` are correct only if every parent precedes its
    children, and *neither raises on its own when that is violated* -- they return a
    plausible, quietly truncated answer. Hence one guard, called by both, rather
    than two copies that could drift apart.
    """
    parents = np.asarray(parent_index, dtype="int64")
    if len(parents) != n_nodes:
        raise ValueError(
            f"parent_index ({len(parents)}) must have length n_nodes ({n_nodes})"
        )
    positions = np.arange(n_nodes)
    out_of_order = (parents >= 0) & (parents >= positions)
    if out_of_order.any():
        bad = np.flatnonzero(out_of_order)[:5].tolist()
        raise ValueError(
            f"parent_index must list every parent before its children (a "
            f"depth-first pre-order); {int(out_of_order.sum())} node(s) violate "
            f"this, first at indices {bad} with parents "
            f"{parents[bad].tolist()}. The single-pass accumulators rely on this "
            f"ordering, and silently under-propagate without it rather than "
            f"failing -- so it is checked instead of assumed."
        )
    return parents


def root_to_tip_divergence(parent_index, branch_mutations, n_nodes):
    """Mutations on the path from the root to every node, inclusive of its own branch.

    ``branch_mutations[i]`` is the count on the branch *entering* node ``i``. A
    parentless node has no entering branch, so its divergence is **0 regardless of
    ``branch_mutations``** there: on a MAT the root's mutation list is its
    difference from the reference, which is not elapsed evolution within the tree
    and would otherwise offset every tip in the plot by a constant.

    One forward pass, which is the mirror of :func:`propagate_min_descendant`'s
    reverse pass and carries the same ordering precondition (checked, via
    :func:`require_parent_before_child`) and the same list-not-array loop for the
    same measured reason: at 21 M nodes NumPy scalar boxing dominates.
    """
    muts = np.asarray(branch_mutations)
    if len(muts) != n_nodes:
        raise ValueError(
            f"branch_mutations ({len(muts)}) must have length n_nodes ({n_nodes})"
        )
    parents = require_parent_before_child(parent_index, n_nodes).tolist()
    counts = muts.astype("int64").tolist()
    divergence = [0] * n_nodes
    for i in range(n_nodes):
        parent = parents[i]
        if parent >= 0:
            divergence[i] = divergence[parent] + counts[i]
    return np.array(divergence, dtype="int64")


def root_to_tip_histogram(day_ordinals, divergence):
    """Exact joint histogram of (observed day, root-to-tip divergence).

    Returns ``(days, rtt, counts)`` as three aligned ``int64`` arrays over the
    non-empty cells, sorted by ``(day, rtt)``.

    **Binning loses nothing here.** Both axes are already integers -- a collection
    date at day resolution and a mutation count -- so the histogram is the joint
    distribution, not a summary of it. Every statistic drawn from it (the regression
    in :func:`weighted_linregress`, the rendered scatter) is therefore identical to
    the one over the per-tip cloud, which is what lets a 17 M-tip figure ship as a
    few tens of kB alongside a 1.6 GB per-tip residuals file that already exists.

    Cells are keyed as ``day * span + rtt`` and passed to a 1-D ``np.unique``
    rather than ``np.unique(..., axis=0)``: the 2-D form lexsorts a 17 M x 2 array
    and measured far slower for an identical result.

    Callers pass **day-resolution tips only** (see ``chronumental_qc.py`` check 4 on
    why imputed month/year dates must not enter a residual or regression), whereas
    chronumental regresses over every tip it dated. The two tip sets differ, so the
    two slopes are close but not equal -- :data:`RTT_SLOPE_REL_TOL` is sized for
    that, and the report states both numbers rather than implying one is the other.
    """
    days = np.asarray(day_ordinals, dtype="float64")
    rtt = np.asarray(divergence)
    if days.shape != rtt.shape:
        raise ValueError(
            f"day_ordinals {days.shape} and divergence {rtt.shape} must align"
        )
    keep = np.isfinite(days)
    days = days[keep].astype("int64")
    rtt = rtt[keep].astype("int64")
    if not days.size:
        empty = np.empty(0, dtype="int64")
        return empty, empty.copy(), empty.copy()
    # Encode with an explicit span so decoding is exact for any rtt >= 0.
    span = int(rtt.max()) + 1
    keys, counts = np.unique(days * span + rtt, return_counts=True)
    return keys // span, keys % span, counts.astype("int64")


def weighted_linregress(x, y, weights):
    """Least squares of ``y`` on ``x`` where each point stands for ``weights`` copies.

    Returns ``{"slope", "intercept", "r", "n"}`` with ``n`` the total weight.

    This is *exactly* the unweighted regression over the expanded point cloud, which
    is the property that makes :func:`root_to_tip_histogram` a lossless
    representation rather than a downsample -- and it is asserted in the unit test
    against an explicitly expanded array, not just documented here.
    """
    x = np.asarray(x, dtype="float64")
    y = np.asarray(y, dtype="float64")
    w = np.asarray(weights, dtype="float64")
    if not (x.shape == y.shape == w.shape):
        raise ValueError(
            f"x {x.shape}, y {y.shape} and weights {w.shape} must all align"
        )
    n = float(w.sum())
    if n <= 0:
        return {"slope": None, "intercept": None, "r": None, "n": 0.0}
    mean_x = float((w * x).sum() / n)
    mean_y = float((w * y).sum() / n)
    dx = x - mean_x
    dy = y - mean_y
    var_x = float((w * dx * dx).sum())
    var_y = float((w * dy * dy).sum())
    cov = float((w * dx * dy).sum())
    if var_x <= 0:
        # Every tip shares one date: a slope is undefined rather than infinite.
        return {"slope": None, "intercept": None, "r": None, "n": n}
    slope = cov / var_x
    r = cov / np.sqrt(var_x * var_y) if var_y > 0 else None
    return {
        "slope": slope,
        "intercept": mean_y - slope * mean_x,
        "r": float(r) if r is not None else None,
        "n": n,
    }


def root_to_tip_note(regression, n_excluded, logged_rate, stdout_path):
    """The report fragment for the recomputed root-to-tip regression.

    Formatted here, next to the math, for the same reason as
    :func:`implied_clock_note`: the one property worth pinning is that the
    comparison against chronumental's own logged slope is either **made or its
    absence is stated**. A line that silently omits the ratio when the log has been
    rotated away is indistinguishable from one where the two agreed.
    """
    slope = regression["slope"]
    if slope is None:
        # Two distinct ways to have no slope, and they mean different things: no
        # day-resolution tips at all (normal on a tree whose metadata is coarse) vs
        # tips that all share one date (a degenerate cloud). Reporting the second
        # wording for the first case would send someone looking for a date problem
        # that is really a resolution problem.
        cause = (
            "no day-resolution tips"
            if regression["n"] == 0
            else "all tips share one collection date"
        )
        return (
            f"root-to-tip regression undefined over {regression['n']:.0f} "
            f"day-resolution tip(s), {n_excluded} excluded: {cause}"
        )
    per_year = slope * DAYS_PER_YEAR
    r = regression["r"]
    r_text = f"{r:.4f}" if r is not None else "n/a"
    note = (
        f"recomputed slope {per_year:.2f} subs/yr (r={r_text}) over "
        f"{regression['n']:.0f} day-resolution tips, {n_excluded} excluded"
    )
    if logged_rate:
        ratio = per_year / logged_rate
        note += (
            f"; chronumental logged {logged_rate:.2f} subs/yr over all its dated "
            f"tips, ratio {ratio:.3f}"
        )
        if abs(ratio - 1.0) > RTT_SLOPE_REL_TOL:
            note += (
                f" -- OUTSIDE the {RTT_SLOPE_REL_TOL:.0%} tripwire; the two "
                f"regressions differ by more than differing tip sets explains"
            )
    else:
        note += (
            f"; chronumental's own slope unavailable for comparison (no parsable "
            f"'{ROOT_TO_TIP_PREFIX}' line at {stdout_path})"
        )
    return note
