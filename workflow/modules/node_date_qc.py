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

**Everything is integer days in histograms, never per-branch rows.** The rule
already costs 349 s at 38.2 GB max RSS over 20,961,614 nodes, and
``node_dates_residuals.csv`` is already 1.6 GB. Δt in whole days is an integer, so a
histogram over it gives **exact** medians and percentiles at a few tens of kB --
there is no accuracy/memory tradeoff being made here.
"""

import numpy as np

#: Mutation-count bins for check A. Upper bounds inclusive; the last is open.
MUT_COUNT_BINS = [0, 1, 2, 3, 4, 5, 10, 20, 50]

#: Rough clock expectation, days per mutation. 8e-4 subs/site/yr x 29903 nt
#: ~= 23.9 subs/yr -> ~15.3 days/substitution. An order-of-magnitude reference
#: only: the MAT is masked and branch-filtered, so the realized rate on counted
#: mutations is lower and the observed slope is expected to exceed this.
DAYS_PER_MUTATION_NOMINAL = 365.25 / (8e-4 * 29903)

#: Plausible window for the SARS-CoV-2 MRCA. Deliberately wide -- this is a
#: tripwire for a mis-scaled time axis, not an estimate.
ROOT_DATE_PLAUSIBLE = ("2019-10-01", "2020-02-29")


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
                "median": None,
                "p95": None,
                "max": None,
                "mean": None,
                "frac_negative": None,
            }
        negatives = sum(c for v, c in self._counts.items() if v < 0)
        return {
            "n": self.n,
            "min": min(self._counts),
            "p05": self.quantile(0.05),
            "median": self.quantile(0.5),
            "p95": self.quantile(0.95),
            "max": max(self._counts),
            "mean": self.total / self.n,
            "frac_negative": negatives / self.n,
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
    """
    if len(parent_index) != n_nodes or len(leaf_values) != n_nodes:
        raise ValueError(
            f"parent_index ({len(parent_index)}) and leaf_values "
            f"({len(leaf_values)}) must both have length n_nodes ({n_nodes})"
        )
    best = np.array(leaf_values, dtype="float64", copy=True)
    for i in range(n_nodes - 1, 0, -1):
        parent = parent_index[i]
        if parent < 0:
            continue
        value = best[i]
        if not np.isnan(value) and (np.isnan(best[parent]) or value < best[parent]):
            best[parent] = value
    return best


def format_stats(stats, unit="d"):
    """Compact one-line rendering for the report's ``detail`` column."""
    if not stats["n"]:
        return "no observations"
    return (
        f"n={stats['n']}, median={stats['median']}{unit}, "
        f"mean={stats['mean']:.1f}{unit}, p05={stats['p05']}{unit}, "
        f"p95={stats['p95']}{unit}, min={stats['min']}{unit}, "
        f"max={stats['max']}{unit}, frac<0={stats['frac_negative']:.4f}"
    )
