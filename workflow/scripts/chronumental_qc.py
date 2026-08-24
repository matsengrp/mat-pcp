#!/usr/bin/env python3
"""Sanity-check chronumental node dates for one tree.

Runs on every input tree (wired into ``rule all`` via ``chronumental_qc``), so
it gates production runs and CI alike. Three hard invariants fail the run on any
violation; the tip-residual check reports every run and fails only when the
median absolute residual exceeds the configured threshold.

Checks:
  1. Node-ID correspondence (hard): every ``bte`` node id -- root included --
     has a date in ``node_dates.tsv``. This is the load-bearing assumption for
     windowed fitness (#117): matUtils Newick node names (which chronumental keys
     on) must equal the ``bte`` node ids the windowing will traverse on the same
     ``.pb``, and windowing relies on all nodes being dated (no contingency).
  2. Time monotonicity (hard): along every edge, parent date <= child date (+eps).
  3. Coverage / no missing dates (hard): every node dated, no unparseable dates.
  4. Tip residuals (report + configurable gate): predicted - observed (days) for
     day-resolution tips; fail if median(|residual|) exceeds the threshold.

Checks 1-4 validate *tips* and *structure*. Checks 5-8 (#139) validate the
**internal** node dates, which are what ``modules/windowing.py`` actually consumes
to bucket every branch into a time window. They are **report-only and gate nothing**
-- see ``modules/node_date_qc.py`` for why, and for why check 4 passing at
``p95|res|=1.0d`` over 17 M tips is close to circular rather than reassuring.

  5. Branch Delta-t vs branch mutation count (report): binned; the clock should make
     these roughly proportional.
  6. Zero-mutation leaf/parent dates (report): both the predicted-date form (a
     one-sided right tail, since check 2 forces Delta-t >= 0) and the non-circular
     form -- parent predicted date vs the *earliest observed* date among its
     zero-mutation leaf children.
  7. Root date plausibility (report): a tripwire for a mis-scaled time axis.
  8. Internal date vs minimum descendant observed date (report): an ancestor cannot
     postdate the earliest sample beneath it.
"""

import os
from importlib import util as importlib_util

import numpy as np
import pandas as pd

import bte

# --- snakemake I/O -----------------------------------------------------------
tree_path = snakemake.input.tree
node_dates_path = snakemake.input.node_dates
metadata_path = snakemake.input.metadata
report_path = snakemake.output.report
residuals_path = snakemake.output.residuals
branch_dt_path = snakemake.output.branch_dt
internal_checks_path = snakemake.output.internal_checks
branch_dt_summary_path = snakemake.output.branch_dt_summary
max_median_abs_tip_residual_days = float(
    snakemake.params.max_median_abs_tip_residual_days
)
# A params path, deliberately NOT an input: chronumental's stdout is a Snakemake
# ``log:`` entry rather than a declared output, so making it an input would put the
# whole (multi-hour) dating rule back in the DAG whenever a log was rotated away.
# Absent is therefore a legitimate state, and check 5 reports it as such.
chronumental_stdout_path = snakemake.params.chronumental_stdout
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None

# The by-path module loader is defined once in workflow/modules/workflow_loader.py.
# Because that loader is itself the mechanism for importing shared workflow code, it
# cannot be imported via itself -- it is wired up here by absolute path (#111).
_loader_spec = importlib_util.spec_from_file_location(
    "workflow_loader",
    os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "modules",
        "workflow_loader.py",
    ),
)
_loader = importlib_util.module_from_spec(_loader_spec)
_loader_spec.loader.exec_module(_loader)
node_date_qc = _loader.load_workflow_module("node_date_qc")

EPS_DAYS = 1  # float/rounding slack for the monotonicity check


def log(msg):
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


# --- load node dates ---------------------------------------------------------
# chronumental --dates_out: header row then two columns (node name, predicted
# date). Column names are "strain"/"predicted_date" but we read positionally so
# the check does not depend on those labels.
nd = pd.read_csv(node_dates_path, sep="\t")
assert len(nd.columns) >= 2, f"node_dates.tsv has < 2 columns: {list(nd.columns)}"
id_col, date_col = nd.columns[0], nd.columns[1]
nd[date_col] = pd.to_datetime(nd[date_col], errors="coerce")
n_unparseable = int(nd[date_col].isna().sum())
node_date = dict(zip(nd[id_col].astype(str), nd[date_col]))

# --- load tree ---------------------------------------------------------------
tree = bte.MATree(tree_path)
nodes = tree.depth_first_expansion()
node_ids = {str(n.id) for n in nodes}
n_nodes = len(node_ids)
assert n_nodes > 0, "tree has no nodes"

# Edges (parent -> child) and root detection, both from the child lists so we
# never rely on a ``.parent`` accessor.
edges = []
child_ids = set()
for n in nodes:
    for c in n.children:
        edges.append((str(n.id), str(c.id)))
        child_ids.add(str(c.id))
root_ids = node_ids - child_ids
assert len(root_ids) == 1, f"expected exactly one root, found {len(root_ids)}"

results = []  # (check, status, detail)
failures = []


def record(check, ok, detail):
    results.append((check, "PASS" if ok else "FAIL", detail))
    if not ok:
        failures.append(f"{check}: {detail}")


# --- check 1: node-ID correspondence ----------------------------------------
# Membership requires a *parseable* date, not just a row in node_dates.tsv, so a
# node with a garbled date string counts as missing here (not silently "dated").
# chronumental dates the FULL time-tree -- every node including the root -- so we
# require every bte node id (root included) to have a date. This is the invariant
# downstream windowing (#117) relies on to key node ids into dates without any
# undated-node contingency; an undated node here means something is broken.
dated_ids = {i for i, d in node_date.items() if pd.notna(d)}
missing = node_ids - dated_ids
record(
    "node_id_correspondence",
    not missing,
    f"{len(missing)} node(s) without a date "
    f"(roots={len(root_ids)}, total nodes={n_nodes})",
)

# --- check 2: time monotonicity ----------------------------------------------
violations = 0
worst_days = 0
for p, c in edges:
    dp, dc = node_date.get(p), node_date.get(c)
    if pd.isna(dp) or pd.isna(dc):
        continue
    delta_days = (dc - dp).days
    if delta_days < -EPS_DAYS:
        violations += 1
        worst_days = min(worst_days, delta_days)
record(
    "time_monotonicity",
    violations == 0,
    f"{violations}/{len(edges)} edges with parent date > child date "
    f"(worst = {worst_days} days)",
)

# --- check 3: coverage / no missing dates ------------------------------------
record(
    "date_coverage",
    n_unparseable == 0 and not missing,
    f"{n_unparseable} unparseable date(s); "
    f"{len(dated_ids)} dated vs {n_nodes} tree nodes",
)

# --- check 4: tip residuals --------------------------------------------------
meta = pd.read_csv(
    metadata_path, sep="\t", usecols=["strain", "date", "date_resolution_str"]
)
# Day resolution only: exclude MONTH/YEAR/MISSING. filter_metadata imputes
# month/year dates to a fabricated day (mid-month / mid-year) and rewrites every
# surviving value to YYYY-MM-DD, so we must filter on its recorded resolution
# rather than re-deriving it from the (already-imputed) date string -- otherwise
# the residual gate would silently compare predictions against fabricated dates
# if the upstream resolution filter is ever loosened.
meta = meta[meta["date_resolution_str"] == "DAY"].copy()
meta["observed"] = pd.to_datetime(meta["date"], errors="coerce")
n_dropped_obs = int(meta["observed"].isna().sum())
if n_dropped_obs:
    log(f"tip_residual: dropped {n_dropped_obs} day-resolution row(s) with unparseable date")
meta = meta.dropna(subset=["observed"])
meta["predicted"] = meta["strain"].astype(str).map(node_date)
res = meta.dropna(subset=["predicted"]).copy()
res["residual_days"] = (res["predicted"] - res["observed"]).dt.days
res[["strain", "predicted", "observed", "residual_days"]].to_csv(
    residuals_path, index=False
)

if len(res):
    abs_res = res["residual_days"].abs()
    median_abs = float(abs_res.median())
    mean_abs = float(abs_res.mean())
    p95_abs = float(abs_res.quantile(0.95))
    n_tips = len(res)
    # The tail, alongside the median the gate uses. The gate is a median over ~17 M
    # tips and is insensitive to a tail of any size: on
    # gisaid-and-public_2026-04-22 it reads median|res|=1.0d against a 30 d
    # threshold while 208,995 tips (1.23%) have predicted LATER than observed, up
    # to +1663 d. That tail is what check 8 then re-reports, amplified ~2x, as
    # "3.14% of internal nodes postdate their earliest sample" -- so reporting only
    # the median here is what made the tail look like a separate discovery about
    # internal nodes. Report-only: the gate itself is unchanged.
    n_late = int((res["residual_days"] > 0).sum())
    n_late_30 = int((res["residual_days"] > 30).sum())
    max_late = int(res["residual_days"].max())
    detail = (
        f"median|res|={median_abs:.1f}d, mean|res|={mean_abs:.1f}d, "
        f"p95|res|={p95_abs:.1f}d, N={n_tips} "
        f"(threshold={max_median_abs_tip_residual_days:.1f}d); "
        f"positive tail (predicted later than observed, which bounds check 8): "
        f"{n_late} tips ({n_late / n_tips:.4%}), {n_late_30} by >30d, "
        f"max={max_late}d"
    )
    record("tip_residual", median_abs <= max_median_abs_tip_residual_days, detail)
else:
    # No day-resolution tips to compare: report, do not gate.
    record(
        "tip_residual",
        True,
        "no day-resolution tips available; residual gate skipped",
    )

# =============================================================================
# Checks 5-8 (#139): the INTERNAL node dates, which windowing actually consumes.
#
# All report-only. ``record_report`` never appends to ``failures``, so nothing
# below can fail the run -- thresholds are chosen in a follow-up once these
# distributions have been looked at.
#
# One pass, integer-day histograms. ``nodes`` is already in memory from check 1
# and is not re-expanded; Delta-t in whole days is an integer, so the histograms
# give exact medians and percentiles rather than approximations.
# =============================================================================
DAY = np.timedelta64(1, "D")


def record_report(check, detail):
    """Record a report-only row: status REPORT, never a failure."""
    results.append((check, "REPORT", detail))


def to_day_ordinal(series):
    """Datetime64 -> float days since epoch, ``nan`` where missing."""
    values = pd.to_datetime(series).to_numpy(dtype="datetime64[D]", copy=False)
    out = values.astype("float64")
    out[pd.isna(series).to_numpy()] = np.nan
    return out


# Index every node once. Parent links come from the child lists, matching how
# ``edges`` is built above, so we never depend on a ``.parent`` accessor.
node_list = nodes
index_of = {str(n.id): i for i, n in enumerate(node_list)}
n_total = len(node_list)
parent_index = np.full(n_total, -1, dtype="int64")
n_branch_mutations = np.zeros(n_total, dtype="int32")
for i, n in enumerate(node_list):
    n_branch_mutations[i] = len(n.mutations)
    for c in n.children:
        parent_index[index_of[str(c.id)]] = i

predicted_days = np.array(
    [
        node_date[nid].to_datetime64().astype("datetime64[D]").astype("float64")
        if node_date.get(nid) is not None and pd.notna(node_date.get(nid))
        else np.nan
        for nid in (str(n.id) for n in node_list)
    ],
    dtype="float64",
)

is_leaf = np.ones(n_total, dtype=bool)
is_leaf[parent_index[parent_index >= 0]] = False

# --- check 5: branch Delta-t vs branch mutation count -------------------------
# Vectorized: a per-node Python loop cost minutes on the 21 M-node production tree,
# against a 349 s baseline for the whole rule. Every branch is one (child, parent)
# index pair, so the whole check is array arithmetic plus one np.unique per bin.
has_parent = parent_index >= 0
child_idx = np.flatnonzero(has_parent)
parent_of_child = parent_index[child_idx]
branch_dt = predicted_days[child_idx] - predicted_days[parent_of_child]
branch_muts = n_branch_mutations[child_idx]

# Bin edges from MUT_COUNT_BINS, so the labels stay the module's single definition.
bin_labels = np.array(
    [node_date_qc.mutation_count_bin(b) for b in node_date_qc.MUT_COUNT_BINS]
    + [node_date_qc.mutation_count_bin(node_date_qc.MUT_COUNT_BINS[-1] + 1)]
)
bin_of_branch = np.searchsorted(
    np.array(node_date_qc.MUT_COUNT_BINS), branch_muts, side="left"
)

dt_by_bin = {}
for bin_pos in np.unique(bin_of_branch):
    label = bin_labels[bin_pos]
    values = branch_dt[bin_of_branch == bin_pos]
    hist = node_date_qc.DayHistogram()
    hist.update_from_array(values)
    if hist.n:
        dt_by_bin[label] = hist

def bin_lower_bound(label):
    """Smallest mutation count in a bin label, for sorting and the slope estimate."""
    return int(str(label).split("-")[0].rstrip("+"))


# Long form: the actual per-bin day distribution, which is what a follow-up needs
# to pick a threshold. Summary statistics go in the report line rather than here so
# there is one source of truth for the distribution.
branch_rows = []
for label in sorted(dt_by_bin, key=bin_lower_bound):
    for days, count in dt_by_bin[label].as_rows():
        branch_rows.append({"mut_count_bin": label, "dt_days": days, "count": count})
pd.DataFrame(
    branch_rows, columns=["mut_count_bin", "dt_days", "count"]
).to_csv(branch_dt_path, index=False)

# Wide form: the per-bin summary. Derived from the same DayHistogram objects that
# produce the long form above, so the two cannot disagree -- and every statistic is
# EXACT rather than interpolated, because Delta-t is whole days and the histogram
# holds every observation (see the DayHistogram docstring). ``iqr`` is p75 - p25 on
# those same exact order statistics, not a separate estimate.
branch_summary = []
for label in sorted(dt_by_bin, key=bin_lower_bound):
    s = dt_by_bin[label].stats()
    branch_summary.append(
        {
            "mut_count_bin": label,
            "n_branches": s["n"],
            "min": s["min"],
            "p05": s["p05"],
            "p25": s["p25"],
            "median": s["median"],
            "p75": s["p75"],
            "iqr": s["iqr"],
            "p95": s["p95"],
            "max": s["max"],
            "mean": round(s["mean"], 2),
            "frac_dt_zero": round(
                dict(dt_by_bin[label].as_rows()).get(0, 0) / s["n"], 6
            ),
        }
    )
# An "all" row, so the marginal distribution does not have to be recombined by hand
# from the bins (and cannot be got wrong by averaging per-bin medians).
all_hist = node_date_qc.DayHistogram()
for hist in dt_by_bin.values():
    all_hist.merge(hist)
s = all_hist.stats()
branch_summary.append(
    {
        "mut_count_bin": "all",
        "n_branches": s["n"],
        "min": s["min"],
        "p05": s["p05"],
        "p25": s["p25"],
        "median": s["median"],
        "p75": s["p75"],
        "iqr": s["iqr"],
        "p95": s["p95"],
        "max": s["max"],
        "mean": round(s["mean"], 2),
        "frac_dt_zero": round(dict(all_hist.as_rows()).get(0, 0) / s["n"], 6),
    }
)
n_dt_defined = int((~np.isnan(branch_dt)).sum())
if s["n"] != n_dt_defined:
    raise ValueError(
        f"the branch Delta-t summary's 'all' row covers {s['n']} branches but "
        f"{n_dt_defined} branches have a defined Delta-t; the mutation-count bins "
        f"are not a partition of them"
    )
pd.DataFrame(branch_summary).to_csv(branch_dt_summary_path, index=False)

# The clock estimate is the INCREMENTAL slope between adjacent single-count bins,
# not median(Delta-t)/n. On the production tree the zero-mutation bin already has a
# median Delta-t of 16 days, so median/n carries that intercept and decays from 23.0
# at n=1 to 2.6 at n>=51 -- a 9x spread that is mostly intercept dilution plus
# saturation on the ~650 branches with >20 mutations, not a clock failure. Reporting
# median/n as "days per mutation" invites reading a defect off an artifact of the
# arithmetic, so the median per bin is reported instead and the slope is taken from
# differences. The nominal ~15 d/mut stays an order-of-magnitude reference: the MAT
# is masked and branch-filtered, so the realized figure is not expected to match.
total_branches = sum(hist.stats()["n"] for hist in dt_by_bin.values())
medians = {}
median_labels = {}
for label, hist in dt_by_bin.items():
    stats = hist.stats()
    if stats["n"]:
        lower = bin_lower_bound(label)
        medians[lower] = stats["median"]
        # Keep the bin's own label for the report: printing the lower bound would
        # render the "6-10" bin as "6", implying a single mutation count.
        median_labels[lower] = label

# Only bins that are a single mutation count (the leading MUT_COUNT_BINS entries)
# have a well-defined increment; the ranged bins at the top do not.
single_counts = sorted(n for n in medians if n in node_date_qc.SINGLE_COUNT_BINS)
steps = sorted(
    medians[b] - medians[a]
    for a, b in zip(single_counts, single_counts[1:])
    if b - a == 1
)
clock_note = ""
if steps:
    slope = steps[len(steps) // 2]
    intercept = medians.get(0)
    clock_note = (
        f", median incremental slope {slope} d/mutation across the single-count "
        f"bins"
        + (f" on a {intercept} d zero-mutation intercept" if intercept is not None else "")
        + f" (nominal ~{node_date_qc.DAYS_PER_MUTATION_NOMINAL:.1f} d/mutation)"
    )
# The aggregate the incremental slope above cannot see (#146). The slope cancels
# the zero-mutation intercept by construction, and on the v11 trees that intercept
# is where all of the duration inflation confirmed in #145 lives -- so check 5 read
# 11 d/mutation against a ~15.3 d nominal and looked healthy while the tree's own
# implied rate sat 4.3x below chronumental's pre-fit root-to-tip estimate. Computed
# from the per-branch arrays, so it is exact rather than reconstructed from the
# ranged bins of the summary CSV (which moves the answer ~+/-7%).
clock_summary = node_date_qc.implied_clock_summary(branch_dt, branch_muts)
aggregate_note = node_date_qc.implied_clock_note(
    clock_summary,
    node_date_qc.parse_root_to_tip_rate(chronumental_stdout_path),
    chronumental_stdout_path,
)
record_report(
    "branch_dt_vs_mutations",
    f"{len(dt_by_bin)} mutation-count bins over {total_branches} dated branches"
    f"{clock_note}{aggregate_note}; median Delta-t by bin: "
    + ", ".join(f"{median_labels[n]}mut:{medians[n]}d" for n in sorted(medians))
    + f"; full per-bin day distributions in {os.path.basename(branch_dt_path)}",
)

# --- check 6: zero-mutation leaf/parent pairs --------------------------------
# The predicted-date form is one-sided: check 2 forces Delta-t >= 0, so this is a
# right tail and a median of 0 is the expected result, not a coincidence.
# ``res`` carries only DAY-resolution tips, inherited from check 4's filter, and that
# inheritance is deliberate: ``filter_metadata`` imputes month/year dates to a
# fabricated day (mid-month / mid-year) and rewrites every surviving value to
# YYYY-MM-DD, so admitting them here would compare predicted dates against invented
# ones and quietly widen both distributions. Checks 6 and 8 therefore cover fewer
# nodes than the tree has, which the reported ``n`` makes visible.
own_observed = np.full(n_total, np.nan, dtype="float64")
observed_days = to_day_ordinal(res["observed"])
observed_pos = np.array(
    [index_of.get(s, -1) for s in res["strain"].astype(str)], dtype="int64"
)
known = observed_pos >= 0
# Counted, not silently dropped. ``res`` is filtered on having a predicted date,
# which comes from node_dates.tsv; ``index_of`` comes from the tree. Those are
# different universes and can diverge (a dating run over a superset tree), and
# the comment above argues the OTHER exclusion here is acceptable precisely
# because the reported n makes it visible. Same standard applies to this one.
n_strains_not_in_tree = int((~known).sum())
own_observed[observed_pos[known]] = observed_days[known]

zero_mut_leaf = is_leaf & (n_branch_mutations == 0) & has_parent
zero_pred = node_date_qc.DayHistogram()
zero_pred.update_from_array(
    predicted_days[zero_mut_leaf] - predicted_days[parent_index[zero_mut_leaf]]
)

# The non-circular form: the parent's predicted date against the EARLIEST observed
# date among its zero-mutation leaf children. Per-leaf equality is the wrong bar --
# identical sequences in the data genuinely carry different collection dates.
# np.minimum.at does the per-parent minimum in one pass; inf is the identity so a
# parent with no observed zero-mutation child stays inf and drops out below.
observed_zero_leaf = zero_mut_leaf & ~np.isnan(own_observed)
min_obs_by_parent = np.full(n_total, np.inf, dtype="float64")
np.minimum.at(
    min_obs_by_parent,
    parent_index[observed_zero_leaf],
    own_observed[observed_zero_leaf],
)
has_min_obs = np.isfinite(min_obs_by_parent) & ~np.isnan(predicted_days)
n_zero_mut_parents = int(has_min_obs.sum())
zero_obs = node_date_qc.DayHistogram()
zero_obs.update_from_array(
    min_obs_by_parent[has_min_obs] - predicted_days[has_min_obs]
)

record_report(
    "zero_mut_leaf_predicted_dt",
    # frac<0 here is forced to 0 by check 2, which is a hard gate -- so quoting it
    # as a confirmed prediction would be quoting an upstream invariant. Two
    # loopholes keep it from being fully forced, and ``min`` is what settles them:
    # EPS_DAYS allows a -1 day edge past check 2, and a 4-decimal fraction over
    # 12 M observations hides up to ~600 negatives.
    "one-sided by construction (check 2 is a hard gate forcing >=0, so frac<0 is "
    "an invariant, not a measurement -- read min); "
    + node_date_qc.format_stats(zero_pred.stats()),
)
record_report(
    "zero_mut_leaf_vs_earliest_observed",
    f"earliest observed date among a parent's zero-mutation leaf children MINUS "
    f"the parent's predicted date, over {n_zero_mut_parents} parents (note this "
    f"is the opposite orientation to check 8: here a NEGATIVE value is the "
    f"anomaly, meaning the parent postdates a leaf identical to it); "
    + node_date_qc.format_stats(zero_obs.stats(), violation="negative")
    + f"; {n_strains_not_in_tree} dated strain(s) had no node in the tree and are "
    f"excluded from checks 6 and 8",
)

# --- check 7: root date plausibility ------------------------------------------
root_id = next(iter(root_ids))
root_date = node_date.get(root_id)
lo, hi = (pd.Timestamp(d) for d in node_date_qc.ROOT_DATE_PLAUSIBLE)
if root_date is None or pd.isna(root_date):
    record_report("root_date_plausible", f"root {root_id} has no parseable date")
else:
    within = lo <= root_date <= hi
    record_report(
        "root_date_plausible",
        f"root {root_id} dated {root_date.date()}; "
        f"{'inside' if within else 'OUTSIDE'} the plausible window "
        f"[{lo.date()}, {hi.date()}]",
    )

# --- check 8: internal date vs minimum descendant observed date ---------------
# own_observed was filled for check 6 above; reuse it rather than rebuilding.
min_desc_observed = node_date_qc.propagate_min_descendant(
    parent_index, own_observed, n_total
)

internal = ~is_leaf & ~np.isnan(min_desc_observed) & ~np.isnan(predicted_days)
n_internal_checked = int(internal.sum())
internal_gap = node_date_qc.DayHistogram()
internal_gap.update_from_array(
    predicted_days[internal] - min_desc_observed[internal]
)

stats = internal_gap.stats()
record_report(
    "internal_date_vs_min_descendant_observed",
    f"predicted internal date minus earliest observed descendant date over "
    f"{n_internal_checked} internal nodes (a POSITIVE value is the anomaly: the "
    f"ancestor postdates its earliest sample); "
    + node_date_qc.format_stats(stats, violation="positive")
    + "; NOTE this is bounded by the tip residuals -- check 2 forces "
    "pred(ancestor) <= pred(tip), so a positive gap requires that tip's own "
    "residual to be at least as large. It is an amplified view of check 4's "
    "positive tail, not independent evidence about internal placement.",
)

# --- check 9: internal node to its CLOSEST descendant, inferred vs inferred ---
# Checks 6 and 8 both compare an inferred date against an *observed* one, so a
# violation there is a three-way inconsistency between the inferred date, the
# topology and a collection date -- and check 8 in particular is bounded by the tip
# residuals (see its note), which makes it a metadata-quality measure wearing an
# algorithm-test label. This check compares two chronumental outputs, so it tests
# the algorithm: given the topology it was handed, how much time does it put
# between an internal node and its nearest descendant?
#
# The minimum over ALL descendants is attained at an immediate CHILD. Check 2
# forces predicted dates to be nondecreasing along every path, so for any deeper
# descendant d reached via child c, pred(d) >= pred(c). That collapses what looks
# like a propagation problem into one np.minimum.at over the edges already built
# for check 5 -- no second traversal, and no reliance on node ordering.
min_child_dt = node_date_qc.min_child_delta(parent_index, predicted_days, n_total)
internal_dated = ~np.isnan(min_child_dt) & ~is_leaf
n_internal_dt = int(internal_dated.sum())
closest_dt = node_date_qc.DayHistogram()
closest_dt.update_from_array(min_child_dt[internal_dated])

# Zero-mutation branches: on the v9 production tree all 12,103,637 of them are
# terminal, so check 6 covers the whole population rather than a leaf subset. That
# identity is reported rather than assumed -- if internal zero-mutation branches
# ever appear, check 6 stops being the complete picture of them.
n_zero_mut_all = int((n_branch_mutations[child_idx] == 0).sum())
n_zero_mut_internal = int(
    ((n_branch_mutations[child_idx] == 0) & ~is_leaf[child_idx]).sum()
)
record_report(
    "internal_dt_to_closest_descendant",
    f"predicted Delta-t from an internal node to its NEAREST descendant, both "
    f"inferred (so this tests the algorithm, not the metadata), over "
    f"{n_internal_dt} internal nodes; "
    + node_date_qc.format_stats(closest_dt.stats())
    + f"; equals min-over-children of the branch Delta-t because check 2 makes "
    f"predicted dates nondecreasing along every path. Of {n_zero_mut_all} "
    f"zero-mutation branches, {n_zero_mut_internal} are internal (the rest are "
    f"terminal and so are covered by check 6)",
)

# Binned artifact for checks 6, 8 and 9, so the follow-up can pick thresholds from a
# distribution rather than a summary line.
pd.concat(
    [
        pd.DataFrame(
            {"check": name, "days": [d for d, _ in hist.as_rows()],
             "count": [c for _, c in hist.as_rows()]}
        )
        for name, hist in (
            ("zero_mut_leaf_predicted_dt", zero_pred),
            ("zero_mut_leaf_vs_earliest_observed", zero_obs),
            ("internal_date_vs_min_descendant_observed", internal_gap),
            ("internal_dt_to_closest_descendant", closest_dt),
        )
    ],
    ignore_index=True,
).to_csv(internal_checks_path, index=False)

# --- write report ------------------------------------------------------------
report = pd.DataFrame(results, columns=["check", "status", "detail"])
report.to_csv(report_path, sep="\t", index=False)
for _, row in report.iterrows():
    log(f"[{row['status']}] {row['check']}: {row['detail']}")

if failures:
    raise SystemExit(
        "chronumental_qc failed:\n  " + "\n  ".join(failures)
    )
log("chronumental_qc: all checks passed")
