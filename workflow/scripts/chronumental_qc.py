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
max_median_abs_tip_residual_days = float(
    snakemake.params.max_median_abs_tip_residual_days
)
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
    detail = (
        f"median|res|={median_abs:.1f}d, mean|res|={mean_abs:.1f}d, "
        f"p95|res|={p95_abs:.1f}d, N={n_tips} "
        f"(threshold={max_median_abs_tip_residual_days:.1f}d)"
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
dt_by_bin = {}
for i in range(n_total):
    p = parent_index[i]
    if p < 0:
        continue
    if np.isnan(predicted_days[i]) or np.isnan(predicted_days[p]):
        continue
    label = node_date_qc.mutation_count_bin(int(n_branch_mutations[i]))
    hist = dt_by_bin.get(label)
    if hist is None:
        hist = dt_by_bin[label] = node_date_qc.DayHistogram()
    hist.add(int(predicted_days[i] - predicted_days[p]))

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

# Days per mutation, over bins with at least one mutation. The nominal ~15 d/mut is
# an order-of-magnitude reference only: the MAT is masked and branch-filtered, so
# the realized rate on counted mutations is lower and the observed slope is expected
# to come out HIGHER than nominal, not equal to it.
ratios = []
total_branches = 0
for label, hist in dt_by_bin.items():
    stats = hist.stats()
    total_branches += stats["n"]
    lower = bin_lower_bound(label)
    if lower >= 1 and stats["n"]:
        ratios.append(stats["median"] / lower)
slope_note = ""
if ratios:
    slope_note = (
        f", median days/mutation across bins with >=1 mutation: "
        f"{min(ratios):.1f}-{max(ratios):.1f} "
        f"(nominal ~{node_date_qc.DAYS_PER_MUTATION_NOMINAL:.1f})"
    )
record_report(
    "branch_dt_vs_mutations",
    f"{len(dt_by_bin)} mutation-count bins over {total_branches} dated branches"
    f"{slope_note}; per-bin day distributions in "
    f"{os.path.basename(branch_dt_path)}",
)

# --- check 6: zero-mutation leaf/parent pairs --------------------------------
# The predicted-date form is one-sided: check 2 forces Delta-t >= 0, so this is a
# right tail and a median of 0 is the expected result, not a coincidence.
zero_pred = node_date_qc.DayHistogram()
# The non-circular form: the parent's predicted date against the EARLIEST observed
# date among its zero-mutation leaf children. Per-leaf equality is the wrong bar --
# identical sequences in the data genuinely carry different collection dates.
observed_by_id = dict(
    zip(res["strain"].astype(str), to_day_ordinal(res["observed"]))
)
min_obs_by_parent = {}
for i in range(n_total):
    p = parent_index[i]
    if p < 0 or not is_leaf[i] or n_branch_mutations[i] != 0:
        continue
    if not np.isnan(predicted_days[i]) and not np.isnan(predicted_days[p]):
        zero_pred.add(int(predicted_days[i] - predicted_days[p]))
    obs = observed_by_id.get(str(node_list[i].id))
    if obs is not None and not np.isnan(obs):
        prev = min_obs_by_parent.get(p)
        if prev is None or obs < prev:
            min_obs_by_parent[p] = obs

zero_obs = node_date_qc.DayHistogram()
for p, obs in min_obs_by_parent.items():
    if not np.isnan(predicted_days[p]):
        zero_obs.add(int(obs - predicted_days[p]))

record_report(
    "zero_mut_leaf_predicted_dt",
    "one-sided by construction (check 2 forces >=0); "
    + node_date_qc.format_stats(zero_pred.stats()),
)
record_report(
    "zero_mut_leaf_vs_earliest_observed",
    f"parent predicted date vs earliest observed date among its zero-mutation "
    f"leaf children over {len(min_obs_by_parent)} parents; "
    + node_date_qc.format_stats(zero_obs.stats()),
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
own_observed = np.full(n_total, np.nan, dtype="float64")
for nid, obs in observed_by_id.items():
    i = index_of.get(nid)
    if i is not None:
        own_observed[i] = obs
min_desc_observed = node_date_qc.propagate_min_descendant(
    parent_index, own_observed, n_total
)

internal_gap = node_date_qc.DayHistogram()
n_internal_checked = 0
for i in range(n_total):
    if is_leaf[i] or np.isnan(min_desc_observed[i]) or np.isnan(predicted_days[i]):
        continue
    n_internal_checked += 1
    internal_gap.add(int(predicted_days[i] - min_desc_observed[i]))

stats = internal_gap.stats()
record_report(
    "internal_date_vs_min_descendant_observed",
    f"predicted internal date minus earliest observed descendant date over "
    f"{n_internal_checked} internal nodes (a positive value means the ancestor "
    f"postdates its earliest sample); " + node_date_qc.format_stats(stats),
)

# Binned artifact for checks 6 and 8, so the follow-up can pick thresholds from a
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
