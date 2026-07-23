#!/usr/bin/env python3
"""Sanity-check chronumental node dates for one tree.

Runs on every input tree (wired into ``rule all`` via ``chronumental_qc``), so
it gates production runs and CI alike. Three hard invariants fail the run on any
violation; the tip-residual check reports every run and fails only when the
median absolute residual exceeds the configured threshold.

Checks:
  1. Node-ID correspondence (hard): every non-root ``bte`` node id has a date in
     ``node_dates.tsv``. This is the load-bearing assumption for windowed fitness
     (#117): matUtils Newick node names (which chronumental keys on) must equal
     the ``bte`` node ids the windowing will traverse on the same ``.pb``.
  2. Time monotonicity (hard): along every edge, parent date <= child date (+eps).
  3. Coverage / no missing dates (hard): every node dated, no unparseable dates.
  4. Tip residuals (report + configurable gate): predicted - observed (days) for
     day-resolution tips; fail if median(|residual|) exceeds the threshold.
"""

import pandas as pd

import bte

# --- snakemake I/O -----------------------------------------------------------
tree_path = snakemake.input.tree
node_dates_path = snakemake.input.node_dates
metadata_path = snakemake.input.metadata
report_path = snakemake.output.report
residuals_path = snakemake.output.residuals
max_median_abs_tip_residual_days = float(
    snakemake.params.max_median_abs_tip_residual_days
)
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None

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
id_col, date_col = nd.columns[0], nd.columns[1]
nd[date_col] = pd.to_datetime(nd[date_col], errors="coerce")
n_unparseable = int(nd[date_col].isna().sum())
node_date = dict(zip(nd[id_col].astype(str), nd[date_col]))

# --- load tree ---------------------------------------------------------------
tree = bte.MATree(tree_path)
nodes = tree.depth_first_expansion()
node_ids = {str(n.id) for n in nodes}
n_nodes = len(node_ids)

# Edges (parent -> child) and root detection, both from the child lists so we
# never rely on a ``.parent`` accessor.
edges = []
child_ids = set()
for n in nodes:
    for c in n.children:
        edges.append((str(n.id), str(c.id)))
        child_ids.add(str(c.id))
root_ids = node_ids - child_ids

results = []  # (check, status, detail)
failures = []


def record(check, ok, detail):
    results.append((check, "PASS" if ok else "FAIL", detail))
    if not ok:
        failures.append(f"{check}: {detail}")


# --- check 1: node-ID correspondence ----------------------------------------
dated_ids = set(node_date)
missing = node_ids - dated_ids
# Allow the root(s) to be undated; anything else is a real mismatch.
missing_non_root = missing - root_ids
record(
    "node_id_correspondence",
    not missing_non_root,
    f"{len(missing_non_root)} non-root node(s) without a date "
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
    n_unparseable == 0 and not missing_non_root,
    f"{n_unparseable} unparseable date(s); "
    f"{len(dated_ids)} dated vs {n_nodes} tree nodes",
)

# --- check 4: tip residuals --------------------------------------------------
meta = pd.read_csv(metadata_path, sep="\t", usecols=["strain", "date"])
# Day resolution only (YYYY-MM-DD): exclude MONTH/YEAR/MISSING, matching
# filter_metadata's treatment of date resolution.
meta = meta[meta["date"].astype(str).str.count("-") == 2].copy()
meta["observed"] = pd.to_datetime(meta["date"], errors="coerce")
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
