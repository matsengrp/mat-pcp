#!/usr/bin/env python3
"""Scatter the root-to-tip regression behind chronumental's clock reference (#146).

``chronumental_qc``'s check 5 divides the timetree's total mutations by its total
branch-years and compares that against the slope chronumental prints *before*
fitting. That pre-fit slope is the only independent reference in the QC -- comparing
against the fitted rate is near-identity under a learned clock -- and until now it
was a single number with nothing behind it. chronumental computes the cloud, prints
the slope and discards the cloud, and one slope describes a tight line and a fan of
recombinants and metadata errors equally well.

This plots that cloud from the exact joint histogram ``chronumental_qc`` writes.
Both axes are integers (a day-resolution collection date and a mutation count), so
the histogram is the joint distribution rather than a summary of it: the figure and
its regression are identical to what a 17 M-point scatter would give, at a few tens
of kB of input.

**Three lines, and the third one is the finding.** The recomputed regression comes from
the same ``modules/node_date_qc.weighted_linregress`` that produced the number in
the report, so the figure cannot disagree with the report. chronumental's own logged
pre-fit slope is drawn alongside it -- anchored at the cloud's centroid, because
chronumental logs its slope but not its intercept -- and the two regressions run
over *different tip sets* (day-resolution tips here; every dated tip there), so they
are expected to be close rather than equal.

The third line is chronumental's **fitted** clock, and unlike the other two it is
not a slope placed near the cloud for comparison: it is the model's own prediction,
``rate * (date - root_fitted_date)``, drawn from zero divergence on the root's
fitted date. That is exact under chronumental's parameterization rather than an
approximation -- see ``node_date_qc.model_implied_divergence`` for why the durations
telescope -- which is what makes the gap between it and the cloud readable as a
quantity. On ``public_2026-08-06_v11`` the fitted clock predicts 62 mutations at the
last collection date, where the cloud's own regression reads 247: the #145 duration
inflation, to scale.

When either logged rate is unavailable the figure **says so on its face** rather
than quietly dropping the line: a plot missing a comparison it never announces
looks exactly like one where the comparison agreed.
"""

import os
from importlib import util as importlib_util

import matplotlib

matplotlib.use("Agg")  # No display on a cluster node; must precede pyplot.

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

# --- snakemake I/O -----------------------------------------------------------
rtt_path = snakemake.input.rtt
rtt_anchor_path = snakemake.input.rtt_anchor
plot_path = snakemake.output.plot
tree_label = snakemake.params.tree_label
# A params path, not an input, for the same reason as in chronumental_qc: this is a
# Snakemake ``log:`` entry rather than a declared output, so requiring it would drag
# the multi-hour dating rule back into the DAG whenever a log was rotated away.
chronumental_stdout_path = snakemake.params.chronumental_stdout
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None

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


def log(msg):
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


cloud = pd.read_csv(rtt_path, parse_dates=["observed_date"])
if cloud.empty:
    # A tree with no day-resolution tips is a legitimate state (the residual gate in
    # chronumental_qc skips for the same reason), so write a figure that says that
    # rather than failing the run or leaving a stale file behind.
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.text(
        0.5,
        0.5,
        "no day-resolution tips:\nno root-to-tip regression to plot",
        ha="center",
        va="center",
        transform=ax.transAxes,
    )
    ax.set_axis_off()
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    log(f"root_to_tip_scatter: {rtt_path} is empty; wrote a placeholder figure")
    raise SystemExit(0)

days = cloud["observed_date"].to_numpy(dtype="datetime64[D]").astype("float64")
muts = cloud["rtt_mutations"].to_numpy(dtype="float64")
counts = cloud["n_tips"].to_numpy(dtype="float64")

regression = node_date_qc.weighted_linregress(days, muts, counts)
logged_rate = node_date_qc.parse_root_to_tip_rate(chronumental_stdout_path)
fitted_rate, fitted_step = node_date_qc.parse_fitted_mutation_rate(
    chronumental_stdout_path
)
# An empty cell means chronumental_qc ran and found no parseable root date (its
# check 7), which is a different state from a missing sidecar and is reported as
# such below rather than collapsing into "no line drawn".
anchor = pd.read_csv(rtt_anchor_path)
root_day = None
if len(anchor) and pd.notna(anchor["root_day_ordinal"].iloc[0]):
    root_day = float(anchor["root_day_ordinal"].iloc[0])

dates = cloud["observed_date"].to_numpy()
fig, ax = plt.subplots(figsize=(10, 6.5))
# Colour is how many tips share a cell, on a log scale: the cloud spans ~6 orders of
# magnitude in density, so a linear scale renders everything but the mode as blank.
points = ax.scatter(
    dates,
    muts,
    c=counts,
    s=6,
    norm=LogNorm(vmin=1, vmax=max(counts.max(), 2)),
    cmap="viridis",
    linewidths=0,
)
fig.colorbar(points, ax=ax, label="tips per (date, divergence) cell")

span = np.array([days.min(), days.max()], dtype="float64")
span_dates = span.astype("datetime64[D]")
lines = []
missing = []  # comparisons that could not be drawn; stated on the figure, not dropped
if regression["slope"] is not None:
    per_year = regression["slope"] * node_date_qc.DAYS_PER_YEAR
    r_text = "" if regression["r"] is None else f" (r={regression['r']:.3f})"
    ax.plot(
        span_dates,
        regression["intercept"] + regression["slope"] * span,
        color="crimson",
        lw=2,
        label=f"recomputed: {per_year:.2f} subs/yr{r_text}",
    )
    lines.append("recomputed")

if logged_rate:
    # chronumental logs its slope but not its intercept, so anchor its line at the
    # cloud's own centroid. That makes this a comparison of SLOPES only, which is
    # what the reference is used for; labelled as such so the line is not read as
    # chronumental's actual fitted line.
    mean_x = float((counts * days).sum() / counts.sum())
    mean_y = float((counts * muts).sum() / counts.sum())
    slope_per_day = logged_rate / node_date_qc.DAYS_PER_YEAR
    ax.plot(
        span_dates,
        mean_y + slope_per_day * (span - mean_x),
        color="black",
        lw=1.6,
        ls="--",
        label=(
            f"chronumental pre-fit: {logged_rate:.2f} subs/yr\n"
            f"(slope only, anchored at the centroid)"
        ),
    )
    lines.append("chronumental")
else:
    missing.append(
        f"chronumental's own pre-fit slope unavailable\n"
        f"(no '{node_date_qc.ROOT_TO_TIP_PREFIX}' line at\n{chronumental_stdout_path})"
    )

# The fitted clock's own predicted line. Anchored at zero divergence on the root's
# fitted date -- not at the centroid like the pre-fit slope above -- because under
# chronumental's parameterization this IS the model's prediction rather than a slope
# placed for comparison, so it has a real intercept and moving it would be a
# fabrication. The two anchoring conventions are therefore different on purpose, and
# both are spelled out in the legend.
# Its span starts at the ROOT, not at the cloud's first date, so the anchor point is
# on the figure: the root predates the earliest sample by a few weeks, and a line
# whose stated origin is cropped out invites the reader to extrapolate it themselves.
implied_span = (
    None
    if root_day is None
    else np.array([min(root_day, days.min()), days.max()], dtype="float64")
)
implied = node_date_qc.model_implied_divergence(implied_span, root_day, fitted_rate)
if implied is not None:
    step_text = "" if fitted_step is None else f", step {fitted_step}"
    root_date = np.array([root_day], dtype="float64").astype("datetime64[D]")[0]
    ax.plot(
        implied_span.astype("datetime64[D]"),
        implied,
        color="darkorange",
        lw=2.2,
        ls="-.",
        label=(
            f"chronumental fitted: {fitted_rate:.2f} subs/yr{step_text}\n"
            f"(model-implied, from 0 at the fitted root {root_date})"
        ),
    )
    # The anchor itself, so the line reads as passing through a specific point
    # rather than as a free-floating slope.
    ax.plot([root_date], [0.0], marker="o", ms=6, color="darkorange", zorder=5)
    lines.append("fitted")
elif fitted_rate is None:
    missing.append(
        f"chronumental's fitted clock unavailable\n"
        f"(no '{node_date_qc.FITTED_RATE_STEP_PREFIX}' progress line at\n"
        f"{chronumental_stdout_path})"
    )
else:
    missing.append(
        f"fitted clock {fitted_rate:.2f} subs/yr present, but the root's fitted "
        f"date is\nblank in {os.path.basename(rtt_anchor_path)}: no anchor, so no "
        f"model-implied line"
    )

if missing:
    # Stated on the figure, never silently omitted.
    ax.text(
        0.02,
        0.96,
        "\n".join(missing),
        transform=ax.transAxes,
        va="top",
        fontsize=8,
        color="darkred",
    )

n_tips = int(counts.sum())
ax.set_xlabel("collection date (day-resolution tips only)")
ax.set_ylabel("root-to-tip divergence (mutations)")
ax.set_title(
    f"{tree_label}: root-to-tip divergence vs date\n"
    f"{n_tips:,} tips in {len(cloud):,} exact (date, divergence) cells"
)
ax.legend(loc="upper left", fontsize=8, framealpha=0.9)
ax.grid(alpha=0.25, lw=0.5)
fig.autofmt_xdate()
fig.savefig(plot_path, dpi=150, bbox_inches="tight")
log(
    f"root_to_tip_scatter: {n_tips} tips, {len(cloud)} cells, "
    f"lines drawn: {', '.join(lines) if lines else 'none'} -> {plot_path}"
)
