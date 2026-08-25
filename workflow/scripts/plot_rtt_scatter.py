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

**Two lines, and the second one is the point.** The recomputed regression comes from
the same ``modules/node_date_qc.weighted_linregress`` that produced the number in
the report, so the figure cannot disagree with the report. chronumental's own logged
slope is drawn alongside it -- anchored at the cloud's centroid, because
chronumental logs its slope but not its intercept -- and the two regressions run
over *different tip sets* (day-resolution tips here; every dated tip there), so they
are expected to be close rather than equal. When the logged slope is unavailable the
figure **says so on its face** rather than quietly dropping the line: a plot missing
a comparison it never announces looks exactly like one where the comparison agreed.
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
    # Stated on the figure, never silently omitted.
    ax.text(
        0.02,
        0.96,
        f"chronumental's own pre-fit slope unavailable\n"
        f"(no '{node_date_qc.ROOT_TO_TIP_PREFIX}' line at\n{chronumental_stdout_path})",
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
ax.legend(loc="upper left", fontsize=9, framealpha=0.9)
ax.grid(alpha=0.25, lw=0.5)
fig.autofmt_xdate()
fig.savefig(plot_path, dpi=150, bbox_inches="tight")
log(
    f"root_to_tip_scatter: {n_tips} tips, {len(cloud)} cells, "
    f"lines drawn: {', '.join(lines) if lines else 'none'} -> {plot_path}"
)
