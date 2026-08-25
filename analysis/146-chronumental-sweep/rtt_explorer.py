"""Interactive explorer for the root-to-tip cloud behind chronumental's clock (#146).

Run with:
    marimo edit --headless --port <port> analysis/146-chronumental-sweep/rtt_explorer.py

Reads whatever ``mat/node_dates_rtt.csv`` files exist under ``results/``, so it works
on the test tree and on a production tree without editing anything.

Three lines, matching ``mat/node_dates_rtt_scatter.png``: the recomputed regression,
chronumental's pre-fit slope (anchored at the centroid, since its intercept is not
logged), and chronumental's **fitted** clock drawn as the model's own prediction from
zero divergence on the root's fitted date. The last one is the interesting one -- the
vertical gap between it and the cloud is the #145 duration inflation in mutations.

**The regression is not re-implemented here.** Both the slope and the parse of
chronumental's logged pre-fit rate come from ``workflow/modules/node_date_qc.py`` --
the same functions ``chronumental_qc`` and ``plot_rtt_scatter.py`` call. An explorer
that computed its own slope could disagree with the pipeline's report and the
committed figure while looking authoritative, which is the failure this whole line of
work is about.
"""

import marimo

__generated_with = "0.20.4"
app = marimo.App(width="medium")


@app.cell
def _():
    import glob
    import importlib.util
    import os

    import altair as alt
    import marimo as mo
    import numpy as np
    import pandas as pd

    REPO = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "..")
    ) if "__file__" in globals() else os.getcwd()

    def load_node_date_qc(repo):
        """Import the pipeline's own module by path, so nothing is re-implemented."""
        path = os.path.join(repo, "workflow", "modules", "node_date_qc.py")
        spec = importlib.util.spec_from_file_location("node_date_qc", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    qc = load_node_date_qc(REPO)
    # Altair caps at 5000 rows by default; a production cloud is far bigger and the
    # cap raises rather than sampling, so it is lifted here and density is controlled
    # by the min-tips slider below instead.
    alt.data_transformers.disable_max_rows()
    return REPO, alt, glob, mo, np, os, pd, qc


@app.cell
def _(REPO, glob, mo, os):
    # Two sources, deliberately in one dropdown. Pipeline trees carry the SHIPPED
    # prior (--variance_dates 0.3); the sweep arms under _ignore/ carry the candidate
    # priors, and their node_dates.tsv files are real dated trees, so their clouds are
    # directly comparable. Keeping them in separate notebooks is how you end up
    # comparing a production cloud against a candidate one by memory.
    clouds = sorted(glob.glob(os.path.join(REPO, "results", "*", "mat", "node_dates_rtt.csv")))
    labels = [p.split(os.sep)[-3] for p in clouds]
    arm_clouds = sorted(
        glob.glob(os.path.join(REPO, "_ignore", "*", "*", "qc", "node_dates_rtt.csv"))
    )
    clouds += arm_clouds
    labels += [f"sweep:{p.split(os.sep)[-3]}" for p in arm_clouds]
    tree = mo.ui.dropdown(
        options=dict(zip(labels, clouds)),
        value=labels[-1] if labels else None,
        label="tree",
    )
    mo.md(
        f"### Root-to-tip cloud explorer\n\n"
        f"{len(clouds)} tree(s) with a cloud on disk. A tree missing from this list "
        f"has not had `chronumental_qc` run since #146 landed; refresh it with "
        f"`--rerun-triggers mtime --forcerun chronumental_qc`."
    )
    return clouds, labels, tree


@app.cell
def _(mo, tree):
    mo.hstack([tree], justify="start")
    return


@app.cell
def _(os, pd, tree):
    cloud = pd.read_csv(tree.value, parse_dates=["observed_date"])
    # chronumental's own pre-fit slope lives in the run's log, which is a Snakemake
    # `log:` entry rather than a declared output -- so it may legitimately be absent,
    # and that is reported rather than silently skipped. The two layouts put it in
    # different places: results/<tree>/log/ for a pipeline tree, and the arm directory
    # itself for a sweep arm (the cloud sits one level deeper, in qc/).
    holder = os.path.dirname(tree.value)
    if os.path.basename(holder) == "qc":
        stdout_path = os.path.join(os.path.dirname(holder), "chronumental.stdout")
    else:
        stdout_path = os.path.join(
            os.path.dirname(holder), "log", "chronumental.stdout"
        )
    # The anchor sidecar always sits beside the cloud, in both layouts.
    anchor_path = os.path.join(holder, "node_dates_rtt_anchor.csv")
    return anchor_path, cloud, stdout_path


@app.cell
def _(cloud, mo):
    min_tips = mo.ui.slider(
        1,
        max(2, int(cloud["n_tips"].max())),
        value=1,
        label="min tips per cell (thins the sparse tail for rendering)",
    )
    date_range = mo.ui.range_slider(
        start=0,
        stop=len(sorted(cloud["observed_date"].unique())) - 1,
        value=(0, len(sorted(cloud["observed_date"].unique())) - 1),
        label="date range (index into the observed dates present)",
    )
    log_color = mo.ui.checkbox(True, label="log colour scale")
    mo.vstack([min_tips, date_range, log_color])
    return date_range, log_color, min_tips


@app.cell
def _(anchor_path, cloud, date_range, min_tips, np, os, pd, qc, stdout_path):
    dates_present = np.array(sorted(cloud["observed_date"].unique()))
    lo, hi = dates_present[date_range.value[0]], dates_present[date_range.value[1]]
    view = cloud[
        (cloud["n_tips"] >= min_tips.value)
        & (cloud["observed_date"] >= lo)
        & (cloud["observed_date"] <= hi)
    ].copy()

    def regress(frame):
        if frame.empty:
            return None
        return qc.weighted_linregress(
            frame["observed_date"].to_numpy(dtype="datetime64[D]").astype("float64"),
            frame["rtt_mutations"].to_numpy(dtype="float64"),
            frame["n_tips"].to_numpy(dtype="float64"),
        )

    # Two regressions on purpose: the full cloud (which is what the pipeline reports)
    # and the current view. Filtering changes the slope, and quoting the filtered
    # number as "the" slope is how an explorer misleads.
    reg_all = regress(cloud)
    reg_view = regress(view)
    logged_rate = qc.parse_root_to_tip_rate(stdout_path)
    # The fitted clock and its anchor. A tree whose chronumental_qc predates this
    # sidecar has no anchor file, which is a different state from a blank root date
    # and is reported as such rather than collapsing into "no line".
    fitted_rate, fitted_step = qc.parse_fitted_mutation_rate(stdout_path)
    root_day, anchor_state = None, "no anchor sidecar (re-run chronumental_qc)"
    if os.path.exists(anchor_path):
        anchor = pd.read_csv(anchor_path)
        if len(anchor) and pd.notna(anchor["root_day_ordinal"].iloc[0]):
            root_day = float(anchor["root_day_ordinal"].iloc[0])
            anchor_state = str(anchor["root_fitted_date"].iloc[0])[:10]
        else:
            anchor_state = "sidecar present, root date blank"
    return (
        anchor_state,
        fitted_rate,
        fitted_step,
        logged_rate,
        reg_all,
        reg_view,
        root_day,
        view,
    )


@app.cell
def _(
    anchor_state,
    fitted_rate,
    fitted_step,
    logged_rate,
    mo,
    qc,
    reg_all,
    reg_view,
    stdout_path,
    view,
):
    def per_year(reg):
        if reg is None or reg["slope"] is None:
            return None
        return reg["slope"] * qc.DAYS_PER_YEAR

    full, shown = per_year(reg_all), per_year(reg_view)
    lines = [
        f"- **full cloud**: {full:.2f} subs/yr (r={reg_all['r']:.4f}) over "
        f"{reg_all['n']:.0f} tips" if full else "- full cloud: no slope",
        f"- **current view**: {shown:.2f} subs/yr over {reg_view['n']:.0f} tips in "
        f"{len(view)} cells" if shown else "- current view: no slope (empty or single-date)",
    ]
    if logged_rate:
        lines.append(
            f"- **chronumental's pre-fit**: {logged_rate:.2f} subs/yr; full-cloud "
            f"ratio {full / logged_rate:.3f}"
            + ("  ⚠️ outside the tripwire" if full and abs(full / logged_rate - 1) > qc.RTT_SLOPE_REL_TOL else "")
        )
    else:
        lines.append(
            f"- **chronumental's pre-fit**: unavailable — no parsable "
            f"`{qc.ROOT_TO_TIP_PREFIX}` line at `{stdout_path}`. Stated rather than "
            f"omitted: a missing comparison must not read like a passing one."
        )
    if fitted_rate is None:
        lines.append(
            f"- **chronumental's fitted clock**: unavailable — no parsable "
            f"`{qc.FITTED_RATE_STEP_PREFIX}` progress line at `{stdout_path}`."
        )
    else:
        note = (
            f"- **chronumental's fitted clock**: {fitted_rate:.2f} subs/yr at its "
            f"last logged step {fitted_step}"
        )
        if logged_rate:
            note += f", **{fitted_rate / logged_rate:.2f}x** its own pre-fit slope"
        note += f"; model line anchored at root {anchor_state}"
        lines.append(note)
    mo.md("\n".join(lines))
    return full, shown


@app.cell
def _(
    alt,
    fitted_rate,
    fitted_step,
    full,
    logged_rate,
    mo,
    np,
    qc,
    reg_view,
    root_day,
    view,
):
    points = (
        alt.Chart(view)
        .mark_circle()
        .encode(
            x=alt.X("observed_date:T", title="collection date (day-resolution tips)"),
            y=alt.Y("rtt_mutations:Q", title="root-to-tip divergence (mutations)"),
            color=alt.Color(
                "n_tips:Q",
                scale=alt.Scale(type="log", scheme="viridis"),
                title="tips in cell",
            ),
            size=alt.Size("n_tips:Q", legend=None, scale=alt.Scale(range=[10, 200])),
            tooltip=["observed_date:T", "rtt_mutations:Q", "n_tips:Q"],
        )
        .properties(height=460)
    )

    layers = [points]
    if not view.empty and reg_view["slope"] is not None:
        days = view["observed_date"].to_numpy(dtype="datetime64[D]").astype("float64")
        span = np.array([days.min(), days.max()])
        fit = alt.Chart(
            __import__("pandas").DataFrame(
                {
                    "observed_date": span.astype("datetime64[D]"),
                    "rtt_mutations": reg_view["intercept"] + reg_view["slope"] * span,
                }
            )
        ).mark_line(color="crimson", size=2).encode(x="observed_date:T", y="rtt_mutations:Q")
        layers.append(fit)

        if logged_rate:
            # Slope only, anchored at the view's centroid: chronumental logs its rate
            # but not its intercept, so this compares slopes and nothing else.
            w = view["n_tips"].to_numpy(dtype="float64")
            mean_x = float((w * days).sum() / w.sum())
            mean_y = float((w * view["rtt_mutations"].to_numpy(dtype="float64")).sum() / w.sum())
            slope_day = logged_rate / qc.DAYS_PER_YEAR
            ref = alt.Chart(
                __import__("pandas").DataFrame(
                    {
                        "observed_date": span.astype("datetime64[D]"),
                        "rtt_mutations": mean_y + slope_day * (span - mean_x),
                    }
                )
            ).mark_line(color="black", strokeDash=[6, 4], size=2).encode(
                x="observed_date:T", y="rtt_mutations:Q"
            )
            layers.append(ref)

        # The fitted clock as the model's own prediction. Its span starts at the
        # ROOT rather than at the view's first date, so the point it passes through
        # is on screen; and it is NOT re-anchored when the view is filtered, because
        # unlike the two regressions above this line does not depend on the view.
        implied_span = (
            None
            if root_day is None
            else np.array([min(root_day, span.min()), span.max()], dtype="float64")
        )
        implied = qc.model_implied_divergence(implied_span, root_day, fitted_rate)
        if implied is not None:
            fit_line = alt.Chart(
                __import__("pandas").DataFrame(
                    {
                        "observed_date": implied_span.astype("datetime64[D]"),
                        "rtt_mutations": implied,
                    }
                )
            ).mark_line(color="darkorange", strokeDash=[8, 3], size=2.5).encode(
                x="observed_date:T", y="rtt_mutations:Q"
            )
            layers.append(fit_line)

    mo.vstack(
        [
            mo.md(
                "Red = regression of the **current view**. Black dashed = "
                "chronumental's pre-fit **slope only**, anchored at the view's "
                "centroid (its intercept is not logged). Orange dashed = the "
                "**fitted** clock's own predicted line, from zero divergence on the "
                "root's fitted date — the only one of the three with a real "
                "intercept, and the only one that does not move when you filter. "
                "Drag to pan, scroll to zoom, hover for cell counts."
            ),
            mo.ui.altair_chart(
                alt.layer(*layers).interactive(), chart_selection=False, legend_selection=False
            ),
        ]
    )
    return


if __name__ == "__main__":
    app.run()
