"""Clade subsampling: how far do fitness effects drift as depth falls? (#151)

Reads `correlation_summary.csv` and the per-clade fitness tables an arm produced,
and draws the two figures the analysis turns on: the per-(clade, fraction)
scatter of arm effect against full depth, and the rho-vs-fraction curves.

Run with:

    pixi run marimo edit notebooks/subsample_correlations.py

Point `OUTPUT_DIR` at any tree's output directory. It defaults to the test
tree, where -- deliberately -- the answer is that there is no answer: every
scored mutation carries exactly one observation, so the coefficients are
arithmetic rather than agreement. The notebook says so rather than drawing a
reassuring line through it.
"""

import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import os

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    return mo, np, os, pd, plt


@app.cell
def _(mo):
    mo.md("""
    # Clade subsampling: fitness effects vs sampling depth

    Each **arm** re-estimates a clade's per-mutation fitness effects from a
    thinned sample list. This notebook asks how far those estimates drift
    from the full-depth ones, and — just as importantly — how many
    mutations survive to be compared at all.

    Two numbers are needed to read any of it:

    - **`n_mutations`** — the count filter raises rho partly by dropping
      hard cases, so a coefficient without its n is not interpretable.
    - **`coverage`** — rho is computed on the intersection of *scored*
      mutations, so it says nothing about mutations that vanish entirely.
    """)
    return


@app.cell
def _(mo):
    output_dir = mo.ui.text(
        value="results/test_data_v12",
        label="Tree output directory",
        full_width=True,
    )
    output_dir
    return (output_dir,)


@app.cell
def _(mo, os, output_dir, pd):
    summary_path = os.path.join(output_dir.value, "correlation_summary.csv")
    mo.stop(
        not os.path.isfile(summary_path),
        mo.md(
            f"**No `correlation_summary.csv` at `{summary_path}`.** "
            "Build it with `snakemake <outdir>/correlation_summary.csv`, or "
            "point the box above at a tree that has arms."
        ),
    )
    corr = pd.read_csv(summary_path)
    subsample_path = os.path.join(output_dir.value, "subsample_summary.csv")
    retention = (
        pd.read_csv(subsample_path) if os.path.isfile(subsample_path) else None
    )
    return corr, retention


@app.cell
def _(corr, mo):
    # The single most important thing this notebook can say. When every scored
    # mutation carries one observation, delta_fitness_ratio is a deterministic
    # function of the expected count, the two columns are collinear by
    # construction, and a scatter of them is a straight line that means nothing.
    # Saying so up front is the point -- a reader who scrolls to the figures
    # first should not be able to miss it.
    max_count = int(corr["max_actual_count_full"].max())
    banner = (
        mo.callout(
            mo.md(
                f"**These coefficients do not measure subsampling noise.** The "
                f"highest full-depth `actual_count` anywhere in this table is "
                f"**{max_count}**. With one observation per scored mutation, "
                f"`delta_fitness_ratio` is a deterministic function of the "
                f"expected count, so any correlation it shows is arithmetic. "
                f"Use a production tree, where a clade has thousands of scored "
                f"mutations at counts much greater than 1."
            ),
            kind="danger",
        )
        if max_count <= 1
        else mo.callout(
            mo.md(
                f"Highest full-depth `actual_count`: **{max_count}**. "
                f"Coefficients below are interpretable to the extent that this "
                f"is comfortably above 1."
            ),
            kind="success",
        )
    )
    banner
    return


@app.cell
def _(corr, mo):
    thinned = corr[corr["thinned"]]
    control = corr[~corr["thinned"]]
    # Control clades are not padding: an arm re-aggregates every clade in its
    # grid but thins only its unit, so a non-unit clade must reproduce full
    # depth exactly. A control rmse above 0 means the arm reached a clade the
    # config never named -- the fan-out bug this feature had to fix once.
    control_bad = control[control["rmse"].fillna(0) > 0]
    mo.md(
        f"""
        **{corr['arm'].nunique()} arm(s)**, {len(thinned)} thinned row(s) and
        {len(control)} control row(s).

        Control clades (not in the arm's unit) must reproduce full depth
        exactly: **{len(control_bad)}** of {len(control)} disagree
        {'— a leak into a clade the config never named' if len(control_bad) else '✓'}.
        """
    )
    return (thinned,)


@app.cell
def _(mo):
    mo.md("""
    ## rho vs sampling fraction, faceted by estimator
    """)
    return


@app.cell
def _(mo, np, plt, thinned):
    mo.stop(
        thinned.empty,
        mo.md("**No thinned rows** — no arm subsampled a surviving clade."),
    )

    estimators = sorted(thinned["estimator"].unique())
    fig, axes = plt.subplots(
        2, len(estimators), figsize=(5.5 * len(estimators), 8), squeeze=False
    )

    for col, est in enumerate(estimators):
        top, bottom = axes[0][col], axes[1][col]
        sub = thinned[thinned["estimator"] == est]
        for filt, grp in sub.groupby("filter"):
            for clade, g in grp.groupby("clade"):
                g = g.sort_values("fraction")
                # Rows with no mutations are plotted as gaps, not as zeros: an
                # empty population is an absent measurement, and joining
                # through it would draw a curve descending to zero correlation
                # that was never measured.
                ok = g["n_mutations"] > 0
                top.plot(
                    g.loc[ok, "fraction"],
                    g.loc[ok, "spearman_rho"],
                    marker="o",
                    label=f"{clade} · {filt}",
                )
                for _, r in g[ok].iterrows():
                    top.annotate(
                        f"n={int(r['n_mutations'])}",
                        (r["fraction"], r["spearman_rho"]),
                        fontsize=7,
                        xytext=(3, 3),
                        textcoords="offset points",
                    )
        top.set_title(f"{est}\nSpearman rho vs full depth")
        top.set_xlabel("nominal sampling fraction")
        top.set_ylabel("spearman rho")
        top.set_xscale("log")
        top.grid(alpha=0.3)
        top.legend(fontsize=7)

        # Coverage is per (clade, fraction, filter) and does not depend on the
        # estimator, so the bottom row repeats across columns by design -- it
        # sits under each estimator so a reader comparing two panels never has
        # to look away to find how many mutations the coefficient rests on.
        cov = sub[sub["filter"] == "unfiltered"]
        for clade, g in cov.groupby("clade"):
            g = g.sort_values("fraction")
            bottom.plot(g["fraction"], g["coverage"], marker="s", label=clade)
        bottom.plot([0, 1], [0, 1], ls=":", c="grey", lw=1, label="y = x")
        bottom.set_title("coverage: scored mutations the arm still scores")
        bottom.set_xlabel("nominal sampling fraction")
        bottom.set_ylabel("coverage")
        bottom.set_ylim(0, 1.05)
        bottom.grid(alpha=0.3)
        bottom.legend(fontsize=7)

    fig.tight_layout()
    _ = np  # keeps the import honest if the annotate loop is edited away
    fig
    return


@app.cell
def _(mo):
    mo.md("""
    ## Per-(clade, fraction) scatter: arm effect vs full depth

    The diagonal is perfect agreement. Scatter *off* the diagonal is the
    drift; mutations missing from the plot entirely are the coverage loss,
    which a scatter cannot show — read it off the panel above.
    """)
    return


@app.cell
def _(mo, os, output_dir, pd, thinned):
    pairs = (
        thinned[["arm", "clade", "fraction"]]
        .drop_duplicates()
        .sort_values(["clade", "fraction"])
    )

    def load_pair(arm, clade):
        """Full-depth and arm fitness tables joined on the invariant key.

        ``clade_founder_aa`` is invariant to sampling fraction, so the same
        ``aa_mutation`` string means the same substitution on both sides.
        """
        d = os.path.join(output_dir.value, "aamut_fitness")
        full = pd.read_csv(os.path.join(d, f"{clade}_aamut_fitness.csv"))
        arm_df = pd.read_csv(os.path.join(d, arm, f"{clade}_aamut_fitness.csv"))
        m = full.merge(
            arm_df,
            on=["cluster", "gene", "aa_mutation"],
            suffixes=("_full", "_arm"),
        )
        # Scored at both depths. Every grid row carries a fitness value even
        # when nothing was observed, and those rows are invariant to sampling
        # fraction -- plotting them would fill the diagonal with points that
        # agree by construction and hide the ones that do not.
        return m[(m["actual_count_full"] > 0) & (m["actual_count_arm"] > 0)]

    mo.stop(pairs.empty, mo.md("**Nothing thinned to scatter.**"))
    return load_pair, pairs


@app.cell
def _(load_pair, pairs, plt):
    n = len(pairs)
    ncol = min(3, n)
    nrow = (n + ncol - 1) // ncol
    fig2, axes2 = plt.subplots(
        nrow, ncol, figsize=(4.6 * ncol, 4.4 * nrow), squeeze=False
    )
    flat = [ax for row in axes2 for ax in row]

    for ax, (_, row) in zip(flat, pairs.iterrows()):
        both = load_pair(row["arm"], row["clade"])
        # ``_est`` rather than ``est``: marimo requires a variable to be
        # defined in exactly one cell, and the rho-vs-fraction cell above
        # already binds ``est``. The underscore marks it cell-private.
        for _est, _colour in [
            ("delta_fitness_ratio", "tab:blue"),
            ("delta_fitness_bayesian", "tab:orange"),
        ]:
            ax.scatter(
                both[f"{_est}_full"],
                both[f"{_est}_arm"],
                s=18,
                alpha=0.6,
                c=_colour,
                label=_est.replace("delta_fitness_", ""),
            )
        lo = min(
            both[["delta_fitness_ratio_full", "delta_fitness_bayesian_full"]]
            .min()
            .min(),
            both[["delta_fitness_ratio_arm", "delta_fitness_bayesian_arm"]]
            .min()
            .min(),
        )
        hi = max(
            both[["delta_fitness_ratio_full", "delta_fitness_bayesian_full"]]
            .max()
            .max(),
            both[["delta_fitness_ratio_arm", "delta_fitness_bayesian_arm"]]
            .max()
            .max(),
        )
        ax.plot([lo, hi], [lo, hi], ls="--", c="grey", lw=1)
        ax.set_title(
            f"{row['clade']} @ {row['fraction']:g}\nn={len(both)} scored in both"
        )
        ax.set_xlabel("full-depth effect")
        ax.set_ylabel("arm effect")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)

    for ax in flat[len(pairs) :]:
        ax.axis("off")
    fig2.tight_layout()
    fig2
    return


@app.cell
def _(mo):
    mo.md("""
    ## The table itself, and the evidence it is keyed to
    """)
    return


@app.cell
def _(corr, mo):
    mo.ui.table(
        corr[
            [
                "arm",
                "clade",
                "fraction",
                "thinned",
                "estimator",
                "filter",
                "n_mutations",
                "coverage",
                "pearson_r",
                "spearman_rho",
                "rmse",
            ]
        ],
        page_size=20,
    )
    return


@app.cell
def _(mo, retention):
    mo.stop(retention is None, mo.md("_No `subsample_summary.csv` found._"))
    mo.md(
        """
        **Realized retention.** A nominal fraction is not an evidence cut: a
        branch survives whenever *any* sample below it is retained, so mutation
        retention runs above the sample fraction and the curve is concave. Key
        every coefficient above to the measured columns here, never to the
        configured percent.
        """
    )
    return


@app.cell
def _(mo, retention):
    mo.stop(retention is None)
    mo.ui.table(
        retention[retention["gene"].isna()][
            [
                "arm",
                "clade",
                "sampling_fraction",
                "included",
                "n_samples_full",
                "n_samples_arm",
                "n_nt_mutations_full",
                "n_nt_mutations_arm",
                "nt_mutation_retention",
                "predicted_nt_mutation_retention",
                "nt_mutation_deviation_sd",
            ]
        ],
        page_size=20,
    )
    return


if __name__ == "__main__":
    app.run()
