"""Clade subsampling: how far do fitness effects drift as depth falls? (#151)

Reads each arm's per-clade fitness tables and draws the three figures the
analysis turns on: r vs sampling fraction, the count-survival curves that say
how much evidence is left, and the per-(clade, fraction) scatter.

Run with:

    pixi run marimo edit notebooks/subsample_correlations.py

**The population is a fixed ground-truth set defined by full depth alone**:
every mutation whose full-depth row has ``max(actual_count, expected_count) >=
GT_THRESHOLD``. Three properties make that the right target and each of the
alternatives wrong.

*It is identical across fractions*, so the coefficients are comparable. Keying
on the arm's own scored rows instead makes the denominator move with the
sampling fraction, and it selects on the outcome -- keeping exactly the
mutations the arm happened to observe well, which is survivorship.

*Every arm has an estimate for every member.* The grid is enumerated from the
clade founder and does not depend on sampling depth, so a mutation the arm never
observed still carries a ``delta_fitness_ratio``. Dropping those rows discards
the largest part of the cost of thinning.

*It is not the whole grid.* Rows with no evidence on either side are invariant to
the sampling fraction by construction, so correlating over them measures
``build_count_table``'s enumeration and returns a confident number describing
nothing.

**Evidence is reported as a count distribution over the whole grid, not as an
observed/unobserved flag on the scored set.** Whether ``actual_count > 0`` is a
poor summary: a strongly deleterious mutation is *expected* to carry zero counts
at any depth, so the flag conflates lost evidence with a correct prediction of
absence. The survival curves give the whole distribution instead, over **every
single-nucleotide mutation away from the clade founder** -- read at nucleotide
level from ``fitness/ntmut_fitness/``. Taking them over the ground-truth subset
would truncate the curve by construction, since the full-depth max curve passes
through 1.0 at the threshold that defined the subset. Sampling shifts the curve
left, and how far it shifts is the cost of thinning.

So the notebook uses two populations on purpose, and labels which is which: the
**fixed ground-truth set** for the correlations, because a comparison needs a
target that does not move; the **whole grid** for the survival curves, because
that is the space a model must predict over.
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

    # Ground truth: a full-depth row needs this many actual *or* expected counts.
    GT_THRESHOLD = 10
    # Reported points on the count-survival curve.
    COUNT_THRESHOLDS = [1, 3, 10, 30, 100]
    JOIN_KEY = ["cluster", "gene", "aa_mutation"]
    ESTIMATORS = ["delta_fitness_ratio", "delta_fitness_bayesian"]
    return (
        COUNT_THRESHOLDS,
        ESTIMATORS,
        GT_THRESHOLD,
        JOIN_KEY,
        mo,
        np,
        os,
        pd,
        plt,
    )


@app.cell
def _(mo):
    mo.md("""
    # Clade subsampling: fitness effects vs sampling depth

    Each **arm** re-estimates a clade's per-mutation fitness effects from a
    thinned sample list. Ground truth is a **fixed** set of mutations, chosen
    by full depth only, so every arm is scored against the same target and the
    coefficients are comparable across fractions.

    Two things have to be read together:

    - **agreement** — Pearson r and RMSE against the full-depth estimate.
    - **evidence** — the count-survival curve. A nominal sampling fraction is
      not an evidence cut, and the count distribution is what actually moves.
      Read every coefficient against it.
    """)
    return


@app.cell
def _(mo):
    output_dir = mo.ui.text(
        value="results/gisaid-and-public_2026-08-06_v12",
        label="Tree output directory",
        full_width=True,
    )
    output_dir
    return (output_dir,)


@app.cell
def _(mo, os, output_dir, pd):
    subsample_path = os.path.join(output_dir.value, "subsample_summary.csv")
    retention = (
        pd.read_csv(subsample_path) if os.path.isfile(subsample_path) else None
    )
    fitness_dir = os.path.join(output_dir.value, "aamut_fitness")
    mo.stop(
        not os.path.isdir(fitness_dir),
        mo.md(f"**No `aamut_fitness/` under `{output_dir.value}`.**"),
    )
    # Two kinds of arm. ``frac_<unit>_<seed>_<bp>bp`` thins a named unit by a
    # sampling fraction; ``cap_<max_samples>_<seed>`` thins every clade above the
    # cap down to it and leaves smaller clades untouched.
    arms = sorted(
        d
        for d in os.listdir(fitness_dir)
        if os.path.isdir(os.path.join(fitness_dir, d))
        and d.startswith(("frac_", "cap_"))
    )
    mo.stop(not arms, mo.md("**No `frac_*` or `cap_*` arm directories found.**"))

    def arm_kind(arm):
        return "cap" if arm.startswith("cap_") else "frac"

    def arm_depth(arm):
        """The configured depth: a sampling fraction, or a per-clade cap."""
        if arm_kind(arm) == "cap":
            return float(arm.split("_")[1])
        return int(arm.rsplit("_", 1)[-1][:-2]) / 10000.0

    def arm_label(arm):
        if arm_kind(arm) == "cap":
            return f"cap {int(arm_depth(arm)):,}"
        return f"frac {arm_depth(arm):g}"

    def thinned_clades(arm):
        """The clades this arm actually thinned, read from the file tree.

        Deliberately not parsed from the arm name. A cap arm's unit is derived
        from the tree's clade sizes and does not appear in its name at all, and
        ``mat/clades/{arm}/`` holds exactly the clades the sampler drew -- the
        same predicate the pipeline used, rather than a second derivation of it.
        A clade absent here is a control: the arm re-aggregated it but read its
        per-clade files from the full tree unchanged.
        """
        d = os.path.join(output_dir.value, "mat", "clades", arm)
        if not os.path.isdir(d):
            return set()
        return {
            n[: -len("_samples.txt")]
            for n in os.listdir(d)
            if n.endswith("_samples.txt")
        }
    return (
        arm_depth,
        arm_kind,
        arm_label,
        arms,
        fitness_dir,
        retention,
        thinned_clades,
    )


@app.cell
def _(
    COUNT_THRESHOLDS,
    ESTIMATORS,
    GT_THRESHOLD,
    JOIN_KEY,
    arm_depth,
    arm_kind,
    arm_label,
    arms,
    fitness_dir,
    np,
    os,
    pd,
    thinned_clades,
):
    def pearson(x, y):
        """Pearson r, or None where it is undefined rather than merely small."""
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
            return None
        return float(np.corrcoef(x, y)[0, 1])

    def rmse(a, b):
        a = np.asarray(a, dtype=float)
        b = np.asarray(b, dtype=float)
        return float(np.sqrt(((a - b) ** 2).mean()))

    def bias(a, b):
        """Signed mean(arm - full). The direction matters and |.| would hide it."""
        return float((np.asarray(b, dtype=float) - np.asarray(a, dtype=float)).mean())

    def slope(a, b):
        """Least-squares slope of arm on full: 1.0 if unbiased in scale.

        r is invariant to any linear transform of the arm's estimates, so it
        cannot see the shrinkage toward zero that thin evidence produces. This
        can: on the completed fraction arms it falls to 0.477 (21J at p=0.01)
        and 0.217 (21I) while r still reads 0.765 and 0.601.
        """
        a = np.asarray(a, dtype=float)
        b = np.asarray(b, dtype=float)
        if len(a) < 2:
            return None
        var = ((a - a.mean()) ** 2).sum()
        if var == 0:
            return None
        return float(((a - a.mean()) * (b - b.mean())).sum() / var)

    def ground_truth(clade):
        """The fixed target set for one clade, from its full-depth table only."""
        full = pd.read_csv(os.path.join(fitness_dir, f"{clade}_aamut_fitness.csv"))
        keep = full[["actual_count", "expected_count"]].max(axis=1) >= GT_THRESHOLD
        return full[keep]

    def joined(arm, clade):
        """Ground truth joined to the arm's estimate for every one of its rows.

        The row count must not change: the grid does not depend on sampling
        depth, so a shrinking join would mean two different mutation sets were
        being compared and every coefficient would be meaningless.
        """
        gt = ground_truth(clade)
        armd = pd.read_csv(
            os.path.join(fitness_dir, arm, f"{clade}_aamut_fitness.csv")
        )
        m = gt.merge(
            armd, on=JOIN_KEY, suffixes=("_full", "_arm"), validate="one_to_one"
        )
        if len(m) != len(gt):
            raise ValueError(
                f"{arm}/{clade}: ground truth has {len(gt)} rows but the join "
                f"yields {len(m)}; the arm's grid does not cover it"
            )
        return m

    def survival(counts):
        """``{x: P(count >= x)}`` -- the complement of the count eCDF."""
        return {t: float((counts >= t).mean()) for t in COUNT_THRESHOLDS}

    rows = []
    for _arm in arms:
        _unit = thinned_clades(_arm)
        for _name in sorted(os.listdir(os.path.join(fitness_dir, _arm))):
            if not _name.endswith("_aamut_fitness.csv") or _name.startswith("all_"):
                continue
            _clade = _name[: -len("_aamut_fitness.csv")]
            _m = joined(_arm, _clade)
            _row = {
                "arm": _arm,
                "arm_label": arm_label(_arm),
                "kind": arm_kind(_arm),
                "clade": _clade,
                "depth": arm_depth(_arm),
                "thinned": _clade in _unit,
                "n_ground_truth": len(_m),
                "median_count_full": float(_m["actual_count_full"].median()),
                "median_count_arm": float(_m["actual_count_arm"].median()),
                "total_counts_arm": int(_m["actual_count_arm"].sum()),
            }
            for _est in ESTIMATORS:
                _row[f"r_{_est}"] = pearson(_m[f"{_est}_full"], _m[f"{_est}_arm"])
                _row[f"rmse_{_est}"] = (
                    rmse(_m[f"{_est}_full"], _m[f"{_est}_arm"]) if len(_m) else None
                )
                _row[f"bias_{_est}"] = (
                    bias(_m[f"{_est}_full"], _m[f"{_est}_arm"]) if len(_m) else None
                )
                _row[f"slope_{_est}"] = (
                    slope(_m[f"{_est}_full"], _m[f"{_est}_arm"]) if len(_m) else None
                )
            for _t, _v in survival(_m["actual_count_arm"]).items():
                _row[f"frac_count_ge_{_t}"] = _v
            for _t, _v in survival(_m["actual_count_full"]).items():
                _row[f"full_frac_count_ge_{_t}"] = _v
            rows.append(_row)

    summary = pd.DataFrame(rows)
    thinned = summary[summary["thinned"]]
    control = summary[~summary["thinned"]]
    return control, ground_truth, joined, summary, thinned


@app.cell
def _(control, mo, summary, thinned):
    # An arm re-aggregates every clade in its grid but thins only its unit, so a
    # non-unit clade must reproduce full depth exactly. A control r below 1 means
    # the arm reached a clade the config never named.
    bad = control[
        (control["r_delta_fitness_ratio"].fillna(1.0) < 1.0)
        | (control["rmse_delta_fitness_ratio"].fillna(0.0) > 0.0)
    ]
    mo.callout(
        mo.md(
            f"**{summary['arm'].nunique()} arm(s)**, {len(thinned)} thinned "
            f"(clade, depth) pair(s), {len(control)} control pair(s). "
            f"Controls must reproduce full depth exactly: **{len(bad)}** of "
            f"{len(control)} disagree"
            + (
                " — an arm reached a clade the config never named."
                if len(bad)
                else " ✓ the grid-derived and subtree-derived fan-outs stayed separate."
            )
        ),
        kind="danger" if len(bad) else "success",
    )
    return


@app.cell
def _(mo):
    mo.md("""
    ## Pearson r vs sampling depth, over the fixed ground-truth set

    n is constant along each curve by construction — the same mutations are
    scored at every depth, so only the arm-side estimate moves.

    The two kinds of arm are plotted separately because their x-axes are
    different quantities: a `frac` arm's depth is a fraction of *each* clade's
    samples, a `cap` arm's is an absolute per-clade sample count. Reading a
    cap's x-value as a fraction would be wrong by two orders of magnitude, and
    the effective fraction it implies differs per clade — that is the whole
    point of a cap, and it is why the per-clade curves of a cap arm are not
    comparable to one another at a fixed x either. Use
    `subsample_summary.csv`'s measured retention for that comparison.

    RMSE is the **RMSE of the `delta_fitness_ratio` estimator** — the
    root-mean-square difference between the arm's estimate and full depth's, in
    log-fitness units. It is not a ratio of anything. It is reported beside r
    because r is scale-free: it cannot see a systematic shift or the shrinkage
    toward zero that thin evidence produces.
    """)
    return


@app.cell
def _(mo, plt, thinned):
    mo.stop(
        thinned.empty,
        mo.md("**No thinned rows** — no arm subsampled a clade in its unit."),
    )

    _panels = [
        (_kind, _clade)
        for _kind in sorted(thinned["kind"].unique())
        for _clade in sorted(thinned[thinned["kind"] == _kind]["clade"].unique())
    ]
    fig, _axes = plt.subplots(
        1, len(_panels), figsize=(5.6 * len(_panels), 4.4), squeeze=False
    )

    for _col, (_kind, _clade) in enumerate(_panels):
        _g = thinned[
            (thinned["kind"] == _kind) & (thinned["clade"] == _clade)
        ].sort_values("depth")
        _ax = _axes[0][_col]
        for _est, _marker in [
            ("delta_fitness_ratio", "o"),
            ("delta_fitness_bayesian", "s"),
        ]:
            _ax.plot(
                _g["depth"],
                _g[f"r_{_est}"],
                marker=_marker,
                label=_est.replace("delta_fitness_", ""),
            )
        _ax2 = _ax.twinx()
        _ax2.plot(
            _g["depth"],
            _g["rmse_delta_fitness_ratio"],
            marker="x",
            ls="--",
            c="grey",
            label="RMSE of delta_fitness_ratio",
        )
        _ax2.set_ylabel("RMSE, log-fitness units")
        _ax2.set_ylim(0, None)
        _ax.set_title(
            f"{_clade} ({_kind}) — fixed ground truth, "
            f"n = {int(_g['n_ground_truth'].iloc[0]):,}"
        )
        _ax.set_xlabel(
            "per-clade cap (samples)" if _kind == "cap"
            else "nominal sampling fraction"
        )
        _ax.set_ylabel("Pearson r vs full depth")
        _ax.set_xscale("log")
        _ax.set_ylim(0, 1.02)
        _ax.grid(alpha=0.3)
        _ax.legend(fontsize=8, loc="lower right")
        _ax2.legend(fontsize=8, loc="upper left")

    fig.tight_layout()
    fig
    return


@app.cell
def _(mo):
    mo.md("""
    ## How much evidence is there: count survival over the whole nucleotide grid

    Each curve is the fraction of mutations carrying at least *x* counts, at one
    sampling depth. The population here is **every single-nucleotide mutation
    away from the clade founder** — the whole grid, not the ground-truth subset
    the correlations are scored on. Two reasons that is the right population for
    this figure:

    - Restricting to the ground-truth set truncates the curve *by construction*:
      the full-depth max panel would pass through 1.0 at the threshold and the
      interesting left tail would be cut away.
    - The grid is the space a model has to predict over. How much of it carries
      any evidence at all is the question, and most of it does not.

    This is read at nucleotide level from `fitness/ntmut_fitness/`, which is one
    row per single-nt change from the founder. (The amino-acid table used for
    the correlations has *more* rows, because a nucleotide in a gene-overlap
    region contributes an amino-acid mutation to each gene.)

    Three panels, because the ground-truth cut is taken on the max of the two
    and the halves behave differently:

    - **actual counts** — what was observed. Thinning drives mass toward zero.
    - **expected counts** — the neutral-model prediction under each arm's own
      `evol_time` normalization. It scales down smoothly with depth and never
      reaches zero, so it degrades without vanishing.
    - **max(actual, expected)** — the quantity the ground-truth cut uses, so the
      dotted line at *x* = 10 is exactly that threshold.
    """)
    return


@app.cell
def _(arm_depth, arm_label, arms, mo, np, os, output_dir, pd, plt, thinned):
    mo.stop(thinned.empty)

    _NT_DIR = os.path.join(output_dir.value, "fitness", "ntmut_fitness")
    mo.stop(
        not os.path.isdir(_NT_DIR),
        mo.md(f"**No `fitness/ntmut_fitness/` under `{output_dir.value}`.**"),
    )

    _PANELS = [
        ("actual_count", "actual counts"),
        ("expected_count", "expected counts"),
        ("max_count", "max(actual, expected)"),
    ]
    _USE = ["actual_count", "expected_count"]

    def _nt_counts(clade, arm=None):
        """The clade's nucleotide grid at one depth: two count columns only.

        No join is needed -- a survival curve is a marginal distribution, and
        the grid is enumerated from the clade founder so it does not depend on
        sampling depth. The row counts are asserted equal below, which is the
        part a marginal cannot check for itself.
        """
        _d = _NT_DIR if arm is None else os.path.join(_NT_DIR, arm)
        return pd.read_csv(
            os.path.join(_d, f"{clade}_ntmut_fitness.csv"), usecols=_USE
        )

    def _quantity(frame, which):
        if which == "max_count":
            return frame[["actual_count", "expected_count"]].max(axis=1)
        return frame[which]

    def _surv(values, xs):
        _v = np.asarray(values, dtype=float)
        return [float((_v >= _x).mean()) for _x in xs]

    _clades = sorted(thinned["clade"].unique())
    fig3, _axes = plt.subplots(
        len(_clades), 3, figsize=(5.4 * 3, 4.3 * len(_clades)), squeeze=False
    )

    for _row, _clade in enumerate(_clades):
        _full = _nt_counts(_clade)
        _by_depth = [(None, _full)]
        for _arm in sorted(arms, key=arm_depth, reverse=True):
            _a = _nt_counts(_clade, _arm)
            if len(_a) != len(_full):
                raise ValueError(
                    f"{_arm}/{_clade}: nucleotide grid has {len(_a)} rows against "
                    f"{len(_full)} at full depth; the grid must not depend on "
                    f"sampling depth, so these survival curves would have "
                    f"different denominators"
                )
            _by_depth.append((arm_label(_arm), _a))

        for _col, (_which, _label) in enumerate(_PANELS):
            _ax = _axes[_row][_col]
            _hi = max(float(np.nanmax(_quantity(_full, _which))), 2.0)
            _xs = np.geomspace(0.05, _hi, 90)
            for _depth, _frame in _by_depth:
                _y = _surv(_quantity(_frame, _which), _xs)
                if _depth is None:
                    _ax.plot(_xs, _y, lw=2.5, c="black", label="full depth")
                else:
                    _ax.plot(_xs, _y, label=_depth)
            _ax.axvline(10, ls=":", c="grey", lw=1)
            _ax.set_xscale("log")
            _ax.set_xlabel(f"x  ({_label})")
            _ax.set_ylabel("fraction of grid with ≥ x")
            _ax.set_title(f"{_clade} — {_label}, grid n = {len(_full):,}")
            _ax.set_ylim(0, 1.02)
            _ax.grid(alpha=0.3)
            if _col == 0:
                _ax.legend(fontsize=8, title="sampling depth")

    fig3.tight_layout()
    fig3
    return


@app.cell
def _(mo):
    mo.md("""
    ## Per-(clade, fraction) scatter: arm estimate vs full depth

    Every panel holds the **same** mutations — only the arm-side estimate
    changes. Points are coloured by how many observations the arm actually has
    for that mutation, so the shrinking-evidence gradient is visible rather
    than collapsed into a yes/no flag. The band along the bottom is where the
    arm's count has fallen far enough that its estimate is driven by the
    expected count rather than by data.
    """)
    return


@app.cell
def _(joined, np, plt, thinned):
    pairs = (
        thinned[["arm", "arm_label", "clade", "depth"]]
        .drop_duplicates()
        .sort_values(["clade", "depth"])
    )
    _ncol = min(4, len(pairs))
    _nrow = (len(pairs) + _ncol - 1) // _ncol
    fig2, _axes = plt.subplots(
        _nrow, _ncol, figsize=(4.7 * _ncol, 4.3 * _nrow), squeeze=False
    )
    _flat = [_a for _row in _axes for _a in _row]

    # Bins rather than a continuous colourmap: the interesting structure is at
    # 0 and at single digits, which a linear norm renders as one indistinct
    # colour, and a log norm cannot show 0 at all.
    _BINS = [
        (0, 0, "0", "tab:red"),
        (1, 2, "1–2", "tab:orange"),
        (3, 9, "3–9", "tab:olive"),
        (10, None, "≥10", "tab:blue"),
    ]

    for _ax, (_, _r) in zip(_flat, pairs.iterrows()):
        _both = joined(_r["arm"], _r["clade"])
        _counts = _both["actual_count_arm"]
        for _lo, _hi, _label, _colour in _BINS:
            _sel = (_counts >= _lo) if _hi is None else ((_counts >= _lo) & (_counts <= _hi))
            _sub = _both[_sel]
            if _sub.empty:
                continue
            _ax.scatter(
                _sub["delta_fitness_ratio_full"],
                _sub["delta_fitness_ratio_arm"],
                s=4,
                alpha=0.3,
                c=_colour,
                label=f"{_label} ({len(_sub):,})",
                linewidths=0,
            )
        _lo_v = float(
            np.nanmin(
                [
                    _both["delta_fitness_ratio_full"].min(),
                    _both["delta_fitness_ratio_arm"].min(),
                ]
            )
        )
        _hi_v = float(
            np.nanmax(
                [
                    _both["delta_fitness_ratio_full"].max(),
                    _both["delta_fitness_ratio_arm"].max(),
                ]
            )
        )
        _ax.plot([_lo_v, _hi_v], [_lo_v, _hi_v], ls="--", c="black", lw=1)
        _ax.set_title(f"{_r['clade']} @ {_r['arm_label']}   n={len(_both):,}")
        _ax.set_xlabel("full-depth delta_fitness_ratio")
        _ax.set_ylabel("arm delta_fitness_ratio")
        _ax.grid(alpha=0.3)
        _ax.legend(fontsize=7, markerscale=2, title="arm counts", title_fontsize=7)

    for _ax in _flat[len(pairs) :]:
        _ax.axis("off")
    fig2.tight_layout()
    fig2
    return


@app.cell
def _(mo):
    mo.md("""
    ## The numbers behind the figures
    """)
    return


@app.cell
def _(COUNT_THRESHOLDS, mo, thinned):
    cols = (
        ["clade", "kind", "arm_label", "n_ground_truth"]
        + ["r_delta_fitness_ratio", "r_delta_fitness_bayesian"]
        + ["rmse_delta_fitness_ratio", "bias_delta_fitness_ratio"]
        + ["slope_delta_fitness_ratio"]
        + [f"frac_count_ge_{t}" for t in COUNT_THRESHOLDS]
        + ["median_count_arm", "median_count_full"]
    )
    mo.ui.table(
        thinned[cols].sort_values(["clade", "kind", "arm_label"]),
        page_size=20,
    )
    return


@app.cell
def _(mo, retention):
    mo.stop(retention is None, mo.md("_No `subsample_summary.csv` found._"))
    mo.md(
        """
        **Realized retention.** A branch survives whenever *any* sample below it
        is retained, so nucleotide-mutation retention runs above the sample
        fraction and the curve is concave.
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
            ]
        ],
        page_size=20,
    )
    return


if __name__ == "__main__":
    app.run()
