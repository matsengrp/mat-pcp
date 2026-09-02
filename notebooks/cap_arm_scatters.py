"""Per-clade cap arms: the statistics that say what a cap costs (#151)

A **cap** arm thins every clade above ``max_samples`` down to exactly that many
samples and leaves smaller clades untouched, so one arm spans the whole depth
curve: on the v12 gisaid tree ``max_samples: 50000`` puts 21J at an effective
1.15% of its samples and 26 smaller clades at 100%.

The notebook opens tree-wide, since what a cap does depends entirely on how
unevenly the clades are sampled: **sequences per clade** against the cap lines,
coloured by when each clade appeared, and then **the total sequences and
mutation counts across all clades** at each cap -- the size of the training set
a DASM would actually see. Controls count in that total: a cap leaves every
clade at or below it untouched, so the tree total falls far less steeply than
its largest clade's does.

Everything after that is one clade at a time, three things per clade:

1. **The count fraction** -- ``sum(actual_count_arm) / sum(actual_count_full)``.
   How much of the clade's observed mutational evidence the arm still has.
2. **The eCDF of per-mutation counts** at each depth, over the clade's whole
   nucleotide grid. How that evidence is distributed, which a single fraction
   cannot show.
3. **The correlation with the ground-truth effects**, drawn as a scatter with
   its fitted line, not only reported as an r.

Run with:

    pixi run marimo edit notebooks/cap_arm_scatters.py

**The count fraction sums observations, not distinct mutations.** A fraction of
*mutations still scored at all* counts each one once however much evidence it
lost, so a mutation falling from 900 observations to 3 reads as fully retained.
That reading is not reported anywhere: it diverges from the count fraction by up
to 30x on real data (21J at cap 5,000, 0.0033 against 0.103) and predicts the
fitted slope worse (Spearman 0.940 against 0.985 over the 120 completed cap
fits). ``n_scored_full`` / ``n_scored_arm`` are still in
``correlation_summary.csv`` as the raw deletion counts, since a mutation that
loses its last observation is *deleted* rather than degraded -- but there is one
retention fraction, not two.

**Read the scatter, not just the r.** Each panel's axes are square, symmetric
and shared across the row, and the dashed line is y = x. A well-recovered clade
lies along it; a thinned one lies along a *shallower* line through the origin,
which is shrinkage rather than noise. Pearson r cannot see that -- it is
invariant to any linear transform of the arm's estimates -- so the fitted slope
is annotated on every panel beside its r. 21J at cap 5,000 reads
r = 0.43 with a slope of 0.12: full depth spans [-7.8, +7.3] where the arm
spans [-2.3, +3.6].

**The correlation population is the fixed ground-truth set**: every mutation
whose **full-depth** row has ``max(actual_count, expected_count) >=
GT_THRESHOLD``. It is therefore the same set at every cap, which is what makes
the panels comparable down a column. Two alternatives are both wrong. Selecting
on the *arm's* counts drops exactly the mutations the thinning destroyed, so the
population gets easier as the cut gets harsher and a degrading estimator draws a
stable-looking panel. And correlating over the whole grid measures
``build_count_table``'s enumeration: most grid rows carry a
``delta_fitness_ratio`` with no evidence behind it on either side, invariant to
sampling depth by construction.

**The eCDFs deliberately use the whole grid, and the nucleotide one.** The
ground-truth subset would truncate the curve by construction -- the full-depth
max curve passes through 1.0 at the very threshold that defined the subset -- so
the distribution is taken over every single-nucleotide mutation away from the
clade founder, read from ``fitness/ntmut_fitness/``. That is the finer grid: the
count fraction measured there agrees with the amino-acid one to better than
0.002 on every (arm, clade) here, the two differing only by double-counting at
gene overlaps.

**Density, not 76,000 overplotted dots.** 21J's ground-truth set is 76,776
mutations, where a raw scatter is an ink blob whose shape is set by draw order.
The panels are hexbin on a log count scale, so the mass of the distribution is
visible rather than only its outline; small clades (under ``SCATTER_MAX_N``) are
drawn as translucent points instead, since a hexbin of 200 points reads as
noise.
"""

import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import os
    import re
    from datetime import date

    import marimo as mo
    import matplotlib as mpl
    import matplotlib.patheffects as path_effects
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from matplotlib.colors import LinearSegmentedColormap, LogNorm

    # Ground truth: a full-depth row needs this many actual *or* expected counts.
    GT_THRESHOLD = 10
    AA_JOIN_KEY = ["cluster", "gene", "aa_mutation"]
    ESTIMATORS = ["delta_fitness_ratio", "delta_fitness_bayesian"]
    # Below this many points a hexbin reads as noise; draw translucent dots.
    SCATTER_MAX_N = 2000
    # A coefficient from fewer than this many mutations has no resolution --
    # 20F at cap 5,000 reads r = 0.999 from six of them. Flagged, not hidden.
    SMALL_N = 500

    # --- Palette -----------------------------------------------------------
    # The caps are an *ordered* scale, not a set of identities, so they take one
    # hue light->dark rather than several categorical hues: bigger cap (more
    # data retained) = darker. Steps 250/450/650 of the blue ramp, validated for
    # the light surface as an ordinal ramp (monotone lightness, every adjacent
    # gap >= 0.06, light end 2.06:1 against the surface, hue spread 3 degrees).
    # The wider 100->700 range is used for the hexbin, where the lightest step
    # means "near zero" and may recede toward the surface.
    BLUE_ORDINAL = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
                    "#2a78d6", "#256abf", "#1c5cab", "#184f95",
                    "#104281", "#0d366b"]
    BLUE_SEQUENTIAL = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5",
                       "#256abf", "#184f95", "#0d366b"]
    ACCENT = "#eb6834"       # categorical slot 2: the fitted line
    SURFACE = "#fcfcfb"
    INK = "#0b0b0b"
    INK_2 = "#52514e"        # full depth, and the y = x reference
    MUTED = "#898781"
    GRID = "#e1e0d9"
    BASELINE = "#c3c2b7"

    DENSITY_CMAP = LinearSegmentedColormap.from_list(
        "blue_seq", BLUE_SEQUENTIAL
    )
    # Clade emergence dates are ordered, so the ramp is monotone in lightness
    # light->dark; but unlike the caps it is *multi-hue* (amber -> pink ->
    # magenta -> violet -> navy, five steps of plasma read from t=0.75 down to
    # 0). A single blue hue has only so much perceptual range to spend, and
    # spread over six years it left 2022, 2023 and 2024 reading as the same
    # mid-blue -- a scale that is ordered and unreadable. Hue does the
    # separating here and lightness still does the ordering, so this is a
    # sequential ramp rather than a rainbow: it passes the ordinal checks
    # (monotone L, every adjacent gap >= 0.06, light end #f89540 at 2.18:1
    # against the surface) and plasma is CVD-safe by construction. Truncated
    # at t=0.75 rather than running to plasma's yellow, which reads 1.5:1
    # against a near-white surface -- a 34-point dot needs 2:1 to be seen at
    # all. Where the caps are colour-coded (the eCDFs, the tree totals) they
    # keep the single-hue BLUE_ORDINAL, so the two scales cannot be confused;
    # on the figure that carries the dates the caps are neutral lines.
    DATE_CMAP = LinearSegmentedColormap.from_list(
        "date_seq",
        ["#f89540", "#da5a6a", "#aa2395", "#6600a7", "#0d0887"],
    )

    mpl.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "axes.edgecolor": BASELINE,
        "axes.labelcolor": INK_2,
        "axes.titlecolor": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelcolor": MUTED,
        "ytick.labelcolor": MUTED,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "text.color": INK,
        "font.size": 9,
        "axes.titlesize": 10,
        "legend.frameon": False,
    })
    return (
        AA_JOIN_KEY,
        ACCENT,
        BLUE_ORDINAL,
        DATE_CMAP,
        DENSITY_CMAP,
        ESTIMATORS,
        GRID,
        GT_THRESHOLD,
        INK,
        INK_2,
        LogNorm,
        MUTED,
        SCATTER_MAX_N,
        SMALL_N,
        SURFACE,
        date,
        mo,
        np,
        os,
        path_effects,
        pd,
        plt,
        re,
    )


@app.cell
def _(mo):
    mo.md("""
    # Per-clade cap arms

    A `cap` arm thins **every** clade above `max_samples` down to it. Membership
    is derived from the tree's clade sizes, not named in the config, so the arm
    directory itself is the roster: `mat/clades/{arm}/` holds exactly the clades
    the sampler drew, and a clade absent from it is a **control** that the arm
    re-aggregated while reading its per-clade files from the full tree unchanged.
    Controls must reproduce full depth exactly, so they are excluded from the
    per-clade figures — their scatter is the identity line by construction and
    their count fraction is 1. They are *included* in the tree-level totals,
    where they are most of what survives a cap.
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
def _(mo, os, output_dir):
    aa_dir = os.path.join(output_dir.value, "aamut_fitness")
    nt_dir = os.path.join(output_dir.value, "fitness", "ntmut_fitness")
    mo.stop(
        not os.path.isdir(aa_dir),
        mo.md(f"**No `aamut_fitness/` under `{output_dir.value}`.**"),
    )
    cap_arms = sorted(
        (d for d in os.listdir(aa_dir)
         if d.startswith("cap_")
         and os.path.isdir(os.path.join(aa_dir, d))),
        key=lambda d: int(d.split("_")[1]),
    )
    mo.stop(
        not cap_arms,
        mo.md(
            f"**No `cap_*` arm directories under `{aa_dir}`.** This notebook "
            "covers the cap kind only; `notebooks/subsample_correlations.py` "
            "covers both kinds."
        ),
    )

    def cap_of(arm):
        """The configured per-clade cap, read off the arm name."""
        return int(arm.split("_")[1])

    def thinned_clades(arm):
        """The clades this arm thinned, read from the file tree.

        Deliberately not derived a second time from the clade sizes. A cap
        arm's unit is resolved by ``subsample.resolve_unit`` against the
        processed roster, and ``mat/clades/{arm}/`` is the record of what that
        resolution actually drew -- so reading it here cannot disagree with the
        pipeline the way a reimplemented ``n > cap`` predicate could.
        """
        d = os.path.join(output_dir.value, "mat", "clades", arm)
        if not os.path.isdir(d):
            return set()
        return {
            n[: -len("_samples.txt")]
            for n in os.listdir(d)
            if n.endswith("_samples.txt")
        }

    def full_depth_size(clade):
        p = os.path.join(
            output_dir.value, "mat", "clades", f"{clade}_samples.txt"
        )
        if not os.path.isfile(p):
            return None
        with open(p) as fh:
            return sum(1 for _ in fh)
    return aa_dir, cap_arms, cap_of, full_depth_size, nt_dir, thinned_clades


@app.cell
def _(mo):
    mo.md("""
    ## Setting the stage: sequences per clade

    A cap intervenes on the tree's **composition**, so what it does depends
    entirely on how unevenly the clades are sampled — and they are very
    unevenly sampled. Each horizontal line below is one cap: every clade
    plotted *above* it is cut back **down to** that line, and every clade below
    it is left alone as a byte-identical control.

    This is the whole difference from a `frac` arm. A fraction takes the same
    share from every clade in its unit, so the imbalance below survives the
    thinning and one arm is one point on a depth curve. A cap flattens the
    right-hand tail, which is what the pooled-vs-per-clade question is actually
    about — and it means a single arm delivers effective sampling fractions
    spanning two orders of magnitude.

    **Dot colour is when the clade appeared** — the 1st percentile of its
    samples' reported collection dates, which is a proxy and is deliberately
    neither the minimum (some records carry a submission date in the
    collection field: 21J's earliest reads 2020-03-02, more than a year before
    Delta) nor the founder node's fitted date (the clade subtrees renumber
    their nodes, so recovering it means walking the whole tree). It is a
    context channel, not an input to any statistic below. Read together with
    the heights it says which cohorts a cap actually spends: the clades a cap
    cuts hardest are the ones that were dominant when sequencing was densest,
    not the oldest ones.
    """)
    return


@app.cell
def _(
    DATE_CMAP,
    GRID,
    INK,
    INK_2,
    MUTED,
    SURFACE,
    cap_arms,
    cap_of,
    date,
    mo,
    np,
    os,
    output_dir,
    path_effects,
    plt,
    re,
):
    _clade_dir0 = os.path.join(output_dir.value, "mat", "clades")
    _date_re0 = re.compile(r"\d{4}-\d{2}-\d{2}")

    def _clade_size_and_emergence(clade):
        """Sample count, and a robust first-appearance date, in one pass.

        The date is the **1st percentile** of the reported collection dates
        already embedded in the sample names (``.../2021|EPI_ISL_...|
        2021-07-20``), not the minimum and not the clade founder node's fitted
        date.

        Not the minimum, because a small number of records carry a submission
        date in the collection field: 21J's earliest reported date here is
        2020-03-02 and 24A's is 2020-01-10, for lineages that did not exist
        until mid-2021 and late 2023 -- the same misdating
        ``pcp_max_predicted_days_past_cutoff`` exists to catch (#149). At 1% of
        4.3 million samples the quantile sits behind ~43,000 records, so a
        handful of placeholders cannot move it, and it lands where the lineage
        is known to have risen (21J 2021-06, 21L 2022-01, 24A 2023-11).

        Not the founder node's own chronumental date, which would be the
        principled quantity: the clade subtrees renumber their nodes
        (``mutations/{clade}_sample_paths.tsv`` starts at ``node_2`` for every
        clade), so recovering the founder's full-tree name means walking the
        8-million-tip tree, which is a pipeline step rather than a notebook
        one. This is a colour channel setting context for the cut, not an
        input to any statistic.

        Counted by value rather than sorted: dates take ~2,000 distinct values
        against 17 million lines, so the quantile comes off a cumulative sum
        over the keys.
        """
        p = os.path.join(_clade_dir0, f"{clade}_samples.txt")
        n = 0
        seen = {}
        with open(p) as fh:
            for line in fh:
                n += 1
                d = line.rstrip("\n")[-10:]
                if _date_re0.fullmatch(d):
                    seen[d] = seen.get(d, 0) + 1
        total = sum(seen.values())
        if not total:
            return n, None
        want = max(1, int(0.01 * total))
        running = 0
        for d in sorted(seen):
            running += seen[d]
            if running >= want:
                return n, d
        return n, None

    _measured0 = {
        _n[: -len("_samples.txt")]: _clade_size_and_emergence(
            _n[: -len("_samples.txt")]
        )
        for _n in os.listdir(_clade_dir0)
        if _n.endswith("_samples.txt")
    }
    _sizes0 = {_k: _v[0] for _k, _v in _measured0.items() if _v[0]}
    _emerged0 = {_k: _v[1] for _k, _v in _measured0.items() if _v[0]}
    mo.stop(not _sizes0, mo.md("**No `mat/clades/*_samples.txt` found.**"))

    # Largest first, left to right.
    _order0 = sorted(_sizes0, key=lambda _c: -_sizes0[_c])
    _x0 = np.arange(len(_order0))
    _vals0 = [_sizes0[_c] for _c in _order0]
    _caps0 = sorted(cap_of(_a) for _a in cap_arms)

    # Colour is the emergence date -- an ordered quantity, so DATE_CMAP is
    # monotone in lightness; see its definition for why it is multi-hue where
    # the caps are single-hue. Encoded as day ordinals so the scale is linear
    # in time, with the colourbar ticked at calendar years. An undated clade is
    # drawn in muted ink rather than at either end of the ramp, which would
    # assert a date.
    _ords0 = [
        date.fromisoformat(_emerged0[_c]).toordinal() if _emerged0[_c] else None
        for _c in _order0
    ]
    _known0 = [_i for _i, _o in enumerate(_ords0) if _o is not None]
    _unknown0 = [_i for _i, _o in enumerate(_ords0) if _o is None]
    _lo0 = min(_ords0[_i] for _i in _known0)
    _hi0 = max(_ords0[_i] for _i in _known0)

    # A *column* dot plot. Horizontal was the obvious orientation and it needed
    # 0.23 in per clade -- 15 inches of scroll for 56 of them, most of it the
    # whitespace a horizontal tick label demands. Clade names are 3 characters,
    # so rotating them 90 degrees lets the points sit ~0.16 in apart and the
    # whole figure fit in 4.6 in of height while every label is still drawn.
    # (Tooltips would compress it further, but they need an interactive backend
    # and would lose the labels in the HTML export, which is how this gets
    # shared.) Dots rather than bars because a log axis has no zero to anchor a
    # bar's length to; the hairline stem is a guide to the axis, not a
    # magnitude.
    _fig0, _ax0 = plt.subplots(figsize=(0.17 * len(_order0) + 2.8, 4.6))
    _ax0.vlines(
        _x0, [min(_vals0) * 0.5] * len(_x0), _vals0,
        color=GRID, linewidth=1.0, zorder=1,
    )
    _sc0 = _ax0.scatter(
        [_x0[_i] for _i in _known0],
        [_vals0[_i] for _i in _known0],
        c=[_ords0[_i] for _i in _known0],
        cmap=DATE_CMAP, vmin=_lo0, vmax=_hi0,
        s=34, edgecolors=SURFACE, linewidths=1.1, zorder=4,
    )
    if _unknown0:
        _ax0.scatter(
            [_x0[_i] for _i in _unknown0],
            [_vals0[_i] for _i in _unknown0],
            s=34, color=MUTED, edgecolors=SURFACE, linewidths=1.1, zorder=4,
        )
    # The caps are reference lines, not data, so they are neutral and directly
    # annotated -- a legend would spend a colour channel that now carries the
    # dates, and three entries in a box are further from the line than the
    # label on it.
    for _cap0 in _caps0:
        _ax0.axhline(
            _cap0, linestyle="--", linewidth=1.4, color=INK_2, zorder=2,
        )
        _ax0.annotate(
            f"cap {_cap0:,}",
            (len(_order0) - 0.4, _cap0),
            textcoords="offset points", xytext=(0, 3),
            ha="right", va="bottom", fontsize=8, color=INK_2,
            path_effects=[
                path_effects.withStroke(linewidth=2.6, foreground=SURFACE)
            ],
            zorder=6,
        )
    # Selective labels only: the extremes, which set the range the caption
    # quotes. A number on all 56 points is unreadable at this width, and the
    # tree-total table below carries the per-cap figures that matter.
    for _i0 in (0, len(_order0) - 1):
        _ax0.annotate(
            f"{_vals0[_i0]:,}",
            (_x0[_i0], _vals0[_i0]),
            textcoords="offset points", xytext=(0, 8),
            ha="center", fontsize=7.5, color=INK_2,
            path_effects=[
                path_effects.withStroke(linewidth=2.4, foreground=SURFACE)
            ],
        )
    _ax0.set_xticks(_x0)
    _ax0.set_xticklabels(
        _order0, fontsize=7, color=INK_2, rotation=90,
    )
    _ax0.set_xlim(-0.8, len(_order0) - 0.2)
    _ax0.set_xlabel("clade")
    _ax0.tick_params(axis="x", length=0, pad=2)
    _ax0.set_yscale("log")
    _ax0.set_ylim(
        bottom=min(_sizes0.values()) * 0.45, top=max(_sizes0.values()) * 3.2
    )
    _ax0.set_ylabel("full-depth samples (log)")
    _ax0.grid(True, axis="y", color=GRID, linewidth=0.6, zorder=0)
    _ax0.set_axisbelow(True)
    _ax0.set_title("Sequences per clade", color=INK, pad=10)

    _cb0 = _fig0.colorbar(_sc0, ax=_ax0, pad=0.02, fraction=0.05)
    _years0 = range(
        date.fromordinal(_lo0).year, date.fromordinal(_hi0).year + 1
    )
    _ticks0 = [
        date(_y, 1, 1).toordinal()
        for _y in _years0
        if _lo0 <= date(_y, 1, 1).toordinal() <= _hi0
    ]
    _cb0.set_ticks(_ticks0)
    _cb0.set_ticklabels([str(date.fromordinal(_t).year) for _t in _ticks0])
    _cb0.set_label("clade emergence\n(1st pct of collection dates)")
    _cb0.outline.set_visible(False)
    _fig0.tight_layout()
    _fig0
    return


@app.cell
def _(mo):
    mo.md("""
    ## The whole training set at each cap

    The plot above is per clade at full depth; this is what the caps add up to.
    What a DASM actually trains on is the tree, so this is the tree-level total
    at each cap: how many sequences survive across **all** clades, and how many
    mutation counts they carry.

    Controls are included and matter here — a cap leaves every clade at or below
    it untouched, so those clades contribute their *full* sequences and counts to
    the total. That is the difference between this and the per-clade statistics
    that follow: 21J losing 98.85% of its samples at cap 50,000 does not mean
    the tree lost 98.85%, because 26 clades lost nothing. Both halves are broken
    out.

    The sequence totals are read from the sample lists — the arm's own
    `mat/clades/{arm}/{clade}_samples.txt` where it thinned, the full tree's
    otherwise — which is the same source the sampler wrote, so this cannot
    disagree with what was actually extracted.
    """)
    return


@app.cell
def _(
    cap_arms,
    cap_of,
    full_depth_size,
    load_nt,
    mo,
    os,
    output_dir,
    pd,
    thinned_clades,
):
    def _arm_sample_count(arm, clade):
        _p = os.path.join(
            output_dir.value, "mat", "clades", arm, f"{clade}_samples.txt"
        )
        with open(_p) as _fh:
            return sum(1 for _ in _fh)

    # Every clade the tree processed, from the full-depth sample lists.
    _clade_dir = os.path.join(output_dir.value, "mat", "clades")
    _all_clades = sorted(
        _n[: -len("_samples.txt")]
        for _n in os.listdir(_clade_dir)
        if _n.endswith("_samples.txt")
    )
    _full_counts = {
        _c: int(load_nt(_c)["actual_count"].sum()) for _c in _all_clades
    }
    _full_samples = {_c: full_depth_size(_c) or 0 for _c in _all_clades}

    _rows = [{
        "arm": "(full depth)",
        "cap": None,
        "clades_thinned": 0,
        "clades_untouched": len(_all_clades),
        "sequences": sum(_full_samples.values()),
        "mutation_counts": sum(_full_counts.values()),
        "sequences_from_thinned": 0,
        "counts_from_thinned": 0,
    }]
    for _arm in cap_arms:
        _thin = thinned_clades(_arm)
        _seq = _cnt = _seq_t = _cnt_t = 0
        for _c in _all_clades:
            if _c in _thin:
                _s = _arm_sample_count(_arm, _c)
                _k = int(load_nt(_c, _arm)["actual_count"].sum())
                _seq_t += _s
                _cnt_t += _k
            else:
                # A control: the arm reads the full tree's files unchanged.
                _s, _k = _full_samples[_c], _full_counts[_c]
            _seq += _s
            _cnt += _k
        _rows.append({
            "arm": _arm,
            "cap": cap_of(_arm),
            "clades_thinned": len(_thin),
            "clades_untouched": len(_all_clades) - len(_thin),
            "sequences": _seq,
            "mutation_counts": _cnt,
            "sequences_from_thinned": _seq_t,
            "counts_from_thinned": _cnt_t,
        })

    tree_totals = pd.DataFrame(_rows)
    _base_seq = tree_totals.sequences.iloc[0]
    _base_cnt = tree_totals.mutation_counts.iloc[0]
    tree_totals["sequence_fraction"] = tree_totals.sequences / _base_seq
    tree_totals["count_fraction"] = tree_totals.mutation_counts / _base_cnt
    # How much of what survives comes from clades the cap actually touched --
    # the rest is controls sitting at full depth, which is why the tree total
    # falls so much less steeply than any single large clade's.
    tree_totals["share_of_seqs_from_thinned"] = (
        tree_totals.sequences_from_thinned / tree_totals.sequences
    )
    tree_totals["share_of_counts_from_thinned"] = (
        tree_totals.counts_from_thinned / tree_totals.mutation_counts
    )
    mo.ui.table(
        tree_totals[[
            "arm", "cap", "clades_thinned", "clades_untouched",
            "sequences", "sequence_fraction",
            "mutation_counts", "count_fraction",
            "share_of_seqs_from_thinned", "share_of_counts_from_thinned",
        ]].round(5),
        selection=None,
    )
    return (tree_totals,)


@app.cell
def _(
    ACCENT,
    BLUE_ORDINAL,
    GRID,
    INK,
    INK_2,
    MUTED,
    SURFACE,
    mo,
    np,
    plt,
    tree_totals,
):
    _t6 = tree_totals[tree_totals.cap.notna()].sort_values("cap")
    mo.stop(not len(_t6), mo.md("**No cap arms.**"))

    _idx6 = (
        np.linspace(0, len(BLUE_ORDINAL) - 1, len(_t6)).round().astype(int)
        if len(_t6) > 1
        else [len(BLUE_ORDINAL) // 2]
    )
    _colors6 = [BLUE_ORDINAL[_i] for _i in _idx6]
    _full_seq = tree_totals.sequences.iloc[0]
    _full_cnt = tree_totals.mutation_counts.iloc[0]

    _fig6, (_b1, _b2) = plt.subplots(1, 2, figsize=(10.4, 4.0))
    for _ax6, _col6, _base6, _title6, _unit6 in (
        (_b1, "sequences", _full_seq, "Sequences in the tree", "sequences"),
        (_b2, "mutation_counts", _full_cnt,
         "Mutation counts in the tree", "counts"),
    ):
        _x6 = np.arange(len(_t6))
        _ax6.bar(
            _x6, _t6[_col6], width=0.55,
            color=_colors6, edgecolor=SURFACE, linewidth=2, zorder=2,
        )
        _ax6.axhline(
            _base6, linestyle="--", linewidth=1.6, color=INK_2, zorder=3,
            label=f"full depth ({_base6:,.0f} {_unit6})",
        )
        for _i6, _r6 in enumerate(_t6.itertuples()):
            _v6 = getattr(_r6, _col6)
            _ax6.annotate(
                f"{_v6:,.0f}\n({_v6 / _base6:.1%})",
                (_i6, _v6),
                textcoords="offset points", xytext=(0, 5),
                ha="center", fontsize=8, color=INK, linespacing=1.3,
            )
        _ax6.set_xticks(_x6)
        _ax6.set_xticklabels([f"cap {int(_c):,}" for _c in _t6.cap])
        # Linear, deliberately: a bar's length is its value, which a log axis
        # would break. These bars *should* look small next to full depth.
        _ax6.set_ylabel(_unit6)
        _ax6.set_ylim(0, _base6 * 1.22)
        _ax6.grid(True, axis="y", color=GRID, linewidth=0.6, zorder=0)
        _ax6.set_axisbelow(True)
        _ax6.set_title(_title6)
        _ax6.legend(loc="upper left", fontsize=8)
    _fig6.suptitle(
        "Tree-level totals across all clades — controls included, so a cap "
        "costs the tree far less than it costs its largest clade",
        y=1.03,
    )
    _fig6.tight_layout()
    _fig6
    return


@app.cell
def _(cap_arms, cap_of, full_depth_size, mo, thinned_clades):
    _all_thinned = sorted(
        {_c for _a in cap_arms for _c in thinned_clades(_a)},
        key=lambda _c: -(full_depth_size(_c) or 0),
    )
    mo.stop(
        not _all_thinned,
        mo.md("**No cap arm thinned any clade** — every clade sits at or below "
              "every cap, so every output byte-equals the full tree's."),
    )
    clade_picker = mo.ui.dropdown(
        options=_all_thinned,
        value=_all_thinned[0],
        label="Clade",
    )
    estimator_picker = mo.ui.dropdown(
        options=["delta_fitness_ratio", "delta_fitness_bayesian"],
        value="delta_fitness_ratio",
        label="Estimator",
    )
    mo.vstack([
        mo.hstack([clade_picker, estimator_picker], justify="start", gap=2),
        mo.md(
            f"{len(_all_thinned)} clade(s) thinned by at least one of "
            + ", ".join(f"`{_a}` (cap {cap_of(_a):,})" for _a in cap_arms)
            + ". Ordered by full-depth size, largest first — the largest clades "
            "are the ones a cap cuts hardest."
        ),
    ])
    return clade_picker, estimator_picker


@app.cell
def _(
    AA_JOIN_KEY,
    ESTIMATORS,
    GT_THRESHOLD,
    aa_dir,
    cap_arms,
    cap_of,
    full_depth_size,
    nt_dir,
    np,
    os,
    pd,
    thinned_clades,
):
    _NT_USE = ["actual_count", "expected_count"]

    def load_nt(clade, arm=None):
        """One clade's nucleotide grid at one depth: the two count columns.

        No join is needed. An eCDF is a marginal distribution and the grid is
        enumerated from the clade founder, so it does not depend on sampling
        depth -- but that is exactly what a marginal cannot check for itself,
        so the row counts are asserted equal by the caller.
        """
        _d = nt_dir if arm is None else os.path.join(nt_dir, arm)
        return pd.read_csv(
            os.path.join(_d, f"{clade}_ntmut_fitness.csv"), usecols=_NT_USE
        )

    def load_aa(clade, arm):
        """Full-depth and arm amino-acid rows, joined, plus the ground-truth mask."""
        _cols = AA_JOIN_KEY + ["actual_count", "expected_count"] + ESTIMATORS
        _full = pd.read_csv(
            os.path.join(aa_dir, f"{clade}_aamut_fitness.csv"), usecols=_cols
        )
        _arm = pd.read_csv(
            os.path.join(aa_dir, arm, f"{clade}_aamut_fitness.csv"),
            usecols=_cols,
        )
        _m = _full.merge(
            _arm,
            on=AA_JOIN_KEY,
            suffixes=("_full", "_arm"),
            validate="one_to_one",
        )
        # The grid does not depend on sampling depth, so a join that lost rows
        # means the founder moved -- which would silently make every panel a
        # comparison of two different mutation sets.
        if not (len(_m) == len(_full) == len(_arm)):
            raise AssertionError(
                f"{arm}/{clade}: joined {len(_full)} full-depth rows to "
                f"{len(_arm)} arm rows on {AA_JOIN_KEY} and got {len(_m)}"
            )
        _gt = (
            _m[["actual_count_full", "expected_count_full"]].max(axis=1)
            >= GT_THRESHOLD
        )
        return _m, _gt

    def fit_line(x, y):
        """OLS slope, intercept and Pearson r of the arm on full depth."""
        if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
            return np.nan, np.nan, np.nan
        _slope, _intercept = np.polyfit(x, y, 1)
        return _slope, _intercept, float(np.corrcoef(x, y)[0, 1])

    def arms_for(clade):
        return [
            _a for _a in cap_arms
            if clade in thinned_clades(_a)
            and os.path.isfile(
                os.path.join(aa_dir, _a, f"{clade}_aamut_fitness.csv")
            )
        ]

    def clade_stats(clade):
        """The three statistics, one row per (arm, estimator)."""
        _nt_full = load_nt(clade)
        _rows = []
        for _arm in arms_for(clade):
            _nt_arm = load_nt(clade, _arm)
            if len(_nt_arm) != len(_nt_full):
                raise ValueError(
                    f"{_arm}/{clade}: nucleotide grid has {len(_nt_arm)} rows "
                    f"against {len(_nt_full)} at full depth; the grid must not "
                    f"depend on sampling depth, so the eCDFs below would have "
                    f"different denominators"
                )
            # Statistic 1, on the nucleotide grid -- the same population the
            # eCDFs are drawn over, so the fraction and the curves agree.
            _sum_full = float(_nt_full["actual_count"].sum())
            _sum_arm = float(_nt_arm["actual_count"].sum())
            _count_fraction = _sum_arm / _sum_full if _sum_full else np.nan

            _m, _gt_mask = load_aa(clade, _arm)
            _gt = _m[_gt_mask]
            # The raw deletion counts, not a second retention fraction: a
            # mutation that loses its last observation leaves the scored set
            # entirely, which is worth seeing, but dividing these gives a
            # number that competes with count_fraction and flatters the arm.
            _scored_full = int((_m["actual_count_full"] > 0).sum())
            _scored_arm = int((_m["actual_count_arm"] > 0).sum())
            for _est in ESTIMATORS:
                _x = _gt[f"{_est}_full"].to_numpy()
                _y = _gt[f"{_est}_arm"].to_numpy()
                _sl, _ic, _r = fit_line(_x, _y)
                _rows.append({
                    "clade": clade,
                    "arm": _arm,
                    "cap": cap_of(_arm),
                    "full_depth_samples": full_depth_size(clade),
                    "estimator": _est,
                    "count_fraction": _count_fraction,
                    "sum_actual_full": int(_sum_full),
                    "sum_actual_arm": int(_sum_arm),
                    "n_ground_truth": len(_gt),
                    "pearson_r": _r,
                    "slope": _sl,
                    "intercept": _ic,
                    "sd_full": float(np.std(_x)) if len(_x) else np.nan,
                    "sd_arm": float(np.std(_y)) if len(_y) else np.nan,
                    "n_scored_full": _scored_full,
                    "n_scored_arm": _scored_arm,
                })
        _df = pd.DataFrame(_rows)
        if len(_df):
            _df["effective_fraction"] = (
                _df["cap"] / _df["full_depth_samples"]
            )
        return _df
    return arms_for, clade_stats, fit_line, load_aa, load_nt


@app.cell
def _(clade_picker, clade_stats, estimator_picker, mo):
    stats = clade_stats(clade_picker.value)
    mo.stop(
        stats.empty,
        mo.md(f"**No cap arm thinned {clade_picker.value}.**"),
    )
    this = stats[stats.estimator == estimator_picker.value].sort_values("cap")
    return stats, this


@app.cell
def _(clade_picker, mo, this):
    mo.md(
        f"## {clade_picker.value} — statistic 1: the count fraction\n\n"
        + "\n".join(
            f"- **cap {_r.cap:,}** — effective sampling fraction "
            f"{_r.effective_fraction:.2%} of "
            f"{int(_r.full_depth_samples):,} samples, and it keeps "
            f"**{_r.count_fraction:.2%}** of the clade's "
            f"{_r.sum_actual_full:,} observed nucleotide mutation counts "
            f"({_r.sum_actual_arm:,}), a "
            f"{_r.count_fraction / _r.effective_fraction:.2f}× "
            f"retention against the nominal cut."
            for _r in this.itertuples()
        )
        + "\n\nRetention runs **above** the nominal fraction because a branch "
        "survives whenever any sample below it is retained — the curve is "
        "concave, so a 100× cut in samples is never a 100× cut in evidence. "
        "That is why every number here keys to the measured fraction rather "
        "than to the configured cap."
    )
    return


@app.cell
def _(mo):
    mo.md("""
    ## Statistic 2: the eCDF of per-mutation counts

    The fraction above is one number; this is the distribution behind it. Each
    curve is the share of the clade's **whole nucleotide grid** — every
    single-nucleotide mutation away from the clade founder — carrying at least
    *x* counts. Thinning shifts the curve left, and how far it shifts is the
    cost.

    The grid is the right population here, not the ground-truth subset: the
    subset was *defined* by a threshold on these counts, so its full-depth
    `max(actual, expected)` curve would pass through 1.0 at exactly that
    threshold and the curve would be truncated by construction. The dotted
    vertical marks the ground-truth threshold, so you can see what share of the
    grid each depth still puts above it.
    """)
    return


@app.cell
def _(
    ACCENT,
    BLUE_ORDINAL,
    GRID,
    GT_THRESHOLD,
    INK_2,
    MUTED,
    arms_for,
    cap_of,
    clade_picker,
    load_nt,
    mo,
    np,
    plt,
):
    _clade3 = clade_picker.value
    _arms3 = arms_for(_clade3)
    mo.stop(not _arms3, mo.md("**No arm thinned this clade.**"))

    _full3 = load_nt(_clade3)
    _idx3 = (
        np.linspace(0, len(BLUE_ORDINAL) - 1, len(_arms3)).round().astype(int)
        if len(_arms3) > 1
        else [len(BLUE_ORDINAL) // 2]
    )
    _color3 = {_a: BLUE_ORDINAL[_i] for _a, _i in zip(_arms3, _idx3)}

    _PANELS = [
        ("actual_count", "actual counts"),
        ("expected_count", "expected counts"),
        ("max_count", "max(actual, expected)"),
    ]

    def _quantity(frame, which):
        if which == "max_count":
            return frame[["actual_count", "expected_count"]].max(axis=1)
        return frame[which]

    def _surv(values, xs):
        _v = np.asarray(values, dtype=float)
        return [float((_v >= _x).mean()) for _x in xs]

    _fig3, _axes3 = plt.subplots(1, 3, figsize=(14.4, 4.2))
    for _col3, (_which3, _label3) in enumerate(_PANELS):
        _ax3 = _axes3[_col3]
        _hi3 = max(float(np.nanmax(_quantity(_full3, _which3))), 2.0)
        _xs3 = np.geomspace(0.05, _hi3, 90)
        _ax3.plot(
            _xs3, _surv(_quantity(_full3, _which3), _xs3),
            linewidth=2.6, color=INK_2, zorder=4, label="full depth",
        )
        for _arm3 in reversed(_arms3):
            _a3 = load_nt(_clade3, _arm3)
            if len(_a3) != len(_full3):
                raise ValueError(
                    f"{_arm3}/{_clade3}: nucleotide grid has {len(_a3)} rows "
                    f"against {len(_full3)} at full depth"
                )
            _ax3.plot(
                _xs3, _surv(_quantity(_a3, _which3), _xs3),
                linewidth=2.0, color=_color3[_arm3], zorder=3,
                label=f"cap {cap_of(_arm3):,}",
            )
        _ax3.axvline(
            GT_THRESHOLD, linestyle=":", linewidth=1.4, color=ACCENT, zorder=2,
        )
        _ax3.annotate(
            f"ground-truth\nthreshold ({GT_THRESHOLD})",
            (GT_THRESHOLD, 0.94),
            textcoords="offset points", xytext=(6, 0),
            fontsize=8, color=MUTED, va="top",
        )
        _ax3.set_xscale("log")
        _ax3.set_xlabel(f"x   ({_label3})")
        _ax3.set_ylim(0, 1.02)
        _ax3.grid(True, color=GRID, linewidth=0.6, zorder=0)
        _ax3.set_title(_label3)
        if _col3 == 0:
            _ax3.set_ylabel("share of the grid with ≥ x")
    _h3, _l3 = _axes3[0].get_legend_handles_labels()
    _fig3.legend(
        _h3, _l3, loc="lower center", ncol=len(_l3),
        bbox_to_anchor=(0.5, -0.07), title="sampling depth",
    )
    _fig3.suptitle(
        f"{_clade3} — count eCDFs over the whole nucleotide grid "
        f"(n = {len(_full3):,})",
        y=1.02,
    )
    _fig3.tight_layout()
    _fig3
    return


@app.cell
def _(mo):
    mo.md("""
    ## Statistic 3: correlation with the ground-truth effects

    One panel per cap, in increasing order (more data to the right). Every panel
    shares one square, symmetric axis pair, so the arm's range collapsing toward
    zero is visible as *geometry* rather than only as a number. Dashed is
    **y = x** — perfect recovery; solid is the **fitted line**, whose slope is
    the shrinkage.
    """)
    return


@app.cell
def _(
    ACCENT,
    DENSITY_CMAP,
    GRID,
    GT_THRESHOLD,
    INK,
    INK_2,
    LogNorm,
    MUTED,
    SCATTER_MAX_N,
    SMALL_N,
    arms_for,
    cap_of,
    clade_picker,
    estimator_picker,
    load_aa,
    mo,
    np,
    plt,
    this,
):
    _clade4 = clade_picker.value
    _est4 = estimator_picker.value
    _arms4 = arms_for(_clade4)
    mo.stop(not _arms4, mo.md("**No arm thinned this clade.**"))

    _panels4 = {}
    for _arm4 in _arms4:
        _m4, _gt4 = load_aa(_clade4, _arm4)
        _panels4[_arm4] = _m4[_gt4]

    # One symmetric limit for the whole row, set by full depth: the reference
    # frame must not shrink with the arm, or the compression is invisible.
    _lim4 = 1.0
    for _d4 in _panels4.values():
        if len(_d4):
            _lim4 = max(
                _lim4,
                float(np.nanmax(np.abs(np.concatenate([
                    _d4[f"{_est4}_full"].to_numpy(),
                    _d4[f"{_est4}_arm"].to_numpy(),
                ])))),
            )
    _lim4 *= 1.05

    _fig4, _axes4 = plt.subplots(
        1, len(_arms4),
        figsize=(3.3 * len(_arms4) + 0.6, 3.9),
        squeeze=False,
    )
    for _j4, _arm4 in enumerate(_arms4):
        _ax4 = _axes4[0][_j4]
        _ax4.set_xlim(-_lim4, _lim4)
        _ax4.set_ylim(-_lim4, _lim4)
        _ax4.set_aspect("equal")
        _ax4.grid(True, color=GRID, linewidth=0.6, zorder=0)
        _ax4.axhline(0, color=GRID, linewidth=0.6, zorder=0)
        _ax4.axvline(0, color=GRID, linewidth=0.6, zorder=0)

        _d4 = _panels4[_arm4]
        _x4 = _d4[f"{_est4}_full"].to_numpy()
        _y4 = _d4[f"{_est4}_arm"].to_numpy()

        if len(_x4) > SCATTER_MAX_N:
            _ax4.hexbin(
                _x4, _y4,
                gridsize=44,
                extent=(-_lim4, _lim4, -_lim4, _lim4),
                cmap=DENSITY_CMAP, norm=LogNorm(), mincnt=1,
                linewidths=0, zorder=1,
            )
        elif len(_x4):
            _ax4.scatter(
                _x4, _y4,
                s=14, alpha=0.5, linewidths=0,
                color=DENSITY_CMAP(0.75), zorder=1,
            )

        _ax4.plot(
            [-_lim4, _lim4], [-_lim4, _lim4],
            linestyle="--", linewidth=1.4, color=INK_2, zorder=3,
            label="y = x (perfect recovery)",
        )
        _row4 = this[this.arm == _arm4]
        if len(_row4) and np.isfinite(_row4.slope.iloc[0]):
            _sl4 = float(_row4.slope.iloc[0])
            _ic4 = float(_row4.intercept.iloc[0])
            _xs4 = np.array([-_lim4, _lim4])
            _ax4.plot(
                _xs4, _sl4 * _xs4 + _ic4,
                linewidth=2.0, color=ACCENT, zorder=4, label="fitted line",
            )
            _n4 = int(_row4.n_ground_truth.iloc[0])
            _note4 = (
                f"r = {_row4.pearson_r.iloc[0]:.3f}\n"
                f"slope = {_sl4:.3f}\n"
                f"n = {_n4:,}"
            )
            if _n4 < SMALL_N:
                _note4 += "  ⚠"
            _ax4.text(
                0.04, 0.96, _note4,
                transform=_ax4.transAxes, ha="left", va="top",
                fontsize=8, color=INK, linespacing=1.35,
            )
            _ax4.text(
                0.96, 0.06,
                f"count fraction {_row4.count_fraction.iloc[0]:.3%}",
                transform=_ax4.transAxes, ha="right", va="bottom",
                fontsize=8, color=MUTED,
            )
        _ax4.set_title(f"cap {cap_of(_arm4):,}", pad=8)
        _ax4.set_xlabel("full-depth estimate")
        if _j4 == 0:
            _ax4.set_ylabel("arm estimate")

    _h4, _l4 = _axes4[0][0].get_legend_handles_labels()
    if _h4:
        _fig4.legend(
            _h4, _l4, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.13)
        )
    _fig4.suptitle(
        f"{_clade4} — {_est4} over the fixed ground-truth set "
        f"(max(actual, expected) ≥ {GT_THRESHOLD:g} at full depth)",
        y=1.02,
    )
    _fig4.tight_layout()
    _fig4
    return


@app.cell
def _(SMALL_N, clade_picker, mo, this):
    _small = this[this.n_ground_truth < SMALL_N]
    mo.md(
        "⚠ **"
        + ", ".join(
            f"cap {_r.cap:,} ({_r.n_ground_truth:,} mutations)"
            for _r in _small.itertuples()
        )
        + f"** rest on fewer than {SMALL_N:,} ground-truth mutations, so their "
        "coefficients have little resolution — 20F at cap 5,000 reads r = 0.999 "
        "from six mutations. Read the panel, not the number."
        if len(_small)
        else f"Every panel above rests on at least {SMALL_N:,} ground-truth "
        f"mutations ({int(this.n_ground_truth.min()):,} at the thinnest cap), "
        "so the coefficients have resolution."
    )
    return


if __name__ == "__main__":
    app.run()
