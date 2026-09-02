#!/usr/bin/env python3
"""Correlate each arm's per-clade fitness effects against full depth (#151).

This is the headline result the arms exist to produce: how far a clade's
estimated fitness effects drift as its sampling depth falls. One row per
``(arm, clade, estimator, population)``, keyed to measured evidence.

Six things about the shape are load-bearing.

*The population is a **fixed** ground-truth set, defined on the full-depth row.*
A mutation is in it when ``max(actual_count, expected_count) >= count_threshold``
**at full depth**, so the set is identical for every arm and the coefficients are
comparable across sampling depths. Filtering on the *arm's* row instead -- the
obvious reading, and what this script did first -- is survivorship bias: it drops
exactly the mutations the thinning destroyed, so the denominator shrinks with the
cut and a degrading estimator reports a flat correlation. Every arm has an
estimate for every mutation in the set, because the grid is enumerated from the
clade founder by ``build_count_table`` and does not depend on sampling depth;
the join therefore asserts ``one_to_one`` **and** that no ground-truth row was
lost, rather than trusting the merge.

*Correlating over the whole grid would measure the grid.* Every row carries a
``delta_fitness_ratio`` even where nothing was observed -- 121,348 of 121,427
rows for 20A on the test data -- and those rows are invariant to sampling depth
by construction, so a grid-wide coefficient is a confident number describing
``build_count_table``'s enumeration: 0.619 for the ratio estimator where the 29
genuinely scored mutations give 1.000. The ``>= count_threshold`` floor is what
makes the population a set of mutations there is evidence *about*.

*r cannot see bias or compression, so it is never reported alone.* Pearson r is
invariant to any linear transform of the arm's estimates, so a systematic shift
and a systematic shrinkage are both invisible to it. At ``cap_5000_0`` on the
v12 gisaid tree 21J's regression slope is **0.117** while r still reads 0.431,
and its estimates' sd falls from 1.77 to 0.48 -- an 8.5x compression that r
prices at a halving. The arms shrink fitness effects toward zero as evidence
thins, which is what the pseudocount does to a mutation whose observations
vanished. ``rmse``, ``bias``, ``regression_slope`` and the two standard
deviations are reported beside every coefficient for that reason.

*What thinning costs is measured as a fraction of observed counts, not of
distinct mutations.* ``count_fraction = sum(actual_count_arm) /
sum(actual_count_full)`` over the whole grid is the headline evidence
statistic. The rejected alternative was a *coverage* fraction -- the share of
full-depth-scored mutations the arm still scores at all -- which counts each
mutation once however much evidence it lost, so a mutation falling from 900
observations to 3 reads as fully retained. It is not the denominator of the
cut, it diverges from one by up to 30x on real data (21J at cap 5,000: 0.0033
against 0.103), and it is the *worse* predictor of the thing the arms exist to
measure: Spearman against the fitted regression slope over the 120 completed cap
fits is 0.985 for ``count_fraction`` and 0.940 for coverage. So it is not
reported at all -- two retention fractions differing threefold with no stated
reason to prefer one is a standing invitation to quote whichever is kinder.

``n_scored_full`` and ``n_scored_arm`` remain, because the deletion channel is
real and separable -- a mutation that loses its last observation leaves the
scored set entirely rather than getting noisier -- and those are its raw
counts. A consumer who wants the ratio can divide; what is not shipped is a
second *fraction* competing with the headline one.

It is taken over the **whole grid** rather than the ground-truth subset,
because it describes the thinning and not the comparison population, and it is
computed on the amino-acid grid this script already holds. That agrees with the
nucleotide-level ratio to better than 0.002 on every measured (arm, clade) --
the two differ only by the double-counting at gene overlaps, which is nearly
identical in numerator and denominator.

*The two estimators are never pooled.* ``delta_fitness_ratio`` and
``delta_fitness_bayesian`` respond differently to thin evidence, so a single
coefficient would hide which one degrades and where.

*A coefficient is meaningless without the n beside it, and sometimes with it.*
A row is **kept with an empty coefficient** rather than dropped when the
population is too small or one side is constant: a silently absent row reads as
"not applicable", where an explicit ``n_mutations`` with a blank ``pearson_r``
reads as "asked, and there was nothing to answer with".
``max_actual_count_full`` is reported for the same reason -- when every scored
mutation carries exactly one observation, ``delta_fitness_ratio`` is a
deterministic function of the expected count, the two columns are collinear by
construction, and the resulting 1.000 is arithmetic rather than agreement. That
is the whole regime on the test data, and no coefficient alone can show it.

*Non-unit clades are included as a control.* An arm re-aggregates every clade in
its grid but thins only its unit, so a clade the arm does not thin must come out
**identical** to full depth. Reporting those rows with ``thinned=False`` turns
the table into its own check: a control clade with rmse > 0 means the arm leaked
into a clade it was never meant to touch, which is the fan-out bug this feature
had to fix once already. The controls are the clades at or below the cap, so
they cost nothing to obtain and cannot be chosen -- they fall out of the arm's
definition rather than being clades someone remembered to leave alone.

Spearman is computed as Pearson-on-ranks rather than via ``pandas``' ``method=
"spearman"``, which imports ``scipy`` -- absent from this environment.
``Series.rank()`` uses average ranks for ties, which is the standard definition.
"""

import os
from importlib import util as importlib_util

import numpy as np
import pandas as pd

_loader_spec = importlib_util.spec_from_file_location(
    "workflow_loader",
    os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "modules",
        "workflow_loader.py",
    ),
)
workflow_loader = importlib_util.module_from_spec(_loader_spec)
_loader_spec.loader.exec_module(workflow_loader)
subsample = workflow_loader.load_workflow_module("subsample")

roster = snakemake.params.roster
outdir = snakemake.params.outdir
count_threshold = float(snakemake.params.count_threshold)
output_path = snakemake.output.summary
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None

JOIN_KEY = ["cluster", "gene", "aa_mutation"]
ESTIMATORS = ["delta_fitness_ratio", "delta_fitness_bayesian"]

COLUMNS = [
    "arm",
    "seed",
    "clade",
    "max_samples",
    "thinned",
    "estimator",
    "population",
    "n_mutations",
    "n_ground_truth",
    "sum_actual_full",
    "sum_actual_arm",
    "count_fraction",
    "n_scored_full",
    "n_scored_arm",
    "max_actual_count_full",
    "pearson_r",
    "spearman_rho",
    "rmse",
    "median_abs_diff",
    "bias",
    "regression_slope",
    "sd_full",
    "sd_arm",
]


def log(msg):
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


def pearson(x, y):
    """Pearson r, or ``None`` when it is undefined rather than merely small.

    Fewer than three points, or a constant on either side, makes the coefficient
    undefined (the denominator is zero). Returning ``None`` puts a blank in the
    cell next to a populated ``n_mutations``, which is distinguishable from a
    genuine 0.0; ``numpy`` would return ``nan`` with a RuntimeWarning and the two
    would print identically.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def spearman(x, y):
    """Pearson on average ranks -- the definition, without the scipy import."""
    if len(x) < 3:
        return None
    return pearson(pd.Series(x).rank(), pd.Series(y).rank())


def stats_for(sub, estimator):
    """Agreement statistics for one population and estimator.

    Returns a dict rather than a tuple because it carries eight numbers whose
    order nothing else depends on, and a positional swap between ``bias`` and
    ``median_abs_diff`` would be invisible in the output.

    ``bias`` is signed ``mean(arm - full)``: the direction matters, because on
    the completed arms it is positive everywhere -- thinned arms estimate
    mutations as *more fit* than full depth does -- and an absolute value would
    hide that. ``regression_slope`` is the least-squares slope of arm on full,
    which is 1.0 for an unbiased-in-scale estimate and falls toward 0 as the
    estimates shrink toward the pseudocount's bound.
    """
    a = sub[f"{estimator}_full"]
    b = sub[f"{estimator}_arm"]
    keep = a.notna() & b.notna()
    a, b = a[keep].to_numpy(dtype=float), b[keep].to_numpy(dtype=float)
    n = int(len(a))
    out = {
        "n_mutations": n,
        "pearson_r": None,
        "spearman_rho": None,
        "rmse": None,
        "median_abs_diff": None,
        "bias": None,
        "regression_slope": None,
        "sd_full": None,
        "sd_arm": None,
    }
    if n == 0:
        return out
    diff = b - a
    out["rmse"] = float(np.sqrt((diff ** 2).mean()))
    out["median_abs_diff"] = float(np.median(np.abs(diff)))
    out["bias"] = float(diff.mean())
    out["pearson_r"] = pearson(a, b)
    out["spearman_rho"] = spearman(a, b)
    if n >= 2:
        out["sd_full"] = float(a.std(ddof=1))
        out["sd_arm"] = float(b.std(ddof=1))
        var_full = float(((a - a.mean()) ** 2).sum())
        if var_full > 0:
            out["regression_slope"] = float(
                ((a - a.mean()) * (b - b.mean())).sum() / var_full
            )
    return out


def fmt(value, places=6):
    return "" if value is None else round(value, places)


rows = []
problems = []
n_thinned_compared = 0

for arm in sorted(roster):
    spec = roster[arm]
    unit = set(spec["unit"])
    included = set(spec["included"])
    seed = int(spec["seed"])
    cap = int(spec["cap"])

    arm_fitness_dir = os.path.join(outdir, "aamut_fitness", arm)
    if not os.path.isdir(arm_fitness_dir):
        problems.append(
            f"arm {arm}: {arm_fitness_dir} does not exist, so its fitness "
            f"effects cannot be compared against full depth"
        )
        continue

    for name in sorted(os.listdir(arm_fitness_dir)):
        if not name.endswith("_aamut_fitness.csv") or name.startswith("all_"):
            continue
        clade = name[: -len("_aamut_fitness.csv")]
        full_path = os.path.join(outdir, "aamut_fitness", name)
        arm_path = os.path.join(arm_fitness_dir, name)
        if not os.path.isfile(full_path):
            problems.append(
                f"arm {arm}, clade {clade}: the arm has a fitness table but the "
                f"full tree does not ({full_path}), so there is nothing to "
                f"correlate against"
            )
            continue

        full = pd.read_csv(full_path)
        armd = pd.read_csv(arm_path)
        merged = full.merge(
            armd, on=JOIN_KEY, suffixes=("_full", "_arm"), validate="one_to_one"
        )
        # The grid is invariant to sampling depth, so an inner join that lost
        # rows means the founder moved or a gene dropped out -- which would make
        # every coefficient below a comparison of two different mutation sets.
        if len(merged) != len(full) or len(merged) != len(armd):
            problems.append(
                f"arm {arm}, clade {clade}: joining {len(full)} full-depth rows "
                f"to {len(armd)} arm rows on {JOIN_KEY} yields {len(merged)}. "
                f"The grid is enumerated from the clade founder and does not "
                f"depend on sampling depth, so these must agree."
            )
            continue

        # The fixed ground-truth set: read off the **full-depth** row, so it is
        # the same set of mutations at every sampling depth and the coefficients
        # below are comparable across arms.
        gt_mask = (
            merged[["actual_count_full", "expected_count_full"]].max(axis=1)
            >= count_threshold
        )
        ground_truth = merged[gt_mask]
        n_ground_truth = int(gt_mask.sum())

        # The evidence statistic: what share of the clade's observed mutation
        # *counts* the arm still has. Taken over the whole grid, not the
        # ground-truth subset, because it is the denominator of the thinning
        # rather than a property of the comparison population -- and summed
        # rather than counted, so a mutation dropping from 900 observations to
        # 3 registers as the loss it is.
        sum_actual_full = float(merged["actual_count_full"].sum())
        sum_actual_arm = float(merged["actual_count_arm"].sum())
        count_fraction = (
            sum_actual_arm / sum_actual_full if sum_actual_full else None
        )

        scored_full = merged["actual_count_full"] > 0
        scored_arm = merged["actual_count_arm"] > 0
        both = merged[scored_full & scored_arm]
        n_scored_full = int(scored_full.sum())
        n_scored_arm = int(scored_arm.sum())
        thinned = clade in unit and clade in included
        if thinned and n_ground_truth:
            n_thinned_compared += 1

        populations = [
            (f"ground_truth_max_full>={count_threshold:g}", ground_truth),
            ("scored_at_both_depths", both),
        ]

        base = {
            "arm": arm,
            "seed": seed,
            "clade": clade,
            "max_samples": cap,
            "thinned": thinned,
            "n_ground_truth": n_ground_truth,
            "sum_actual_full": int(sum_actual_full),
            "sum_actual_arm": int(sum_actual_arm),
            "count_fraction": fmt(count_fraction),
            "n_scored_full": n_scored_full,
            "n_scored_arm": n_scored_arm,
            "max_actual_count_full": (
                int(merged["actual_count_full"].max()) if len(merged) else 0
            ),
        }

        for population_name, population in populations:
            for estimator in ESTIMATORS:
                stats = stats_for(population, estimator)
                rows.append({
                    **base,
                    "estimator": estimator,
                    "population": population_name,
                    "n_mutations": stats["n_mutations"],
                    "pearson_r": fmt(stats["pearson_r"]),
                    "spearman_rho": fmt(stats["spearman_rho"]),
                    "rmse": fmt(stats["rmse"]),
                    "median_abs_diff": fmt(stats["median_abs_diff"]),
                    "bias": fmt(stats["bias"]),
                    "regression_slope": fmt(stats["regression_slope"]),
                    "sd_full": fmt(stats["sd_full"]),
                    "sd_arm": fmt(stats["sd_arm"]),
                })

                # A control clade must reproduce full depth exactly. Checked on
                # rmse rather than on the coefficient, since r is undefined when
                # a population is constant and would pass by being blank. Run on
                # the ground-truth population because it is non-empty wherever
                # the clade has any evidence at all.
                if (
                    not thinned
                    and population_name.startswith("ground_truth")
                    and stats["rmse"] is not None
                    and stats["rmse"] > 0
                ):
                    problems.append(
                        f"arm {arm}, clade {clade} is not thinned by this arm "
                        f"({subsample.depth_label(spec)}), so its fitness "
                        f"effects must equal full depth's, but {estimator} "
                        f"differs (rmse {stats['rmse']:.6g} over "
                        f"{stats['n_mutations']} mutations). The arm reached a "
                        f"clade it was never meant to touch."
                    )

        log(
            f"{arm}/{clade}{'' if thinned else ' [control]'}: "
            f"{n_ground_truth} ground-truth mutations "
            f"(max(actual,expected) >= {count_threshold:g} at full depth), "
            f"count fraction {fmt(count_fraction)} "
            f"({int(sum_actual_arm):,}/{int(sum_actual_full):,} observations), "
            f"{n_scored_arm}/{n_scored_full} distinct scored, "
            f"max full-depth count {base['max_actual_count_full']}"
        )

if problems:
    raise ValueError(
        "subsample_correlations found "
        + str(len(problems))
        + " problem(s):\n  "
        + "\n  ".join(problems)
    )

summary = pd.DataFrame(rows, columns=COLUMNS)
summary.to_csv(output_path, index=False)

if not n_thinned_compared:
    log(
        "VACUOUS: no arm thinned a clade with any ground-truth mutation, so no "
        "row of this table compares a subsampled estimate against full depth."
    )
else:
    thinned_rows = summary[summary["thinned"] & (summary["n_mutations"] > 0)]
    log(
        f"Wrote {len(summary)} rows for {summary['arm'].nunique()} arm(s); "
        f"{n_thinned_compared} thinned (arm, clade) pair(s) compared, "
        f"{len(thinned_rows)} populated row(s) among them -> {output_path}"
    )
    maxcount = int(summary["max_actual_count_full"].max())
    if maxcount <= 1:
        log(
            f"CAUTION: the highest full-depth actual_count anywhere in this "
            f"table is {maxcount}. With one observation per scored mutation, "
            f"delta_fitness_ratio is a deterministic function of the expected "
            f"count and any correlation it shows is arithmetic, not agreement. "
            f"NOT A NOISE MEASUREMENT."
        )
