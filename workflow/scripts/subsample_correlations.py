#!/usr/bin/env python3
"""Correlate each arm's per-clade fitness effects against full depth (#151).

This is the headline result the arms exist to produce: how far a clade's
estimated fitness effects drift as its sampling depth falls. One row per
``(arm, clade, estimator, filter)``, keyed to measured evidence.

Five things about the shape are load-bearing.

*The correlation is computed over **scored** mutations, not the grid.* Every row
of ``{clade}_aamut_fitness.csv`` carries a ``delta_fitness_ratio`` -- the grid is
enumerated from the clade founder and filled with zero counts, so a mutation that
was never observed still gets a number, driven entirely by its expected count.
On the test data that is 121,348 of 121,427 rows for 20A. Correlating over all of
them measures the grid's construction, which is invariant to sampling fraction by
design, and returns a confident coefficient describing nothing: 0.619 for the
ratio estimator where the 29 genuinely scored mutations give 1.000. ``scored``
means ``actual_count > 0``, and the population is the **intersection** -- scored
at both depths -- because a mutation the arm never observed has no estimate to
compare, only an absence, which is what ``coverage`` reports instead.

*The count filter is applied to the arm's row, never the full-depth row.* The arm
is the side whose evidence is being cut; filtering on the target would move the
mutation set for two reasons at once and make the two filters incomparable.

*The two estimators are never pooled.* ``delta_fitness_ratio`` and
``delta_fitness_bayesian`` respond differently to thin evidence -- the Bayesian
estimator shrinks toward zero as counts fall while the ratio does not -- so a
single coefficient would hide which one degrades and where.

*A coefficient is meaningless without the n beside it, and sometimes with it.*
The count filter raises rho partly by dropping hard cases, so ``n_mutations`` is
reported on every row and a row is **kept with an empty coefficient** rather than
dropped when the population is too small or one side is constant. A silently
absent row reads as "not applicable"; an explicit ``n_mutations`` with a blank
``pearson_r`` reads as "asked, and there was nothing to answer with".
``max_actual_count_full`` is reported for the same reason: when every scored
mutation carries exactly one observation, ``delta_fitness_ratio`` is a
deterministic function of the expected count, the two columns are collinear by
construction, and the resulting 1.000 is arithmetic rather than agreement. That
is the whole regime on the test data, and no coefficient alone can show it.

*Non-unit clades are included as a control.* An arm re-aggregates every clade in
its grid but thins only its unit, so a non-unit clade's fitness must come out
**identical** to full depth. Reporting those rows with ``thinned=False`` turns the
table into its own check: a control clade with rho < 1 means the arm leaked into a
clade the config never named, which is the fan-out bug this feature had to fix
once already. Dropping them would remove the evidence that it stayed fixed.

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
    "unit",
    "seed",
    "clade",
    "fraction",
    "thinned",
    "estimator",
    "filter",
    "n_mutations",
    "n_scored_full",
    "n_scored_arm",
    "coverage",
    "max_actual_count_full",
    "pearson_r",
    "spearman_rho",
    "rmse",
    "median_abs_diff",
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
    """``(n, r, rho, rmse, median_abs_diff)`` for one population and estimator."""
    a = sub[f"{estimator}_full"]
    b = sub[f"{estimator}_arm"]
    keep = a.notna() & b.notna()
    a, b = a[keep], b[keep]
    n = int(len(a))
    if n == 0:
        return n, None, None, None, None
    diff = (a - b).abs()
    return (
        n,
        pearson(a, b),
        spearman(a, b),
        float(np.sqrt(((a - b) ** 2).mean())),
        float(diff.median()),
    )


def fmt(value, places=6):
    return "" if value is None else round(value, places)


rows = []
problems = []
n_thinned_compared = 0

for arm in sorted(roster):
    spec = roster[arm]
    unit = list(spec["unit"])
    frac = float(spec["frac"])
    seed = int(spec["seed"])
    included = set(spec["included"])
    unit_label = "-".join(unit)

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
        # The grid is invariant to sampling fraction, so an inner join that lost
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

        scored_full = merged["actual_count_full"] > 0
        scored_arm = merged["actual_count_arm"] > 0
        both = merged[scored_full & scored_arm]
        n_scored_full = int(scored_full.sum())
        n_scored_arm = int(scored_arm.sum())
        coverage = len(both) / n_scored_full if n_scored_full else None
        thinned = clade in unit and clade in included
        if thinned and len(both):
            n_thinned_compared += 1

        # The threshold is read off the *arm's* row: it is the side whose
        # evidence is being cut, and filtering on the full-depth row would move
        # the mutation set for two reasons at once.
        strong = both[
            both[["actual_count_arm", "expected_count_arm"]].max(axis=1)
            >= count_threshold
        ]
        populations = [
            ("unfiltered", both),
            (f"max(actual,expected)_arm>={count_threshold:g}", strong),
        ]

        base = {
            "arm": arm,
            "unit": unit_label,
            "seed": seed,
            "clade": clade,
            "fraction": frac,
            "thinned": thinned,
            "n_scored_full": n_scored_full,
            "n_scored_arm": n_scored_arm,
            "coverage": fmt(coverage),
            "max_actual_count_full": (
                int(merged["actual_count_full"].max()) if len(merged) else 0
            ),
        }

        for filter_name, population in populations:
            for estimator in ESTIMATORS:
                n, r, rho, rmse, mad = stats_for(population, estimator)
                rows.append({
                    **base,
                    "estimator": estimator,
                    "filter": filter_name,
                    "n_mutations": n,
                    "pearson_r": fmt(r),
                    "spearman_rho": fmt(rho),
                    "rmse": fmt(rmse),
                    "median_abs_diff": fmt(mad),
                })

                # A control clade must reproduce full depth exactly. Checked on
                # the numbers rather than on the coefficient, since r is
                # undefined when a population is constant and would pass by
                # being blank.
                if (
                    not thinned
                    and filter_name == "unfiltered"
                    and rmse is not None
                    and rmse > 0
                ):
                    problems.append(
                        f"arm {arm}, clade {clade} is not in the arm's unit "
                        f"{unit}, so its fitness effects must equal full "
                        f"depth's, but {estimator} differs (rmse {rmse:.6g} over "
                        f"{n} mutations). The arm reached a clade the config "
                        f"never named."
                    )

        log(
            f"{arm}/{clade}{'' if thinned else ' [control]'}: "
            f"{n_scored_arm}/{n_scored_full} scored mutations "
            f"(coverage {fmt(coverage)}), max full-depth count "
            f"{base['max_actual_count_full']}"
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
        "VACUOUS: no arm thinned a clade that retained any scored mutation, so "
        "no row of this table compares a subsampled estimate against full depth."
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
            f"These coefficients do not measure subsampling noise."
        )
