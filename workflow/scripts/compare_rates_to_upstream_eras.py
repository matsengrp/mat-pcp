#!/usr/bin/env python
"""Compare mat-pcp GLM rate tables against ``SARS2-mut-fitness-v2``'s per-era ones.

Read-only diagnostic over a completed run. Not a pipeline rule::

    # before #137 Task 2 -- one pooled table against both eras
    pixi run python workflow/scripts/compare_rates_to_upstream_eras.py \\
        --pooled results/gisaid-and-public_2026-04-22_v9/glm/master_table.csv \\
        --outdir analysis/137-era-offsets/upstream

    # after #137 Task 2 -- era against matching era
    pixi run python workflow/scripts/compare_rates_to_upstream_eras.py \\
        --pre-omicron results/<tree>/glm/master_table_pre_omicron.csv \\
        --omicron     results/<tree>/glm/master_table_omicron.csv \\
        --outdir analysis/137-era-offsets/upstream

**This is a different kind of check from
``tests/test_glm_reference_equivalence.py``, and conflating them would be a
mistake in both directions.** That test runs *upstream's code* on *mat-pcp's
fixture input* and compares at floating-point-noise tolerance -- it asks "is the
port faithful?", and any real difference fails it. This script compares
*upstream's published rates, fit on upstream's own data* against mat-pcp's, fit
on a different tree with a different date cutoff and different filtering. The
two can only agree *approximately*, and demanding FP tolerance here would be
meaningless. So this reports agreement statistics and asserts nothing about
absolute equality.

What it does assert: that the 512-cell grids join **1:1** on
``(mut_type, motif, unpaired, nt_site_before_boundary)``. A silently partial join
would compute every statistic below over a subset and report healthy-looking
correlations, which is the #133 ``gene == {gene}`` vacuity shape.

Note ``nt_site_boundary`` is deliberately **not** part of the join key even
though both sides carry it: it is the configured light-switch boundary, and every
shipped config sets it to 21562 (see CLAUDE.md pattern 8). Joining on a
configuration constant would make the key look more specific than it is while
adding nothing, and would break the join outright if upstream ever retuned it.

The headline numbers this produced on v9 (pooled) are in the README beside the
vendored upstream tables. Two are worth restating because they are what justifies
the era split at all:

* the two upstream **eras differ from each other more** (median |log ratio|
  0.1317) than mat-pcp's pooled fit differs from either of them (0.087, 0.070),
  so the era effect is larger than the implementation-plus-data noise floor; and
* upstream's own per-era rates rank **G->T first by 2.3x** (+0.616 mean
  log(pre_O/omicron) against +0.271 for the next), independently reproducing
  #137 Task 1's conclusion from a different quantity -- Task 1 measured
  synonymous residuals against a pooled model, this compares two fitted rate
  tables.
"""

import argparse
import os

import numpy as np
import pandas as pd

#: Join key. See the module docstring on why nt_site_boundary is excluded.
KEY = ["mut_type", "motif", "unpaired", "nt_site_before_boundary"]
GRID_CELLS = 512
ERAS = ("pre_omicron", "omicron")


def load_grid(path, label):
    """A rate table, validated as a full one-row-per-cell grid."""
    table = pd.read_csv(path)
    missing = [c for c in KEY + ["rate"] if c not in table.columns]
    if missing:
        raise ValueError(f"{label} ({path}) lacks required column(s) {missing}")
    if len(table) != GRID_CELLS:
        raise ValueError(
            f"{label} ({path}) has {len(table)} rows, expected the full "
            f"{GRID_CELLS}-cell grid. A partial table would make every statistic "
            f"below an average over an unstated subset."
        )
    if table.duplicated(KEY).any():
        raise ValueError(f"{label} ({path}) is not unique on {KEY}")
    if not (table["rate"] > 0).all():
        raise ValueError(
            f"{label} ({path}) has non-positive rates; the log-ratio statistics "
            f"below would be undefined"
        )
    return table[KEY + ["rate"]].rename(columns={"rate": label})


def spearman(a, b):
    """Rank correlation without scipy (not in the default env)."""
    return float(np.corrcoef(pd.Series(a).rank(), pd.Series(b).rank())[0, 1])


def agreement(frame, left, right):
    log_ratio = np.log(frame[left] / frame[right])
    return {
        "left": left,
        "right": right,
        "n_cells": len(frame),
        "pearson_log": round(
            float(np.corrcoef(np.log(frame[left]), np.log(frame[right]))[0, 1]), 6
        ),
        "spearman": round(spearman(frame[left], frame[right]), 6),
        "median_abs_log_ratio": round(float(log_ratio.abs().median()), 6),
        "p95_abs_log_ratio": round(float(log_ratio.abs().quantile(0.95)), 6),
        "max_abs_log_ratio": round(float(log_ratio.abs().max()), 6),
    }


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--pooled",
        help="mat-pcp glm/master_table.csv from a pooled (pre-Task-2) run. "
        "Compared against BOTH upstream eras.",
    )
    parser.add_argument("--pre-omicron", help="mat-pcp pre-Omicron master table")
    parser.add_argument("--omicron", help="mat-pcp Omicron master table")
    parser.add_argument(
        "--upstream-dir",
        default=os.path.join("analysis", "137-era-offsets", "upstream"),
        help="directory holding master_table_{pre_omicron,omicron}.csv from "
        "neherlab/SARS2-mut-fitness-v2 (vendored; see that directory's README)",
    )
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    if not args.pooled and not (args.pre_omicron and args.omicron):
        parser.error(
            "give either --pooled, or both --pre-omicron and --omicron. Before "
            "#137 Task 2 only a pooled table exists; after it, compare era to "
            "matching era."
        )

    upstream = {
        era: load_grid(
            os.path.join(args.upstream_dir, f"master_table_{era}.csv"),
            f"upstream_{era}",
        )
        for era in ERAS
    }

    merged = upstream["pre_omicron"].merge(
        upstream["omicron"], on=KEY, how="inner", validate="1:1"
    )
    ours = {}
    if args.pooled:
        ours["matpcp_pooled"] = load_grid(args.pooled, "matpcp_pooled")
    else:
        ours["matpcp_pre_omicron"] = load_grid(args.pre_omicron, "matpcp_pre_omicron")
        ours["matpcp_omicron"] = load_grid(args.omicron, "matpcp_omicron")
    for label, table in ours.items():
        merged = merged.merge(table, on=KEY, how="inner", validate="1:1")

    # The join must not have quietly shrunk. Every statistic below is an average
    # over these rows, so a subset would report healthy correlations over whatever
    # happened to match -- the vacuity shape CLAUDE.md pattern 10 records.
    if len(merged) != GRID_CELLS:
        raise ValueError(
            f"the grids joined to {len(merged)} cells, not {GRID_CELLS}. mat-pcp "
            f"and upstream disagree on the {KEY} grid, so these tables are not "
            f"comparable cell-by-cell and the statistics would be computed over an "
            f"unstated subset."
        )

    pairs = []
    if args.pooled:
        pairs += [("matpcp_pooled", "upstream_pre_omicron"),
                  ("matpcp_pooled", "upstream_omicron")]
    else:
        pairs += [("matpcp_pre_omicron", "upstream_pre_omicron"),
                  ("matpcp_omicron", "upstream_omicron"),
                  ("matpcp_pre_omicron", "matpcp_omicron")]
    pairs.append(("upstream_pre_omicron", "upstream_omicron"))

    rows = [agreement(merged, a, b) for a, b in pairs]
    summary = pd.DataFrame(rows)
    os.makedirs(args.outdir, exist_ok=True)
    summary.to_csv(os.path.join(args.outdir, "rate_agreement.csv"), index=False)

    print(f"grid joined 1:1 on all {GRID_CELLS} cells\n")
    print(summary.to_string(index=False))

    # The era difference per mutation type, from upstream's fitted rates. This is
    # independent of Task 1, which measured synonymous residuals against a pooled
    # model rather than comparing two fitted rate tables.
    merged["upstream_log_ratio"] = np.log(
        merged["upstream_pre_omicron"] / merged["upstream_omicron"]
    )
    by_type = (
        merged.groupby("mut_type")["upstream_log_ratio"]
        .agg(["mean", "median", "count"])
        .sort_values("mean", ascending=False)
        .reset_index()
    )
    if args.pooled:
        between = (
            merged["matpcp_pooled"].between(
                merged[["upstream_pre_omicron", "upstream_omicron"]].min(axis=1),
                merged[["upstream_pre_omicron", "upstream_omicron"]].max(axis=1),
            )
        ).mean()
        print(
            f"\nmat-pcp pooled lies between the two upstream era rates in "
            f"{100 * between:.1f}% of cells"
        )
    else:
        merged["matpcp_log_ratio"] = np.log(
            merged["matpcp_pre_omicron"] / merged["matpcp_omicron"]
        )
        by_type = by_type.merge(
            merged.groupby("mut_type")["matpcp_log_ratio"]
            .mean()
            .rename("matpcp_mean")
            .reset_index(),
            on="mut_type",
        )
    by_type.to_csv(
        os.path.join(args.outdir, "era_ratio_by_mut_type.csv"), index=False
    )
    print("\n--- era difference per mut_type: log(pre_omicron rate / omicron rate) ---")
    print(by_type.to_string(index=False, float_format=lambda x: f"{x:+.4f}"))
    print(f"\nwrote tables to {args.outdir}/")


if __name__ == "__main__":
    main()
