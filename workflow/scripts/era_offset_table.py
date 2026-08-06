#!/usr/bin/env python
"""Task 1 of #137: the all-clade x mut_type synonymous offset table.

Read-only diagnostic over a completed run's ``fitness/ntmut_fitness_aggregated.csv``.
Not a pipeline rule — run it by hand:

    pixi run python workflow/scripts/era_offset_table.py \
        --fitness results/gisaid-and-public_2026-04-22_v9/fitness/ntmut_fitness_aggregated.csv \
        --outdir analysis/137-era-offsets

Synonymous mutations are a free control: they are excluded from the amino-acid
fitness analyses, and under a well-calibrated neutral model the residual should be
~0. So for each (clade, mutation type) cell define the **synonymous offset**

    offset(clade, mut_type) = log( sum actual_count / sum expected_count )

over synonymous nucleotide mutations of that type in that clade. A model that
needs an era split will show that offset separating by era.

**The decision rule is fixed before the table is read** (#137, #138 step 2): for
each mutation type, require *complete separation*,
``min(offset over pre-Omicron cells) > max(offset over Omicron cells)``. G->T is
expected to pass — it is already measured at a gap of 0.393 with no overlap — so
this is a confirmation. C->G and G->C are expected to *fail*, which is evidence
that an era split is necessary but not sufficient, not a reason to widen #137.
If G->T fails, the premise of the whole v10 release is wrong.

**The ``source_clade != "all"`` filter is load-bearing.**
``aggregate_ntmut_fitness`` is fed ``extra_clades=["all"]`` (``Snakefile:1542``),
so the file carries an extra block of cross-clade pooled counts in which ``clade``
is NaN, because ``ntmut_counts_by_clade.py:107`` drops ``clade`` from that
aggregation. Grouping on ``clade`` gets the right answer only because pandas drops
NaN group keys by default. Re-implementing this on ``source_clade``, or with
``dropna=False``, double-counts every mutation into a bogus extra row. This script
therefore filters explicitly *and* asserts that no NaN clade survives, so the
correctness does not rest on a pandas default.

Cells below ``--min-count`` observed counts are excluded from the separation
verdict but **listed, not silently dropped** — a truncated table reads as
"covered everything". Cells with zero counts yield ``-inf`` and are reported
separately again rather than propagated.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "modules"))

import era as era_mod  # noqa: E402

USECOLS = [
    "clade",
    "source_clade",
    "mut_type",
    "actual_count",
    "expected_count",
    "synonymous",
]


def load_offsets(fitness_csv, chunksize):
    """Sum synonymous actual/expected counts per (clade, mut_type).

    Chunked because the production file is 1.26 GB / 4.6 M rows.
    """
    partials = []
    n_rows = n_kept = n_pooled = 0
    for chunk in pd.read_csv(
        fitness_csv, usecols=USECOLS, chunksize=chunksize, low_memory=False
    ):
        n_rows += len(chunk)
        pooled = chunk["source_clade"] == "all"
        n_pooled += int(pooled.sum())
        chunk = chunk[(chunk["synonymous"] == True) & ~pooled]  # noqa: E712
        n_kept += len(chunk)
        if chunk.empty:
            continue
        partials.append(
            chunk.groupby(["clade", "mut_type"], as_index=False)[
                ["actual_count", "expected_count"]
            ].sum()
        )
    if not partials:
        raise ValueError(f"no synonymous non-pooled rows in {fitness_csv}")

    grouped = (
        pd.concat(partials)
        .groupby(["clade", "mut_type"], as_index=False)[
            ["actual_count", "expected_count"]
        ]
        .sum()
    )
    # Not a pandas default we rely on: the pooled block carries clade == NaN, and
    # dropping it above must be what removes it, not groupby's dropna.
    if grouped["clade"].isna().any():
        raise ValueError(
            "NaN clade survived the source_clade != 'all' filter -- the pooled "
            "cross-clade block is being counted as a clade"
        )
    with np.errstate(divide="ignore", invalid="ignore"):
        grouped["offset"] = np.log(
            grouped["actual_count"] / grouped["expected_count"]
        )
    stats = {
        "rows_read": n_rows,
        "rows_kept": n_kept,
        "rows_source_clade_all": n_pooled,
    }
    return grouped, stats


def separation_verdict(cells, min_count):
    """Complete-separation verdict per mutation type over qualifying cells."""
    rows = []
    for mut_type, sub in cells.groupby("mut_type"):
        qualifying = sub[
            (sub["actual_count"] >= min_count) & np.isfinite(sub["offset"])
        ]
        pre = qualifying[qualifying["era"] == era_mod.PRE_OMICRON]["offset"]
        omi = qualifying[qualifying["era"] == era_mod.OMICRON]["offset"]
        if pre.empty or omi.empty:
            rows.append(
                {
                    "mut_type": mut_type,
                    "n_pre_omicron": len(pre),
                    "n_omicron": len(omi),
                    "pre_min": np.nan,
                    "pre_max": np.nan,
                    "omicron_min": np.nan,
                    "omicron_max": np.nan,
                    "gap": np.nan,
                    "separates": False,
                    "note": "one era has no qualifying cell -- cannot separate",
                }
            )
            continue
        gap = pre.min() - omi.max()
        rows.append(
            {
                "mut_type": mut_type,
                "n_pre_omicron": len(pre),
                "n_omicron": len(omi),
                "pre_min": pre.min(),
                "pre_max": pre.max(),
                "omicron_min": omi.min(),
                "omicron_max": omi.max(),
                "gap": gap,
                "separates": bool(gap > 0),
                "note": "",
            }
        )
    return pd.DataFrame(rows).sort_values("gap", ascending=False)


def weighted_era_means(cells):
    """Count-weighted mean offset per (era, mut_type).

    Weighted by ``actual_count`` so a 588-sample clade does not carry the same
    weight as a 4.3 M-sample one.
    """
    finite = cells[np.isfinite(cells["offset"])]
    rows = []
    for (era, mut_type), sub in finite.groupby(["era", "mut_type"]):
        rows.append(
            {
                "era": era,
                "mut_type": mut_type,
                "n_cells": len(sub),
                "total_actual": sub["actual_count"].sum(),
                "weighted_mean_offset": np.average(
                    sub["offset"], weights=sub["actual_count"]
                ),
                "unweighted_mean_offset": sub["offset"].mean(),
            }
        )
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fitness", required=True, help="ntmut_fitness_aggregated.csv")
    parser.add_argument("--outdir", required=True, help="directory for output tables")
    parser.add_argument(
        "--min-count",
        type=int,
        default=100,
        help="minimum summed actual_count for a cell to enter the verdict",
    )
    parser.add_argument("--chunksize", type=int, default=2_000_000)
    args = parser.parse_args()

    n_audited = era_mod.audit_clade_era_mapping()
    print(f"era mapping: {n_audited} clades, year-prefix audit passed")

    cells, stats = load_offsets(args.fitness, args.chunksize)
    print(
        f"read {stats['rows_read']:,} rows; kept {stats['rows_kept']:,} synonymous "
        f"non-pooled; excluded {stats['rows_source_clade_all']:,} "
        f"source_clade == 'all'"
    )

    cells["era"] = cells["clade"].map(era_mod.era_for_clade)
    unmapped = sorted(set(cells["clade"]) - set(era_mod.CLADE_ERA))
    if unmapped:
        raise KeyError(f"clades in the fitness table with no era: {unmapped}")

    cells = cells.sort_values(["clade", "mut_type"]).reset_index(drop=True)
    os.makedirs(args.outdir, exist_ok=True)

    n_clades = cells["clade"].nunique()
    n_types = cells["mut_type"].nunique()
    zero_count = cells[cells["actual_count"] == 0]
    sub_threshold = cells[
        (cells["actual_count"] < args.min_count) & (cells["actual_count"] > 0)
    ]
    qualifying = cells[cells["actual_count"] >= args.min_count]

    print(
        f"\n{n_clades} clades x {n_types} mutation types = {len(cells)} cells\n"
        f"  {len(qualifying)} cells with actual_count >= {args.min_count}\n"
        f"  {len(sub_threshold)} cells below threshold (listed in "
        f"sub_threshold_cells.csv, not dropped silently)\n"
        f"  {len(zero_count)} cells with zero counts (offset -inf, excluded)"
    )

    cells.to_csv(os.path.join(args.outdir, "offset_cells.csv"), index=False)
    sub_threshold.to_csv(
        os.path.join(args.outdir, "sub_threshold_cells.csv"), index=False
    )
    zero_count.to_csv(os.path.join(args.outdir, "zero_count_cells.csv"), index=False)

    # Clade x mut_type matrix in Nextstrain clade order, era marked.
    matrix = cells.pivot(index="clade", columns="mut_type", values="offset")
    matrix.insert(0, "era", [era_mod.era_for_clade(c) for c in matrix.index])
    matrix.round(4).to_csv(os.path.join(args.outdir, "offset_matrix.csv"))

    verdict = separation_verdict(cells, args.min_count)
    verdict.round(4).to_csv(os.path.join(args.outdir, "separation_verdict.csv"), index=False)

    means = weighted_era_means(cells)
    means.round(4).to_csv(os.path.join(args.outdir, "era_weighted_means.csv"), index=False)

    print("\n=== complete-separation verdict "
          "(min(pre-Omicron offset) > max(Omicron offset)) ===")
    print(
        verdict[
            [
                "mut_type",
                "n_pre_omicron",
                "n_omicron",
                "pre_min",
                "pre_max",
                "omicron_min",
                "omicron_max",
                "gap",
                "separates",
            ]
        ].to_string(index=False, float_format=lambda v: f"{v:7.3f}")
    )

    passing = verdict[verdict["separates"]]["mut_type"].tolist()
    failing = verdict[~verdict["separates"]]["mut_type"].tolist()
    print(f"\nseparates ({len(passing)}/{len(verdict)}): {', '.join(passing)}")
    print(f"does not separate ({len(failing)}): {', '.join(failing)}")

    gt = verdict[verdict["mut_type"] == "GT"]
    if gt.empty:
        raise ValueError("no G->T row in the verdict -- cannot evaluate the gate")
    gt_separates = bool(gt.iloc[0]["separates"])
    gt_gap = float(gt.iloc[0]["gap"])
    print(f"\nG->T: separates={gt_separates}, gap={gt_gap:.3f}")
    if not gt_separates:
        print(
            "\n*** GATE FAILED ***\n"
            "G->T does not separate completely by era. This is the pre-stated kill "
            "criterion for #137 and therefore for the whole v10 era split: the "
            "release's premise is that G->T's miscalibration is an era effect. Stop "
            "and escalate before any fit code is written.",
            file=sys.stderr,
        )
        return 1

    print(
        f"\nGATE PASSED: G->T separates with a gap of {gt_gap:.3f}. "
        f"Proceeding to Task 2 is justified."
    )
    print(f"\nwrote 6 tables to {args.outdir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
