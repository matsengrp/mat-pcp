#!/usr/bin/env python
"""#135: decompose ``build_count_table``'s "absent from the founder grid" drop.

Read-only diagnostic over a completed run. Not a pipeline rule:

    pixi run python workflow/scripts/audit_off_grid_drop.py \
        --run-dir results/gisaid-and-public_2026-04-22_v9 \
        --outdir analysis/135-off-grid-drop/v9

``build_count_table.py`` joins observed counts onto a grid enumerating every
possible mutation *away from each clade's founder base*, and logs two drops:

    dropped 7910 rows mutated codon background (14113 counts)       <- #130's filter
    dropped 46284 rows absent from the founder grid (210111 counts) <- pre-existing

The second is 5.9x the rows and 14.9x the counts of the filter #130 was written to
add, and has never been validated. This script reconstructs it from the same two
inputs and classifies every dropped row.

**Grid membership is a three-part predicate**, which is what makes the
decomposition tractable. ``(clade, nt_mutation)`` is on the grid iff

1. ``clade_founder_nts.csv`` has a row for ``(clade, nt_site)``;
2. the mutation's *from* base equals that clade founder base; and
3. the mutation's *to* base differs from it (the grid enumerates mutations away
   from the founder, so founder->founder is not a row it can have).

Each off-grid row therefore fails at least one part, and the cascade below assigns
it to exactly the first one it fails. The categories are asserted to partition the
drop, and the totals are reconciled against the run's own log line -- a
decomposition that does not add up is not a decomposition.

**Three causes named in #135 are structurally impossible here**, and the script
asserts the first two rather than just claiming it:

* *Outside ``site_include_range``* and *masked site*: ``compute_exclusions`` runs in
  step 3 of ``build_count_table.py``, **after** the join, and sets the ``exclude`` /
  ``masked_in_usher`` flags on grid rows. It never removes a row, so it cannot move
  an observation off the grid.
* *Non-coding*: the founder grid covers all 29,903 sites including non-coding ones
  (``gene == "noncoding"``), so being non-coding does not withhold a grid row.
* *No founder codon*: a codon-level notion. Grid membership here is decided purely
  at the nucleotide level; the codon-level filter is #130's, reported separately.
"""

import argparse
import os
import sys

import pandas as pd

ACGT = frozenset("ACGT")

# Ordered cascade. Each row is assigned the first category whose predicate holds,
# which is what makes the categories mutually exclusive by construction rather
# than by hope.
CATEGORY_ORDER = [
    "no_founder_row",
    "founder_nt_not_acgt",
    "mutation_nt_not_acgt",
    "self_mutation",
    "reversion_to_clade_founder",
    "off_founder_nt_background",
    "unexplained",
]

CATEGORY_NOTES = {
    "no_founder_row": (
        "clade_founder_nts.csv has no row for this (clade, nt_site), so no grid "
        "row could exist at any base"
    ),
    "founder_nt_not_acgt": (
        "the clade founder base is not A/C/G/T, so the grid's rows at this site "
        "carry a from-base no observation can match"
    ),
    "mutation_nt_not_acgt": (
        "the observed mutation's from- or to-base is not A/C/G/T (indel or "
        "ambiguity notation)"
    ),
    "self_mutation": "the observed from- and to-base are identical (a no-op mutation)",
    "reversion_to_clade_founder": (
        "the observed to-base IS the clade founder base: a reversion *toward* the "
        "founder. The grid enumerates mutations *away from* the founder, so it has "
        "no such row by construction -- this is excluded by design, not missing"
    ),
    "off_founder_nt_background": (
        "the from-base is neither the founder base nor reverting to it: a lineage "
        "mutated founder->X and then X->Y. The nucleotide-level analogue of #130's "
        "codon-level mutated-background case"
    ),
    "unexplained": "fails grid membership for none of the reasons above",
}


def classify(observed):
    """Assign each off-grid row exactly one category, first match wins."""
    from_nt = observed["from_nt"]
    to_nt = observed["to_nt"]
    founder = observed["grid_founder_nt"]

    category = pd.Series("unexplained", index=observed.index, dtype=object)
    unassigned = pd.Series(True, index=observed.index)

    def assign(name, mask):
        nonlocal unassigned
        hit = unassigned & mask.fillna(False)
        category[hit] = name
        unassigned &= ~hit

    assign("no_founder_row", founder.isna())
    assign("founder_nt_not_acgt", ~founder.isin(ACGT))
    assign("mutation_nt_not_acgt", ~from_nt.isin(ACGT) | ~to_nt.isin(ACGT))
    assign("self_mutation", from_nt == to_nt)
    assign("reversion_to_clade_founder", to_nt == founder)
    assign("off_founder_nt_background", from_nt != founder)
    return category


def load_run(run_dir):
    """Observed counts joined to the founder grid's base, plus the off-grid mask."""
    founder_csv = os.path.join(run_dir, "clade_founder_nts.csv")
    actual_csv = os.path.join(run_dir, "mutation_counts", "aggregated.csv")
    for path in (founder_csv, actual_csv):
        if not os.path.exists(path):
            raise FileNotFoundError(path)

    founder = pd.read_csv(founder_csv, usecols=["clade", "site", "nt", "gene"]).rename(
        columns={"site": "nt_site", "nt": "grid_founder_nt", "gene": "grid_gene"}
    )
    usecols = ["clade", "nt_mutation", "nt_site", "count"]
    header = pd.read_csv(actual_csv, nrows=0).columns
    for extra in ("clade_founder_nt", "reference_nt", "noncoding", "protein"):
        if extra in header:
            usecols.append(extra)
    actual = pd.read_csv(actual_csv, usecols=usecols)

    actual["from_nt"] = actual["nt_mutation"].str[0]
    actual["to_nt"] = actual["nt_mutation"].str[-1]

    merged = actual.merge(founder, on=["clade", "nt_site"], how="left")
    if len(merged) != len(actual):
        raise ValueError(
            f"joining founder bases changed the observation count from "
            f"{len(actual)} to {len(merged)}; clade_founder_nts.csv is not unique "
            f"on (clade, site)"
        )

    on_grid = (
        merged["grid_founder_nt"].notna()
        & (merged["from_nt"] == merged["grid_founder_nt"])
        & (merged["to_nt"] != merged["grid_founder_nt"])
    )
    return merged, ~on_grid


def logged_totals(run_dir):
    """The run's own drop totals, for reconciliation. ``None`` if unavailable.

    v8 predates #130, so its log has no drop lines at all and the off-grid drop was
    entirely silent there -- which is the point of #135's criterion 4.
    """
    log_path = os.path.join(run_dir, "log", "build_count_table.log")
    if not os.path.exists(log_path):
        return None
    marker = "rows absent from the founder grid ("
    with open(log_path) as handle:
        for line in handle:
            if marker in line:
                head, tail = line.split(marker)
                rows = int(head.strip().split()[-1])
                counts = int(tail.split()[0].lstrip("("))
                return rows, counts
    return None


def summarize(frame, keys, count_col="count"):
    out = frame.groupby(keys, dropna=False).agg(
        rows=("nt_mutation", "size"), counts=(count_col, "sum")
    )
    return out.reset_index().sort_values("counts", ascending=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, help="results/<tree>_<version>")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    merged, off_grid_mask = load_run(args.run_dir)
    off = merged[off_grid_mask].copy()
    n_rows = len(off)
    n_counts = int(off["count"].sum())

    print(f"run: {args.run_dir}")
    print(f"observed rows: {len(merged):,}")
    print(f"off-grid rows: {n_rows:,}  counts: {n_counts:,}")

    logged = logged_totals(args.run_dir)
    if logged is None:
        print(
            "  (this run's log has no 'absent from the founder grid' line -- the "
            "drop was silent here, so there is nothing to reconcile against)"
        )
    else:
        if logged != (n_rows, n_counts):
            raise ValueError(
                f"reconciliation FAILED: this script computes {n_rows} rows / "
                f"{n_counts} counts but the run's own log reports {logged[0]} rows "
                f"/ {logged[1]} counts. The reconstruction does not match the "
                f"pipeline, so the decomposition below describes something else."
            )
        print(f"  reconciles exactly against the run's log line: {logged[0]:,} rows / {logged[1]:,} counts")

    off["category"] = classify(off)

    # Exhaustive and mutually exclusive, asserted rather than assumed.
    unknown = set(off["category"]) - set(CATEGORY_ORDER)
    if unknown:
        raise ValueError(f"classify() produced categories outside the cascade: {unknown}")
    if off["category"].isna().any():
        raise ValueError("some off-grid rows were left unclassified")
    by_cat = summarize(off, ["category"])
    if int(by_cat["rows"].sum()) != n_rows or int(by_cat["counts"].sum()) != n_counts:
        raise ValueError(
            f"categories sum to {int(by_cat['rows'].sum())} rows / "
            f"{int(by_cat['counts'].sum())} counts, expected {n_rows} / {n_counts}"
        )

    os.makedirs(args.outdir, exist_ok=True)
    by_cat["pct_rows"] = (100 * by_cat["rows"] / n_rows).round(2)
    by_cat["pct_counts"] = (100 * by_cat["counts"] / n_counts).round(2)
    by_cat["note"] = by_cat["category"].map(CATEGORY_NOTES)
    by_cat.to_csv(os.path.join(args.outdir, "by_category.csv"), index=False)

    print("\n=== decomposition (every row in exactly one category) ===")
    print(
        by_cat[["category", "rows", "pct_rows", "counts", "pct_counts"]].to_string(
            index=False
        )
    )

    summarize(off, ["clade", "category"]).to_csv(
        os.path.join(args.outdir, "by_clade_category.csv"), index=False
    )
    gene_key = "grid_gene" if "grid_gene" in off.columns else "protein"
    summarize(off, [gene_key, "category"]).to_csv(
        os.path.join(args.outdir, "by_gene_category.csv"), index=False
    )
    summarize(off, ["category", "from_nt", "to_nt", "grid_founder_nt"]).to_csv(
        os.path.join(args.outdir, "by_base_change.csv"), index=False
    )
    summarize(off, ["clade"]).to_csv(
        os.path.join(args.outdir, "by_clade.csv"), index=False
    )

    # Highest-count individual drops, for eyeballing whether any single site
    # dominates -- binned artifacts elsewhere, but the tail matters here.
    cols = [
        c
        for c in ["clade", "nt_mutation", "nt_site", gene_key, "grid_founder_nt", "count", "category"]
        if c in off.columns
    ]
    off.nlargest(200, "count")[cols].to_csv(
        os.path.join(args.outdir, "top_200_by_count.csv"), index=False
    )

    unexplained = off[off["category"] == "unexplained"]
    print(
        f"\nunexplained: {len(unexplained):,} rows / "
        f"{int(unexplained['count'].sum()):,} counts"
    )
    if not unexplained.empty:
        unexplained.to_csv(
            os.path.join(args.outdir, "unexplained_rows.csv"), index=False
        )
        print("  written to unexplained_rows.csv -- these need their own issue")

    print(f"\nwrote tables to {args.outdir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
