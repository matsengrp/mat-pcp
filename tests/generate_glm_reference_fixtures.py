#!/usr/bin/env python3
"""Regenerate the GLM reference-equivalence fixtures for #126.

Builds two files under ``tests/data``:

* ``glm_curated_counts_fixture.csv`` -- a small stratified subsample of a real
  ``glm/curated_counts.csv``, retaining all four ``mut_class`` values so that a
  test training on it can distinguish synonymous-only from pooled training.
* ``glm_reference_rates.csv`` -- the rate table that the *reference*
  implementation (``SARS2-mut-fitness-v2``) produces from that fixture.

This script is the only place the reference repository is read. It is run by
hand when the reference implementation changes; CI only ever consumes the
committed outputs via ``tests/test_glm_reference_equivalence.py``.

Usage:
    python tests/generate_glm_reference_fixtures.py \
        --reference-repo ../SARS2-mut-fitness-v2 \
        --curated-counts results/test_data_v7/glm/curated_counts.csv
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Columns the fixture must carry: everything the reference load/GLM path reads,
# plus mut_class so tests can assert on the training-set composition.
FIXTURE_COLS = [
    "nt_mutation",
    "nt_site",
    "mut_type",
    "motif",
    "unpaired",
    "nt_site_before_boundary",
    "actual_count",
    "synonymous",
    "noncoding",
    "four_fold_degenerate",
    "mut_class",
    "gene",
    "clade_founder_aa",
    "mutant_aa",
]

# Keep every Nth row within each (mut_type, mut_class) stratum. Chosen so the
# fixture stays a few thousand rows while leaving every mutation type more
# synonymous observations than the GLM has parameters (8 without a light-switch
# boundary, 9 with one). The thinnest strata at this stride are GC and GT, at 45
# synonymous rows each.
STRIDE = 6

# Non-synonymous rows are thinned by a further factor on top of STRIDE. They are
# never trained on, so they do not affect the golden table at all -- they exist
# only so the tests can show that pooling them in *would* change the fit. A few
# thousand is ample for that, and dropping the rest takes the committed fixture
# from ~1.1 MB to ~0.4 MB.
NON_SYNONYMOUS_EXTRA_STRIDE = 4

# The reference implementation hard-codes these light-switch boundaries in
# modules/glm.py (13467 = end of ORF1ab, 21562 = end of S), and mat-pcp mirrors
# them in rates.DEFAULT_LIGHTSWITCH_BOUNDARIES. Used here only to check that the
# fixture spans them.
REFERENCE_BOUNDARIES = {"CT": 13467, "AT": 21562, "CG": 21562, "GC": 21562}


def build_fixture(curated_counts_path):
    """Stratified-subsample a curated counts table down to a test fixture."""
    df = pd.read_csv(curated_counts_path)

    missing = [c for c in FIXTURE_COLS if c not in df.columns]
    if missing:
        sys.exit(f"ERROR: {curated_counts_path} lacks columns {missing}")

    df = df[FIXTURE_COLS].copy()

    # Stratify by (mut_type, mut_class) so every mutation type keeps its
    # synonymous/nonsynonymous/nonsense/noncoding mix.
    # (cumcount rather than groupby.apply: the latter's include_groups
    # semantics differ across the pandas versions this repo supports.)
    df = df.sort_values("nt_mutation", kind="mergesort").reset_index(drop=True)
    within_stratum = df.groupby(["mut_type", "mut_class"], sort=False).cumcount()
    stride = np.where(
        df["synonymous"] == True,  # noqa: E712
        STRIDE,
        STRIDE * NON_SYNONYMOUS_EXTRA_STRIDE,
    )
    fixture = df[within_stratum % stride == 0].reset_index(drop=True)

    fixture = fixture.sort_values("nt_mutation", kind="mergesort").reset_index(drop=True)

    # ``Rates.populate_rates`` evaluates the light-switch factor at
    # ``df_train["nt_site"].max()``, and df_train is the *synonymous* subset. So
    # the check has to be on the synonymous rows, not on the fixture as a whole:
    # if the synonymous rows all sat below the largest boundary, s_max would too,
    # every grid row would get global_context = 0, and the True/False halves of
    # the AT/CG/GC/CT grid would collapse into duplicate rates -- silently
    # halving what the golden table actually pins. Appending a high-site row of
    # some other class would not help, since the loader filters it back out.
    max_boundary = max(REFERENCE_BOUNDARIES.values())
    syn_max_site = fixture.loc[fixture["synonymous"] == True, "nt_site"].max()  # noqa: E712
    if not syn_max_site > max_boundary:
        sys.exit(
            f"ERROR: highest synonymous site in the fixture is {syn_max_site}, "
            f"which does not exceed the largest light-switch boundary "
            f"({max_boundary}). The post-boundary half of the rate grid would "
            "not be exercised. Lower STRIDE or use a source table with "
            "synonymous mutations past that site."
        )

    # Both directions matter. Without non-synonymous rows a dropped filter would
    # satisfy every assertion; without synonymous rows the training set is empty
    # and the checks compare 0 == 0.
    classes = set(fixture["mut_class"].dropna().unique())
    if "synonymous" not in classes:
        sys.exit("ERROR: fixture has no synonymous rows; nothing to train on")
    if not classes - {"synonymous"}:
        sys.exit(
            "ERROR: fixture is entirely synonymous; it must contain other "
            "mut_class values or the synonymous-only test is vacuous"
        )
    return fixture


def reference_rates(reference_repo, fixture_path):
    """Run the reference implementation's GLM on the fixture."""
    reference_repo = Path(reference_repo).resolve()
    if not (reference_repo / "modules" / "rates.py").is_file():
        sys.exit(f"ERROR: no modules/rates.py under {reference_repo}")

    # The reference modules import each other as ``from modules import glm``,
    # so the repo root -- not modules/ -- goes on the path.
    sys.path.insert(0, str(reference_repo))
    from modules import load as ref_load
    from modules import rates as ref_rates

    # This is the reference contract under test: master_tables.py.ipynb trains
    # on load_synonymous_muts(...), NOT on the curated table as a whole.
    train = ref_load.load_synonymous_muts(str(fixture_path))
    print(f"reference trained on {len(train)} synonymous rows")

    rate = ref_rates.Rates()
    rate.populate_rates(train)
    return rate.rates.drop(columns=["condition"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference-repo",
        required=True,
        help="Path to a SARS2-mut-fitness-v2 clone",
    )
    parser.add_argument(
        "--curated-counts",
        required=True,
        help="A pipeline glm/curated_counts.csv to subsample",
    )
    parser.add_argument(
        "--outdir",
        default=str(Path(__file__).parent / "data"),
        help="Where to write the fixtures (default: tests/data)",
    )
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fixture_path = outdir / "glm_curated_counts_fixture.csv"
    rates_path = outdir / "glm_reference_rates.csv"

    fixture = build_fixture(args.curated_counts)
    fixture.to_csv(fixture_path, index=False)
    print(f"wrote {fixture_path} ({len(fixture)} rows)")
    print(fixture["mut_class"].value_counts().to_string())

    rates = reference_rates(args.reference_repo, fixture_path)
    rates.to_csv(rates_path, index=False)
    print(f"wrote {rates_path} ({len(rates)} rows)")


if __name__ == "__main__":
    main()
