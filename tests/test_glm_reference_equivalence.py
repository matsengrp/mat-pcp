#!/usr/bin/env python3
"""Pin the GLM neutral-rate model to the SARS2-mut-fitness-v2 reference (#126).

Two properties are checked against frozen fixtures in ``tests/data`` -- no
pipeline run and no reference clone required, so this is safe to run early in CI:

1. **The GLM trains on synonymous mutations only.** ``glm/curated_counts.csv``
   intentionally holds every mutation class, and the synonymous restriction
   lives in the *loader* (``modules/load.py::load_curated_muts``) rather than in
   ``curate_counts.py``. That is easy to lose by "simplifying"
   ``master_table.py`` to a bare ``pd.read_csv``, which would silently train the
   neutral model on mutations under purifying selection. See #126.

2. **The rates are numerically identical to the reference implementation.**
   mat-pcp's ``modules/rates.py`` + ``modules/glm.py`` are a port of
   ``SARS2-mut-fitness-v2``; this asserts the port still agrees to floating-point
   noise on a fixed input.

Two further checks keep those two honest: one confirms the fixture could actually
detect a dropped filter, so 1 and 2 cannot pass vacuously; the other asserts the
shipped configs still carry the reference light-switch boundaries, since the
pipeline passes those explicitly rather than using the module default, and
equivalence holds only at the reference's values.

Why a frozen fixture instead of pipeline output: ``Rates.populate_rates``
rescales rates by ``n_training_rows / total_mutations``, so every rate moves with
tree size. Comparing the two implementations on *the same* curated counts removes
that variable and makes the comparison exact rather than approximate. Regenerate
the fixtures with ``tests/generate_glm_reference_fixtures.py``.

Expect a ``WARNING: N conditions have negative rates`` block on stdout even when
this test passes. The fixture is a subsample, so the thinnest strata (45
synonymous rows against 8-9 GLM parameters) extrapolate to negative rates for
some conditions, and the frozen reference table contains them too. That is fine
here: these checks pin *agreement with the reference*, not the plausibility of
the fit, and the reference produces the same values.

Usage:
    python tests/test_glm_reference_equivalence.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
FIXTURE = REPO_ROOT / "tests" / "data" / "glm_curated_counts_fixture.csv"
GOLDEN = REPO_ROOT / "tests" / "data" / "glm_reference_rates.csv"

# The reference implementation hard-codes these light-switch boundaries in
# modules/glm.py (13467 = end of ORF1ab, 21562 = end of S). Reference
# equivalence only holds at these values, so the golden table is tied to them.
REFERENCE_BOUNDARIES = {"CT": 13467, "AT": 21562, "CG": 21562, "GC": 21562}

# Rate-grid identity columns, matching the reference's `condition` tuple.
GRID_COLS = ["mut_type", "motif", "unpaired", "nt_site_before_boundary"]

# Both implementations solve the same normal equations in float64, but not in a
# bit-identical operation order, so allow float noise and nothing more.
RATE_RTOL = 1e-9
RATE_ATOL = 1e-12

# Floor for check_guard_is_sensitive. Currently 433 of 512 conditions move when
# the synonymous filter is dropped; this leaves headroom for fixture churn while
# still failing if the fixture stops discriminating.
MIN_SENSITIVITY_CONDITIONS = 400


def _load_workflow_modules():
    """Import workflow/modules/{load,rates}.py the way pipeline scripts do."""
    sys.path.insert(0, str(REPO_ROOT / "workflow" / "modules"))
    from workflow_loader import load_workflow_module

    workflow_dir = str(REPO_ROOT / "workflow")
    return (
        load_workflow_module("load", workflow_dir),
        load_workflow_module("rates", workflow_dir),
    )


def _normalize(df):
    """Put grid columns in comparable dtypes before joining.

    Currently no-ops: ``read_csv`` infers ``bool`` for these columns and
    ``Rates`` already emits them that way. Kept so the join key cannot silently
    become dtype-sensitive. Do not reuse this idiom on a nullable column --
    ``astype(bool)`` on object-dtype strings makes ``"False"`` truthy.
    """
    df = df.copy()
    df["nt_site_before_boundary"] = df["nt_site_before_boundary"].astype(bool)
    df["unpaired"] = df["unpaired"].astype(int)
    df["mut_type"] = df["mut_type"].astype(str)
    df["motif"] = df["motif"].astype(str)
    return df


def check_training_set_is_synonymous_only(load_module):
    """The loader used by master_table.py must yield synonymous rows only."""
    failures = []
    fixture = pd.read_csv(FIXTURE)

    # Guard against a vacuous test: if the fixture were synonymous-only, a
    # dropped filter would still pass every assertion below.
    classes = set(fixture["mut_class"].dropna().unique())
    non_syn = classes - {"synonymous"}
    if not non_syn:
        failures.append(
            "fixture contains only synonymous mutations, so it cannot detect a "
            "dropped synonymous filter -- regenerate it"
        )
        return failures
    # The converse also has to hold: with zero synonymous rows the checks below
    # would compare 0 == 0 and pass for the wrong reason.
    if "synonymous" not in classes:
        failures.append(
            "fixture contains no synonymous mutations, so the checks below would "
            "pass vacuously on an empty training set -- regenerate it"
        )
        return failures
    print(f"  fixture mut_class values: {sorted(classes)} ({len(fixture)} rows)")

    train = load_module.load_curated_muts(str(FIXTURE))
    n_expected = int((fixture["synonymous"] == True).sum())  # noqa: E712

    if len(train) != n_expected:
        failures.append(
            f"load_curated_muts returned {len(train)} rows, expected the "
            f"{n_expected} rows with synonymous == True. The GLM neutral model "
            "must not be trained on non-synonymous mutations (#126)."
        )

    bad_flag = int((train["synonymous"] != True).sum())  # noqa: E712
    if bad_flag:
        failures.append(
            f"{bad_flag} training rows have synonymous != True (#126)"
        )

    bad_class = sorted(set(train["mut_class"].dropna().unique()) - {"synonymous"})
    if bad_class:
        failures.append(
            f"training set contains mut_class values {bad_class}; only "
            "synonymous mutations may inform a neutral rate (#126)"
        )

    if not failures:
        print(f"  trained on {len(train)}/{len(fixture)} rows, all synonymous")
    return failures


def _fit_rates(rates_module, train):
    """Fit the rate grid on a given training frame."""
    rate = rates_module.Rates(REFERENCE_BOUNDARIES)
    rate.populate_rates(train)
    return _normalize(rate.rates)


def check_rates_match_reference(load_module, rates_module):
    """mat-pcp's GLM must reproduce the reference rate table on the fixture."""
    failures = []
    golden = _normalize(pd.read_csv(GOLDEN))

    # Exercise the module defaults as well as the values: the golden table is
    # only meaningful at the reference's boundaries.
    defaults = getattr(rates_module, "DEFAULT_LIGHTSWITCH_BOUNDARIES", None)
    if defaults != REFERENCE_BOUNDARIES:
        failures.append(
            f"DEFAULT_LIGHTSWITCH_BOUNDARIES is {defaults}, reference uses "
            f"{REFERENCE_BOUNDARIES}"
        )

    actual = _fit_rates(rates_module, load_module.load_curated_muts(str(FIXTURE)))

    if len(actual) != len(golden):
        failures.append(
            f"rate grid has {len(actual)} rows, reference has {len(golden)}"
        )
        return failures

    merged = golden[GRID_COLS + ["rate"]].merge(
        actual[GRID_COLS + ["rate"]],
        on=GRID_COLS,
        how="outer",
        suffixes=("_ref", "_actual"),
        indicator=True,
    )
    unmatched = merged[merged["_merge"] != "both"]
    if len(unmatched):
        failures.append(
            f"{len(unmatched)} rate-grid conditions differ between "
            f"implementations, e.g.\n{unmatched.head().to_string()}"
        )
        return failures

    close = np.isclose(
        merged["rate_actual"], merged["rate_ref"], rtol=RATE_RTOL, atol=RATE_ATOL
    )
    if not close.all():
        worst = merged[~close].assign(
            absdiff=(merged["rate_actual"] - merged["rate_ref"]).abs()
        ).nlargest(5, "absdiff")
        failures.append(
            f"{(~close).sum()}/{len(merged)} rates deviate from the reference "
            f"implementation, worst offenders:\n{worst.to_string()}"
        )
    else:
        max_diff = np.abs(merged["rate_actual"] - merged["rate_ref"]).max()
        print(
            f"  all {len(merged)} rates match the reference "
            f"(max abs diff {max_diff:.3e})"
        )
    return failures


def check_configs_use_reference_boundaries():
    """Shipped configs must carry the reference's light-switch boundaries.

    The pipeline does not use ``DEFAULT_LIGHTSWITCH_BOUNDARIES``: ``master_table``
    and ``expected_counts_by_clade`` are passed ``config["glm_lightswitch_boundaries"]``.
    So pinning the module default is not enough -- a config-only edit would move
    every production rate away from the reference while the rest of this test,
    which fits at ``REFERENCE_BOUNDARIES`` directly, stayed green.
    """
    failures = []
    config_dir = REPO_ROOT / "config"
    configs = sorted(config_dir.glob("*.yaml"))
    if not configs:
        return [f"no configs found under {config_dir}"]

    for path in configs:
        with open(path) as f:
            cfg = yaml.safe_load(f) or {}
        if "glm_lightswitch_boundaries" not in cfg:
            continue
        boundaries = {str(k): int(v) for k, v in cfg["glm_lightswitch_boundaries"].items()}
        if boundaries != REFERENCE_BOUNDARIES:
            failures.append(
                f"{path.name} sets glm_lightswitch_boundaries to {boundaries}, "
                f"but reference equivalence holds only at {REFERENCE_BOUNDARIES}. "
                "Changing these is a deliberate deviation from "
                "SARS2-mut-fitness-v2 (see #125)."
            )
    if not failures:
        print(f"  all {len(configs)} configs use the reference boundaries")
    return failures


def check_guard_is_sensitive(load_module, rates_module):
    """The committed fixture must be able to tell the two training sets apart.

    Checks A and B both pass vacuously if the fixture's synonymous subset is
    close enough to the whole fixture that pooling barely moves the fit. Rather
    than trust that, fit the grid both ways on the *same* fixture -- once via
    ``load_curated_muts``, once via a bare ``pd.read_csv`` (the #126 scenario) --
    and require that most conditions actually differ.

    Note this asserts *divergence*, not that the pooled fit looks broken: on
    this sparse fixture the pooled rates are not uniformly lower and have fewer
    negative values than the correct ones, because the
    ``n_rows / total_mutations`` rescaling changes along with the training set.
    Direction is not the signal; disagreement is.
    """
    failures = []
    syn = _fit_rates(rates_module, load_module.load_curated_muts(str(FIXTURE)))
    pooled = _fit_rates(rates_module, pd.read_csv(str(FIXTURE)))

    merged = syn[GRID_COLS + ["rate"]].merge(
        pooled[GRID_COLS + ["rate"]], on=GRID_COLS, suffixes=("_syn", "_pooled")
    )
    differing = int(
        (
            ~np.isclose(
                merged["rate_pooled"],
                merged["rate_syn"],
                rtol=RATE_RTOL,
                atol=RATE_ATOL,
            )
        ).sum()
    )
    if differing < MIN_SENSITIVITY_CONDITIONS:
        failures.append(
            f"pooling non-synonymous mutations into the training set changes "
            f"only {differing}/{len(merged)} rate-grid conditions (require at "
            f"least {MIN_SENSITIVITY_CONDITIONS}), so checks A and B can no "
            "longer demonstrate that they would catch a dropped synonymous "
            "filter. Either the loader's filter was removed (check A's result "
            "first) or the fixture needs regenerating (#126)."
        )
    else:
        print(
            f"  pooling would move {differing}/{len(merged)} conditions, so the "
            "fixture can detect a dropped filter"
        )
    return failures


def run_all_checks():
    """Run every check and return the accumulated failure messages."""
    for path in (FIXTURE, GOLDEN):
        if not path.is_file():
            return [f"missing fixture {path}"]

    load_module, rates_module = _load_workflow_modules()

    failures = []
    print("Checking GLM training set is synonymous-only (#126)...")
    failures += check_training_set_is_synonymous_only(load_module)
    print("Checking GLM rates match SARS2-mut-fitness-v2...")
    failures += check_rates_match_reference(load_module, rates_module)
    print("Checking the fixture can detect a dropped synonymous filter...")
    failures += check_guard_is_sensitive(load_module, rates_module)
    print("Checking configs use the reference light-switch boundaries...")
    failures += check_configs_use_reference_boundaries()
    return failures


def test_glm_reference_equivalence():
    """Pytest entry point.

    This repo runs its tests as standalone scripts (``pixi run python
    tests/<file>.py``) and does not depend on pytest. This wrapper exists so that
    anyone who does reach for ``pytest tests/`` gets a real result: without a
    ``test_``-prefixed *function*, pytest would collect zero tests here and exit
    green, which reads as "the guard passed" when nothing ran.
    """
    failures = run_all_checks()
    assert not failures, "\n".join(failures)


def main():
    failures = run_all_checks()
    if failures:
        print(f"\nFAILED ({len(failures)} problem(s)):")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)

    print("\nPASSED: GLM is synonymous-only and matches the reference.")


if __name__ == "__main__":
    main()
