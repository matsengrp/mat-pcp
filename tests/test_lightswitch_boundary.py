#!/usr/bin/env python3
"""Pin the light-switch boundary convention across all three codepaths (#125).

The GLM splits each mutation type into two regimes at a configured genome
position. Through v8 the *training* predicate and the *serving* label disagreed
at the boundary site itself:

* ``modules/glm.py`` trains with ``global_context = (nt_site > boundary)``, so
  ``site == boundary`` lands in the before-bucket;
* ``scripts/annotate_counts.py`` labelled ``nt_site_before_boundary = (site <
  boundary)``, so ``site == boundary`` was served the *after*-bucket's rate.

``not (site < b)`` is ``site >= b``, not ``site > b`` -- the gap is exactly one
site. Both label producers now use ``modules/lightswitch.py::before_boundary``
(``site <= boundary``).

What this test pins, and why each check exists:

1. **Complementarity against the real training predicate.** Calls
   ``glm.GeneralLinearModel._compute_global_context_feature`` -- the actual
   shipped code, not a restatement of it -- and asserts
   ``before_boundary == (global_context == 0)`` at every site around the
   boundary. This is the invariant the bug violated, and it is the one check
   that would fail if someone "fixed" #125 by flipping ``glm.py`` instead.
2. **The boundary site specifically.** Checks 1 and 3 both pass for a
   consistently-wrong-by-one convention, so the boundary site is asserted
   directly.
3. **The two label producers agree.** ``annotate_counts.py``'s per-row scalar
   path and ``expected_counts.py``'s vectorized path are separate call sites
   (EPIC #131 F2). Fixing one alone leaves the clade and windowed fitness paths
   disagreeing at the boundary site with nothing to catch it, since #127's
   fixture pins *stored* labels rather than computed ones.
4. **Per-mutation-type boundary lookup.** ``CG`` and ``GC`` used to read the
   ``AT`` key, so a config setting them differently was honored in training and
   ignored at serving. Every shipped config sets AT/CG/GC to 21562, which masks
   it -- and ``test_glm_reference_equivalence.py``'s config check pins those
   values, so that guard *hid* this. This check therefore uses deliberately
   *distinct* per-type boundaries, which no config does.
5. **Unbounded mutation types stay False.** ``glm.py`` adds no
   ``global_context`` factor and ``rates.py`` builds only
   ``nt_site_before_boundary=False`` grid rows for them, so a ``True`` label
   would have no rate row to match. This is the one behavior that could have
   regressed when the hardcoded 21562/13467 fallbacks were removed from
   ``light_switch``.

``annotate_counts.py`` is a Snakemake script that reads a module-level
``snakemake`` object at import time, and exercising its ``light_switch`` would
mean building count/founder/RNA-structure fixtures for the whole script. Check 3
instead extracts that one function by ``exec``-ing the shipped file with the
module-level statements stripped, so the assertion is against the shipped
predicate rather than a copy of it. See ``_light_switch_from_script``.

Usage:
    python tests/test_lightswitch_boundary.py
"""

import ast
import sys
import traceback
from importlib import util as importlib_util
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
WORKFLOW_DIR = REPO_ROOT / "workflow"

# Deliberately distinct per-type values, unlike every shipped config (check 4).
DISTINCT_BOUNDARIES = {"CT": 5000, "AT": 10000, "CG": 15000, "GC": 20000}

# What the shipped configs actually use, so the real values are exercised too.
SHIPPED_BOUNDARIES = {"CT": 13467, "AT": 21562, "CG": 21562, "GC": 21562}

# Mutation types with no boundary feature: rates.mut_types minus
# rates.DEFAULT_LIGHTSWITCH_BOUNDARIES. Listed explicitly rather than derived so
# that a change to either set has to be acknowledged here; the cross-check in
# check_unbounded_mut_types_are_false enforces that this stays in sync.
UNBOUNDED_MUT_TYPES = ["AC", "AG", "CA", "GA", "GT", "TA", "TC", "TG"]


def _load_module(name, path):
    spec = importlib_util.spec_from_file_location(name, path)
    module = importlib_util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _load_lightswitch():
    return _load_module("lightswitch", WORKFLOW_DIR / "modules" / "lightswitch.py")


def _load_glm():
    return _load_module("glm", WORKFLOW_DIR / "modules" / "glm.py")


def _light_switch_from_script():
    """Extract ``light_switch`` from the shipped ``annotate_counts.py``.

    The script runs top-level statements against a module-level ``snakemake``
    object, so it cannot simply be imported. Rather than restate the predicate
    here -- which would pass even if the shipped one regressed -- this executes
    the file's *prologue* (everything before the first reference to ``snakemake``,
    i.e. the imports and the by-path module loader bootstrap that ``light_switch``
    depends on) plus every function definition after it.

    Both failure modes are loud rather than silent: if the loader bootstrap moves
    below the ``snakemake`` block the prologue stops providing ``lightswitch`` and
    the exec raises, and if ``light_switch`` stops being a module-level function
    the explicit check below raises.
    """
    path = WORKFLOW_DIR / "scripts" / "annotate_counts.py"
    tree = ast.parse(path.read_text())

    def _reads_snakemake(node):
        return any(
            isinstance(sub, ast.Name) and sub.id == "snakemake"
            for sub in ast.walk(node)
        )

    first_snakemake = next(
        (i for i, node in enumerate(tree.body) if _reads_snakemake(node)), None
    )
    if first_snakemake is None:
        raise AssertionError(
            f"{path.name} no longer reads a module-level snakemake object; "
            "this extraction helper needs revisiting"
        )

    keep = list(tree.body[:first_snakemake]) + [
        node
        for node in tree.body[first_snakemake:]
        if isinstance(node, ast.FunctionDef)
    ]
    namespace = {"__file__": str(path)}
    exec(
        compile(ast.Module(body=keep, type_ignores=[]), str(path), "exec"),
        namespace,
    )
    if "light_switch" not in namespace:
        raise AssertionError(
            "light_switch is no longer a module-level function in "
            "annotate_counts.py; this test would otherwise pass vacuously"
        )
    return namespace["light_switch"]


def check_complement_of_training_predicate():
    """before_boundary must be the exact complement of glm.py's global_context."""
    failures = []
    lightswitch = _load_lightswitch()
    glm = _load_glm()
    import pandas as pd

    model = glm.GeneralLinearModel(included_factors=["global_context"])

    for mut_type, boundary in SHIPPED_BOUNDARIES.items():
        sites = [boundary - 2, boundary - 1, boundary, boundary + 1, boundary + 2]
        df = pd.DataFrame({"nt_site": sites})
        global_context = model._compute_global_context_feature(
            df, mut_type, SHIPPED_BOUNDARIES
        )
        if global_context is None:
            failures.append(
                f"glm returned no global_context for {mut_type}, which has a "
                f"configured boundary at {boundary}"
            )
            continue
        if len(global_context) != len(sites):
            failures.append(
                f"{mut_type}: glm returned {len(global_context)} labels for "
                f"{len(sites)} sites; zip would silently compare a prefix"
            )
            continue
        for site, gc in zip(sites, global_context):
            label = lightswitch.before_boundary(site, boundary)
            expected = gc == 0
            if bool(label) != bool(expected):
                failures.append(
                    f"{mut_type} site {site} (boundary {boundary}): "
                    f"before_boundary={label} but glm global_context={gc} "
                    f"(bucket 'before' means global_context==0)"
                )
    if not failures:
        print("  before_boundary is the exact complement of glm's global_context")
    return failures


def check_boundary_site_is_before():
    """The boundary site itself must be labelled True (the fitted bucket).

    Asserted directly because an off-by-one convention that is *consistently*
    wrong would satisfy the complementarity check above.
    """
    failures = []
    lightswitch = _load_lightswitch()
    for mut_type, boundary in SHIPPED_BOUNDARIES.items():
        if not lightswitch.before_boundary(boundary, boundary):
            failures.append(
                f"{mut_type}: site {boundary} == boundary is labelled "
                f"after-boundary, but glm.py fits it into the before-bucket "
                f"(this is the v8 bug)"
            )
        if lightswitch.before_boundary(boundary + 1, boundary):
            failures.append(
                f"{mut_type}: site {boundary + 1} is labelled before-boundary"
            )
        if not lightswitch.before_boundary(boundary - 1, boundary):
            failures.append(
                f"{mut_type}: site {boundary - 1} is labelled after-boundary"
            )
    if not failures:
        print("  boundary site is labelled before-boundary; its neighbours flip")
    return failures


def _windowed_labels_from_module(mut_type, boundary, sites):
    """Run ``expected_counts.py``'s *own* label computation over ``sites``.

    The windowed producer is ``ExpectedCounts._add_context_features``, which
    merges a per-mut_type ``nt_site_boundary`` column and then calls
    ``lightswitch.before_boundary``. Calling ``before_boundary`` directly instead
    -- which this check did originally -- compares ``annotate_counts.py`` against
    the *shared module* and never executes ``expected_counts.py`` at all: a
    two-producer test wearing a one-producer implementation. Reverting the call
    site to ``the_df.site < the_df.nt_site_boundary`` (the v8 bug, in the path
    CLAUDE.md pattern 8 says is uncovered) left all five test files green.

    The instance is built with ``object.__new__`` and given only the two
    attributes the method reads, because ``__init__`` demands a full coding-sites
    table and reference genome that have nothing to do with the label. The method
    under test is the real one, unpatched.
    """
    import pandas as pd

    expected_counts = _load_module(
        "expected_counts", WORKFLOW_DIR / "modules" / "expected_counts.py"
    )
    engine = object.__new__(expected_counts.ExpectedCounts)
    engine.basepair_df = pd.DataFrame({"site": sites, "unpaired": [0] * len(sites)})
    boundary_series = pd.Series(
        {mut_type: boundary}, name="nt_site_boundary", dtype="int64"
    )
    boundary_series.index.name = "mut_type"
    engine.boundary_df = boundary_series

    the_df = pd.DataFrame({
        "site": sites,
        "wt_nt": [mut_type[0]] * len(sites),
        "mut_nt": [mut_type[1]] * len(sites),
    })
    # Long enough for the motif slice at site 29903; the motif is not asserted.
    annotated = engine._add_context_features(the_df, "A" * 29910)
    return annotated["nt_site_before_boundary"].tolist()


def check_both_label_producers_agree():
    """The clade path and the windowed path must compute the same label.

    Both sides are extracted from the shipped files -- ``annotate_counts.py`` for
    the clade path, ``expected_counts.py`` for the windowed one -- so that fixing
    the convention in one and not the other fails here. That asymmetry is the
    whole point of #125: the two paths must label the boundary site identically or
    the ``evol_time`` normalization stops making clade and windowed fitness
    comparable.
    """
    failures = []
    light_switch = _light_switch_from_script()

    for boundaries in (SHIPPED_BOUNDARIES, DISTINCT_BOUNDARIES):
        for mut_type, boundary in boundaries.items():
            sites = [1, boundary - 1, boundary, boundary + 1, 29903]
            windowed = _windowed_labels_from_module(mut_type, boundary, sites)
            if len(windowed) != len(sites):
                failures.append(
                    f"{mut_type}: windowed producer returned {len(windowed)} "
                    f"labels for {len(sites)} sites, so the comparison below "
                    f"would silently skip rows"
                )
                continue
            for site, window_label in zip(sites, windowed):
                clade_label = light_switch(mut_type, site, boundaries)
                if bool(clade_label) != bool(window_label):
                    failures.append(
                        f"{mut_type} site {site} (boundary {boundary}): clade path "
                        f"says {clade_label}, windowed path says {window_label}"
                    )
    if not failures:
        print("  annotate_counts.py and expected_counts.py agree at every site "
              "tested, including the boundary site itself")
    return failures


def check_per_mut_type_boundary_lookup():
    """Each mutation type must read its own configured boundary, not AT's."""
    failures = []
    light_switch = _light_switch_from_script()

    at_boundary = DISTINCT_BOUNDARIES["AT"]
    for mut_type, boundary in DISTINCT_BOUNDARIES.items():
        if not light_switch(mut_type, boundary, DISTINCT_BOUNDARIES):
            failures.append(
                f"{mut_type}: site {boundary} (its own boundary) not labelled "
                f"before-boundary"
            )
        if light_switch(mut_type, boundary + 1, DISTINCT_BOUNDARIES):
            failures.append(
                f"{mut_type}: site {boundary + 1} labelled before-boundary"
            )

    # The specific v8 defect: CG and GC read the AT key. With AT=10000 and
    # CG=15000, a site between them would be mislabelled under the old code.
    for mut_type in ("CG", "GC"):
        probe = at_boundary + 1
        own_boundary = DISTINCT_BOUNDARIES[mut_type]
        if not light_switch(mut_type, probe, DISTINCT_BOUNDARIES):
            failures.append(
                f"{mut_type}: site {probe} is before its own boundary "
                f"({own_boundary}) but was labelled after-boundary -- the "
                f"AT-key lookup bug is back"
            )
    if not failures:
        print("  each mutation type reads its own boundary (CG/GC no longer read AT)")
    return failures


def check_unbounded_mut_types_are_false():
    """Mutation types with no configured boundary must never be labelled True.

    ``rates.py`` builds only ``nt_site_before_boundary=False`` grid rows for
    these, so a True label has no rate row to match.
    """
    failures = []
    light_switch = _light_switch_from_script()

    for mut_type in UNBOUNDED_MUT_TYPES:
        for site in (1, 13467, 21562, 29903):
            if light_switch(mut_type, site, SHIPPED_BOUNDARIES):
                failures.append(
                    f"{mut_type} has no configured boundary but site {site} was "
                    f"labelled before-boundary"
                )

    # rates.DEFAULT_LIGHTSWITCH_BOUNDARIES must stay consistent with the list
    # above, or this check silently stops covering the real unbounded types.
    rates = _load_module("rates", WORKFLOW_DIR / "modules" / "rates.py")
    defaults = set(rates.DEFAULT_LIGHTSWITCH_BOUNDARIES)
    unbounded = set(rates.mut_types) - defaults
    if unbounded != set(UNBOUNDED_MUT_TYPES):
        failures.append(
            f"UNBOUNDED_MUT_TYPES is {sorted(UNBOUNDED_MUT_TYPES)} but rates.py "
            f"implies {sorted(unbounded)}; update this test"
        )
    if not failures:
        print("  unbounded mutation types are labelled False everywhere")
    return failures


CHECKS = [
    ("complement of glm's training predicate", check_complement_of_training_predicate),
    ("boundary site belongs to the before-bucket", check_boundary_site_is_before),
    ("both label producers agree", check_both_label_producers_agree),
    ("per-mutation-type boundary lookup", check_per_mut_type_boundary_lookup),
    ("unbounded mutation types stay False", check_unbounded_mut_types_are_false),
]


def run_all_checks():
    failures = []
    for label, check in CHECKS:
        print(f"Checking {label}...")
        try:
            failures.extend(check())
        except Exception:
            failures.append(f"{label} raised:\n{traceback.format_exc()}")
    return failures


def test_lightswitch_boundary():
    """Pytest entry point.

    This repo runs tests as standalone scripts (``pixi run python tests/<f>.py``)
    and does not depend on pytest. This wrapper exists so that anyone who does
    reach for ``pytest tests/`` gets a real result rather than a green run that
    collected nothing.
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
    print(
        "\nPASSED: the light-switch boundary site is labelled into the same "
        "bucket it is fitted into, in both label producers."
    )


if __name__ == "__main__":
    main()
