#!/usr/bin/env python3
"""Unit tests for the clade-subsampling arms (#151).

Runs the real ``workflow/scripts/subsample_clade_samples.py`` and
``subsample_correlations.py`` against a stub ``snakemake`` global, and the
mutation tests exec *deliberately corrupted* copies -- with every corruption applied to
``workflow/modules/subsample.py`` rather than to the script text. That
distinction is the point: an operator implementing the sampler works in the
module, so a test that patched script text would report the guard as covered
while the code it guards sat one file away (CLAUDE.md pattern 14's observation,
one axis over).

**What a pipeline run cannot check, and this file therefore must.** No shipped
production config sets ``clade_subsamples`` at all, and the CI config's clades
(19A at 2 samples, 19B at 5, 20A at 46) are far too small for the shipped caps
of 50000/20000/5000 to touch anything, so every property below is unreachable
from a production run and most are unreachable from CI too.

*Order preservation is invisible end-to-end.*
``validate_subsample_arms`` asserts set containment and the retained count, both
of which a reordering leaves intact, so a sampler that sorted its output -- or
emitted its picks in draw order -- passes every pipeline check. The fixture's
seven samples are therefore deliberately unsorted and the draw is partial, so
file order and draw order differ observably.

*The membership boundary is unreachable in CI.* Membership is ``n > cap``
strictly, and the difference from ``n >= cap`` shows only for a clade of exactly
``cap`` samples. No CI clade has 1000 or 20 samples, so a widened predicate
moves nothing there; the fixture caps at exactly ``len(FULL_SAMPLES)``.

*An empty draw is a hard failure and a below-``min_clade_size`` draw is an
exclusion, and CI reaches neither.* ``validate_entries`` forces ``cap >= 1`` and
``resolve_unit`` admits only ``n > cap``, so no accepted config can produce an
empty draw -- the guard is defence in depth against the failure that would
otherwise surface as a crash inside ``matUtils`` with nothing able to attribute
it to clade size (#149's SIGSEGV, one axis over). The exclusion needs
``min_clade_size > cap``, which no config sets -- and since the arms dropped
``subsample_summary.csv``, the exclusion is no longer recorded anywhere a run
could reconcile, so ``check_min_clade_size_excludes_rather_than_raises`` below
is its only coverage of any kind.

*The removed fraction keys cannot be rejected by any shipped config, because
none carries them.* ``sampling_fractions`` and ``unit`` belonged to the fraction
arms and are refused **by name** so that a config written against them says
what happened rather than failing as an unknown-key schema error. Nothing
end-to-end exercises that path.

Hermetic: synthetic sample lists and a synthetic output tree. Nothing reads a
pipeline output.
"""

import hashlib
import os
import shutil
import sys
import tempfile
import traceback
from importlib import util as importlib_util
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
DRAW_SCRIPT = REPO_ROOT / "workflow" / "scripts" / "subsample_clade_samples.py"
CORR_SCRIPT = REPO_ROOT / "workflow" / "scripts" / "subsample_correlations.py"
MODULE_REL = "modules/subsample.py"


def _load(module_name):
    spec = importlib_util.spec_from_file_location(
        module_name, REPO_ROOT / "workflow" / "modules" / f"{module_name}.py"
    )
    module = importlib_util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


subsample = _load("subsample")


class _Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def get(self, key, default=None):
        return self.__dict__.get(key, default)

    def __getitem__(self, key):
        return self.__dict__[key]


def _exec_script(path, script_globals):
    text = path.read_text()
    script_globals.setdefault("__name__", "__main__")
    script_globals.setdefault("__file__", str(path))
    exec(compile(text, str(path), "exec"), script_globals)


def _load_mutated(tmp, edits, module_name):
    """Import ``modules/subsample.py`` from a corrupted copy of ``workflow/``.

    Used by the checks that exercise the module directly rather than through a
    script. The corruption still lands in the *module*, which is the file an
    operator implementing the sampler edits.
    """
    root = Path(_mutated_workflow(tmp, edits)).parent
    spec = importlib_util.spec_from_file_location(
        module_name, root / "modules" / "subsample.py"
    )
    module = importlib_util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _mutated_workflow(tmp, edits):
    """A throwaway copy of ``workflow/`` with text edits applied.

    Copied rather than edited in place: the scripts resolve
    ``modules/<name>.py`` relative to their own location, so a mutation to a
    *module* can only be reached by relocating the pair -- and mutating the real
    file would corrupt the repo for any concurrent run and leave it broken if
    the process were killed.
    """
    root = Path(tmp) / "wf"
    shutil.copytree(REPO_ROOT / "workflow" / "modules", root / "modules")
    shutil.copytree(REPO_ROOT / "workflow" / "scripts", root / "scripts")
    for relative, (old, new) in edits.items():
        path = root / relative
        original = path.read_text()
        assert original.count(old) == 1, (
            f"could not locate the text to mutate in {relative}: {old[:70]!r}"
        )
        path.write_text(original.replace(old, new))
    return root / "scripts"


# ---------------------------------------------------------------------------
# Fixture
# ---------------------------------------------------------------------------
#
# Seven samples, deliberately NOT in sorted order: ``select_cap`` must return
# them in *this* order, so a sampler that sorted its output -- or emitted picks
# in draw order -- is caught. Seven is also one more than any cap used below, so
# every draw is partial and the ordering is observable.
FULL_SAMPLES = [
    "zeta/2021",
    "alpha/2020",
    "mu/2022",
    "beta/2020",
    "omega/2023",
    "gamma/2021",
    "delta/2021",
]
CLADE = "20A"
SEED = 0
# The cap used by every fixture below. Four of seven, so the draw is partial --
# which is what makes order preservation and subset containment observable at
# all -- and it is what ``_write_summary_fixture``'s ``n_arm`` default matches.
CAP = 4


def _write_samples(path, samples, trailing_newline=True):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(samples) + ("\n" if trailing_newline else ""))


def cap_spec(cap, seed=SEED):
    return {"seed": seed, "cap": cap}


def _run_draw_spec(tmp, spec, arm, clade=CLADE, scripts_dir=None,
                   samples=None, trailing_newline=True):
    samples = FULL_SAMPLES if samples is None else samples
    full = Path(tmp) / "mat" / "clades" / f"{clade}_samples.txt"
    _write_samples(full, samples, trailing_newline)
    out = Path(tmp) / "mat" / "clades" / arm / f"{clade}_samples.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    snake = _Namespace(
        input=_Namespace(samples=str(full)),
        output=_Namespace(samples=str(out)),
        params=_Namespace(spec=spec),
        wildcards=_Namespace(clade=clade, arm=arm),
        log=[],
    )
    script = DRAW_SCRIPT if scripts_dir is None else (
        Path(scripts_dir) / "subsample_clade_samples.py"
    )
    _exec_script(script, {"snakemake": snake})
    return out


def run_cap_draw(tmp, cap, clade=CLADE, scripts_dir=None, samples=None,
                 trailing_newline=True, seed=SEED):
    return _run_draw_spec(
        tmp,
        cap_spec(cap, seed),
        subsample.arm_name(cap, seed),
        clade=clade,
        scripts_dir=scripts_dir,
        samples=samples,
        trailing_newline=trailing_newline,
    )


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------
def check_a_cap_at_or_above_the_size_returns_the_list_verbatim():
    """``select_cap`` at or above ``len(samples)`` returns the list unchanged.

    Unreachable through the pipeline -- ``resolve_unit`` admits a clade only
    when ``n > cap`` and the sampler raises if an untouched one arrives anyway
    -- so this exercises the library invariant directly. It is what makes
    "a clade at or below the cap is not thinned" a property of the draw rather
    than only of the fan-out, and it is the one place a future sampler that
    drew every index instead of returning early would be caught.
    """
    failures = []
    n = len(FULL_SAMPLES)
    for cap in (n, n + 1, n * 10):
        got = subsample.select_cap(FULL_SAMPLES, cap, SEED, CLADE)
        if got != FULL_SAMPLES:
            failures.append(
                f"cap {cap} over {n} samples returned {got}, expected the "
                f"input list verbatim"
            )
    return failures


def check_draw_preserves_file_order():
    """A partial draw keeps the input's order, not the order it was drawn in.

    The fixture's seven samples are deliberately unsorted, so a sampler that
    sorted its output -- or emitted its picks in draw order -- lands in a
    different order than the input's and is caught. No pipeline run checks
    this: ``validate_subsample_arms`` asserts set containment and the retained
    count, both of which a reordering leaves intact.
    """
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        drawn = run_cap_draw(tmp, CAP).read_text().split("\n")
        drawn = [d for d in drawn if d]
        positions = [FULL_SAMPLES.index(d) for d in drawn]
        if positions != sorted(positions):
            failures.append(
                f"draw is not in file order: {drawn} sits at positions {positions}"
            )
        if not set(drawn) <= set(FULL_SAMPLES):
            failures.append(f"draw {drawn} is not a subset of the full list")
        if len(set(drawn)) != len(drawn):
            failures.append(f"draw {drawn} contains duplicates")
        if len(drawn) != CAP:
            failures.append(
                f"draw kept {len(drawn)} samples, expected min(7, {CAP}) = {CAP}"
            )
    return failures


def check_emitting_draw_order_is_caught():
    """Mutation test: drop the sort and the emitted order stops matching the file.

    Checked over several seeds rather than one. A four-of-seven draw whose
    indices happen to come out ascending is in file order by luck, so a single
    seed could leave the mutant looking correct -- the check needs at least one
    seed where the unsorted draw is observably out of order.
    """
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        scripts = _mutated_workflow(
            tmp,
            {MODULE_REL: (
                "idx = np.sort(rng.choice(n, size=cap, replace=False))",
                "idx = rng.choice(n, size=cap, replace=False)",
            )},
        )
        reordered = False
        for seed in range(8):
            out = run_cap_draw(tmp, CAP, scripts_dir=scripts, seed=seed)
            drawn = [d for d in out.read_text().split("\n") if d]
            positions = [FULL_SAMPLES.index(d) for d in drawn]
            if positions != sorted(positions):
                reordered = True
                break
        if not reordered:
            failures.append(
                "removing the sort from select_cap still emitted every draw in "
                "file order across 8 seeds, so the ordering guarantee is untested"
            )
    return failures
def check_seeding_from_iteration_order_is_caught():
    """Mutation test: collapse the RNG key and every clade draws alike.

    ``cap_rng``'s key is a string literal, so a refactor can rewrite it without
    any assertion noticing -- and a key that dropped the clade would give every
    clade of a given size the same sample list, which is well-formed and wrong.
    """
    failures = []
    scripts_edits = {MODULE_REL: (
        'key = "|".join([ARM_PREFIX, cap_token(cap), str(seed), clade])',
        'key = "|".join([ARM_PREFIX, str(seed)])',
    )}
    with tempfile.TemporaryDirectory() as tmp:
        module = _load_mutated(tmp, scripts_edits, "subsample_seedkey")
        wide = [f"s{i}" for i in range(200)]
        a = module.select_cap(wide, 100, SEED, CLADE)
        b = module.select_cap(wide, 100, SEED, "21J")
        # Under the mutation the key no longer mentions the clade or the cap, so
        # every arm draws the same indices. If the two draws still differ,
        # something *other* than ``cap_rng``'s key is keying the draw and the
        # purity check above is proving less than it claims.
        if a != b:
            failures.append(
                "collapsing the RNG key to the seed alone still gave different "
                "draws per clade, so something else is seeding the draw"
            )
    return failures
def check_min_clade_size_excludes_rather_than_raises():
    """The two thresholds are different in kind, and must stay so.

    ``matUtils extract -s`` on an empty sample list is undefined, so a draw of
    zero is a **hard failure**. ``min_clade_size`` says how big a clade the user
    wants, so a clade thinned below it is merely **excluded** -- and recorded,
    since an absent directory is otherwise indistinguishable from a wiring
    mistake. Reachable exactly when ``min_clade_size > cap``, which no shipped
    config sets.
    """
    failures = []
    arms = subsample.validate_entries([{"max_samples": [2], "seed": 0}])
    arm = next(iter(arms))
    kept = subsample.check_arm_sizes(arms, {"A": 10}, min_clade_size=5)
    if kept[arm]:
        failures.append(
            f"cap 2 below min_clade_size 5 was kept: {kept[arm]}"
        )
    kept = subsample.check_arm_sizes(arms, {"A": 10}, min_clade_size=2)
    if kept[arm] != {"A": 2}:
        failures.append(
            f"cap 2 at min_clade_size 2 should be kept, got {kept[arm]}"
        )
    # A clade at or below the cap is never in the unit, so it is neither kept
    # nor excluded -- it simply is not this arm's business.
    kept = subsample.check_arm_sizes(arms, {"A": 2, "B": 1}, min_clade_size=0)
    if kept[arm] != {}:
        failures.append(
            f"clades at or below the cap appeared in the arm: {kept[arm]}"
        )
    return failures


def check_an_empty_draw_is_a_hard_failure_not_an_exclusion():
    """A draw of zero raises where a small draw is excluded.

    ``validate_entries`` forces ``cap >= 1`` and ``resolve_unit`` admits only
    ``n > cap``, so ``check_arm_sizes`` cannot see an empty draw from any
    accepted config -- the spec here is hand-built to reach the branch. That is
    the point of keeping it: it is defence in depth for the one failure that
    would otherwise surface as a crash inside ``matUtils`` with nothing able to
    attribute it to clade size (#149's SIGSEGV, one axis over).
    """
    failures = []
    arms = {"cap_0_0": {"seed": 0, "cap": 0}}
    try:
        subsample.check_arm_sizes(arms, {"A": 10}, min_clade_size=0)
        failures.append(
            "a cap of 0 did not raise; matUtils extract -s on an empty list is "
            "undefined"
        )
    except ValueError:
        pass
    # And the same via the draw itself, which is the other route to an empty
    # list reaching matUtils.
    try:
        subsample.select_cap(FULL_SAMPLES, 0, SEED, CLADE)
        failures.append("select_cap at cap 0 did not raise")
    except ValueError:
        pass
    return failures


def check_treating_a_zero_draw_as_an_exclusion_is_caught():
    """Mutation test: demote the empty-draw failure and the distinction dies."""
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        module = _load_mutated(
            tmp,
            {MODULE_REL: (
                "            if k < 1:\n                empty.append((clade, n, k))\n"
                "            elif k >= min_clade_size:",
                "            if k >= min_clade_size and k >= 1:",
            )},
            "subsample_demoted",
        )
        try:
            module.check_arm_sizes(
                {"cap_0_0": {"seed": 0, "cap": 0}}, {"A": 10}, min_clade_size=0
            )
        except ValueError:
            failures.append(
                "demoting the empty-draw check still raised, so this check "
                "proves nothing"
            )
    return failures
def check_an_empty_full_list_raises_about_emptiness():
    """An empty full list raises, and about the *right* thing.

    Order matters here rather than merely the fact of raising. Every cap is
    ``>= 1``, so an empty list is also at or below the cap: if the membership
    check ran first, an empty list would be reported as "this clade is not one
    of the arm's" -- true, and the wrong bug. An empty list means the clade
    extraction disagrees with the checkpoint, which has a different owner.
    """
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        try:
            # ``trailing_newline=False`` so the file is genuinely 0 bytes,
            # which is what ``samples_by_clade``'s awk writes for a clade with
            # no rows. With a trailing newline it would hold one blank line --
            # a different, malformed input.
            run_cap_draw(tmp, CAP, samples=[], trailing_newline=False)
            failures.append("an empty full sample list did not raise")
        except ValueError as exc:
            if "is empty" not in str(exc):
                failures.append(
                    f"an empty list raised, but about membership rather than "
                    f"emptiness: {exc}"
                )
    return failures


# ---------------------------------------------------------------------------
# The roster param
# ---------------------------------------------------------------------------

def _roster(arm, cap, included, unit=None):
    """The param the Snakefile passes, for one arm over one clade.

    ``unit`` is the *resolved* list the fan-out used. There is no configured
    unit to fall back on -- membership is derived from the tree's clade sizes --
    so ``thinned`` has to come from the function that drove the fan-out rather
    than from a second evaluation of the sizing arithmetic.
    """
    return {
        arm: {
            "unit": list([CLADE] if unit is None else unit),
            "seed": SEED,
            "cap": cap,
            "included": included,
        }
    }


# ---------------------------------------------------------------------------
# correlation_summary.csv (#151): the headline result
# ---------------------------------------------------------------------------

def _write_fitness(path, rows):
    """One clade's aamut_fitness table. ``rows`` is (aa_mutation, actual, ratio, bayes)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            {
                "cluster": CLADE,
                "gene": "S",
                "aa_site": i + 1,
                "aa_mutation": mut,
                "clade_founder_aa": "A",
                "mutant_aa": "V",
                "expected_count": 0.5,
                "actual_count": actual,
                "delta_fitness_ratio": ratio,
                "subset_of_ORF1ab": False,
                "delta_fitness_bayesian": bayes,
                "uncertainty": 1.0,
            }
            for i, (mut, actual, ratio, bayes) in enumerate(rows)
        ]
    ).to_csv(path, index=False)


def _corr_fixture(root, arm, scored_disagree=True, arm_actual=12):
    """A grid where the ground-truth and whole-grid correlations MUST disagree.

    Eighty unscored rows (``actual_count`` 0 on both sides, ``expected_count``
    0.5) carry values that agree perfectly, as real unscored grid rows do --
    they are functions of the expected count, which is invariant to sampling
    depth. Six rows carry ``actual_count`` 12 at full depth, clearing the
    ground-truth threshold of 10, and values that are *anti*-correlated between
    the two depths.

    So the whole-grid Pearson r is near +1 and the ground-truth r is near -1. A
    script that correlated over the grid would report a confident positive
    number, which is precisely the mistake this fixture exists to catch.

    ``arm_actual`` sets the *arm's* count on those six rows. Lowering it below
    the threshold leaves the ground-truth set unchanged (it is defined on the
    full-depth row) while making an arm-side filter drop them -- which is how
    the survivorship-bias mutant is distinguished from the correct one.
    """
    root = Path(root)
    d = root / "aamut_fitness"
    unscored = [(f"A{i}V", 0, 0.1 * i, 0.2 * i) for i in range(80)]
    scored_full = [(f"S{i}V", 12, float(i), float(i)) for i in range(6)]
    scored_arm = (
        [(f"S{i}V", arm_actual, float(5 - i), float(5 - i)) for i in range(6)]
        if scored_disagree
        else [(f"S{i}V", arm_actual, float(i), float(i)) for i in range(6)]
    )
    _write_fitness(d / f"{CLADE}_aamut_fitness.csv", unscored + scored_full)
    _write_fitness(d / arm / f"{CLADE}_aamut_fitness.csv", unscored + scored_arm)
    return root


GT_POPULATION = "ground_truth_max_full>=10"


def run_correlations(root, roster, count_threshold=10):
    out = Path(root) / "correlation_summary.csv"
    snake = _Namespace(
        input=_Namespace(fitness=[]),
        output=_Namespace(summary=str(out)),
        params=_Namespace(
            roster=roster, outdir=str(root), count_threshold=count_threshold
        ),
        log=[],
    )
    _exec_script(CORR_SCRIPT, {"snakemake": snake})
    return pd.read_csv(out)


def check_correlation_uses_the_ground_truth_set_not_the_grid():
    """The population is the fixed ground-truth set, never the whole grid.

    Every grid row carries a fitness value even when nothing was observed, and
    those rows are invariant to sampling depth, so correlating over them
    measures the grid's construction and returns a confident number describing
    nothing. Measured on the real test tree: 0.619 over 121,427 grid rows where
    the 29 genuinely scored mutations give 1.000.
    """
    failures = []
    arm = subsample.arm_name(CAP, SEED)
    with tempfile.TemporaryDirectory() as tmp:
        _corr_fixture(tmp, arm)
        got = run_correlations(tmp, _roster(arm, CAP, [CLADE]))
        row = got[
            (got["estimator"] == "delta_fitness_ratio")
            & (got["population"] == GT_POPULATION)
        ].iloc[0]
        if int(row["n_mutations"]) != 6:
            failures.append(
                f"correlated over {row['n_mutations']} mutations, expected the "
                f"6 ground-truth ones (the fixture has 80 sub-threshold grid "
                f"rows that must be excluded)"
            )
        if row["pearson_r"] > -0.99:
            failures.append(
                f"pearson_r is {row['pearson_r']}, expected about -1.0. The "
                f"fixture's ground-truth rows are anti-correlated and its "
                f"sub-threshold rows agree perfectly, so a positive value means "
                f"the whole grid was correlated instead."
            )
        if int(row["n_ground_truth"]) != 6:
            failures.append(
                f"n_ground_truth is {row['n_ground_truth']}, expected 6"
            )
        # r is scale-free; the fixture's arm values are the full-depth ones
        # reversed, so the slope must be about -1 and both sds equal.
        if abs(row["regression_slope"] + 1.0) > 1e-6:
            failures.append(
                f"regression_slope is {row['regression_slope']}, expected -1.0 "
                f"for a fixture whose arm values are the full-depth values "
                f"reversed"
            )
        if int(row["n_scored_full"]) != 6 or int(row["n_scored_arm"]) != 6:
            failures.append(
                f"n_scored_full/n_scored_arm are {row['n_scored_full']}/"
                f"{row['n_scored_arm']}, expected 6/6 (all six scored "
                f"mutations are scored at both depths)"
            )
    return failures


def check_correlating_the_whole_grid_is_caught():
    """Mutating the script to drop the scored filter must move the answer."""
    failures = []
    arm = subsample.arm_name(CAP, SEED)
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp) / "tree"
        _corr_fixture(root, arm)
        scripts = _mutated_workflow(
            tmp,
            {
                "scripts/subsample_correlations.py": (
                    "        ground_truth = merged[gt_mask]",
                    "        ground_truth = merged",
                )
            },
        )
        out = root / "correlation_summary.csv"
        snake = _Namespace(
            input=_Namespace(fitness=[]),
            output=_Namespace(summary=str(out)),
            params=_Namespace(
                roster=_roster(arm, CAP, [CLADE]),
                outdir=str(root),
                count_threshold=10,
            ),
            log=[],
        )
        try:
            _exec_script(
                scripts / "subsample_correlations.py", {"snakemake": snake}
            )
        except Exception:
            return failures  # raising is an acceptable way to be caught
        got = pd.read_csv(out)
        row = got[
            (got["estimator"] == "delta_fitness_ratio")
            & (got["population"] == GT_POPULATION)
        ].iloc[0]
        if int(row["n_mutations"]) == 6 and row["pearson_r"] < -0.99:
            failures.append(
                "correlating the whole grid produced the same answer as "
                "correlating the ground-truth set, so the test cannot tell "
                "them apart"
            )
    return failures


def check_count_fraction_counts_observations_not_mutations():
    """``count_fraction`` sums observations; the rejected alternative counted mutations.

    A mutation falling from 900 observations to 3 has lost almost everything,
    and a *coverage* fraction -- the share of full-depth-scored mutations the
    arm still scores at all -- calls it fully retained. That alternative is not
    shipped (it diverges from ``count_fraction`` by up to 30x on real data and
    predicts the fitted slope worse: Spearman 0.940 against 0.985), but the
    mistake it represents is one line away from the correct implementation, so
    the fixture is built to catch it.

    All six scored mutations are still scored in the arm and each falls from 12
    counts to 3, so summing observations gives exactly 0.25 while counting
    distinct mutations gives exactly 1.0. The two readings are therefore
    separated by the widest possible margin on this fixture, rather than by a
    few percent.

    The denominator is the **whole grid**, not the ground-truth subset: the
    fraction describes the thinning, not the comparison population. The 80
    unscored rows contribute 0 to both sums, so they do not move it, which is
    what makes the expected value exactly 0.25 rather than approximately so.
    """
    failures = []
    arm = subsample.arm_name(CAP, SEED)
    with tempfile.TemporaryDirectory() as tmp:
        # Still scored on both sides, but at a quarter of the counts.
        _corr_fixture(tmp, arm, arm_actual=3)
        got = run_correlations(tmp, _roster(arm, CAP, [CLADE]))
        row = got[
            (got["estimator"] == "delta_fitness_ratio")
            & (got["population"] == GT_POPULATION)
        ].iloc[0]
        if int(row["sum_actual_full"]) != 72:
            failures.append(
                f"sum_actual_full is {row['sum_actual_full']}, expected 72 "
                f"(six scored mutations at 12 counts each, over the whole grid)"
            )
        if int(row["sum_actual_arm"]) != 18:
            failures.append(
                f"sum_actual_arm is {row['sum_actual_arm']}, expected 18 "
                f"(six scored mutations at 3 counts each)"
            )
        if abs(row["count_fraction"] - 0.25) > 1e-9:
            failures.append(
                f"count_fraction is {row['count_fraction']}, expected 0.25. "
                f"The arm keeps every scored mutation but only a quarter of "
                f"the observations, so a fraction that counted mutations "
                f"instead of counts would read 1.0"
            )
        if int(row["n_scored_arm"]) != 6:
            failures.append(
                f"n_scored_arm is {row['n_scored_arm']}, expected 6 -- the "
                f"fixture keeps every mutation scored so that summing "
                f"observations and counting mutations must disagree, so if "
                f"this has moved the two readings are no longer separated"
            )
    return failures


def check_ground_truth_set_is_fixed_across_sampling_depths():
    """The threshold is read off the **full-depth** row, never the arm's.

    Filtering on the arm's row is survivorship bias: it drops exactly the
    mutations the thinning destroyed, so the denominator shrinks with the cut
    and a degrading estimator reports a flat correlation over an ever-easier
    population. The two readings are indistinguishable on any fixture where the
    arm's counts match the full tree's, which is why this one drops the arm's
    counts below the threshold while leaving the full-depth counts above it: the
    correct implementation still compares 6 mutations, an arm-side filter
    compares 0.

    Nothing end-to-end can see this. Every coefficient is still well-formed
    under the wrong reading, and on a real tree the population it produces looks
    entirely plausible -- it is simply a different set at every depth.
    """
    failures = []
    arm = subsample.arm_name(CAP, SEED)
    with tempfile.TemporaryDirectory() as tmp:
        # arm_actual=1: still scored, but far below the threshold of 10.
        _corr_fixture(tmp, arm, arm_actual=1)
        got = run_correlations(tmp, _roster(arm, CAP, [CLADE]))
        row = got[
            (got["estimator"] == "delta_fitness_ratio")
            & (got["population"] == GT_POPULATION)
        ].iloc[0]
        if int(row["n_mutations"]) != 6:
            failures.append(
                f"the arm's counts fell below the threshold and the population "
                f"shrank to {row['n_mutations']}; the ground-truth set is "
                f"defined on the full-depth row and must stay at 6"
            )

    # And the survivorship reading must give a different answer, or this check
    # cannot tell the two apart.
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp) / "tree"
        _corr_fixture(root, arm, arm_actual=1)
        scripts = _mutated_workflow(
            tmp,
            {
                "scripts/subsample_correlations.py": (
                    '            merged[["actual_count_full", '
                    '"expected_count_full"]].max(axis=1)',
                    '            merged[["actual_count_arm", '
                    '"expected_count_arm"]].max(axis=1)',
                )
            },
        )
        out = root / "correlation_summary.csv"
        snake = _Namespace(
            input=_Namespace(fitness=[]),
            output=_Namespace(summary=str(out)),
            params=_Namespace(
                roster=_roster(arm, CAP, [CLADE]),
                outdir=str(root),
                count_threshold=10,
            ),
            log=[],
        )
        try:
            _exec_script(
                scripts / "subsample_correlations.py", {"snakemake": snake}
            )
            got = pd.read_csv(out)
            row = got[
                (got["estimator"] == "delta_fitness_ratio")
                & (got["population"] == GT_POPULATION)
            ].iloc[0]
            if int(row["n_mutations"]) == 6:
                failures.append(
                    "filtering on the arm's row gave the same population as "
                    "filtering on the full-depth row, so the test cannot "
                    "distinguish survivorship bias from the fixed set"
                )
        except Exception:
            pass  # raising is an acceptable way to be caught
    return failures


def check_correlation_keeps_empty_populations_as_rows():
    """An unmeasurable population is a row with n=0, not a missing row.

    A silently absent row reads as "not applicable"; an explicit ``n_mutations``
    of 0 next to a blank coefficient reads as "asked, and there was nothing to
    answer with". On the real test tree the ``>=10`` filter matches nothing, so
    this is the ordinary case rather than a corner.
    """
    failures = []
    arm = subsample.arm_name(CAP, SEED)
    with tempfile.TemporaryDirectory() as tmp:
        _corr_fixture(tmp, arm)
        # Threshold far above any count in the fixture, so the filtered
        # population is empty for both estimators.
        got = run_correlations(
            tmp, _roster(arm, CAP, [CLADE]), count_threshold=1000
        )
        filtered = got[got["population"].str.contains(">=1000")]
        if len(filtered) != 2:
            failures.append(
                f"expected 2 ground-truth rows (one per estimator), got "
                f"{len(filtered)} -- an empty population was dropped instead "
                f"of recorded"
            )
        if not (filtered["n_mutations"] == 0).all():
            failures.append(
                f"filtered rows report n_mutations "
                f"{filtered['n_mutations'].tolist()}, expected 0"
            )
        if not filtered["pearson_r"].isna().all():
            failures.append(
                f"filtered rows report a coefficient "
                f"{filtered['pearson_r'].tolist()} over an empty population"
            )
    return failures


def check_a_control_clade_that_disagrees_raises():
    """A non-unit clade must reproduce full depth, or the arm leaked into it.

    This is the fan-out bug (#151 development) turned into an assertion: an arm
    re-aggregates every clade in its grid but thins only its unit, so a clade
    outside the unit whose effects moved means the sampler reached a clade the
    config never named.
    """
    failures = []
    arm = subsample.arm_name(CAP, SEED)
    with tempfile.TemporaryDirectory() as tmp:
        _corr_fixture(tmp, arm)
        # Same tree, but the roster says this clade is NOT in the arm's unit --
        # so its (deliberately disagreeing) effects must be reported as a leak.
        roster = _roster(arm, CAP, ["21J"], unit=["21J"])
        try:
            run_correlations(tmp, roster)
            failures.append(
                "a control clade whose fitness effects differ from full depth "
                "did not raise; an arm reaching a clade outside its unit would "
                "go unreported"
            )
        except ValueError as exc:
            if "never meant to touch" not in str(exc):
                failures.append(f"raised, but not about the leak: {exc}")
    return failures



# ---------------------------------------------------------------------------
# Membership: which clades an arm thins
# ---------------------------------------------------------------------------


def check_cap_thins_only_clades_above_it():
    """``n > cap`` is drawn; ``n <= cap`` is returned unchanged and undrawn.

    The strictness matters at exactly one clade size and nowhere else: a clade
    of exactly ``cap`` samples is already at the target, so thinning it would
    write an identical file into the arm directory and make the arm claim a
    clade it did not thin. That clade is then no longer a control.
    """
    failures = []
    n = len(FULL_SAMPLES)
    for cap, expect_drawn in ((n - 1, True), (n, False), (n + 1, False)):
        drawn = subsample.select_cap(FULL_SAMPLES, cap, SEED, CLADE)
        if expect_drawn:
            if len(drawn) != cap:
                failures.append(
                    f"cap {cap} over {n} samples retained {len(drawn)}, "
                    f"expected {cap}"
                )
        elif drawn != FULL_SAMPLES:
            failures.append(
                f"cap {cap} over {n} samples altered the list; a clade at or "
                f"below the cap must be returned unchanged"
            )
    # The same predicate, at the level the fan-out uses.
    sizes = {"big": n + 1, "exact": n, "small": n - 1}
    unit = subsample.resolve_unit(cap_spec(n), sizes)
    if unit != ["big"]:
        failures.append(
            f"resolve_unit at cap {n} returned {unit}, expected ['big'] -- a "
            f"clade of exactly cap samples is untouched and must not be in the "
            f"unit"
        )
    return failures


def check_including_a_clade_at_the_cap_is_caught():
    """``n > cap`` mutated to ``n >= cap`` must change the resolved unit.

    Mutation-tested in ``modules/subsample.py``, not in the script: the operator
    implementing membership works in the module, and patching script text would
    report the guard as covered while the code it guards sat one file away.
    """
    failures = []
    n = len(FULL_SAMPLES)
    with tempfile.TemporaryDirectory() as tmp:
        mutated = _load_mutated(
            tmp,
            {
                MODULE_REL: (
                    "    return [c for c, n in clade_sizes.items() if int(n) > cap]",
                    "    return [c for c, n in clade_sizes.items() if int(n) >= cap]",
                )
            },
            "subsample_cap_boundary",
        )
        sizes = {"big": n + 1, "exact": n, "small": n - 1}
        unit = mutated.resolve_unit(cap_spec(n), sizes)
        if unit == ["big"]:
            failures.append(
                "widening the membership test to '>=' produced the same unit, "
                "so the boundary is untested"
            )
    return failures


def check_the_sampler_refuses_a_clade_at_or_below_the_cap():
    """A fan-out routing an untouched clade to the sampler is a hard failure.

    The sampler decides membership from the sample list it has already read,
    rather than from a resolved unit param -- a cap arm's unit depends on the
    ``get_clades`` checkpoint, and resolving one inside ``params`` would tie the
    rule's scheduling to it. The predicate must still be the same one.
    """
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        try:
            run_cap_draw(tmp, len(FULL_SAMPLES))
            failures.append(
                "the sampler thinned a clade at exactly the cap; an arm whose "
                "'unchanged' clades are quietly rewritten would look correct "
                "in every count"
            )
        except ValueError:
            pass
    with tempfile.TemporaryDirectory() as tmp:
        out = run_cap_draw(tmp, len(FULL_SAMPLES) - 1)
        got = out.read_text().rstrip("\n").split("\n")
        if len(got) != len(FULL_SAMPLES) - 1:
            failures.append(
                f"the sampler wrote {len(got)} samples at cap "
                f"{len(FULL_SAMPLES) - 1}"
            )
        if got != [s for s in FULL_SAMPLES if s in got]:
            failures.append(
                f"the cap draw did not preserve the input file's order: {got}"
            )
    return failures


def check_cap_draw_is_a_pure_function_of_cap_seed_clade():
    """An arm's list cannot move when another cap is added or removed.

    The key deliberately omits the **unit**: it is derived from the tree's clade
    sizes, so keying the draw on it would make every clade's sample list depend
    on how many other clades happened to exceed the cap. A tree gaining one
    large clade would then silently redraw every arm, and nothing downstream
    could attribute it.
    """
    failures = []
    samples = [f"s{i}/2021" for i in range(200)]
    baseline = None
    for caps in ([50], [50, 20], [20, 50], [5, 20, 50], [50]):
        arms = subsample.validate_entries([{"max_samples": caps, "seed": SEED}])
        name = subsample.arm_name(50, SEED)
        if name not in arms:
            failures.append(f"cap 50 missing from arms built with {caps}")
            continue
        drawn = subsample.select_cap(samples, arms[name]["cap"], SEED, CLADE)
        digest = hashlib.md5("\n".join(drawn).encode()).hexdigest()
        if baseline is None:
            baseline = digest
        elif digest != baseline:
            failures.append(
                f"the cap-50 arm's draw moved when the config listed {caps}"
            )
    # And two clades of the same size must draw differently, or the seed is not
    # reaching the per-clade key at all.
    arms = subsample.validate_entries([{"max_samples": [50], "seed": SEED}])
    spec = arms[subsample.arm_name(50, SEED)]
    if subsample.select_cap(samples, spec["cap"], SEED, "20A") == (
        subsample.select_cap(samples, spec["cap"], SEED, "21J")
    ):
        failures.append(
            "two different clades drew identically, so the clade is not part "
            "of the RNG key"
        )
    return failures


def check_cap_validation():
    """``max_samples`` must be a non-empty list of integers >= 1."""
    failures = []
    cases = {
        "no depth key": {"seed": 0},
        "cap not a list": {"max_samples": 10},
        "empty cap list": {"max_samples": []},
        "zero cap": {"max_samples": [0]},
        "negative cap": {"max_samples": [-5]},
        "float cap": {"max_samples": [50000.0]},
        "string cap": {"max_samples": ["50000"]},
        "bool cap": {"max_samples": [True]},
        "non-int seed": {"max_samples": [10], "seed": 0.5},
    }
    for label, entry in cases.items():
        try:
            subsample.validate_entries([entry])
            failures.append(f"{label} was accepted: {entry!r}")
        except ValueError:
            pass
    # Duplicate arm names must collide loudly rather than silently share outputs.
    try:
        subsample.validate_entries(
            [{"max_samples": [50000], "seed": 0}, {"max_samples": [50000], "seed": 0}]
        )
        failures.append("two entries yielding cap_50000_0 did not raise")
    except ValueError:
        pass
    arms = subsample.validate_entries([{"max_samples": [50000, 20000], "seed": 0}])
    if sorted(arms) != ["cap_20000_0", "cap_50000_0"]:
        failures.append(f"unexpected cap arm names: {sorted(arms)}")
    return failures


def check_check_arm_sizes_reports_the_cap_and_honours_min_clade_size():
    """``check_arm_sizes`` reports ``min(n, cap)`` and honours min_clade_size."""
    failures = []
    arms = subsample.validate_entries([{"max_samples": [100], "seed": 0}])
    sizes = {"big": 1000, "at_cap": 100, "small": 50}
    kept = subsample.check_arm_sizes(arms, sizes, min_clade_size=0)["cap_100_0"]
    if kept != {"big": 100}:
        failures.append(
            f"expected only the above-cap clade at its cap, got {kept}"
        )
    # A min_clade_size above the cap excludes rather than raises: the cap says
    # what can be built, min_clade_size says how big a clade the user wants.
    kept = subsample.check_arm_sizes(arms, sizes, min_clade_size=500)["cap_100_0"]
    if kept != {}:
        failures.append(
            f"a cap below min_clade_size should exclude the clade, got {kept}"
        )
    return failures


def check_the_removed_fraction_keys_are_rejected_by_name():
    """``sampling_fractions`` and ``unit`` must fail with an explanation.

    They belonged to the fraction arms, which applied one sampling fraction to
    an explicitly named unit of clades. A config written against them would
    otherwise fail as an unknown-key schema error, which says nothing about
    where the feature went -- and the wrong reading is available: an operator
    seeing "unknown key" is as likely to conclude the key is misspelled as that
    the mechanism is gone.
    """
    failures = []
    cases = {
        "fractions alone": {"unit": ["20A"], "sampling_fractions": [0.5]},
        "fractions beside a cap": {"max_samples": [10], "sampling_fractions": [0.5]},
        "unit beside a cap": {"unit": ["20A"], "max_samples": [10]},
    }
    for label, entry in cases.items():
        try:
            subsample.validate_entries([entry])
            failures.append(f"{label} was accepted: {entry!r}")
            continue
        except ValueError as exc:
            if "#151" not in str(exc) or "max_samples" not in str(exc):
                failures.append(
                    f"{label} raised, but without saying what replaced it: {exc}"
                )
    # An ordinary unknown key must still say "unknown", not point at the
    # fraction arms -- otherwise the message above is just the generic one.
    try:
        subsample.validate_entries([{"max_samples": [10], "oops": 1}])
        failures.append("an unknown key was accepted")
    except ValueError as exc:
        if "unknown" not in str(exc):
            failures.append(f"an unknown key raised the wrong message: {exc}")
    return failures


CHECKS = [
    ("a cap at or above the size returns the list verbatim", check_a_cap_at_or_above_the_size_returns_the_list_verbatim),
    ("draw preserves file order", check_draw_preserves_file_order),
    ("emitting draw order is caught", check_emitting_draw_order_is_caught),
    ("seeding from iteration order is caught", check_seeding_from_iteration_order_is_caught),
    ("min_clade_size excludes rather than raises", check_min_clade_size_excludes_rather_than_raises),
    ("an empty draw is a hard failure, not an exclusion", check_an_empty_draw_is_a_hard_failure_not_an_exclusion),
    ("demoting the zero-draw failure is caught", check_treating_a_zero_draw_as_an_exclusion_is_caught),
    ("an empty full list raises about emptiness", check_an_empty_full_list_raises_about_emptiness),
    ("the removed fraction keys are rejected by name", check_the_removed_fraction_keys_are_rejected_by_name),
    ("cap thins only clades above it", check_cap_thins_only_clades_above_it),
    ("including a clade at the cap is caught", check_including_a_clade_at_the_cap_is_caught),
    ("the sampler refuses a clade at or below the cap", check_the_sampler_refuses_a_clade_at_or_below_the_cap),
    ("the draw is a pure function of (cap, seed, clade)", check_cap_draw_is_a_pure_function_of_cap_seed_clade),
    ("cap validation", check_cap_validation),
    ("check_arm_sizes reports the cap and honours min_clade_size", check_check_arm_sizes_reports_the_cap_and_honours_min_clade_size),
    ("correlation uses the ground-truth set, not the grid", check_correlation_uses_the_ground_truth_set_not_the_grid),
    ("correlating the whole grid is caught", check_correlating_the_whole_grid_is_caught),
    ("count_fraction sums observations, not distinct mutations", check_count_fraction_counts_observations_not_mutations),
    ("the ground-truth set is fixed across sampling depths", check_ground_truth_set_is_fixed_across_sampling_depths),
    ("empty populations stay as rows with n=0", check_correlation_keeps_empty_populations_as_rows),
    ("a control clade that disagrees raises", check_a_control_clade_that_disagrees_raises),
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


def test_clade_subsamples():
    """pytest entry point.

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
        "\nPASSED: membership is strictly above the cap, the draw is a pure "
        "function of (cap, seed, clade) in file order, an empty draw is a hard "
        "failure while a below-min_clade_size draw is an exclusion, the "
        "correlations use a ground-truth set fixed at full depth over a "
        "population that does not move with the cut, and every guard fires on "
        "a deliberately broken module."
    )


if __name__ == "__main__":
    main()
