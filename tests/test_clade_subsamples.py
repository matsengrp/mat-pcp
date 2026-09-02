#!/usr/bin/env python3
"""Unit tests for the clade-subsampling arms (#151).

Runs the real ``workflow/scripts/subsample_clade_samples.py`` and
``subsample_summary.py`` against a stub ``snakemake`` global, and the mutation
tests exec *deliberately corrupted* copies -- with every corruption applied to
``workflow/modules/subsample.py`` rather than to the script text. That
distinction is the point: an operator implementing the sampler works in the
module, so a test that patched script text would report the guard as covered
while the code it guards sat one file away (CLAUDE.md pattern 14's observation,
one axis over).

**What a pipeline run cannot check, and this file therefore must.**

*Order preservation is invisible on the test data.* The degenerate arm's
byte-equality is what proves the subsampling machinery left the ordinary path
alone, and it rests on the draw preserving the input file's order. But at
``frac: 1.0`` a sampler that emitted its picks in *draw* order would still emit
the same set, and no count, no roster and no row total would move -- only the
file's bytes would, and CI byte-compares only that one arm. The fixture here
uses a fraction where draw order and file order differ observably.

*Banker's rounding is unobservable at the CI clade sizes.* ``round(frac * n)``
and ``floor(frac * n)`` agree for 19B at 0.2 (1 either way) and for 20A at 0.5
(23 either way), so a sampler using either passes CI. The fixture uses ``n=7``
at ``0.5``, where ``round`` gives 4 and ``floor`` gives 3.

*Trailing-zero arm collisions cannot arise from a shipped config.* YAML parses
``0.10`` as ``0.1``, so ``str(frac)`` would name two arms identically -- but no
config carries both spellings, so the collision is unreachable end-to-end.

*A zero-sample draw is a hard failure and a below-``min_clade_size`` draw is an
exclusion, and CI reaches only the second.* ``config_test_data.yaml``'s smallest
arm retains 1 sample, so nothing exercises the ``matUtils extract -s`` on an
empty list case -- the same latent shape #149 hit with its own floor.

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
SUMMARY_SCRIPT = REPO_ROOT / "workflow" / "scripts" / "subsample_summary.py"
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
# Seven samples, deliberately NOT in sorted order: ``select_subsample`` must
# return them in *this* order, so a sampler that sorted its output -- or emitted
# picks in draw order -- is caught. Seven at 0.5 is also where banker's rounding
# bites: round(3.5) == 4 while floor(3.5) == 3.
FULL_SAMPLES = [
    "zeta/2021",
    "alpha/2020",
    "mu/2022",
    "beta/2020",
    "omega/2023",
    "gamma/2021",
    "delta/2021",
]
UNIT = ["20A"]
CLADE = "20A"
SEED = 0


def _write_samples(path, samples, trailing_newline=True):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(samples) + ("\n" if trailing_newline else ""))


def frac_spec(frac, unit=None):
    return {
        "kind": subsample.KIND_FRAC,
        "unit": list(UNIT if unit is None else unit),
        "seed": SEED,
        "frac": frac,
        "cap": None,
    }


def cap_spec(cap, seed=SEED):
    return {
        "kind": subsample.KIND_CAP,
        "unit": None,
        "seed": seed,
        "frac": None,
        "cap": cap,
    }


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


def run_draw(tmp, frac, clade=CLADE, unit=None, scripts_dir=None,
             samples=None, trailing_newline=True):
    unit = UNIT if unit is None else unit
    return _run_draw_spec(
        tmp,
        frac_spec(frac, unit),
        subsample.arm_name(unit, SEED, frac),
        clade=clade,
        scripts_dir=scripts_dir,
        samples=samples,
        trailing_newline=trailing_newline,
    )


def run_cap_draw(tmp, cap, clade=CLADE, scripts_dir=None, samples=None,
                 trailing_newline=True, seed=SEED):
    return _run_draw_spec(
        tmp,
        cap_spec(cap, seed),
        subsample.cap_arm_name(cap, seed),
        clade=clade,
        scripts_dir=scripts_dir,
        samples=samples,
        trailing_newline=trailing_newline,
    )


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------


def check_degenerate_arm_reproduces_the_input_byte_for_byte():
    """``sampling_fractions: [1.0]`` must emit the input file verbatim.

    This is the only end-to-end evidence that the arm machinery did not perturb
    the ordinary path: every output below a degenerate arm inherits the identity
    from this file.
    """
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        out = run_draw(tmp, 1.0)
        full = Path(tmp) / "mat" / "clades" / f"{CLADE}_samples.txt"
        if out.read_bytes() != full.read_bytes():
            failures.append(
                "degenerate arm (frac=1.0) is not byte-identical to its input:\n"
                f"  full: {full.read_bytes()!r}\n  arm:  {out.read_bytes()!r}"
            )
    # And with no trailing newline, so the round-trip is not merely tolerant of
    # one particular ending.
    with tempfile.TemporaryDirectory() as tmp:
        out = run_draw(tmp, 1.0, trailing_newline=False)
        full = Path(tmp) / "mat" / "clades" / f"{CLADE}_samples.txt"
        if out.read_bytes() != full.read_bytes():
            failures.append(
                "degenerate arm is not byte-identical when the input has no "
                "trailing newline"
            )
    return failures


def check_draw_preserves_file_order():
    """A partial draw keeps the input's order, not the order it was drawn in."""
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        drawn = run_draw(tmp, 0.5).read_text().split("\n")
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
        if len(drawn) != 4:
            failures.append(
                f"draw kept {len(drawn)} samples, expected round(0.5 * 7) = 4"
            )
    return failures


def check_emitting_draw_order_is_caught():
    """Mutation test: drop the sort and the degenerate arm stops matching."""
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        scripts = _mutated_workflow(
            tmp,
            {MODULE_REL: (
                "idx = np.sort(rng.choice(n, size=k, replace=False))",
                "idx = rng.choice(n, size=k, replace=False)",
            )},
        )
        out = run_draw(tmp, 1.0, scripts_dir=scripts)
        full = Path(tmp) / "mat" / "clades" / f"{CLADE}_samples.txt"
        if out.read_bytes() == full.read_bytes():
            failures.append(
                "removing the sort from select_subsample left the degenerate "
                "arm byte-identical, so the ordering guarantee is untested"
            )
    return failures


def check_flooring_instead_of_rounding_is_caught():
    """Mutation test: ``int()`` for ``round()`` changes the retained count."""
    failures = []
    expected = subsample.retained_count(len(FULL_SAMPLES), 0.5)
    if expected != 4:
        failures.append(
            f"fixture no longer distinguishes round from floor: round(0.5 * "
            f"{len(FULL_SAMPLES)}) = {expected}"
        )
    with tempfile.TemporaryDirectory() as tmp:
        scripts = _mutated_workflow(
            tmp,
            {MODULE_REL: (
                "return int(round(float(frac) * int(n_samples)))",
                "return int(float(frac) * int(n_samples))",
            )},
        )
        drawn = [d for d in run_draw(tmp, 0.5, scripts_dir=scripts)
                 .read_text().split("\n") if d]
        if len(drawn) == expected:
            failures.append(
                "flooring instead of rounding produced the same count, so the "
                "size rule is untested"
            )
    return failures


def check_draw_is_a_pure_function_of_unit_seed_clade_fraction():
    """An arm's list cannot move when another fraction is added or removed.

    Otherwise a later config edit silently redefines earlier rungs of the
    ladder, and every number already computed against them becomes wrong with
    nothing to attribute it to.
    """
    failures = []

    def digest(fractions):
        arms = subsample.validate_entries(
            [{"unit": UNIT, "sampling_fractions": fractions, "seed": SEED}]
        )
        name = subsample.arm_name(UNIT, SEED, 0.5)
        assert name in arms
        drawn = subsample.select_subsample(FULL_SAMPLES, 0.5, UNIT, SEED, CLADE)
        return hashlib.md5("\n".join(drawn).encode()).hexdigest()

    baseline = digest([0.5])
    for fractions in ([0.5, 0.25], [0.25, 0.5], [1.0, 0.5, 0.25], [0.5]):
        if digest(fractions) != baseline:
            failures.append(
                f"the 0.5 arm's draw moved when the entry's fractions were "
                f"{fractions}"
            )
    # A different clade, unit or seed *must* give a different draw, or the
    # "pure function of the tuple" claim is vacuous. Checked on a list long
    # enough that a collision is not luck: at 7 samples there are only
    # C(7, 4) = 35 sorted subsets, so two unrelated seeds coincide routinely.
    wide = [f"s{i}" for i in range(200)]
    variants = {
        "clade": subsample.select_subsample(wide, 0.5, UNIT, SEED, "21J"),
        "unit": subsample.select_subsample(wide, 0.5, ["21J"], SEED, CLADE),
        "seed": subsample.select_subsample(wide, 0.5, UNIT, 1, CLADE),
        "fraction": subsample.select_subsample(wide, 0.25, UNIT, SEED, CLADE),
    }
    base_draw = subsample.select_subsample(wide, 0.5, UNIT, SEED, CLADE)
    for key, draw in variants.items():
        if draw == base_draw:
            failures.append(
                f"changing the {key} did not change the draw, so the seed does "
                f"not depend on it"
            )
    return failures


def check_seeding_from_iteration_order_is_caught():
    """Mutation test: seed from position and the draw stops being stable."""
    failures = []
    scripts_edits = {MODULE_REL: (
        'key = "|".join([ARM_PREFIX, "-".join(unit), str(seed), clade, frac_token(frac)])',
        'key = "|".join([ARM_PREFIX, str(seed)])',
    )}
    with tempfile.TemporaryDirectory() as tmp:
        module = _load_mutated(tmp, scripts_edits, "subsample_seedkey")
        wide = [f"s{i}" for i in range(200)]
        a = module.select_subsample(wide, 0.5, UNIT, SEED, CLADE)
        b = module.select_subsample(wide, 0.5, ["21J"], SEED, "21J")
        # Under the mutation the key no longer mentions the clade, unit or
        # fraction, so every arm draws the same indices. If the two draws still
        # differ, something *other* than ``arm_rng``'s key is keying the draw
        # and the purity check above is proving less than it claims.
        if a != b:
            failures.append(
                "collapsing the RNG key to the seed alone still gave different "
                "draws per clade, so something else is seeding the draw"
            )
    return failures


def check_trailing_zero_fractions_name_one_arm():
    """``0.10`` and ``0.1`` are the same fraction and must collide loudly."""
    failures = []
    if subsample.arm_name(UNIT, SEED, 0.10) != subsample.arm_name(UNIT, SEED, 0.1):
        failures.append("0.10 and 0.1 produced different arm names")
    try:
        subsample.validate_entries(
            [
                {"unit": UNIT, "sampling_fractions": [0.1], "seed": SEED},
                {"unit": UNIT, "sampling_fractions": [0.10], "seed": SEED},
            ]
        )
        failures.append(
            "two entries yielding the same arm name did not raise, so one would "
            "silently overwrite the other's outputs"
        )
    except ValueError:
        pass
    return failures


def check_str_frac_naming_is_caught():
    """Mutation test: ``str(frac)`` is not a canonical fraction encoding.

    ``frac_token`` exists so that any two float values standing for the same
    basis-point fraction name the same arm. Float arithmetic does not respect
    that -- ``0.1 * 0.5`` is ``0.05000000000000001`` -- so under ``str`` a
    config that computed a fraction rather than writing it out would silently
    get a *different* arm directory, and its outputs would not join to the ones
    already computed at that depth.
    """
    failures = []
    computed = 0.1 * 0.5
    if computed == 0.05:
        failures.append(
            "the fixture no longer distinguishes the encodings: 0.1 * 0.5 == "
            "0.05 exactly on this platform"
        )
        return failures
    if subsample.frac_token(computed) != subsample.frac_token(0.05):
        failures.append(
            f"frac_token({computed}) != frac_token(0.05), so equal fractions do "
            f"not name one arm"
        )
    with tempfile.TemporaryDirectory() as tmp:
        module = _load_mutated(
            tmp,
            {MODULE_REL: (
                'return f"{round(float(frac) * 1e4):d}bp"',
                "return str(frac)",
            )},
            "subsample_strfrac",
        )
        if module.frac_token(computed) == module.frac_token(0.05):
            failures.append(
                "str(frac) gave the same token for 0.1 * 0.5 and 0.05, so this "
                "check proves nothing about the encoding"
            )
        # And the pinned form must contain no path- or wildcard-hostile
        # characters, which is the other half of why it is pinned.
    token = subsample.frac_token(0.05)
    if not token.isalnum():
        failures.append(f"frac_token produced a non-alphanumeric token: {token!r}")
    return failures


def check_fraction_range_is_enforced():
    """A fraction outside ``0 < frac <= 1.0`` raises rather than being clamped."""
    failures = []
    for bad in (0.0, -0.1, 1.5, 2):
        try:
            subsample.validate_entries(
                [{"unit": UNIT, "sampling_fractions": [bad], "seed": SEED}]
            )
            failures.append(
                f"sampling fraction {bad} was accepted; a clamped fraction "
                f"produces a well-formed arm at the wrong depth"
            )
        except ValueError:
            pass
    for bad in ("0.5", None, True):
        try:
            subsample.validate_entries(
                [{"unit": UNIT, "sampling_fractions": [bad], "seed": SEED}]
            )
            failures.append(f"non-numeric sampling fraction {bad!r} was accepted")
        except ValueError:
            pass
    return failures


def check_unit_validation():
    """Unit names must be path-safe, dash-free, unique and known to the tree."""
    failures = []
    cases = {
        "dashed clade": {"unit": ["21J-1"], "sampling_fractions": [0.5]},
        "empty unit": {"unit": [], "sampling_fractions": [0.5]},
        "repeated clade": {"unit": ["21J", "21J"], "sampling_fractions": [0.5]},
        "slash in clade": {"unit": ["21J/x"], "sampling_fractions": [0.5]},
        "unknown key": {"unit": ["21J"], "sampling_fractions": [0.5], "oops": 1},
        "non-int seed": {"unit": ["21J"], "sampling_fractions": [0.5], "seed": 0.5},
    }
    for label, entry in cases.items():
        try:
            subsample.validate_entries([entry])
            failures.append(f"{label} was accepted: {entry!r}")
        except ValueError:
            pass
    try:
        subsample.validate_entries(
            [{"unit": ["21J"], "sampling_fractions": [0.5]}],
            known_clades={"19A", "20A"},
        )
        failures.append(
            "a unit naming a clade the tree does not process was accepted; that "
            "clade has no subtree to thin"
        )
    except ValueError:
        pass
    return failures


def check_zero_draw_raises_and_small_draw_is_merely_excluded():
    """The two thresholds are different in kind, and must stay so.

    ``matUtils extract -s`` on an empty sample list is undefined, so a clade
    thinned to zero is a **hard failure**. ``min_clade_size`` says how big a
    clade the user wants, so a clade thinned below it is merely **excluded** --
    and recorded, since an absent directory is otherwise indistinguishable from
    a wiring mistake.
    """
    failures = []
    arms = subsample.validate_entries(
        [{"unit": ["A"], "sampling_fractions": [0.01], "seed": 0}]
    )
    try:
        subsample.check_arm_sizes(arms, {"A": 10}, min_clade_size=0)
        failures.append(
            "round(0.01 * 10) = 0 did not raise; matUtils extract -s on an "
            "empty list is undefined"
        )
    except ValueError:
        pass

    arms = subsample.validate_entries(
        [{"unit": ["A"], "sampling_fractions": [0.2], "seed": 0}]
    )
    kept = subsample.check_arm_sizes(arms, {"A": 10}, min_clade_size=5)
    arm = next(iter(arms))
    if kept[arm]:
        failures.append(
            f"round(0.2 * 10) = 2 below min_clade_size 5 was kept: {kept[arm]}"
        )
    kept = subsample.check_arm_sizes(arms, {"A": 10}, min_clade_size=2)
    if kept[arm] != {"A": 2}:
        failures.append(
            f"round(0.2 * 10) = 2 at min_clade_size 2 should be kept, got "
            f"{kept[arm]}"
        )
    try:
        subsample.check_arm_sizes(arms, {"B": 10})
        failures.append("a unit clade absent from the clade counts did not raise")
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
        arms = module.validate_entries(
            [{"unit": ["A"], "sampling_fractions": [0.01], "seed": 0}]
        )
        try:
            module.check_arm_sizes(arms, {"A": 10}, min_clade_size=0)
        except ValueError:
            failures.append(
                "demoting the empty-draw check still raised, so this check "
                "proves nothing"
            )
    return failures


def check_a_non_unit_clade_is_refused_by_the_sampler():
    """Thinning a clade the config never named must raise, not silently draw.

    ``per_unit_clade_inputs`` routes non-unit clades to the full tree's files,
    so this can only be reached by a wiring mistake -- and a wiring mistake here
    produces an arm whose "unchanged" clades are quietly thinned, which nothing
    downstream can detect.
    """
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        try:
            run_draw(tmp, 0.5, clade="19A", unit=["20A"])
            failures.append(
                "the sampler thinned 19A for an arm whose unit is ['20A']"
            )
        except ValueError:
            pass
    return failures


def check_an_empty_full_list_raises():
    """A clade whose full sample list is empty raises rather than writing one."""
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        try:
            run_draw(tmp, 0.5, samples=[])
            failures.append("an empty full sample list did not raise")
        except ValueError:
            pass
    return failures


def check_retention_exceeds_the_sampling_fraction():
    """Realized retention runs above nominal, concavely, and is hand-checkable.

    A branch survives whenever *any* sample below it is retained, which is why
    a nominal fraction is not an evidence cut and why every downstream number
    has to be keyed to the measured retention instead.
    """
    failures = []
    # Three nodes, two samples and one mutation each, at p = 0.5:
    #   per node 1 - 0.5^2 = 0.75, so E = 2.25 of 3 -- 75% of the mutations
    #   from 50% of the samples.
    expected, sd = subsample.expected_mutation_retention([2, 2, 2], [1, 1, 1], 0.5)
    if abs(expected - 2.25) > 1e-12:
        failures.append(f"closed form gave {expected}, expected 2.25")
    if abs(sd - 0.75) > 1e-12:
        failures.append(f"closed-form sd gave {sd}, expected 0.75")

    # Zero-mutation branches contribute nothing, so passing them in is free.
    with_zeros, _ = subsample.expected_mutation_retention(
        [2, 2, 2, 9, 9], [1, 1, 1, 0, 0], 0.5
    )
    if abs(with_zeros - expected) > 1e-12:
        failures.append(
            f"zero-mutation branches moved the closed form: {with_zeros} vs "
            f"{expected}"
        )

    # Concavity: retained/fraction must fall as the fraction rises.
    node_samples = [2] * 500
    node_muts = [1] * 500
    ratios = []
    for p in (0.01, 0.05, 0.10, 0.20):
        e, _ = subsample.expected_mutation_retention(node_samples, node_muts, p)
        ratios.append((e / len(node_muts)) / p)
    if ratios != sorted(ratios, reverse=True):
        failures.append(f"retention/fraction is not decreasing: {ratios}")
    if ratios[-1] <= 1.0:
        failures.append(
            f"retention did not exceed the sampling fraction at p=0.20: "
            f"{ratios[-1]}"
        )
    return failures


# ---------------------------------------------------------------------------
# subsample_summary.py
# ---------------------------------------------------------------------------

QC_ROWS = pd.DataFrame(
    [
        ("outcome", "kept", 6),
        ("outcome", "max_nt_mutations", 2),
        ("outcome", "max_reversions_to_ref", 1),
        ("outcome", "max_reversions_to_clade_founder", 0),
        ("outcome", "no_mutations", 3),
        ("outcome", "duplicate_codon_position", 0),
        ("outcome", "identical_slice", 1),
        ("total", "branches_scored", 13),
    ],
    columns=["category", "key", "n_branches"],
)


def _write_summary_fixture(root, arm, n_full=7, n_arm=4, drop_arm_samples=False):
    """A minimal output tree for one arm over one clade and one gene."""
    root = Path(root)
    clades = root / "mat" / "clades"
    _write_samples(clades / f"{CLADE}_samples.txt", FULL_SAMPLES[:n_full])
    if not drop_arm_samples:
        _write_samples(clades / arm / f"{CLADE}_samples.txt", FULL_SAMPLES[:n_arm])
    for arm_dir in ("", arm):
        muts = root / "mutations" / arm_dir / f"{CLADE}_nt_mutations.csv"
        muts.parent.mkdir(parents=True, exist_ok=True)
        # Two branches, three mutations at full depth; one branch, one at the arm.
        rows = (
            [("n1", "A1T;C2G"), ("n2", "G3A")]
            if arm_dir == ""
            else [("n1", "A1T")]
        )
        pd.DataFrame(rows, columns=["node_id", "nt_mutations"]).to_csv(
            muts, index=False
        )
        # ``leaves_sharing_mutations`` is the ``n_node`` the closed form needs.
        # Two samples under each branch, so at p = 4/7 the closed form predicts
        # 3 * (1 - (3/7)^2) = 2.4490 of the 3 mutations -- deliberately well
        # above the sample fraction, which is the whole point of the statistic.
        tsv = root / "mutations" / arm_dir / f"{CLADE}_mutations.tsv"
        pd.DataFrame(
            [(node, 2) for node, _ in rows],
            columns=["node_id", "leaves_sharing_mutations"],
        ).to_csv(tsv, sep="\t", index=False)
        qc = root / "mutation_counts" / arm_dir / "branch_qc" / f"{CLADE}.csv"
        qc.parent.mkdir(parents=True, exist_ok=True)
        QC_ROWS.to_csv(qc, index=False)
        gqc = root / "pcps" / "S" / arm_dir / "branch_qc" / f"{CLADE}.csv"
        gqc.parent.mkdir(parents=True, exist_ok=True)
        QC_ROWS.to_csv(gqc, index=False)


def run_summary(root, roster, min_clade_size=2, genes=("S",)):
    out = Path(root) / "subsample_summary.csv"
    snake = _Namespace(
        input=_Namespace(samples=[], branch_qc=[]),
        output=_Namespace(summary=str(out)),
        params=_Namespace(
            roster=roster,
            genes=list(genes),
            min_clade_size=min_clade_size,
            outdir=str(root),
        ),
        log=[],
    )
    _exec_script(SUMMARY_SCRIPT, {"snakemake": snake})
    return pd.read_csv(out)


def _roster(arm, frac, included, excluded, n_retained, unit=None, cap=None):
    unit = UNIT if unit is None else unit
    return {
        arm: {
            "kind": subsample.KIND_CAP if cap is not None else subsample.KIND_FRAC,
            "unit": list(unit),
            "unit_label": "all" if cap is not None else "-".join(unit),
            "seed": SEED,
            "frac": frac,
            "cap": cap,
            "included": included,
            "excluded": excluded,
            "n_retained": n_retained,
        }
    }


def check_summary_reports_realized_retention():
    """The happy path: one counts row and one per-gene row, ratios measured."""
    failures = []
    arm = subsample.arm_name(UNIT, SEED, 0.5)
    with tempfile.TemporaryDirectory() as tmp:
        _write_summary_fixture(tmp, arm)
        df = run_summary(tmp, _roster(arm, 0.5, [CLADE], [], {CLADE: 4}))
    counts = df[df["gene"].isna() | (df["gene"] == "")]
    genes = df[df["gene"] == "S"]
    if len(counts) != 1 or len(genes) != 1:
        failures.append(
            f"expected one counts row and one gene row, got {len(counts)} and "
            f"{len(genes)}"
        )
        return failures
    row = counts.iloc[0]
    if int(row["n_samples_arm"]) != 4 or int(row["n_samples_full"]) != 7:
        failures.append(f"sample counts wrong: {row['n_samples_arm']}/{row['n_samples_full']}")
    # 1 of 3 nucleotide mutations retained from 4 of 7 samples.
    if abs(float(row["nt_mutation_retention"]) - 1 / 3) > 1e-6:
        failures.append(f"nt_mutation_retention wrong: {row['nt_mutation_retention']}")
    if abs(float(row["sample_fraction_realized"]) - 4 / 7) > 1e-6:
        failures.append(
            f"sample_fraction_realized wrong: {row['sample_fraction_realized']}"
        )
    if int(row["branches_dropped_max_nt_mutations"]) != 2:
        failures.append(
            f"branch QC attribution wrong: {row['branches_dropped_max_nt_mutations']}"
        )
    # The closed form must predict retention *above* the sample fraction: two
    # samples under each of the two full-depth branches at p = 4/7 gives
    # 3 * (1 - (3/7)^2) = 2.4490 of 3 mutations = 0.8163, against 0.5714 of the
    # samples. A prediction at or below the sample fraction means the survival
    # term was dropped.
    predicted = float(row["predicted_nt_mutations_arm"])
    if abs(predicted - 2.45) > 0.01:
        failures.append(
            f"closed-form prediction wrong: {predicted}, expected 2.45"
        )
    if float(row["predicted_nt_mutation_retention"]) <= float(
        row["sample_fraction_realized"]
    ):
        failures.append(
            "predicted retention did not exceed the sample fraction, so the "
            "any-sample-survives term is missing"
        )
    if pd.isna(row["nt_mutation_deviation_sd"]):
        failures.append("nt_mutation_deviation_sd was not reported")
    # The gene-level drop reasons exist only on the PCP side, since
    # count_mutations.py scores the whole genome.
    grow = genes.iloc[0]
    if int(grow["branches_dropped_identical_slice"]) != 1:
        failures.append(
            f"identical_slice missing from the gene row: "
            f"{grow['branches_dropped_identical_slice']}"
        )
    if not pd.isna(row.get("branches_dropped_identical_slice")):
        failures.append(
            "identical_slice was reported on the counts row, where the limits "
            "are applied genome-wide and there is no gene slice"
        )
    return failures


def check_summary_records_an_excluded_clade():
    """A clade thinned below ``min_clade_size`` is named with ``included=False``."""
    failures = []
    arm = subsample.arm_name(UNIT, SEED, 0.2)
    with tempfile.TemporaryDirectory() as tmp:
        _write_summary_fixture(tmp, arm, drop_arm_samples=True)
        df = run_summary(tmp, _roster(arm, 0.2, [], [CLADE], {}), min_clade_size=2)
    if len(df) != 1:
        failures.append(f"expected exactly one row for an excluded clade, got {len(df)}")
        return failures
    row = df.iloc[0]
    if bool(row["included"]):
        failures.append("an excluded clade was recorded as included")
    if int(row["n_samples_arm"]) != 1:
        failures.append(
            f"expected round(0.2 * 7) = 1 recorded for the excluded clade, got "
            f"{row['n_samples_arm']}"
        )
    return failures


def check_summary_reconciles_in_both_directions():
    """Requested-but-absent and present-but-declined both raise."""
    failures = []
    arm = subsample.arm_name(UNIT, SEED, 0.5)

    with tempfile.TemporaryDirectory() as tmp:
        _write_summary_fixture(tmp, arm, drop_arm_samples=True)
        try:
            run_summary(tmp, _roster(arm, 0.5, [CLADE], [], {CLADE: 4}))
            failures.append(
                "a clade the roster says is included, with no sample list on "
                "disk, did not raise"
            )
        except ValueError:
            pass

    with tempfile.TemporaryDirectory() as tmp:
        _write_summary_fixture(tmp, arm)
        try:
            run_summary(tmp, _roster(arm, 0.5, [], [CLADE], {}))
            failures.append(
                "a clade the roster declined, whose subtree nonetheless exists "
                "on disk, did not raise -- a stale intermediate would be read "
                "as a legitimate exclusion"
            )
        except ValueError:
            pass

    # A clade the arithmetic approves must not be silently declined: that is the
    # broken-fan-out case the roster exists to expose (#137's shape).
    with tempfile.TemporaryDirectory() as tmp:
        _write_summary_fixture(tmp, arm, drop_arm_samples=True)
        try:
            run_summary(
                tmp, _roster(arm, 0.5, [], [CLADE], {}), min_clade_size=2
            )
            failures.append(
                "round(0.5 * 7) = 4 >= min_clade_size 2 was recorded as an "
                "exclusion without raising"
            )
        except ValueError:
            pass
    return failures


def check_summary_catches_a_draw_of_the_wrong_size():
    """The recorded count is asserted against the *file*, not the checkpoint.

    The clade-counts checkpoint keys on the metadata column's first
    whitespace-delimited token while ``samples_by_clade`` matches
    ``^clade( |$)``, so the two can disagree -- and a disagreement would
    otherwise produce an arm at a different depth than the one recorded.
    """
    failures = []
    arm = subsample.arm_name(UNIT, SEED, 0.5)
    with tempfile.TemporaryDirectory() as tmp:
        _write_summary_fixture(tmp, arm, n_arm=3)  # should be round(0.5 * 7) = 4
        try:
            run_summary(tmp, _roster(arm, 0.5, [CLADE], [], {CLADE: 4}))
            failures.append(
                "an arm holding 3 samples where round(0.5 * 7) = 4 did not raise"
            )
        except ValueError:
            pass
    return failures



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
    arm = subsample.arm_name(UNIT, SEED, 0.5)
    with tempfile.TemporaryDirectory() as tmp:
        _corr_fixture(tmp, arm)
        got = run_correlations(tmp, _roster(arm, 0.5, [CLADE], [], {CLADE: 4}))
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
    arm = subsample.arm_name(UNIT, SEED, 0.5)
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
                roster=_roster(arm, 0.5, [CLADE], [], {CLADE: 4}),
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
    arm = subsample.arm_name(UNIT, SEED, 0.5)
    with tempfile.TemporaryDirectory() as tmp:
        # Still scored on both sides, but at a quarter of the counts.
        _corr_fixture(tmp, arm, arm_actual=3)
        got = run_correlations(tmp, _roster(arm, 0.5, [CLADE], [], {CLADE: 4}))
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
    arm = subsample.arm_name(UNIT, SEED, 0.5)
    with tempfile.TemporaryDirectory() as tmp:
        # arm_actual=1: still scored, but far below the threshold of 10.
        _corr_fixture(tmp, arm, arm_actual=1)
        got = run_correlations(tmp, _roster(arm, 0.5, [CLADE], [], {CLADE: 4}))
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
                roster=_roster(arm, 0.5, [CLADE], [], {CLADE: 4}),
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
    arm = subsample.arm_name(UNIT, SEED, 0.5)
    with tempfile.TemporaryDirectory() as tmp:
        _corr_fixture(tmp, arm)
        # Threshold far above any count in the fixture, so the filtered
        # population is empty for both estimators.
        got = run_correlations(
            tmp, _roster(arm, 0.5, [CLADE], [], {CLADE: 4}), count_threshold=1000
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
    arm = subsample.arm_name(UNIT, SEED, 0.5)
    with tempfile.TemporaryDirectory() as tmp:
        _corr_fixture(tmp, arm)
        # Same tree, but the roster says this clade is NOT in the arm's unit --
        # so its (deliberately disagreeing) effects must be reported as a leak.
        roster = {
            arm: {
                "kind": subsample.KIND_FRAC,
                "unit": ["21J"],
                "unit_label": "21J",
                "seed": SEED,
                "frac": 0.5,
                "cap": None,
                "included": ["21J"],
                "excluded": [],
                "n_retained": {"21J": 4},
            }
        }
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
# Cap arms (#151 amendment)
# ---------------------------------------------------------------------------
#
# No shipped production config sets a cap either, and the CI config's clades are
# far too small for the shipped values (50000/20000/5000) to touch anything, so
# every property below is unreachable from any pipeline run.


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
                    "        return [c for c, n in clade_sizes.items() if int(n) > cap]",
                    "        return [c for c, n in clade_sizes.items() if int(n) >= cap]",
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
    """A cap arm's list cannot move when another cap is added or removed.

    Stronger than the fraction arms' guarantee, and deliberately so: a cap arm's
    *unit* is derived from the tree's clade sizes, so keying the draw on the unit
    would make every clade's sample list depend on how many other clades
    happened to exceed the cap. A tree gaining one large clade would then
    silently redraw every arm, and nothing downstream could attribute it.
    """
    failures = []
    samples = [f"s{i}/2021" for i in range(200)]
    baseline = None
    for caps in ([50], [50, 20], [20, 50], [5, 20, 50], [50]):
        arms = subsample.validate_entries([{"max_samples": caps, "seed": SEED}])
        name = subsample.cap_arm_name(50, SEED)
        if name not in arms:
            failures.append(f"cap 50 missing from arms built with {caps}")
            continue
        drawn = subsample.select_for_spec(samples, arms[name], CLADE)
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
    spec = arms[subsample.cap_arm_name(50, SEED)]
    if subsample.select_for_spec(samples, spec, "20A") == subsample.select_for_spec(
        samples, spec, "21J"
    ):
        failures.append(
            "two different clades drew identically, so the clade is not part "
            "of the RNG key"
        )
    return failures


def check_cap_validation():
    """``max_samples`` must be a non-empty list of integers >= 1, and alone."""
    failures = []
    cases = {
        "cap with a unit": {"unit": ["20A"], "max_samples": [10]},
        "cap and fractions": {"max_samples": [10], "sampling_fractions": [0.5]},
        "neither depth key": {"unit": ["20A"]},
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


def check_cap_arms_size_and_exclude_like_fraction_arms():
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


CHECKS = [
    ("degenerate arm is byte-identical", check_degenerate_arm_reproduces_the_input_byte_for_byte),
    ("draw preserves file order", check_draw_preserves_file_order),
    ("emitting draw order is caught", check_emitting_draw_order_is_caught),
    ("flooring instead of rounding is caught", check_flooring_instead_of_rounding_is_caught),
    ("draw is a pure function of (unit, seed, clade, frac)", check_draw_is_a_pure_function_of_unit_seed_clade_fraction),
    ("trailing-zero fractions name one arm", check_trailing_zero_fractions_name_one_arm),
    ("fraction range is enforced", check_fraction_range_is_enforced),
    ("unit validation", check_unit_validation),
    ("zero draw raises, small draw is excluded", check_zero_draw_raises_and_small_draw_is_merely_excluded),
    ("demoting the zero-draw failure is caught", check_treating_a_zero_draw_as_an_exclusion_is_caught),
    ("a non-unit clade is refused by the sampler", check_a_non_unit_clade_is_refused_by_the_sampler),
    ("an empty full list raises", check_an_empty_full_list_raises),
    ("retention exceeds the sampling fraction", check_retention_exceeds_the_sampling_fraction),
    ("summary reports realized retention", check_summary_reports_realized_retention),
    ("summary records an excluded clade", check_summary_records_an_excluded_clade),
    ("summary reconciles in both directions", check_summary_reconciles_in_both_directions),
    ("summary catches a wrong-size draw", check_summary_catches_a_draw_of_the_wrong_size),
    ("cap thins only clades above it", check_cap_thins_only_clades_above_it),
    ("including a clade at the cap is caught", check_including_a_clade_at_the_cap_is_caught),
    ("the sampler refuses a clade at or below the cap", check_the_sampler_refuses_a_clade_at_or_below_the_cap),
    ("cap draw is a pure function of (cap, seed, clade)", check_cap_draw_is_a_pure_function_of_cap_seed_clade),
    ("cap validation", check_cap_validation),
    ("cap arms size and exclude like fraction arms", check_cap_arms_size_and_exclude_like_fraction_arms),
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
        "\nPASSED: the degenerate arm is byte-identical, the draw is a pure "
        "function of (unit, seed, clade, fraction) in file order, a zero draw "
        "is a hard failure while a below-min_clade_size draw is a recorded "
        "exclusion, the summary reconciles its roster against the file tree in "
        "both directions, and every guard fires on a deliberately broken module."
    )


if __name__ == "__main__":
    main()
