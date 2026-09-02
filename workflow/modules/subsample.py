#!/usr/bin/env python3
"""Clade subsampling arms: naming, validation and the draw itself (#151).

Per-clade fitness effects are estimated from one clade's own subtree, so a
thinly sampled clade gets a noisy estimate and cannot borrow strength from the
rest of the tree. Measuring how thin a clade has to get before a pooled model
wins needs the *same* clade re-estimated at several sampling depths, and this
module is the hook: an **arm** is one sampling depth applied to one **unit** of
clades, and everything an arm rebuilds follows from re-extracting each unit
clade's subtree from a subsampled sample list.

There are two kinds of arm, and the difference is what they intervene on.

*A ``frac`` arm* applies one sampling fraction to an explicitly named unit of
clades. It intervenes on **size**: every clade in the unit keeps the same share
of its samples, so the tree's sampling imbalance survives intact and each arm
gives one point on the depth curve.

*A ``cap`` arm* thins every clade with more than ``max_samples`` samples down to
exactly that many and leaves smaller clades untouched. It intervenes on
**composition**, which is the intervention the pooled-vs-per-clade question
actually asks about. Two consequences follow from the unit being *derived* rather
than named, and both are why the cap is not merely a different number:

- **One arm yields the whole depth curve.** At ``max_samples: 50000`` on the v12
  gisaid tree, 30 of 56 clades are capped and their effective fractions span
  1.15% (21J, 4,344,573 samples) to 100%, so a single arm measures
  correlation-vs-depth across two orders of magnitude. The 26 untouched clades
  are byte-identical controls **by construction** rather than by being
  deliberately left unthinned, which is how a ``frac`` arm gets its controls.
- **The unit cannot be named in the config, so it is resolved from the clade
  counts checkpoint.** ``resolve_unit`` is therefore the single definition of
  "which clades does this arm thin", called by the fan-out helpers, the summary
  and the sampler's own defensive check. A second derivation would let the
  fan-out and the sampler disagree about one clade, which is exactly the failure
  ``per_unit_clade_inputs`` exists to prevent.

The cap equalizes *across clades*, not across time: a uniform draw within a
clade preserves that clade's temporal profile and merely scales it down. A
genuinely uniform-in-time draw would cap per ``(clade, time bin)``; that is a
different draw rule and is deliberately not implemented here, so the
"more uniform sampling" framing is not overclaimed.

Four things about the shape are load-bearing, and each is why a function lives
here rather than inline in a script:

*The draw is a pure function of the arm's own identity and the clade.*
``arm_rng`` hashes ``(unit, seed, clade, frac)`` and ``cap_rng`` hashes
``(cap, seed, clade)``, never iteration order or a module-level default, so an
arm's sample list cannot move when another fraction or cap is added to or
removed from the config. Determinism is then structural rather than a property of
Snakemake happening not to re-run the sampler. Note what ``cap_rng``
deliberately omits: the **unit**, which for a cap arm is derived from the cap and
the clade sizes. Keying on it would make every clade's draw depend on how many
*other* clades happened to exceed the cap, so adding a clade to the tree would
silently redraw all of them.

*The degenerate arm must be byte-identical, so the draw preserves order.*
``select_subsample`` sorts the chosen indices back into the input file's order,
which makes ``frac: 1.0`` emit the input list unchanged -- and therefore makes
every downstream output of that arm byte-equal to the full tree's. A sampler
that emitted its picks in draw order would still be a correct uniform sample
and would silently retire the only end-to-end check that the arm machinery left
the ordinary path alone (#149's pattern, one axis over).

*A fraction token is not ``str(frac)``.* YAML parses ``0.10`` as ``0.1``, so two
configs differing only in a trailing zero would name the same arm; and ``0.05``
formats as ``0.05`` while ``0.010000000000000002`` would not. ``frac_token``
pins the encoding to basis points.
"""

import hashlib
import re

import numpy as np

# The prefix a fraction arm's directory carries. It exists so ``arm_dir``'s
# optional alternative (see the Snakefile's ``wildcard_constraints``) cannot be
# confused with a clade name: no Nextstrain clade starts with ``frac_``.
ARM_PREFIX = "frac"

# The prefix a cap arm's directory carries, on the same reasoning. Kept as a
# separate literal rather than derived from ``KIND_CAP`` because it also appears
# inside ``cap_rng``'s hashed key: renaming the kind must not silently redraw
# every cap arm's sample list.
CAP_PREFIX = "cap"

ARM_PREFIXES = (ARM_PREFIX, CAP_PREFIX)

KIND_FRAC = "frac"
KIND_CAP = "cap"

_PREFIX_ALT = "|".join(ARM_PREFIXES)

# Matches an arm directory *including* its trailing slash, or the empty string
# for the full-tree path. Kept here rather than in the Snakefile so the tests
# and the scripts agree with the DAG on what an arm path looks like.
ARM_DIR_PATTERN = rf"|(?:{_PREFIX_ALT})_[A-Za-z0-9._+-]+/"

# Matches an arm name (no trailing slash).
ARM_PATTERN = rf"(?:{_PREFIX_ALT})_[A-Za-z0-9._+-]+"


def frac_token(frac):
    """Encode a sampling fraction as a basis-point token, e.g. ``0.05`` -> ``500bp``.

    Pinned rather than ``str(frac)``: YAML parses ``0.10`` as ``0.1``, so two
    configs differing only in a trailing zero would collide on one arm name and
    silently share its outputs.
    """
    return f"{round(float(frac) * 1e4):d}bp"


def arm_name(unit, seed, frac):
    """``frac_{unit}_{seed}_{frac}``, with ``unit`` joined by ``-``.

    ``unit`` is a list of clades because clade labels are finer than
    backgrounds: 21I's spike founder differs from 21J's by one amino acid while
    carrying ~10% of its samples, so thinning 21J alone leaves a ~10% evidence
    floor. Clade names containing ``-`` are rejected by ``validate_entries`` so
    the join stays invertible.
    """
    return f"{ARM_PREFIX}_{'-'.join(unit)}_{seed}_{frac_token(frac)}"


def cap_token(cap):
    """Encode a per-clade cap as its plain integer, e.g. ``50000`` -> ``"50000"``.

    Unlike ``frac_token`` this needs no pinned encoding: a cap is a sample count,
    so YAML parses it as an ``int`` and there is no trailing-zero spelling that
    means the same number. ``int()`` rather than ``str()`` so ``50000.0`` and
    ``50000`` cannot name two arms that would share every output.
    """
    return f"{int(cap):d}"


def cap_arm_name(cap, seed):
    """``cap_{max_samples}_{seed}``.

    No unit in the name, because a cap arm's unit is *derived* from the cap and
    the tree's clade sizes rather than named -- see ``resolve_unit``. Putting the
    derived unit in the name would make the arm directory change when the tree
    gains a clade, orphaning every output already built under the old name.
    """
    return f"{CAP_PREFIX}_{cap_token(cap)}_{seed}"


def retained_count(n_samples, frac):
    """``round(frac * n_samples)`` -- the arm's exact sample count.

    Python's banker's rounding is used deliberately and unwrapped: the number is
    reported in ``subsample_summary.csv`` and asserted by the membership check,
    so both sides must come from this one function rather than from two
    re-derivations that agree except at the halfway case.
    """
    return int(round(float(frac) * int(n_samples)))


def retained_cap_count(n_samples, cap):
    """``min(n_samples, cap)`` -- a cap arm's exact sample count for one clade.

    The ``min`` is the whole difference between a cap and a fraction: a clade
    already at or below the cap is **not thinned at all**, which is what makes
    the small clades of a cap arm byte-identical controls rather than
    near-identical ones.
    """
    return min(int(n_samples), int(cap))


def retained_for_spec(n_samples, spec):
    """The arm's retained sample count for a clade of ``n_samples``, either kind.

    One function so the DAG's sizing arithmetic, the sampler and the summary's
    reconciliation cannot drift: the summary asserts the *file's* line count
    against this number, and a second re-derivation agreeing except at a
    boundary would report an arm at a depth it was not built at.
    """
    if spec_kind(spec) == KIND_CAP:
        return retained_cap_count(n_samples, spec["cap"])
    return retained_count(n_samples, spec["frac"])


def spec_kind(spec):
    """``"frac"`` or ``"cap"`` for one validated arm spec.

    Read from the explicit ``kind`` key rather than inferred from which of
    ``frac``/``cap`` is non-``None``: an inference would silently pick a kind for
    a malformed spec carrying both, where this raises.
    """
    kind = spec.get("kind")
    if kind not in (KIND_FRAC, KIND_CAP):
        raise ValueError(
            f"arm spec has kind {kind!r}; expected {KIND_FRAC!r} or {KIND_CAP!r}. "
            f"Spec: {spec!r}"
        )
    return kind


def arm_rng(unit, seed, clade, frac):
    """A generator seeded from ``(unit, seed, clade, frac)`` alone.

    ``hashlib`` rather than ``hash()``: Python randomizes string hashing per
    process, so ``hash()`` would make an arm's sample list depend on which
    interpreter drew it -- reproducible within a run and not across runs, which
    is the worst of both.
    """
    key = "|".join([ARM_PREFIX, "-".join(unit), str(seed), clade, frac_token(frac)])
    digest = hashlib.sha256(key.encode()).digest()
    return np.random.default_rng(int.from_bytes(digest[:8], "big"))


def cap_rng(cap, seed, clade):
    """A generator seeded from ``(cap, seed, clade)`` alone.

    Deliberately **not** keyed on the unit, unlike ``arm_rng``. A cap arm's unit
    is derived from the cap and the tree's clade sizes, so keying on it would make
    each clade's draw depend on how many other clades happened to exceed the cap:
    a tree gaining one large clade would silently redraw every arm. Keying on the
    cap alone makes a clade's sample list a property of that clade and that
    number, which is also what lets two cap arms share nothing but still be
    reproducible independently.
    """
    key = "|".join([CAP_PREFIX, cap_token(cap), str(seed), clade])
    digest = hashlib.sha256(key.encode()).digest()
    return np.random.default_rng(int.from_bytes(digest[:8], "big"))


def select_subsample(samples, frac, unit, seed, clade):
    """``round(frac * len(samples))`` of ``samples``, in the input's own order.

    Order preservation is what makes the degenerate arm exact: at ``frac ==
    1.0`` every index is chosen and sorting them recovers the input list
    verbatim, so the arm's subtree -- and every output below it -- is
    byte-identical to the full tree's.
    """
    n = len(samples)
    k = retained_count(n, frac)
    if k < 1:
        raise ValueError(
            f"Subsample of clade {clade!r} at fraction {frac} would retain "
            f"round({frac} * {n}) = {k} samples; `matUtils extract -s` on an "
            f"empty sample list is undefined. Raise the fraction or drop the "
            f"clade from the unit."
        )
    rng = arm_rng(unit, seed, clade, frac)
    idx = np.sort(rng.choice(n, size=k, replace=False))
    return [samples[i] for i in idx]


def select_cap(samples, cap, seed, clade):
    """At most ``cap`` of ``samples``, in the input's own order.

    A clade already at or below the cap is returned **unchanged and undrawn** --
    not drawn with every index selected. The two are equal here because
    ``select_subsample`` sorts its indices, but returning early says the invariant
    out loud and keeps it true of any future sampler: an untouched clade of a cap
    arm must be byte-identical to the full tree's, and that is the property the
    arm's control clades rest on.
    """
    n = len(samples)
    cap = int(cap)
    if cap < 1:
        raise ValueError(
            f"Cap arm would retain {cap} samples for clade {clade!r}; "
            f"`matUtils extract -s` on an empty sample list is undefined. "
            f"'max_samples' must be >= 1."
        )
    if n <= cap:
        return list(samples)
    rng = cap_rng(cap, seed, clade)
    idx = np.sort(rng.choice(n, size=cap, replace=False))
    return [samples[i] for i in idx]


def select_for_spec(samples, spec, clade):
    """Draw one clade's sample list for either kind of arm."""
    if spec_kind(spec) == KIND_CAP:
        return select_cap(samples, spec["cap"], spec["seed"], clade)
    return select_subsample(
        samples, spec["frac"], spec["unit"], spec["seed"], clade
    )


def resolve_unit(spec, clade_sizes):
    """The clades this arm thins: named for a ``frac`` arm, derived for a ``cap``.

    **The single definition of arm membership**, called by the Snakefile's
    fan-out helpers, the summary and the sampler's own defensive check. A cap
    arm's unit cannot be written in the config -- it depends on the tree's clade
    sizes, which only the ``get_clades`` checkpoint knows -- so a second
    derivation anywhere would let the fan-out and the sampler disagree about one
    clade, producing an arm whose "unchanged" clades are quietly thinned or whose
    thinned clades are quietly read from the full tree.

    ``clade_sizes`` order is preserved, so the unit inherits the full roster's
    order and the aggregate rules concatenate their inputs in the same sequence
    the full tree does -- which is what the degenerate arm's byte-equality rests
    on.

    A ``cap`` arm's unit is ``n_clade > cap``, strictly: a clade of exactly
    ``cap`` samples is untouched, so thinning it would write an identical file
    into the arm directory and make the "arm directory holds only thinned clades"
    assertion false for a clade that is not thinned.
    """
    if spec_kind(spec) == KIND_CAP:
        cap = int(spec["cap"])
        return [c for c, n in clade_sizes.items() if int(n) > cap]
    missing = [c for c in spec["unit"] if c not in clade_sizes]
    if missing:
        raise ValueError(
            f"arm unit names {missing}, absent from the clade counts "
            f"checkpoint. Clades in 'clades_to_exclude' or beyond 'max_clades' "
            f"have no subtree to thin. Available: {sorted(clade_sizes)}"
        )
    return list(spec["unit"])


def arm_units(arms, clade_sizes):
    """``{arm: [clades it thins]}`` for every arm, resolved against one tree."""
    return {arm: resolve_unit(spec, clade_sizes) for arm, spec in arms.items()}


def unit_label(spec):
    """The ``unit`` column's value: the joined clade list, or ``all`` for a cap.

    A cap arm's unit is every clade above the cap -- 30 of 56 on the v12 gisaid
    tree -- so joining it would put a 200-character cell on every row and make
    the column change when the tree does. The per-clade rows already name their
    clade, and ``max_samples`` already says which clades qualify.
    """
    if spec_kind(spec) == KIND_CAP:
        return "all"
    return "-".join(spec["unit"])


def depth_label(spec):
    """A short human-readable depth, for logs: ``frac=0.05`` or ``cap=50000``."""
    if spec_kind(spec) == KIND_CAP:
        return f"cap={cap_token(spec['cap'])}"
    return f"frac={spec['frac']}"


def expected_mutation_retention(node_samples, node_muts, frac):
    """``(expected, sd)`` nucleotide mutations retained at sampling fraction ``frac``.

    A branch survives whenever **any** sample below it is retained, so retention
    of mutations runs well *above* the fraction of samples drawn:

        E[retained](p) = sum over nodes of  k_node * (1 - (1 - p) ** n_node)

    with ``n_node`` the clade samples below the node and ``k_node`` the
    nucleotide mutations on its own branch. Both come from
    ``mutations/{clade}_mutations.tsv``, whose ``leaves_sharing_mutations``
    column *is* ``n_node`` -- one row per node, so the sum is cheap.
    ``{clade}_sample_paths.tsv`` carries the same information one row per *tip*
    (multi-GB on 21J) and is deliberately not read. On the v11 gisaid tree this makes 1%
    of 21J's samples a 24x evidence cut rather than a 100x one (4.10% of its
    mutations retained), and the curve is concave: retained/fraction runs 1.79x
    at p=0.20, 2.22x at 0.10, 2.72x at 0.05 and 4.10x at 0.01.

    ``sd`` is returned as a **scale for judging a deviation, not a p-value**.
    Sibling survival is negatively correlated under fixed-size sampling, so the
    Bernoulli variance below is an approximation; and the draw being fixed-size
    rather than Bernoulli makes the mean itself a slight over-prediction. What
    the pair supports is "a same-signed deviation of many sd at two independent
    seeds indicates a bug", which is the question actually being asked.

    Zero-mutation branches contribute exactly zero to both sums, which is why
    they may be passed in freely -- unlike the aggregate-rate statistic of
    CLAUDE.md pattern 13, where dropping them is the whole bug.
    """
    q = 1.0 - float(frac)
    expected = 0.0
    variance = 0.0
    for n, k in zip(node_samples, node_muts):
        if not k:
            continue
        survive = 1.0 - q ** int(n)
        expected += k * survive
        variance += (k ** 2) * (q ** int(n)) * survive
    return expected, variance ** 0.5


def validate_entries(entries, known_clades=None):
    """Validate ``clade_subsamples`` and return the expanded arm table.

    Returns ``{arm_name: spec}`` in config order, where ``spec`` carries
    ``kind``, ``unit``, ``seed``, ``frac`` and ``cap`` -- the same five keys for
    both kinds, with the inapplicable depth key ``None``. A uniform shape rather
    than two: every consumer then reads ``spec_kind`` and dispatches, instead of
    testing for a key's presence and silently treating a malformed spec as the
    other kind.

    An entry sets **exactly one** of ``sampling_fractions`` (a ``frac`` arm) or
    ``max_samples`` (a ``cap`` arm). Neither is a hard failure and so is both:
    an entry setting both would have two defensible readings, and picking one
    would build an arm at a depth the config does not unambiguously ask for.

    Every failure raises at DAG-construction time with the offending numbers in
    the message, rather than being clamped: a silently clamped fraction or cap
    produces a well-formed arm at the wrong depth, which is indistinguishable
    from the depth that was asked for.
    """
    arms = {}
    if not entries:
        return arms
    if not isinstance(entries, list):
        raise ValueError(
            f"config 'clade_subsamples' must be a list of entries; got "
            f"{type(entries).__name__}"
        )
    for i, entry in enumerate(entries):
        if not isinstance(entry, dict) or not set(entry) <= {
            "unit",
            "sampling_fractions",
            "max_samples",
            "seed",
        }:
            raise ValueError(
                f"config 'clade_subsamples[{i}]' must be a mapping with keys "
                f"'unit' + 'sampling_fractions' (a frac arm) or 'max_samples' "
                f"(a cap arm), and optionally 'seed'; got {entry!r}"
            )
        has_fractions = "sampling_fractions" in entry
        has_cap = "max_samples" in entry
        if has_fractions == has_cap:
            raise ValueError(
                f"config 'clade_subsamples[{i}]' must set exactly one of "
                f"'sampling_fractions' (thin a named unit of clades by a "
                f"fraction) or 'max_samples' (thin every clade above the cap "
                f"down to it); got keys {sorted(entry)}."
            )
        seed = entry.get("seed", 0)
        if not isinstance(seed, int) or isinstance(seed, bool):
            raise ValueError(
                f"config 'clade_subsamples[{i}].seed' must be an integer; got "
                f"{seed!r}"
            )
        if has_cap:
            _add_cap_arms(arms, i, entry, seed)
            continue
        unit = entry.get("unit")
        if not isinstance(unit, list) or not unit:
            raise ValueError(
                f"config 'clade_subsamples[{i}].unit' must be a non-empty list "
                f"of clade names; got {unit!r}"
            )
        bad_dash = [c for c in unit if "-" in str(c)]
        if bad_dash:
            raise ValueError(
                f"config 'clade_subsamples[{i}].unit' names {bad_dash}, which "
                f"contain '-'. Arm names join the unit with '-', so such a name "
                f"would make the arm directory ambiguous."
            )
        bad_chars = [c for c in unit if not re.fullmatch(r"[A-Za-z0-9._+]+", str(c))]
        if bad_chars:
            raise ValueError(
                f"config 'clade_subsamples[{i}].unit' names {bad_chars}, which "
                f"are not valid path-safe clade names ([A-Za-z0-9._+])."
            )
        if len(set(unit)) != len(unit):
            raise ValueError(
                f"config 'clade_subsamples[{i}].unit' repeats a clade: {unit!r}"
            )
        if known_clades is not None:
            unknown = [c for c in unit if c not in known_clades]
            if unknown:
                raise ValueError(
                    f"config 'clade_subsamples[{i}].unit' names {unknown}, which "
                    f"this tree does not process. Clades in 'clades_to_exclude' "
                    f"or beyond 'max_clades' have no subtree to thin. Available: "
                    f"{sorted(known_clades)}"
                )
        fractions = entry.get("sampling_fractions")
        if not isinstance(fractions, list) or not fractions:
            raise ValueError(
                f"config 'clade_subsamples[{i}].sampling_fractions' must be a "
                f"non-empty list; got {fractions!r}"
            )
        for frac in fractions:
            if not isinstance(frac, (int, float)) or isinstance(frac, bool):
                raise ValueError(
                    f"config 'clade_subsamples[{i}].sampling_fractions' contains "
                    f"{frac!r}, which is not a number."
                )
            frac = float(frac)
            if not (0.0 < frac <= 1.0):
                raise ValueError(
                    f"config 'clade_subsamples[{i}].sampling_fractions' contains "
                    f"{frac}, outside the allowed range 0 < frac <= 1.0."
                )
            name = arm_name(unit, seed, frac)
            if name in arms:
                raise ValueError(
                    f"config 'clade_subsamples' yields the arm name {name!r} "
                    f"twice. Two entries share a unit, a seed and a sampling "
                    f"fraction (note {frac} and its trailing-zero spellings are "
                    f"the same fraction)."
                )
            arms[name] = {
                "kind": KIND_FRAC,
                "unit": list(unit),
                "seed": seed,
                "frac": frac,
                "cap": None,
            }
    return arms


def _add_cap_arms(arms, i, entry, seed):
    """Expand one ``max_samples`` entry into its cap arms, validating as it goes.

    ``unit`` is rejected rather than intersected: a cap arm's unit is *derived*
    from the cap and the tree's clade sizes (``resolve_unit``), so an entry
    carrying both would ask for two membership rules at once. Allowing it and
    intersecting them would make "which clades did this arm thin" answerable only
    by re-running the intersection, which is the ambiguity the derived unit
    exists to remove.
    """
    if "unit" in entry:
        raise ValueError(
            f"config 'clade_subsamples[{i}]' sets both 'max_samples' and "
            f"'unit'. A cap arm thins every clade above the cap, so its unit is "
            f"derived from the tree's clade sizes rather than named; drop 'unit'."
        )
    caps = entry["max_samples"]
    if not isinstance(caps, list) or not caps:
        raise ValueError(
            f"config 'clade_subsamples[{i}].max_samples' must be a non-empty "
            f"list of sample counts; got {caps!r}"
        )
    for cap in caps:
        if not isinstance(cap, int) or isinstance(cap, bool):
            raise ValueError(
                f"config 'clade_subsamples[{i}].max_samples' contains {cap!r}, "
                f"which is not an integer sample count. A cap is a number of "
                f"sequences, not a fraction."
            )
        if cap < 1:
            raise ValueError(
                f"config 'clade_subsamples[{i}].max_samples' contains {cap}, "
                f"which would retain no samples. `matUtils extract -s` on an "
                f"empty sample list is undefined, so this is a hard failure "
                f"rather than a clamp."
            )
        name = cap_arm_name(cap, seed)
        if name in arms:
            raise ValueError(
                f"config 'clade_subsamples' yields the arm name {name!r} twice. "
                f"Two entries share a cap and a seed."
            )
        arms[name] = {
            "kind": KIND_CAP,
            "unit": None,
            "seed": seed,
            "frac": None,
            "cap": int(cap),
        }


def check_arm_sizes(arms, clade_sizes, min_clade_size=0):
    """Raise on any (arm, clade) whose draw would be empty; report the rest.

    Returns ``{arm: {clade: n_retained}}`` for the clades an arm keeps. A clade
    thinned to zero samples is a **hard failure**, because ``matUtils extract
    -s`` on an empty list is undefined; a clade thinned below ``min_clade_size``
    is merely *excluded* from the arm, since that threshold expresses how big a
    clade the user wants rather than what can be built at all.

    ``get_clades`` applies ``min_clade_size`` to the full metadata and the
    checkpoint's consumers read only its output, so neither can see a clade that
    falls below the threshold only after thinning -- which is why the arithmetic
    happens here, over the per-clade counts that checkpoint already wrote.

    Both branches are unreachable for a ``cap`` arm by construction, and
    deliberately kept rather than special-cased: a cap arm's unit is
    ``n_clade > cap``, so its retained count is exactly ``cap``, which
    ``validate_entries`` has already forced to be ``>= 1``. It can still fall
    below a ``min_clade_size`` larger than the cap, which is a real
    configuration and is excluded-and-recorded like any other.
    """
    kept = {}
    for arm, spec in arms.items():
        empty = []
        kept[arm] = {}
        # ``resolve_unit`` is the single definition of membership and raises on a
        # frac arm naming a clade this tree does not process; a cap arm's unit is
        # derived from these very sizes and so cannot name an absent clade.
        for clade in resolve_unit(spec, clade_sizes):
            n = int(clade_sizes[clade])
            k = retained_for_spec(n, spec)
            if k < 1:
                empty.append((clade, n, k))
            elif k >= min_clade_size:
                kept[arm][clade] = k
        if empty:
            detail = ", ".join(f"{c}: {n} -> {k}" for c, n, k in empty)
            raise ValueError(
                f"arm {arm!r} ({depth_label(spec)}) would draw zero samples for "
                f"{len(empty)} clade(s) ({detail}). `matUtils extract -s` on an "
                f"empty sample list is undefined, so this is a hard failure "
                f"rather than a skip. Raise the sampling depth or drop the clade "
                f"from the unit."
            )
    return kept
