#!/usr/bin/env python3
"""Clade subsampling arms: naming, validation and the draw itself (#151).

Per-clade fitness effects are estimated from one clade's own subtree, so a
thinly sampled clade gets a noisy estimate and cannot borrow strength from the
rest of the tree. Measuring how thin a clade has to get before a pooled model
wins needs the *same* clade re-estimated at several sampling depths, and this
module is the hook: an **arm** is one per-clade cap, and everything an arm
rebuilds follows from re-extracting each thinned clade's subtree from a
subsampled sample list.

*A ``cap`` arm* thins every clade with more than ``max_samples`` samples down to
exactly that many and leaves smaller clades untouched. It intervenes on the
tree's **composition**, which is the intervention the pooled-vs-per-clade
question actually asks about. Two consequences follow from an arm's membership
being *derived* rather than named, and both are why the cap is not merely a
different number:

- **One arm yields the whole depth curve.** At ``max_samples: 50000`` on the v12
  gisaid tree, 30 of 56 clades are capped and their effective fractions span
  1.15% (21J, 4,344,573 samples) to 100%, so a single arm measures
  correlation-vs-depth across two orders of magnitude. The 26 untouched clades
  are byte-identical controls **by construction**, not by being deliberately
  left unthinned.
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

**A fraction arm is deliberately not offered.** Through an earlier revision the
module also carried a ``frac`` kind, which applied one sampling fraction to an
explicitly named unit of clades -- intervening on **size**, so every clade in
the unit kept the same share and the tree's sampling imbalance survived intact.
It was removed rather than kept as a second option, for three reasons that are
properties of the question rather than of the implementation. It answers a
different question (how noisy is *this clade* at depth p) than the one the arms
exist for (does a pooled model beat per-clade counts when the tree's
composition is levelled). Its controls have to be *chosen* -- a clade left
deliberately unthinned -- where a cap's fall out of the definition. And naming a
unit honestly is harder than it looks: 21I's spike founder differs from 21J's by
one amino acid while carrying ~10% of its samples, so a unit of ``[21J]`` alone
leaves a ~10% evidence floor, and on the v11 gisaid tree 21L alone leaves
2,737,183 samples within 20 nt of 21J's founder -- 139% of that unit -- so an
honest ladder for one clade needs a five-clade unit. A cap leaves no large clade
at full depth at all, which dissolves the problem instead of documenting it.
``validate_entries`` therefore rejects ``sampling_fractions`` and ``unit`` by
name, with that message, rather than as unknown keys.

Three things about the shape are load-bearing, and each is why a function lives
here rather than inline in a script:

*The draw is a pure function of the arm's own identity and the clade.*
``cap_rng`` hashes ``(cap, seed, clade)``, never iteration order or a
module-level default, so an arm's sample list cannot move when another cap is
added to or removed from the config. Determinism is then structural rather than a
property of Snakemake happening not to re-run the sampler. Note what the key
deliberately omits: the **unit**, which is derived from the cap and the clade
sizes. Keying on it would make every clade's draw depend on how many *other*
clades happened to exceed the cap, so adding a clade to the tree would silently
redraw all of them.

*An untouched clade's byte-identity comes from not being drawn at all.* Note
where it does **not** come from: ``select_cap`` returns a clade at or below the
cap unchanged and undrawn, but nothing in the pipeline reaches that branch --
``resolve_unit`` admits a clade only when ``n > cap``, so the fan-out never asks
for an untouched clade's arm sample list and ``subsample_clade_samples`` raises
if one arrives anyway. The control's guarantee is that the arm *reads the full
tree's file*, which is stronger than reproducing it. The early return is kept
because it states the invariant for any other caller, and because a future
sampler that drew every index would then be wrong in one place rather than two.

That in turn is what makes an arm whose cap sits above **every** clade's size
the degenerate case: its unit is empty, every clade is a control, and every
output byte-equals the full tree's -- which is the end-to-end check that the arm
machinery left the ordinary path alone (#149's pattern, one axis over).

*The draw still preserves the input file's order*, so an arm's sample list is a
**subsequence** of the full tree's rather than a permutation of a subset. That is
not what the controls rest on, but it keeps the two lists comparable line-for-line
and removes sample order as a variable when an arm's subtree is compared against
full depth.

*A cap token is ``int``, not ``str``.* ``50000.0`` and ``50000`` are the same
depth, so they must not name two arms that would share every output;
``validate_entries`` rejects the float outright and ``cap_token`` pins the
encoding for everything downstream of it.
"""

import hashlib

import numpy as np

# The prefix an arm's directory carries. It exists so ``arm_dir``'s optional
# alternative (see the Snakefile's ``wildcard_constraints``) cannot be confused
# with a clade name: no Nextstrain clade starts with ``cap_``.
#
# **Its value is hashed into ``cap_rng``'s key**, so renaming the constant is
# free but changing the string silently redraws every arm's sample list. The
# four production cap arms built before the frac kind was removed reproduce
# byte-for-byte because this stayed ``"cap"``.
ARM_PREFIX = "cap"

# Matches an arm directory *including* its trailing slash, or the empty string
# for the full-tree path. Kept here rather than in the Snakefile so the tests
# and the scripts agree with the DAG on what an arm path looks like.
ARM_DIR_PATTERN = rf"|{ARM_PREFIX}_[A-Za-z0-9._+-]+/"

# Matches an arm name (no trailing slash).
ARM_PATTERN = rf"{ARM_PREFIX}_[A-Za-z0-9._+-]+"


def cap_token(cap):
    """Encode a per-clade cap as its plain integer, e.g. ``50000`` -> ``"50000"``.

    ``int()`` rather than ``str()`` so ``50000.0`` and ``50000`` cannot name two
    arms that would share every output. A cap needs no further pinning: it is a
    sample count, so YAML parses it as an ``int`` and there is no trailing-zero
    spelling that means the same number.
    """
    return f"{int(cap):d}"


def arm_name(cap, seed):
    """``cap_{max_samples}_{seed}``.

    No unit in the name, because an arm's unit is *derived* from the cap and the
    tree's clade sizes rather than named -- see ``resolve_unit``. Putting the
    derived unit in the name would make the arm directory change when the tree
    gains a clade, orphaning every output already built under the old name.
    """
    return f"{ARM_PREFIX}_{cap_token(cap)}_{seed}"


def retained_count(n_samples, cap):
    """``min(n_samples, cap)`` -- an arm's exact sample count for one clade.

    The ``min`` is the whole point of a cap: a clade already at or below it is
    **not thinned at all**, which is what makes an arm's small clades
    byte-identical controls rather than near-identical ones.

    One function so the DAG's sizing arithmetic, the sampler and the summary's
    reconciliation cannot drift: the summary asserts the *file's* line count
    against this number, and a second re-derivation agreeing except at a
    boundary would report an arm at a depth it was not built at.
    """
    return min(int(n_samples), int(cap))


def cap_rng(cap, seed, clade):
    """A generator seeded from ``(cap, seed, clade)`` alone.

    ``hashlib`` rather than ``hash()``: Python randomizes string hashing per
    process, so ``hash()`` would make an arm's sample list depend on which
    interpreter drew it -- reproducible within a run and not across runs, which
    is the worst of both.

    Deliberately **not** keyed on the unit. An arm's unit is derived from the cap
    and the tree's clade sizes, so keying on it would make each clade's draw
    depend on how many other clades happened to exceed the cap: a tree gaining
    one large clade would silently redraw every arm. Keying on the cap alone
    makes a clade's sample list a property of that clade and that number, which
    is also what lets two arms share nothing but still be reproducible
    independently.
    """
    key = "|".join([ARM_PREFIX, cap_token(cap), str(seed), clade])
    digest = hashlib.sha256(key.encode()).digest()
    return np.random.default_rng(int.from_bytes(digest[:8], "big"))


def select_cap(samples, cap, seed, clade):
    """At most ``cap`` of ``samples``, in the input's own order.

    A clade already at or below the cap is returned **unchanged and undrawn** --
    not drawn with every index selected. The two are equal here because the draw
    sorts its indices, but returning early says the invariant out loud and keeps
    it true of any future sampler: an untouched clade must be byte-identical to
    the full tree's, and that is the property an arm's control clades rest on.
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


def resolve_unit(spec, clade_sizes):
    """The clades this arm thins, derived from the cap and the tree's sizes.

    **The single definition of arm membership**, called by the Snakefile's
    fan-out helpers, the summary and the sampler's own defensive check. It cannot
    be written in the config -- it depends on the tree's clade sizes, which only
    the ``get_clades`` checkpoint knows -- so a second derivation anywhere would
    let the fan-out and the sampler disagree about one clade, producing an arm
    whose "unchanged" clades are quietly thinned or whose thinned clades are
    quietly read from the full tree.

    ``clade_sizes`` order is preserved, so the unit inherits the full roster's
    order and the aggregate rules concatenate their inputs in the same sequence
    the full tree does -- which is what an untouched clade's byte-equality rests
    on.

    Membership is ``n_clade > cap``, **strictly**: a clade of exactly ``cap``
    samples is already at the target, so thinning it would write an identical
    file into the arm directory and make the "arm directory holds only thinned
    clades" assertion false for a clade that is not thinned.
    """
    cap = int(spec["cap"])
    return [c for c, n in clade_sizes.items() if int(n) > cap]


def arm_units(arms, clade_sizes):
    """``{arm: [clades it thins]}`` for every arm, resolved against one tree."""
    return {arm: resolve_unit(spec, clade_sizes) for arm, spec in arms.items()}


def depth_label(spec):
    """A short human-readable depth, for logs and error messages."""
    return f"cap={cap_token(spec['cap'])}"


# Keys a ``clade_subsamples`` entry may set. ``unit`` and ``sampling_fractions``
# are absent deliberately and are rejected *by name* below rather than as
# unknown keys, so a config written against the removed frac kind gets told what
# happened instead of a generic schema error.
_ENTRY_KEYS = {"max_samples", "seed"}
_REMOVED_KEYS = {"unit", "sampling_fractions"}


def validate_entries(entries, known_clades=None):
    """Validate ``clade_subsamples`` and return the expanded arm table.

    Returns ``{arm_name: spec}`` in config order, where ``spec`` carries ``seed``
    and ``cap``. ``known_clades`` is accepted and unused: an arm's membership is
    derived from the clade sizes rather than named, so there is no clade list in
    the config that could name a clade this tree does not process. The parameter
    is kept because the Snakefile passes it at two call sites and dropping it
    would make a future named-membership check a signature change rather than a
    body change.

    Every failure raises at DAG-construction time with the offending numbers in
    the message, rather than being clamped: a silently clamped cap produces a
    well-formed arm at the wrong depth, which is indistinguishable from the depth
    that was asked for.
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
        if not isinstance(entry, dict):
            raise ValueError(
                f"config 'clade_subsamples[{i}]' must be a mapping with key "
                f"'max_samples' and optionally 'seed'; got {entry!r}"
            )
        removed = sorted(set(entry) & _REMOVED_KEYS)
        if removed:
            raise ValueError(
                f"config 'clade_subsamples[{i}]' sets {removed}, which belonged "
                f"to the fraction arms removed in #151. An arm now thins every "
                f"clade above 'max_samples' down to it, so its membership is "
                f"derived from the tree's clade sizes rather than named -- "
                f"replace the entry with {{'max_samples': [...], 'seed': N}}."
            )
        unknown = sorted(set(entry) - _ENTRY_KEYS)
        if unknown:
            raise ValueError(
                f"config 'clade_subsamples[{i}]' sets unknown key(s) {unknown}; "
                f"allowed keys are {sorted(_ENTRY_KEYS)}."
            )
        if "max_samples" not in entry:
            raise ValueError(
                f"config 'clade_subsamples[{i}]' must set 'max_samples', a list "
                f"of per-clade sample caps; got keys {sorted(entry)}."
            )
        seed = entry.get("seed", 0)
        if not isinstance(seed, int) or isinstance(seed, bool):
            raise ValueError(
                f"config 'clade_subsamples[{i}].seed' must be an integer; got "
                f"{seed!r}"
            )
        caps = entry["max_samples"]
        if not isinstance(caps, list) or not caps:
            raise ValueError(
                f"config 'clade_subsamples[{i}].max_samples' must be a "
                f"non-empty list of sample counts; got {caps!r}"
            )
        for cap in caps:
            if not isinstance(cap, int) or isinstance(cap, bool):
                raise ValueError(
                    f"config 'clade_subsamples[{i}].max_samples' contains "
                    f"{cap!r}, which is not an integer sample count. A cap is a "
                    f"number of sequences, not a fraction."
                )
            if cap < 1:
                raise ValueError(
                    f"config 'clade_subsamples[{i}].max_samples' contains {cap}, "
                    f"which would retain no samples. `matUtils extract -s` on an "
                    f"empty sample list is undefined, so this is a hard failure "
                    f"rather than a clamp."
                )
            name = arm_name(cap, seed)
            if name in arms:
                raise ValueError(
                    f"config 'clade_subsamples' yields the arm name {name!r} "
                    f"twice. Two entries share a cap and a seed."
                )
            arms[name] = {"seed": seed, "cap": int(cap)}
    return arms


def check_arm_sizes(arms, clade_sizes, min_clade_size=0):
    """Raise on any (arm, clade) whose draw would be empty; report the rest.

    Returns ``{arm: {clade: n_retained}}`` for the clades an arm keeps. A clade
    thinned below ``min_clade_size`` is merely *excluded* from the arm, since
    that threshold expresses how big a clade the user wants rather than what can
    be built at all; a clade thinned to zero samples is a **hard failure**,
    because ``matUtils extract -s`` on an empty list is undefined.

    ``get_clades`` applies ``min_clade_size`` to the full metadata and the
    checkpoint's consumers read only its output, so neither can see a clade that
    falls below the threshold only after thinning -- which is why the arithmetic
    happens here, over the per-clade counts that checkpoint already wrote. The
    exclusion is reachable exactly when ``min_clade_size > cap``, which is a real
    configuration.

    The empty-draw branch is *unreachable* as the two sides currently stand --
    ``resolve_unit`` admits a clade only when ``n > cap`` and
    ``validate_entries`` has already forced ``cap >= 1``, so the retained count
    is ``cap`` and at least 1. It is kept rather than dropped because it is the
    one place that says so: a change to either side that let an empty list reach
    ``matUtils`` would otherwise fail inside a C++ tool with no attribution.
    """
    kept = {}
    for arm, spec in arms.items():
        empty = []
        kept[arm] = {}
        # ``resolve_unit`` is the single definition of membership.
        for clade in resolve_unit(spec, clade_sizes):
            n = int(clade_sizes[clade])
            k = retained_count(n, spec["cap"])
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
                f"rather than a skip. Raise the cap."
            )
    return kept
