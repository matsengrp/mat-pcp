#!/usr/bin/env python3
"""Clade subsampling arms: naming, validation and the draw itself (#151).

Per-clade fitness effects are estimated from one clade's own subtree, so a
thinly sampled clade gets a noisy estimate and cannot borrow strength from the
rest of the tree. Measuring how thin a clade has to get before a pooled model
wins needs the *same* clade re-estimated at several sampling depths, and this
module is the hook: an **arm** is one sampling fraction applied to one **unit**
of clades, and everything an arm rebuilds follows from re-extracting each unit
clade's subtree from a subsampled sample list.

Three things about the shape are load-bearing, and each is why a function lives
here rather than inline in a script:

*The draw is a pure function of ``(unit, seed, clade, frac)``.* ``arm_rng``
derives its seed by hashing that tuple, never from iteration order or a
module-level default, so an arm's sample list cannot move when another fraction
is added to or removed from the config. Determinism is then structural rather
than a property of Snakemake happening not to re-run the sampler.

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

# The prefix every arm directory carries. It exists so ``arm_dir``'s optional
# alternative (see the Snakefile's ``wildcard_constraints``) cannot be confused
# with a clade name: no Nextstrain clade starts with ``frac_``.
ARM_PREFIX = "frac"

# Matches an arm directory *including* its trailing slash, or the empty string
# for the full-tree path. Kept here rather than in the Snakefile so the tests
# and the scripts agree with the DAG on what an arm path looks like.
ARM_DIR_PATTERN = rf"|(?:{ARM_PREFIX})_[A-Za-z0-9._+-]+/"

# Matches an arm name (no trailing slash).
ARM_PATTERN = rf"(?:{ARM_PREFIX})_[A-Za-z0-9._+-]+"


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


def retained_count(n_samples, frac):
    """``round(frac * n_samples)`` -- the arm's exact sample count.

    Python's banker's rounding is used deliberately and unwrapped: the number is
    reported in ``subsample_summary.csv`` and asserted by the membership check,
    so both sides must come from this one function rather than from two
    re-derivations that agree except at the halfway case.
    """
    return int(round(float(frac) * int(n_samples)))


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

    Returns ``{arm_name: {"unit": [...], "seed": int, "frac": float}}`` in config
    order. Every failure raises at DAG-construction time with the offending
    numbers in the message, rather than being clamped: a silently clamped
    fraction produces a well-formed arm at the wrong depth, which is
    indistinguishable from the depth that was asked for.
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
            "seed",
        }:
            raise ValueError(
                f"config 'clade_subsamples[{i}]' must be a mapping with keys "
                f"'unit', 'sampling_fractions' and optionally 'seed'; got {entry!r}"
            )
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
        seed = entry.get("seed", 0)
        if not isinstance(seed, int) or isinstance(seed, bool):
            raise ValueError(
                f"config 'clade_subsamples[{i}].seed' must be an integer; got "
                f"{seed!r}"
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
            arms[name] = {"unit": list(unit), "seed": seed, "frac": frac}
    return arms


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
    """
    kept = {}
    for arm, spec in arms.items():
        frac = spec["frac"]
        empty = []
        kept[arm] = {}
        for clade in spec["unit"]:
            if clade not in clade_sizes:
                raise ValueError(
                    f"arm {arm!r} names clade {clade!r}, absent from the clade "
                    f"counts checkpoint. Available: {sorted(clade_sizes)}"
                )
            n = int(clade_sizes[clade])
            k = retained_count(n, frac)
            if k < 1:
                empty.append((clade, n, k))
            elif k >= min_clade_size:
                kept[arm][clade] = k
        if empty:
            detail = ", ".join(
                f"{c}: round({frac} * {n}) = {k}" for c, n, k in empty
            )
            raise ValueError(
                f"arm {arm!r} would draw zero samples for {len(empty)} clade(s) "
                f"({detail}). `matUtils extract -s` on an empty sample list is "
                f"undefined, so this is a hard failure rather than a skip. "
                f"Raise the sampling fraction or drop the clade from the unit."
            )
    return kept
