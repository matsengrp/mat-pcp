#!/usr/bin/env python3
"""Pin the two-era neutral-rate fit and its per-clade dispatch (#137 Task 2).

``modules/era.py`` (pinned by ``tests/test_era_mapping.py``) says *which* era a
clade belongs to. This file pins what the pipeline then *does* with that: one GLM
fit per era, each clade served by its own era's rate table.

The failure mode throughout is a misrouted clade rather than a mis-specified
model. Both eras' tables are well-formed 512-row grids either way, every clade
gets a finite ``expected_count`` either way, and the two estimators
(``delta_fitness_bayesian`` and ``delta_fitness_ratio``) are corrupted
*consistently*, so no internal cross-check disagrees -- the #129/#130/#133 shape.
That is why these are unit tests over constructed data with a known answer rather
than pipeline assertions.

What each check exists for:

1. **The pooled path is the degenerate per-era path, exactly.**
   ``add_expected_count`` is now a one-era call to
   ``add_expected_count_by_era``, so ISSUE-137 criterion 5 ("single-era
   degeneracy is exact") is structural. The check *demonstrates* that by running
   both entry points on the same data and requiring bit-level equality -- if a
   future edit re-implements the pooled path separately, this fails.

2. **No positional dependence.** The direct regression test against upstream's
   ``clades_pre_o = clades[0:11]``, which indexes ``count_df.clade.unique()`` and
   is correct only because the counts file happens to arrive in Nextstrain clade
   order. Shuffling the row order of both the training and the count frames must
   change nothing. Note the training frame is shuffled too: the GLM's design
   matrix is built row-wise, and a fit that depended on row order would make
   every other check here unreproducible.

3. **Each clade is served by its own era's table.** Two eras are trained on data
   differing *only* in G->T counts by a known factor, and each clade's G->T
   expected counts must reflect its own era. The check also confirms the two
   eras' A->G expectations agree, which is what stops it from passing on a
   dispatch that simply serves everyone the same (wrong) table -- with a single
   table the G->T ratio would be 1.0 for both, and with the eras swapped it
   would be inverted. Both are distinguished.

4. **The eras partition the clades.** Disjoint, covering, and no row left with a
   null ``expected_count``. Disjointness is not free: the dispatch writes per
   era in sequence, so a clade in both lists would be silently overwritten by
   whichever era ran second.

5. **A clade routed to an unfitted era raises**, naming both the clade and the
   fitted eras. This is the failure that must stay *distinct* from the legitimate
   empty-era skip: an absent era is handled by not requesting its table at DAG
   construction, so reaching the dispatch with an unfitted era means the two
   disagree. One ``try/except`` around both would turn a broken fit into a
   silently absent table.

6. **Synthetic two-era recovery.** Counts are generated from two rate tables
   differing only in G->T by a known factor; the *fitted* tables must recover
   that factor for G->T and ~1.0 for every other type. Check 3 pins the
   dispatch; this pins that the fit itself sees the era difference, which is the
   thing the release is for.

7. **The two tables share one grid.** 512 rows, identical keys, so the era
   dimension is the only difference between them. A fit that silently lost a
   condition would still produce a usable-looking table.

8. **The real ``curate_counts.py`` partitions by era without moving the
   conserved-site set.** Executes the shipped script three times against a stub
   ``snakemake`` (pooled, then each era) and requires: the pooled output
   unchanged, the two era outputs partitioning the counts, and -- the load-bearing
   one -- **the same conserved-site set in all three**. Upstream computes
   conservation over all clades and only then queries by era; computing it within
   era would shrink each era's training sites and break criterion 5's guard while
   still producing three plausible files. Re-implementing the script's filter here
   instead of running it would be the CLAUDE.md pattern-10 vacuity, so the script
   is exec'd.

Usage:
    python tests/test_era_rate_dispatch.py
"""

import runpy
import sys
import tempfile
import traceback
from importlib import util as importlib_util
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
MODULES = REPO_ROOT / "workflow" / "modules"
SCRIPTS = REPO_ROOT / "workflow" / "scripts"

# Every shipped config sets these; restated rather than read so the checks below
# do not silently change meaning when a config is retuned.
BOUNDARIES = {"CT": 13467, "AT": 21562, "CG": 21562, "GC": 21562}

# The full condition grid: 4 boundary types x 16 motifs x 2 unpaired x 2 sides,
# plus 8 non-boundary types x 16 motifs x 2 unpaired x 1 side.
GRID_ROWS = 512

#: The known era difference the synthetic checks inject: pre-Omicron G->T counts
#: are this multiple of Omicron's. Chosen well above the 1.3x that ISSUE-137's
#: validation list asks the real tables to clear, so a recovered ratio is
#: unambiguous rather than within fit noise.
GT_ERA_FACTOR = 2.0

LETTERS = ["A", "C", "G", "T"]
MUT_TYPES = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]

# Clades from both eras, deliberately interleaved so that no positional slice of
# this list separates the eras -- upstream's clades[0:11] would mis-split it.
PRE_OMICRON_CLADES = ["20A", "21J"]
OMICRON_CLADES = ["21K", "23A"]
MIXED_CLADES = ["20A", "21K", "21J", "23A"]


def _load(module_name):
    path = MODULES / f"{module_name}.py"
    spec = importlib_util.spec_from_file_location(module_name, path)
    module = importlib_util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def synthetic_training_frame(gt_scale=1.0, seed=0):
    """A curated-counts frame spanning the whole condition grid.

    One row per (mut_type, motif, unpaired, boundary side) condition, with counts
    drawn from a fixed per-condition mean so the GLM has a real signal to fit.
    ``gt_scale`` multiplies the G->T means only, which is the era difference the
    fit is asked to recover.

    ``nt_site`` is placed on the correct side of the type's boundary for the
    condition's ``nt_site_before_boundary`` label, because ``glm.py`` derives the
    boundary feature from ``nt_site`` at train time while ``rates.py`` serves rates
    keyed on the label. Those two disagreeing at the boundary site is #125; here
    they must agree or the fit and the lookup describe different conditions.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for mut_type in MUT_TYPES:
        sides = [False, True] if mut_type in BOUNDARIES else [False]
        for left in LETTERS:
            for right in LETTERS:
                motif = left + mut_type[0] + right
                for unpaired in (0, 1):
                    for before in sides:
                        boundary = BOUNDARIES.get(mut_type, 0)
                        # Strictly below / strictly above, so the label and the
                        # nt_site-derived feature cannot disagree.
                        nt_site = boundary - 100 if before else boundary + 100
                        mean = 20.0 + 8.0 * unpaired + 4.0 * LETTERS.index(left)
                        if mut_type == "GT":
                            mean *= gt_scale
                        rows.append(
                            {
                                "clade": "pooled",
                                "nt_site": int(nt_site),
                                "nt_mutation": f"{mut_type[0]}{int(nt_site)}{mut_type[1]}",
                                "mut_type": mut_type,
                                "motif": motif,
                                "unpaired": unpaired,
                                "nt_site_before_boundary": before,
                                "mut_class": "synonymous",
                                "actual_count": int(rng.poisson(mean)),
                            }
                        )
    return pd.DataFrame(rows)


def synthetic_count_frame(clades):
    """A count frame over ``clades``, one row per clade per condition.

    Sites and conditions match :func:`synthetic_training_frame`, so every row is
    in the training site set and therefore contributes to ``T``.
    """
    template = synthetic_training_frame()
    frames = []
    for i, clade in enumerate(clades):
        frame = template.copy()
        frame["clade"] = clade
        # A per-clade count scale, so T differs between clades and a dispatch
        # that ignored the clade loop could not pass by coincidence.
        frame["actual_count"] = frame["actual_count"] * (i + 1)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True).reset_index(drop=True)


def check_pooled_is_degenerate_two_era():
    """Criterion 5: the pooled entry point equals the one-era per-era call."""
    rates = _load("rates")

    train_a = synthetic_training_frame()
    train_b = synthetic_training_frame()
    counts_a = synthetic_count_frame(PRE_OMICRON_CLADES)
    counts_b = counts_a.copy()

    rates.add_expected_count(
        train_df=train_a,
        count_df=counts_a,
        clades=PRE_OMICRON_CLADES,
        lightswitch_boundaries=BOUNDARIES,
    )
    rates.add_expected_count_by_era(
        train_by_era={"pre_omicron": train_b},
        count_df=counts_b,
        clades=PRE_OMICRON_CLADES,
        era_of_clade=lambda _clade: "pre_omicron",
        lightswitch_boundaries=BOUNDARIES,
    )

    for column in ("expected_count", "tau_squared"):
        left = counts_a[column].to_numpy()
        right = counts_b[column].to_numpy()
        if not np.array_equal(left, right):
            worst = np.nanmax(np.abs(left - right))
            raise AssertionError(
                f"pooled add_expected_count and the one-era "
                f"add_expected_count_by_era disagree on {column} "
                f"(max abs difference {worst:.3g}). The pooled path is supposed "
                f"to BE the degenerate per-era path, not a second implementation "
                f"of it (ISSUE-137 criterion 5)."
            )
    if not counts_a["expected_count"].notna().all():
        raise AssertionError("the pooled fit left null expected_count values")


def check_no_positional_dependence():
    """Row order of the training and count frames must not change any result.

    The regression test against upstream's ``clades[0:11]`` positional dispatch.
    """
    rates = _load("rates")
    era = _load("era")

    train_by_era = {
        "pre_omicron": synthetic_training_frame(seed=1),
        "omicron": synthetic_training_frame(gt_scale=1 / GT_ERA_FACTOR, seed=2),
    }
    counts = synthetic_count_frame(MIXED_CLADES)

    rates.add_expected_count_by_era(
        train_by_era={k: v.copy() for k, v in train_by_era.items()},
        count_df=counts,
        clades=MIXED_CLADES,
        era_of_clade=era.era_for_clade,
        lightswitch_boundaries=BOUNDARIES,
    )
    ordered = counts.set_index(
        ["clade", "mut_type", "motif", "unpaired", "nt_site_before_boundary"]
    )["expected_count"]

    rng = np.random.default_rng(7)
    shuffled_counts = counts.drop(columns=["expected_count", "tau_squared"]).sample(
        frac=1.0, random_state=11
    )
    shuffled_counts = shuffled_counts.reset_index(drop=True)
    shuffled_train = {
        key: frame.sample(frac=1.0, random_state=int(rng.integers(1 << 30))).reset_index(
            drop=True
        )
        for key, frame in train_by_era.items()
    }
    # Reverse the clade order too: the dispatch loops over this list.
    rates.add_expected_count_by_era(
        train_by_era=shuffled_train,
        count_df=shuffled_counts,
        clades=list(reversed(MIXED_CLADES)),
        era_of_clade=era.era_for_clade,
        lightswitch_boundaries=BOUNDARIES,
    )
    reordered = shuffled_counts.set_index(
        ["clade", "mut_type", "motif", "unpaired", "nt_site_before_boundary"]
    )["expected_count"]

    aligned = reordered.reindex(ordered.index)
    if aligned.isna().any():
        raise AssertionError(
            "shuffling the input changed which (clade, condition) keys exist"
        )
    worst = float(np.max(np.abs(aligned.to_numpy() - ordered.to_numpy())))
    if worst > 1e-9:
        raise AssertionError(
            f"expected counts depend on input row order (max abs difference "
            f"{worst:.3g}). This is upstream's clades[0:11] failure mode: a "
            f"dispatch that reads position rather than the clade label."
        )


def check_correct_era_table_dispatched():
    """Each clade's G->T expectation must come from its own era's fit."""
    rates = _load("rates")
    era = _load("era")

    # Two training sets identical except that pre-Omicron G->T counts are
    # GT_ERA_FACTOR times Omicron's. Same seed, so nothing else differs.
    train_by_era = {
        "pre_omicron": synthetic_training_frame(gt_scale=GT_ERA_FACTOR, seed=3),
        "omicron": synthetic_training_frame(gt_scale=1.0, seed=3),
    }
    counts = synthetic_count_frame(MIXED_CLADES)
    rates.add_expected_count_by_era(
        train_by_era=train_by_era,
        count_df=counts,
        clades=MIXED_CLADES,
        era_of_clade=era.era_for_clade,
        lightswitch_boundaries=BOUNDARIES,
    )

    # Per clade, the mean expected count for G->T relative to a type with no
    # injected era difference. Normalizing by A->G divides out the per-clade T,
    # so what is left is the era's own rate structure.
    def ratio(clade):
        sub = counts[counts["clade"] == clade]
        gt = sub.loc[sub["mut_type"] == "GT", "expected_count"].mean()
        ag = sub.loc[sub["mut_type"] == "AG", "expected_count"].mean()
        return gt / ag

    pre_ratios = [ratio(c) for c in PRE_OMICRON_CLADES]
    om_ratios = [ratio(c) for c in OMICRON_CLADES]

    # Within an era the ratio must be identical (T cancels); across eras it must
    # differ by the injected factor.
    for name, values in (("pre_omicron", pre_ratios), ("omicron", om_ratios)):
        if np.ptp(values) > 1e-9:
            raise AssertionError(
                f"clades within era {name} got different G->T/A->G structure "
                f"({values}); they share one rate table, so this can only mean "
                f"the tables were rebuilt per clade."
            )
    observed = pre_ratios[0] / om_ratios[0]
    if not 1.5 < observed < 2.6:
        raise AssertionError(
            f"pre-Omicron/Omicron G->T structure ratio is {observed:.3f}, "
            f"expected near the injected {GT_ERA_FACTOR}. A ratio near 1.0 means "
            f"one table served every clade; near {1 / GT_ERA_FACTOR:.2f} means "
            f"the eras are wired backwards."
        )


def check_eras_partition_the_clades():
    """Disjoint, covering, and no row left unserved."""
    rates = _load("rates")
    era = _load("era")

    partition = era.partition_clades(MIXED_CLADES)
    buckets = list(partition.values())
    for i, left in enumerate(buckets):
        for right in buckets[i + 1:]:
            overlap = set(left) & set(right)
            if overlap:
                raise AssertionError(
                    f"eras are not disjoint: {sorted(overlap)} appear in two "
                    f"buckets. The dispatch writes per era in sequence, so an "
                    f"overlapping clade is silently overwritten by the second."
                )
    covered = sorted(c for bucket in buckets for c in bucket)
    if covered != sorted(MIXED_CLADES):
        raise AssertionError(
            f"the eras do not cover the clade list: {covered} vs "
            f"{sorted(MIXED_CLADES)}"
        )

    counts = synthetic_count_frame(MIXED_CLADES)
    rates.add_expected_count_by_era(
        train_by_era={
            "pre_omicron": synthetic_training_frame(seed=4),
            "omicron": synthetic_training_frame(seed=5),
        },
        count_df=counts,
        clades=MIXED_CLADES,
        era_of_clade=era.era_for_clade,
        lightswitch_boundaries=BOUNDARIES,
    )
    n_null = int(counts["expected_count"].isna().sum())
    if n_null:
        raise AssertionError(
            f"{n_null} count rows have a null expected_count after dispatch; "
            f"every clade present must be served by exactly one era's model"
        )


def check_unfitted_era_raises():
    """A clade routed to an era with no fit must raise, naming clade and eras."""
    rates = _load("rates")
    era = _load("era")

    counts = synthetic_count_frame(MIXED_CLADES)
    try:
        rates.add_expected_count_by_era(
            # Only pre-Omicron fitted, but Omicron clades are present.
            train_by_era={"pre_omicron": synthetic_training_frame(seed=6)},
            count_df=counts,
            clades=MIXED_CLADES,
            era_of_clade=era.era_for_clade,
            lightswitch_boundaries=BOUNDARIES,
        )
    except ValueError as exc:
        message = str(exc)
        for needed in ("21K", "omicron", "pre_omicron"):
            if needed not in message:
                raise AssertionError(
                    f"the raise for an unfitted era does not name {needed!r}: "
                    f"{message}"
                )
    else:
        raise AssertionError(
            "dispatching an Omicron clade against a pre-Omicron-only fit did "
            "not raise. Silently serving the only available table is exactly "
            "the fallback #137 exists to remove."
        )

    # And the empty case, which must not be mistaken for a legitimate skip.
    try:
        rates.add_expected_count_by_era(
            train_by_era={},
            count_df=counts,
            clades=MIXED_CLADES,
            era_of_clade=era.era_for_clade,
            lightswitch_boundaries=BOUNDARIES,
        )
    except ValueError:
        pass
    else:
        raise AssertionError("an empty train_by_era did not raise")


def check_synthetic_two_era_recovery():
    """The fitted tables must recover the injected G->T era factor."""
    rates = _load("rates")

    tables = {}
    for name, scale, seed in (
        ("pre_omicron", GT_ERA_FACTOR, 21),
        ("omicron", 1.0, 21),
    ):
        train = synthetic_training_frame(gt_scale=scale, seed=seed)
        rate = rates.Rates(BOUNDARIES)
        rate.populate_rates(train)
        tables[name] = rate.rates.set_index(
            ["mut_type", "motif", "unpaired", "nt_site_before_boundary"]
        )["rate"]

    # populate_rates rescales by n_sites / total mutations, so raising G->T
    # counts lowers every other type's rate slightly. Normalizing each type by
    # the table's non-G->T median removes that global factor, which is what makes
    # "~1.0 for every other type" a meaningful requirement.
    ratio = tables["pre_omicron"] / tables["omicron"]
    by_type = ratio.groupby("mut_type").mean()
    baseline = by_type.drop("GT").median()
    normalized = by_type / baseline

    if not 1.5 < normalized["GT"] < 2.6:
        raise AssertionError(
            f"the fit did not recover the injected G->T era factor "
            f"{GT_ERA_FACTOR}: normalized ratio {normalized['GT']:.3f}. The two "
            f"tables differ only in G->T counts, so this is the era signal the "
            f"release exists to capture."
        )
    others = normalized.drop("GT")
    if float(np.max(np.abs(others - 1.0))) > 0.15:
        raise AssertionError(
            f"types with no injected era difference moved: "
            f"{others[np.abs(others - 1.0) > 0.15].to_dict()}. The era split "
            f"must not manufacture structure where there is none."
        )


def check_era_tables_share_one_grid():
    """Both eras' tables must be the same 512-row grid."""
    rates = _load("rates")

    keys = {}
    for name, scale in (("pre_omicron", GT_ERA_FACTOR), ("omicron", 1.0)):
        rate = rates.Rates(BOUNDARIES)
        rate.populate_rates(synthetic_training_frame(gt_scale=scale, seed=31))
        table = rate.rates
        if len(table) != GRID_ROWS:
            raise AssertionError(
                f"era {name} table has {len(table)} rows, expected {GRID_ROWS}"
            )
        keys[name] = set(
            map(
                tuple,
                table[
                    ["mut_type", "motif", "unpaired", "nt_site_before_boundary"]
                ].to_numpy(),
            )
        )
    if keys["pre_omicron"] != keys["omicron"]:
        only = keys["pre_omicron"] ^ keys["omicron"]
        raise AssertionError(
            f"the two era tables do not share a grid; {len(only)} conditions "
            f"appear in only one of them. The era must be the ONLY difference "
            f"between the tables, or they cannot be compared cell by cell."
        )


class _StubSnakemake:
    """Minimal stand-in for the ``snakemake`` global a script sees."""

    class _Named(dict):
        def __getattr__(self, name):
            try:
                return self[name]
            except KeyError as exc:
                raise AttributeError(name) from exc

        def get(self, name, default=None):
            return dict.get(self, name, default)

    def __init__(self, inputs, outputs, params, log_path):
        self.input = self._Named(inputs)
        self.output = self._Named(outputs)
        self.params = self._Named(params)
        self.log = [str(log_path)]


def _curation_fixture(tmp):
    """An annotated-counts / founder pair spanning both eras.

    Site 100 is conserved in every clade; site 200 is conserved in the Omicron
    clades only. So *if* conservation were computed within era, the pre-Omicron
    training set would keep site 100 only while Omicron's would keep both -- a
    difference this fixture makes visible.
    """
    clades = PRE_OMICRON_CLADES + OMICRON_CLADES
    founder_rows = []
    count_rows = []
    for clade in clades:
        for site in (100, 200):
            divergent = site == 200 and clade in PRE_OMICRON_CLADES
            founder_rows.append(
                {
                    "clade": clade,
                    "nt_site": site,
                    "codon": "AAT" if divergent else "AAA",
                    "ref_codon": "AAA",
                    "motif": "AAT" if divergent else "AAA",
                    "ref_motif": "AAA",
                }
            )
            count_rows.append(
                {
                    "clade": clade,
                    "nt_site": site,
                    "nt_mutation": f"A{site}G",
                    "mut_type": "AG",
                    "motif": "AAA",
                    "unpaired": 0,
                    "nt_site_before_boundary": False,
                    "synonymous": True,
                    "noncoding": False,
                    "gene": "S",
                    "actual_count": 5,
                    "exclude": False,
                    "masked_in_usher": False,
                }
            )
    counts_path = tmp / "annotated_counts.csv"
    founder_path = tmp / "clade_founder_annotated.csv"
    pd.DataFrame(count_rows).to_csv(counts_path, index=False)
    pd.DataFrame(founder_rows).to_csv(founder_path, index=False)
    return counts_path, founder_path


def _run_curate(tmp, counts_path, founder_path, era_name, tag):
    out = tmp / f"curated_{tag}.csv"
    stub = _StubSnakemake(
        inputs={"counts": str(counts_path), "clade_founder": str(founder_path)},
        outputs={"curated": str(out)},
        params={"era": era_name},
        log_path=tmp / f"curate_{tag}.log",
    )
    saved = sys.modules.get("snakemake")
    sys.modules["snakemake"] = stub
    builtins_snakemake = {"snakemake": stub}
    try:
        runpy.run_path(
            str(SCRIPTS / "curate_counts.py"),
            init_globals=builtins_snakemake,
            run_name="__main__",
        )
    finally:
        if saved is None:
            sys.modules.pop("snakemake", None)
        else:
            sys.modules["snakemake"] = saved
    return pd.read_csv(out)


def check_real_curate_counts_partitions_by_era():
    """Execute the shipped ``curate_counts.py`` pooled and per era."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        counts_path, founder_path = _curation_fixture(tmp)

        pooled = _run_curate(tmp, counts_path, founder_path, None, "pooled")
        pre = _run_curate(tmp, counts_path, founder_path, "pre_omicron", "pre")
        om = _run_curate(tmp, counts_path, founder_path, "omicron", "om")

        # The conserved-site set must be identical in all three. Site 200 is
        # divergent in the pre-Omicron clades, so conservation computed WITHIN
        # era would admit it for Omicron and not for pre-Omicron.
        sites = {
            "pooled": set(pooled["nt_site"]),
            "pre_omicron": set(pre["nt_site"]),
            "omicron": set(om["nt_site"]),
        }
        if sites["pooled"] != {100}:
            raise AssertionError(
                f"the pooled conserved-site set is {sorted(sites['pooled'])}, "
                f"expected just site 100 (site 200 is divergent in two clades)"
            )
        for name in ("pre_omicron", "omicron"):
            if sites[name] != sites["pooled"]:
                raise AssertionError(
                    f"era {name} trains on sites {sorted(sites[name])} but the "
                    f"pooled fit uses {sorted(sites['pooled'])}. Conservation "
                    f"must be computed over ALL clades before the era "
                    f"partition -- computing it within era changes each era's "
                    f"training site set and breaks parity with upstream and "
                    f"with the pooled control arm."
                )

        # The eras must partition the counts: each era carries its own clades'
        # counts only, and together they sum to the pooled total.
        pooled_total = int(pooled["actual_count"].sum())
        split_total = int(pre["actual_count"].sum()) + int(om["actual_count"].sum())
        if split_total != pooled_total:
            raise AssertionError(
                f"the era partition does not conserve counts: "
                f"{split_total} split vs {pooled_total} pooled"
            )
        expected_pre = 5 * len(PRE_OMICRON_CLADES)
        if int(pre["actual_count"].sum()) != expected_pre:
            raise AssertionError(
                f"pre-Omicron curated counts sum to "
                f"{int(pre['actual_count'].sum())}, expected {expected_pre} "
                f"({len(PRE_OMICRON_CLADES)} clades x 5 at the one conserved site)"
            )

        # An era with no clades must raise rather than write an empty table:
        # the legitimate skip happens at DAG construction, not here.
        only_pre = tmp / "only_pre_counts.csv"
        only_pre_founder = tmp / "only_pre_founder.csv"
        pd.read_csv(counts_path).query("clade in @PRE_OMICRON_CLADES").to_csv(
            only_pre, index=False
        )
        pd.read_csv(founder_path).query("clade in @PRE_OMICRON_CLADES").to_csv(
            only_pre_founder, index=False
        )
        try:
            _run_curate(tmp, only_pre, only_pre_founder, "omicron", "empty")
        except ValueError as exc:
            if "omicron" not in str(exc):
                raise AssertionError(
                    f"the empty-era raise does not name the era: {exc}"
                )
        else:
            raise AssertionError(
                "curating an era with no clades present did not raise. Writing "
                "an empty training table would ship a rate table that looks "
                "like a fit and contains nothing."
            )


def check_real_expected_counts_script_dispatches_two_eras():
    """Execute the shipped ``expected_counts_by_clade.py`` against two eras.

    **This is the only two-era coverage of the production script.** Every shipped
    config resolves to a single era in CI (``config_test_data.yaml``'s
    ``max_clades: 3`` gives 19A/19B/20A, all pre-Omicron), on which the script
    reduces to the v9 pooled fit -- so a pipeline run cannot tell whether the
    ``curated_{era}`` inputs are read at all, whether they are matched to the
    right era, or whether ``era_for_clade`` is consulted. Checking
    ``rates.add_expected_count_by_era`` directly and calling that two-era coverage
    would be the CLAUDE.md pattern-8 vacuity: a test named for the dispatch that
    never runs the code doing the dispatching.

    The era labels are recoverable from the output because the two training sets
    differ in G->T by a known factor, so a swapped pair inverts the ratio rather
    than merely perturbing it.
    """
    era = _load("era")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)

        # Curated training tables per era, differing only in G->T. Written with
        # the columns load_curated_muts requires (it filters on synonymous /
        # noncoding / gene -- it is the whole neutral-rate restriction, #126).
        paths = {}
        for era_name, scale, seed in (
            ("pre_omicron", GT_ERA_FACTOR, 41),
            ("omicron", 1.0, 41),
        ):
            frame = synthetic_training_frame(gt_scale=scale, seed=seed)
            frame["synonymous"] = True
            frame["noncoding"] = False
            frame["gene"] = "S"
            path = tmp / f"curated_counts_{era_name}.csv"
            frame.to_csv(path, index=False)
            paths[era_name] = path

        counts = synthetic_count_frame(MIXED_CLADES)
        counts_path = tmp / "annotated_counts.csv"
        counts.to_csv(counts_path, index=False)
        out = tmp / "expected_counts_by_clade.csv"

        stub = _StubSnakemake(
            inputs={
                "counts": str(counts_path),
                "curated_pre_omicron": str(paths["pre_omicron"]),
                "curated_omicron": str(paths["omicron"]),
            },
            outputs={"csv": str(out)},
            params={"lightswitch_boundaries": BOUNDARIES},
            log_path=tmp / "expected_counts_by_clade.log",
        )
        saved = sys.modules.get("snakemake")
        sys.modules["snakemake"] = stub
        try:
            runpy.run_path(
                str(SCRIPTS / "expected_counts_by_clade.py"),
                init_globals={"snakemake": stub},
                run_name="__main__",
            )
        finally:
            if saved is None:
                sys.modules.pop("snakemake", None)
            else:
                sys.modules["snakemake"] = saved

        result = pd.read_csv(out)
        if result["expected_count"].isna().any():
            raise AssertionError(
                "the script left null expected_count values, which downstream "
                "reads as a legitimate near-zero expectation"
            )

        def ratio(clade):
            sub = result[result["clade"] == clade]
            gt = sub.loc[sub["mut_type"] == "GT", "expected_count"].mean()
            ag = sub.loc[sub["mut_type"] == "AG", "expected_count"].mean()
            return gt / ag

        pre = ratio(PRE_OMICRON_CLADES[0])
        om = ratio(OMICRON_CLADES[0])
        observed = pre / om
        if not 1.5 < observed < 2.6:
            raise AssertionError(
                f"the shipped script's G->T era ratio is {observed:.3f}, expected "
                f"near the injected {GT_ERA_FACTOR}. Near 1.0 means one era's "
                f"table served every clade; near {1 / GT_ERA_FACTOR:.2f} means "
                f"the curated_{{era}} inputs are matched to the wrong eras."
            )

        # An era present in the counts but absent from the inputs must raise,
        # not quietly serve the one table that did arrive.
        stub_missing = _StubSnakemake(
            inputs={
                "counts": str(counts_path),
                "curated_pre_omicron": str(paths["pre_omicron"]),
            },
            outputs={"csv": str(tmp / "partial.csv")},
            params={"lightswitch_boundaries": BOUNDARIES},
            log_path=tmp / "partial.log",
        )
        saved = sys.modules.get("snakemake")
        sys.modules["snakemake"] = stub_missing
        try:
            runpy.run_path(
                str(SCRIPTS / "expected_counts_by_clade.py"),
                init_globals={"snakemake": stub_missing},
                run_name="__main__",
            )
        except ValueError as exc:
            if "21K" not in str(exc):
                raise AssertionError(
                    f"the raise does not name the unserved clade: {exc}"
                )
        else:
            raise AssertionError(
                "the script accepted Omicron clades with no Omicron training "
                "table and produced a full output file"
            )
        finally:
            if saved is None:
                sys.modules.pop("snakemake", None)
            else:
                sys.modules["snakemake"] = saved


CHECKS = [
    ("the pooled path is the degenerate two-era path", check_pooled_is_degenerate_two_era),
    ("nothing depends on input row order", check_no_positional_dependence),
    ("each clade is served its own era's table", check_correct_era_table_dispatched),
    ("the eras partition the clades", check_eras_partition_the_clades),
    ("a clade routed to an unfitted era raises", check_unfitted_era_raises),
    ("the fit recovers a synthetic era difference", check_synthetic_two_era_recovery),
    ("both era tables share one grid", check_era_tables_share_one_grid),
    (
        "the real curate_counts.py partitions by era without moving conserved sites",
        check_real_curate_counts_partitions_by_era,
    ),
    (
        "the real expected_counts_by_clade.py dispatches two eras",
        check_real_expected_counts_script_dispatches_two_eras,
    ),
]


def run_all_checks():
    failures = []
    for label, check in CHECKS:
        try:
            check()
        except Exception as exc:  # noqa: BLE001 - report every check
            failures.append(f"{label}: {exc}")
            print(f"FAIL  {label}")
            traceback.print_exc()
        else:
            print(f"ok    {label}")
    return failures


def test_era_rate_dispatch():
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
        for failure in failures:
            print(f"  - {failure}")
        sys.exit(1)
    print(
        "\nPASSED: one fit per era, each clade served by its own era's table, "
        "the pooled path exactly degenerate, and conserved sites computed before "
        "the partition."
    )


if __name__ == "__main__":
    main()
