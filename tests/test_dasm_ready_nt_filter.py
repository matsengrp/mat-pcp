#!/usr/bin/env python3
"""Unit tests for the DASM-ready drop filter's nucleotide mask (#158).

**No pipeline run can see this change, so these tests are the only coverage.**
Applying the old codon predicate and the new nucleotide one to the v12 test
tree's PCPs gives *identical* row sets for every gene -- S 8/8, ORF1ab 23/23,
ORF1b 13/13, E 1/1 kept -- because no test-data PCP happens to mutate an
unmasked nucleotide inside a masked codon. So a full ``config_test_data.yaml``
run plus ``validate_test.py`` is green whether this filter is nucleotide-
granular, codon-granular, or shifted by one, and the superset invariant is
VACUOUS in CI. Four things need a fixture:

* **The recovery case itself** -- a PCP whose only mutation sits at an unmasked
  nucleotide of a masked codon. That is the entire point of the change and it
  occurs nowhere in the test data.
* **The off-by-one.** The retired ``{gene}_codon_sites_to_mask.txt`` was
  1-based and its reader subtracted one; ``masked_nt_sites.csv``'s ``nt_index``
  is 0-based and must be read as-is. Both directions of the shift are checked,
  because a mask shifted one way and one the other are different wrong answers.
* **The one-to-many frameshift index.** ORF1b reads genomic 13468 at nt indices
  26 *and* 27, in different codons, so a mask there is two rows. No
  ProblematicSites VCF or config entry falls in ORF1b's leading block, so this
  is unreachable from any run (the same gap #144 and #153 record).
* **Whether the reader's guard fires.** A mask that fails to parse and a mask
  that is empty produce the same observable outcome -- every masked position
  trained on as live -- so the missing-column case must raise rather than fall
  back.

Runs the real ``workflow/scripts/curate_pcps_for_dasm.py`` against a stub
``snakemake`` global, never a re-implementation of it. The mutation tests exec
deliberately corrupted copies of ``workflow/modules/masking.py``, where
``read_masked_nt_indices`` actually lives, so the corruption lands in the file
an operator would edit rather than in a paraphrase of it (CLAUDE.md patterns 10,
14 and 15).

Hermetic: synthetic sequences, a synthetic mask CSV. Nothing reads a pipeline
output.
"""

import csv
import json
import random
import shutil
import sys
import tempfile
import traceback
from importlib import util as importlib_util
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT_REL = "scripts/curate_pcps_for_dasm.py"
MASKING_REL = "modules/masking.py"

# ORF1b's real CDS blocks, as verified in #133 against the shipped founders.
# They share nucleotide 13468 so the -1 ribosomal frameshift stays in frame,
# which is what makes the genomic->index map one-to-many.
ORF1B_BLOCKS = [(13442, 13468), (13468, 21552)]


def _load(module_name):
    """Load a workflow module by path, as the pipeline scripts do."""
    spec = importlib_util.spec_from_file_location(
        module_name, REPO_ROOT / "workflow" / "modules" / f"{module_name}.py"
    )
    module = importlib_util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


masking = _load("masking")
sequences = _load("sequences")


class _Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class _Params(_Namespace):
    """Snakemake's ``params`` supports ``.get``, and the script uses it."""

    def get(self, key, default=None):
        return self.__dict__.get(key, default)


# ---------------------------------------------------------------------------
# The reference implementation
# ---------------------------------------------------------------------------
#
# Plain Python, deliberately: the shipped predicate is vectorized over encoded
# buffers because it runs over 4.3 M PCPs x 13,203 nt on gisaid ORF1a, and a
# fast implementation and a legible one should be checked against each other
# rather than one of them serving as both.


def reference_keep_nt(parent, child, masked_indices):
    """Keep iff some position differs at an index outside ``masked_indices``."""
    return any(
        p != c and i not in masked_indices
        for i, (p, c) in enumerate(zip(parent, child))
    )


def reference_keep_codon(parent, child, masked_codons):
    """The retired codon rule, kept as the thing the new one must dominate.

    ``masked_codons`` is a set of 0-based codon positions. This is what
    ``_has_unmasked_codon_mutation`` did through v12.
    """
    for codon_index, start in enumerate(range(0, len(parent), 3)):
        if codon_index in masked_codons:
            continue
        if parent[start : start + 3] != child[start : start + 3]:
            return True
    return False


def widen(masked_indices):
    """The codon set the retired file held: every codon touching a masked nt.

    Built from the nucleotide set rather than asserted against it independently,
    because the check that the two agreed on disk was ``masked_nt_sites.py``'s
    guard 2 -- which #158 deletes along with the file it compared against. The
    superset invariant below is only meaningful for a codon set that *is* this
    widening, so it is constructed that way here (#158).
    """
    return {index // 3 for index in masked_indices}


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _write_mask(path, rows, columns=("nt_site", "nt_index", "reason")):
    """Write a ``masked_nt_sites.csv``-shaped file.

    ``rows`` is a list of ``(nt_site, nt_index, reason)`` tuples. The real file
    is long-form -- one row per ``(nt_index, reason)`` -- so a two-reason site
    appears twice with one index, and a frameshift site appears twice with two
    indices. Both multiplicities are exercised below.
    """
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(list(columns))
        writer.writerows(rows)
    return path


def _mask_rows(masked_indices, nt_site_of=None):
    """One ``sites_to_exclude`` row per index, with plausible genomic sites."""
    return [
        ((nt_site_of or (lambda i: 21563 + i))(index), index, "sites_to_exclude")
        for index in sorted(masked_indices)
    ]


def _mutate_at(sequence, index):
    """``sequence`` with position ``index`` changed to a different base.

    Chosen rather than written literally: a fixture that substituted the base
    already there would give a PCP with zero differences, which the predicate
    drops for a reason unrelated to the mask -- and the test would report the
    mask as broken.
    """
    replacement = "A" if sequence[index] != "A" else "C"
    return sequence[:index] + replacement + sequence[index + 1 :]


def _write_pcps(path, pcps):
    """Write a per-clade PCP table: ``pcps`` is a list of (name, parent, child)."""
    frame = pd.DataFrame(
        [
            {
                "parent_name": f"node_{name}",
                "child_name": f"child_{name}",
                "parent": parent,
                "child": child,
                "branch_length": 1,
            }
            for name, parent, child in pcps
        ]
    )
    frame.to_csv(path, index=False, compression="gzip")
    return path


def _mutated_workflow(tmp, edits):
    """A throwaway copy of ``workflow/`` with text edits applied.

    Copied rather than edited in place: the script resolves
    ``modules/<name>.py`` relative to its own location, so a mutation to a
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
    return root


def run_curate(
    tmpdir,
    pcps,
    mask_rows,
    mask_columns=("nt_site", "nt_index", "reason"),
    seed=1729,
    val_fraction=0.10,
    workflow_root=None,
):
    """Execute curate_pcps_for_dasm.py the way Snakemake does.

    Returns ``(kept_names, stats, log_text)``. exec rather than import, because
    the script runs top-level statements against the ``snakemake`` global --
    exactly how Snakemake invokes it. ``workflow_root`` points the run at a
    corrupted copy of ``workflow/`` for the mutation tests.
    """
    tmpdir = Path(tmpdir)
    tmpdir.mkdir(parents=True, exist_ok=True)
    root = Path(workflow_root) if workflow_root else REPO_ROOT / "workflow"
    script = root / SCRIPT_REL

    pcp_path = _write_pcps(tmpdir / "clade_pcps.csv.gz", pcps)
    mask_path = _write_mask(tmpdir / "masked_nt_sites.csv", mask_rows, mask_columns)

    out_path = tmpdir / "clade_pcps_dasm_ready.csv.gz"
    report_path = tmpdir / "clade_dasm_filter_report.txt"
    stats_path = tmpdir / "clade_dasm_filter_stats.json"
    log_path = tmpdir / "curate.log"

    stub = _Namespace(
        input=_Namespace(pcps=str(pcp_path), mask_sites=str(mask_path)),
        output=_Namespace(
            dasm_ready=str(out_path),
            report=str(report_path),
            stats=str(stats_path),
        ),
        params=_Params(seed=seed, val_fraction=val_fraction),
        log=[str(log_path)],
    )
    script_globals = {
        "__name__": "__main__",
        "__file__": str(script),
        "snakemake": stub,
    }
    try:
        exec(compile(script.read_text(), str(script), "exec"), script_globals)
    except SystemExit as exit_signal:
        # The empty-input branch raises SystemExit(0) by design.
        if exit_signal.code not in (0, None):
            raise

    kept = pd.read_csv(out_path)
    names = [] if kept.empty else [
        name.replace("child_", "") for name in kept["child_name"]
    ]
    stats = json.loads(stats_path.read_text())
    return names, stats, log_path.read_text()


# ---------------------------------------------------------------------------
# T1 -- known answer, including the recovery case
# ---------------------------------------------------------------------------

# 9 nt / 3 codons, mask = nt index 3 (the first nucleotide of codon 2). Row A is
# the whole point of #158: its only mutation is at index 5, an UNMASKED
# nucleotide of the MASKED codon 1, so the codon rule drops it and the
# nucleotide rule keeps it.
T1_PARENT = "ATGAAACCC"
T1_PCPS = [
    ("A", T1_PARENT, "ATGAATCCC"),  # differs at 5      -> codon: drop, nt: KEEP
    ("B", T1_PARENT, "ATGCAACCC"),  # differs at 3      -> codon: drop, nt: drop
    ("C", T1_PARENT, "CTGAAACCC"),  # differs at 0      -> codon: keep, nt: keep
    ("D", T1_PARENT, "ATGCATCCC"),  # differs at 3 and 5 -> codon: drop, nt: KEEP
]
T1_MASK = {3}


def check_recovery_case_is_kept():
    """The shipped script must keep {A, C, D}; a codon rule keeps only {C}."""
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        names, stats, log = run_curate(
            Path(tmp) / "t1", T1_PCPS, _mask_rows(T1_MASK)
        )
    if sorted(names) != ["A", "C", "D"]:
        failures.append(
            f"T1: the nucleotide filter kept {sorted(names)}, expected "
            f"['A', 'C', 'D']. A, C and D each differ at an unmasked "
            f"nucleotide; only B's single difference is at masked index 3."
        )

    # The two rules must actually disagree on this fixture, or the check above
    # would pass for a script that never changed.
    codon_kept = sorted(
        name
        for name, parent, child in T1_PCPS
        if reference_keep_codon(parent, child, widen(T1_MASK))
    )
    if codon_kept != ["C"]:
        failures.append(
            f"T1 fixture is not discriminating: the retired codon rule keeps "
            f"{codon_kept}, expected ['C'] -- so this test would pass on a "
            f"reverted script."
        )

    # And the shipped predicate must agree with the plain-Python reference.
    reference = sorted(
        name
        for name, parent, child in T1_PCPS
        if reference_keep_nt(parent, child, T1_MASK)
    )
    if sorted(names) != reference:
        failures.append(
            f"T1: the vectorized predicate kept {sorted(names)} but the "
            f"plain-Python reference keeps {reference}."
        )

    if stats["num_masked_sites"] != 1:
        failures.append(
            f"T1: num_masked_sites is {stats['num_masked_sites']}, expected 1 "
            f"distinct nt_index. It must count distinct indices, not rows."
        )
    if stats["removed_masked_only"] != 1:
        failures.append(
            f"T1: removed_masked_only is {stats['removed_masked_only']}, "
            f"expected 1 (row B)."
        )
    if "nucleotide" not in log:
        failures.append(
            "T1: the log line does not mention nucleotides, so a run's log "
            "still reads as though the mask were codon-granular."
        )
    return failures


def check_multi_reason_site_counts_once():
    """A two-reason site is two rows and one masked position.

    18 of the 121 masked sites on the six v12 production genes carry both a VCF
    and a config reason, so this is the ordinary case rather than an edge one. A
    reader that counted rows would report 139 masked positions where there are
    121, and ``num_masked_sites`` is documented in ``OUTPUT_SCHEMA.md`` as a
    downstream contract.
    """
    failures = []
    rows = [
        (21566, 3, "sites_to_exclude"),
        (21566, 3, "vcf_masked_sites"),
    ]
    with tempfile.TemporaryDirectory() as tmp:
        names, stats, _ = run_curate(Path(tmp) / "dual", T1_PCPS, rows)
    if stats["num_masked_sites"] != 1:
        failures.append(
            f"a site carrying two reasons reported "
            f"{stats['num_masked_sites']} masked positions, expected 1: the "
            f"reader must de-duplicate nt_index."
        )
    if sorted(names) != ["A", "C", "D"]:
        failures.append(
            f"a duplicated mask row changed the kept set to {sorted(names)}, "
            f"expected ['A', 'C', 'D']."
        )
    return failures


# ---------------------------------------------------------------------------
# T2 -- the off-by-one, caught in both directions
# ---------------------------------------------------------------------------

# E differs at index 2, one before the masked index 3. Under the correct
# 0-based reading (mask {3}) E is kept and B dropped. A reader that kept the
# retired codon file's ``- 1`` (mask {2}) inverts BOTH: E dropped, B kept. So
# neither direction of the shift can pass by luck.
T2_PCPS = T1_PCPS + [("E", T1_PARENT, "ATTAAACCC")]


def check_index_is_read_as_zero_based():
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        names, _, _ = run_curate(Path(tmp) / "t2", T2_PCPS, _mask_rows(T1_MASK))
    kept = set(names)
    if "E" not in kept:
        failures.append(
            "T2: E (differs at index 2) was dropped under mask {3}, so the "
            "reader is shifting nt_index down by one -- the retired codon "
            "file's 1-based convention leaking in."
        )
    if "B" in kept:
        failures.append(
            "T2: B (differs at index 3) was kept under mask {3}, so the "
            "reader is shifting nt_index up by one."
        )
    return failures


def check_a_shifted_reader_is_caught():
    """Mutation test: reinstate the ``- 1`` in the module and expect it to move.

    Planted in ``workflow/modules/masking.py`` rather than in the script,
    because that is where ``read_masked_nt_indices`` lives and so where an
    operator would reintroduce the subtraction (patterns 14 and 15).
    """
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        root = _mutated_workflow(
            tmp,
            {
                MASKING_REL: (
                    "return sorted({int(value) for value in "
                    "mask_df[MASKED_NT_INDEX_COLUMN]})",
                    "return sorted({int(value) - 1 for value in "
                    "mask_df[MASKED_NT_INDEX_COLUMN]})",
                )
            },
        )
        names, _, _ = run_curate(
            Path(tmp) / "shifted",
            T2_PCPS,
            _mask_rows(T1_MASK),
            workflow_root=root,
        )
    kept = set(names)
    if "E" in kept or "B" not in kept:
        failures.append(
            f"a reader subtracting 1 from nt_index kept {sorted(kept)}, which "
            f"does not invert E and B -- so this test cannot see the shift it "
            f"exists to catch."
        )
    return failures


# ---------------------------------------------------------------------------
# T3 -- the superset invariant, non-vacuously
# ---------------------------------------------------------------------------


def check_nt_filter_keeps_a_superset():
    """Every PCP the codon rule keeps, the nt rule keeps; strictly for some.

    The warrant for this on real data was ``masked_nt_sites.py``'s guard 2,
    which #158 deletes -- so the codon set here is *constructed* as
    ``widen(mask)`` rather than asserted independently. On the v12 test tree
    the containment is not strict (zero PCPs are recovered), which is why this
    needs a fixture built to make it strict.
    """
    failures = []
    rng = random.Random(20260903)
    length = 30
    # No T, so no random draw can produce TAA/TAG/TGA. Filter 5 drops a PCP
    # with an in-frame stop codon, and those drops are invisible to the codon
    # reference below -- an earlier version of this fixture used all four bases
    # and reported 11 stop-codon drops as superset violations.
    bases = "ACG"
    parent = "".join(rng.choice(bases) for _ in range(length))
    masked = {index for index in range(length) if rng.random() < 0.25}
    masked_codons = widen(masked)

    pcps = []
    for row in range(200):
        child = list(parent)
        for _ in range(rng.randint(1, 3)):
            position = rng.randrange(length)
            child[position] = rng.choice([b for b in bases if b != parent[position]])
        child = "".join(child)
        if child != parent:
            pcps.append((str(row), parent, child))

    with tempfile.TemporaryDirectory() as tmp:
        names, _, _ = run_curate(Path(tmp) / "t3", pcps, _mask_rows(masked))
    kept_nt = set(names)
    kept_codon = {
        name
        for name, p, c in pcps
        if reference_keep_codon(p, c, masked_codons)
    }

    missing = sorted(kept_codon - kept_nt)
    if missing:
        failures.append(
            f"the nucleotide filter dropped {len(missing)} PCPs the codon rule "
            f"keeps ({missing[:5]}), so it is not a superset. Every masked "
            f"nucleotide lies in a masked codon by construction, so an "
            f"unmasked codon can contain no masked nucleotide."
        )
    if not (kept_nt - kept_codon):
        failures.append(
            "no PCP was recovered on a fixture built to recover some, so the "
            "containment is vacuous here and the check proves nothing."
        )
    return failures


# ---------------------------------------------------------------------------
# T4 -- the frameshift junction, unreachable from any pipeline run
# ---------------------------------------------------------------------------


def check_frameshift_junction_masks_both_indices():
    """A one-to-many site is two rows, and both must mask.

    Genomic 13468 occupies ORF1b's nt indices 26 and 27, in different codons.
    The indices are taken from the real ``sequences.genomic_to_pcp_indices``
    rather than written down, so the fixture cannot claim a junction the
    coordinate code does not have. The filter reads pre-expanded rows and never
    re-derives an offset; this pins that a consumer of the file gets both.
    """
    failures = []
    junction = sequences.genomic_to_pcp_indices(ORF1B_BLOCKS, 13468)
    if junction != [26, 27]:
        return [
            f"sequences.genomic_to_pcp_indices puts ORF1b's genomic 13468 at "
            f"{junction}, not [26, 27]; the fixture below is no longer about "
            f"the frameshift junction."
        ]
    if junction[0] // 3 == junction[1] // 3:
        failures.append(
            f"ORF1b's junction indices {junction} fall in the same codon, so "
            f"the case this test exists for is gone."
        )

    length = 30
    parent = "ATG" * (length // 3)
    pcps = [(f"only_{i}", parent, _mutate_at(parent, i)) for i in (25, 26, 27)]
    for name, p, c in pcps:
        assert sum(a != b for a, b in zip(p, c)) == 1, name

    rows = [(13468, index, "vcf_masked_sites") for index in junction]
    with tempfile.TemporaryDirectory() as tmp:
        names, stats, _ = run_curate(Path(tmp) / "t4", pcps, rows)
    kept = set(names)
    if kept != {"only_25"}:
        failures.append(
            f"the junction mask kept {sorted(kept)}, expected ['only_25']: "
            f"both 26 and 27 are masked, and 25 is not."
        )
    if stats["num_masked_sites"] != 2:
        failures.append(
            f"the junction site reported {stats['num_masked_sites']} masked "
            f"positions, expected 2: one genomic site, two indices."
        )
    return failures


def check_deduplicating_by_site_is_caught():
    """Mutation test: a reader keying on nt_site loses one junction index."""
    failures = []
    rows = [(13468, index, "vcf_masked_sites") for index in (26, 27)]
    length = 30
    parent = "ATG" * (length // 3)
    pcps = [("only_27", parent, _mutate_at(parent, 27))]
    with tempfile.TemporaryDirectory() as tmp:
        root = _mutated_workflow(
            tmp,
            {
                MASKING_REL: (
                    "return sorted({int(value) for value in "
                    "mask_df[MASKED_NT_INDEX_COLUMN]})",
                    "return sorted({int(value) for value in "
                    'mask_df.drop_duplicates("nt_site")'
                    "[MASKED_NT_INDEX_COLUMN]})",
                )
            },
        )
        names, _, _ = run_curate(
            Path(tmp) / "dedup", pcps, rows, workflow_root=root
        )
    if names != ["only_27"]:
        failures.append(
            "a reader de-duplicating the mask by nt_site did not keep the PCP "
            "mutating only index 27, so this test cannot see a collapsed "
            "one-to-many site."
        )
    return failures


# ---------------------------------------------------------------------------
# The reader's guards
# ---------------------------------------------------------------------------


def check_missing_column_raises_attributably():
    """A mask that cannot be parsed must raise, not read as an empty mask.

    An unparseable mask and an empty one have the same observable outcome --
    every masked position trained on as live -- so the reader names the path
    and refuses rather than falling back.
    """
    failures = []
    rows = [(21566, 3, "sites_to_exclude")]
    with tempfile.TemporaryDirectory() as tmp:
        try:
            run_curate(
                Path(tmp) / "bad",
                T1_PCPS,
                rows,
                mask_columns=("nt_site", "codon_site", "reason"),
            )
        except ValueError as error:
            message = str(error)
            if "nt_index" not in message or "masked_nt_sites.csv" not in message:
                failures.append(
                    f"the missing-column error does not name both the column "
                    f"and the file: {message!r}"
                )
        except Exception:
            failures.append(
                f"a mask without an nt_index column raised the wrong error "
                f"type:\n{traceback.format_exc()}"
            )
        else:
            failures.append(
                "a mask without an nt_index column did NOT raise; the filter "
                "silently applied no mask."
            )
    return failures


def check_a_silent_fallback_is_caught():
    """Mutation test: a reader falling back to an empty mask keeps everything.

    The mutation returns ``[]`` rather than deleting the ``if`` outright,
    because deleting it yields pandas' bare ``KeyError`` -- a loud failure,
    which is not the outcome the guard prevents. What it prevents is a *quiet*
    unmasked run, so that is what is planted.
    """
    failures = []
    rows = [(21566, 3, "sites_to_exclude")]
    with tempfile.TemporaryDirectory() as tmp:
        root = _mutated_workflow(
            tmp,
            {
                MASKING_REL: (
                    "    if MASKED_NT_INDEX_COLUMN not in mask_df.columns:\n"
                    "        raise ValueError(",
                    "    if MASKED_NT_INDEX_COLUMN not in mask_df.columns:\n"
                    "        return []\n"
                    "    if False:\n"
                    "        raise ValueError(",
                )
            },
        )
        try:
            names, stats, _ = run_curate(
                Path(tmp) / "fallback",
                T1_PCPS,
                rows,
                mask_columns=("nt_site", "codon_site", "reason"),
                workflow_root=root,
            )
        except Exception:
            failures.append(
                "with the guard disabled the run still failed, so the guard "
                "is not what produces the attributable error above:\n"
                + traceback.format_exc()
            )
        else:
            if stats["num_masked_sites"] != 0 or sorted(names) != [
                "A",
                "B",
                "C",
                "D",
            ]:
                failures.append(
                    f"with the guard disabled the run kept {sorted(names)} at "
                    f"{stats['num_masked_sites']} masked positions; expected "
                    f"an unmasked run keeping everything, which is the outcome "
                    f"the guard prevents."
                )
    return failures


def check_empty_mask_keeps_every_valid_pcp():
    """A header-only mask is legal and must skip the filter, not empty the file.

    ``masked_nt_sites.py`` declares the header-only case legitimate for a gene
    with no globally masked site, so the ``num_masked_sites > 0`` gate in the
    script is load-bearing rather than an optimization.
    """
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        names, stats, _ = run_curate(Path(tmp) / "empty", T1_PCPS, [])
    if sorted(names) != ["A", "B", "C", "D"]:
        failures.append(
            f"a header-only mask kept {sorted(names)}, expected every valid "
            f"PCP."
        )
    if stats["num_masked_sites"] != 0 or stats["removed_masked_only"] != 0:
        failures.append(
            f"a header-only mask reported {stats['num_masked_sites']} masked "
            f"positions and {stats['removed_masked_only']} drops, expected 0 "
            f"and 0."
        )
    return failures


def check_train_val_split_is_untouched():
    """Filters 1-5 and the RNG are unchanged, so the split is reproducible.

    The equivalence to ``dasm2 preprocess`` now covers filters 1-5 plus the
    split rather than all six (#158), and this pins the part that still holds:
    the same row set and seed must give the same assignment, and the two seeds
    the pipeline uses must disagree -- the guardrail against producing the
    aggregate with the per-clade seed.
    """
    failures = []
    rng_reference = __import__("numpy").random.default_rng(seed=1729)
    expected = list(~(rng_reference.random(3) < 0.10))

    with tempfile.TemporaryDirectory() as tmp:
        _write_pcps(Path(tmp) / "x.csv.gz", T1_PCPS)  # keep tmp non-empty
        names, _, _ = run_curate(Path(tmp) / "split", T1_PCPS, _mask_rows(T1_MASK))
        kept = pd.read_csv(Path(tmp) / "split" / "clade_pcps_dasm_ready.csv.gz")
        got = list(kept["in_train"])
        if got != expected:
            failures.append(
                f"in_train is {got} for the 3 kept rows but the documented RNG "
                f"call sequence gives {expected}."
            )
        names_42, _, _ = run_curate(
            Path(tmp) / "split42", T1_PCPS, _mask_rows(T1_MASK), seed=42
        )
        kept_42 = pd.read_csv(Path(tmp) / "split42" / "clade_pcps_dasm_ready.csv.gz")
        if sorted(names) != sorted(names_42):
            failures.append(
                f"changing the seed changed the kept set from {sorted(names)} "
                f"to {sorted(names_42)}; the seed must affect only the split."
            )
        if list(kept_42["in_train"]) == got:
            failures.append(
                "seeds 1729 and 42 gave the same in_train assignment, so the "
                "aggregate/per-clade guardrail is not actually distinguishing "
                "them on this fixture."
            )
    return failures


CHECKS = [
    ("T1 recovery case is kept", check_recovery_case_is_kept),
    ("T1 a two-reason site counts once", check_multi_reason_site_counts_once),
    ("T2 nt_index is read 0-based", check_index_is_read_as_zero_based),
    ("T2 mutation test: a shifted reader is caught", check_a_shifted_reader_is_caught),
    ("T3 nt filter keeps a superset", check_nt_filter_keeps_a_superset),
    ("T4 junction masks both indices", check_frameshift_junction_masks_both_indices),
    (
        "T4 mutation test: de-duplicating by site is caught",
        check_deduplicating_by_site_is_caught,
    ),
    ("guard: missing nt_index raises", check_missing_column_raises_attributably),
    ("guard mutation test: a silent fallback is caught", check_a_silent_fallback_is_caught),
    ("empty mask keeps every valid PCP", check_empty_mask_keeps_every_valid_pcp),
    ("train/val split is untouched", check_train_val_split_is_untouched),
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


def test_dasm_ready_nt_filter():
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
        for failure in failures:
            print(f"  - {failure}")
        sys.exit(1)
    print(
        "\nPASSED: the DASM-ready filter keeps a PCP whose only mutation is at "
        "an unmasked nucleotide of a masked codon, reads nt_index 0-based in "
        "both directions, masks both indices of the frameshift junction, and "
        "refuses an unparseable mask rather than applying none."
    )


if __name__ == "__main__":
    main()
