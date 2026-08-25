#!/usr/bin/env python3
"""Unit tests for the date-cutoff PCP subtrees (#149).

Runs the real ``workflow/scripts/samples_upto_cutoff.py`` and
``clade_counts_by_cutoff.py`` against a stub ``snakemake`` global, and the
mutation tests exec *deliberately corrupted* copies of that shipped text --
never a re-implementation. CLAUDE.md pattern 10 records why: a test that
re-implements a script's logic and asserts against its own copy stays green
while three separate breaks of the shipped script ship.

**What a pipeline run cannot check, and this file therefore must.** The CI arms
do cover the two things a run can see: ``config_test_data.yaml``'s ``2020-06-01``
arm drops 19A below ``min_clade_size``, and its ``2022-06-01`` arm is degenerate
and byte-compared against the full tree by
``validate_test.py::validate_cutoff_pcps``. What no run reaches is the
*boundary*: every test-data sample date is interior to both arms, so a cutoff
comparison written ``<`` instead of ``<=`` would drop nothing observable, and the
counts, the roster and the byte-equality would all still agree. Inclusivity is
checked here on a synthetic fixture built so that off-by-one moves a sample.

Nor can a run reach the *ordering* guarantee. The degenerate arm's byte-equality
depends on the kept sample list preserving the input file's order, and the test
data's clade lists happen to be sorted, so a script that sorted its output would
pass CI. The fixture here is deliberately unsorted.

Hermetic: synthetic metadata and synthetic sample lists. Nothing reads a
pipeline output.
"""

import os
import shutil
import sys
import tempfile
import traceback
from importlib import util as importlib_util
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SAMPLES_SCRIPT = REPO_ROOT / "workflow" / "scripts" / "samples_upto_cutoff.py"
COUNTS_SCRIPT = REPO_ROOT / "workflow" / "scripts" / "clade_counts_by_cutoff.py"


def _load(module_name):
    """Load a workflow module by path, as the pipeline scripts do."""
    spec = importlib_util.spec_from_file_location(
        module_name, REPO_ROOT / "workflow" / "modules" / f"{module_name}.py"
    )
    module = importlib_util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


cutoffs = _load("cutoffs")


class _Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def get(self, key, default=None):
        return self.__dict__.get(key, default)

    def __getitem__(self, key):
        return self.__dict__[key]


# ---------------------------------------------------------------------------
# Fixture
# ---------------------------------------------------------------------------
#
# ``boundary`` sits exactly on the 2021-01-01 cutoff, so an exclusive comparison
# drops it and an inclusive one keeps it -- the distinction no pipeline run can
# see. ``early``/``late`` bracket it. The clade list order is deliberately NOT
# sorted (``d_early`` before ``a_boundary``), so a script that sorted its output
# would break the degenerate arm's byte-equality and is caught here.
CUTOFF = "2021-01-01"
CLADE_ROWS = {
    "20A": [
        ("d_early", "2020-06-15"),
        ("a_boundary", "2021-01-01"),
        ("c_late", "2021-06-01"),
        ("b_early", "2020-12-31"),
    ],
    "20B": [
        ("e_late", "2021-03-01"),
        ("f_late", "2022-01-01"),
    ],
    "19A": [
        ("g_early", "2020-02-01"),
        ("h_early", "2020-03-01"),
        ("i_late", "2021-09-09"),
    ],
}
# Strains at or before 2021-01-01, in each clade's own file order.
EXPECTED_UPTO = {
    "20A": ["d_early", "a_boundary", "b_early"],
    "20B": [],
    "19A": ["g_early", "h_early"],
}


def _write_fixture(root):
    """Write a synthetic filtered metadata TSV and per-clade sample lists."""
    root = Path(root)
    meta = root / "filtered_metadata.tsv.gz"
    import pandas as pd

    rows = [
        {"strain": strain, "date": date, "host": "Human", "length": 29903}
        for clade_rows in CLADE_ROWS.values()
        for strain, date in clade_rows
    ]
    pd.DataFrame(rows).to_csv(meta, sep="\t", index=False, compression="gzip")

    sample_files = {}
    clades_dir = root / "clades"
    clades_dir.mkdir(exist_ok=True)
    for clade, clade_rows in CLADE_ROWS.items():
        path = clades_dir / f"{clade}_samples.txt"
        path.write_text("".join(f"{strain}\n" for strain, _ in clade_rows))
        sample_files[clade] = str(path)
    return str(meta), sample_files


def _exec_script(path, script_globals, mutate=None):
    """Exec a shipped script, optionally against a corrupted copy of its text."""
    text = path.read_text()
    if mutate is not None:
        mutated = mutate(text)
        assert mutated != text, f"mutation was a no-op on {path.name}"
        text = mutated
    script_globals.setdefault("__name__", "__main__")
    script_globals.setdefault("__file__", str(path))
    exec(compile(text, str(path), "exec"), script_globals)


def _mutated_workflow(tmp, edits):
    """A throwaway copy of ``workflow/`` with text edits applied.

    Copied rather than edited in place. The scripts resolve
    ``modules/<name>.py`` relative to their own location, so a mutation to a
    *module* can only be reached by relocating the pair -- and mutating the real
    file, even with a ``finally`` that restores it, would corrupt the repo for
    any concurrent run and leave it broken if the process were killed.

    Returns the copied ``scripts/`` directory.
    """
    root = Path(tmp) / "wf"
    shutil.copytree(REPO_ROOT / "workflow" / "modules", root / "modules")
    shutil.copytree(REPO_ROOT / "workflow" / "scripts", root / "scripts")
    for relative, (old, new) in edits.items():
        path = root / relative
        original = path.read_text()
        assert original.count(old) == 1, (
            f"could not locate the text to mutate in {relative}: {old[:60]!r}"
        )
        path.write_text(original.replace(old, new))
    return root / "scripts"


def run_samples_upto(meta, samples_in, out_path, cutoff, clade, mutate=None,
                     scripts_dir=None):
    snake = _Namespace(
        input=_Namespace(samples=samples_in, metadata=meta),
        output=_Namespace(samples=str(out_path)),
        wildcards=_Namespace(cutoff=cutoff, clade=clade),
        log=[],
    )
    script = SAMPLES_SCRIPT if scripts_dir is None else (
        Path(scripts_dir) / "samples_upto_cutoff.py"
    )
    _exec_script(script, {"snakemake": snake}, mutate)
    return Path(out_path).read_text().split()


def run_counts(meta, sample_files, out_path, cutoff_list, min_clade_size, mutate=None):
    snake = _Namespace(
        input=_Namespace(metadata=meta, samples=sorted(sample_files.values())),
        output=_Namespace(counts=str(out_path)),
        params=_Namespace(
            cutoffs=cutoff_list,
            min_clade_size=min_clade_size,
            clade_column="Nextstrain_clade",
        ),
        log=[],
    )
    _exec_script(COUNTS_SCRIPT, {"snakemake": snake}, mutate)
    import pandas as pd

    return pd.read_csv(out_path)


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------


def check_cutoff_is_inclusive_and_order_preserving():
    """A sample collected *on* the cutoff is kept, and file order survives.

    Both properties are invisible to every pipeline run (see the module
    docstring), and both are load-bearing: inclusivity is what "at or before"
    means, and order is what makes the degenerate arm byte-identical to the full
    tree.
    """
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        meta, sample_files = _write_fixture(tmp)
        for clade, expected in EXPECTED_UPTO.items():
            if not expected:
                continue  # empty clades raise; covered separately
            got = run_samples_upto(
                meta, sample_files[clade], Path(tmp) / f"{clade}.txt", CUTOFF, clade
            )
            if got != expected:
                failures.append(
                    f"{clade} at {CUTOFF}: expected {expected}, got {got}"
                )
        # State the boundary claim directly, so a future reader does not have to
        # decode it from the fixture.
        got_20a = run_samples_upto(
            meta, sample_files["20A"], Path(tmp) / "b.txt", CUTOFF, "20A"
        )
        if "a_boundary" not in got_20a:
            failures.append(
                "a sample collected exactly on the cutoff was dropped; 'at or "
                "before' must include the cutoff date itself"
            )
        # Guard the fixture itself: if the expected order ever became sorted,
        # check_sorting_the_output_is_caught would silently stop testing anything.
        if sorted(got_20a) == got_20a:
            failures.append(
                "the kept list came out sorted, but the fixture's input order is "
                "not sorted -- so this check can no longer detect a script that "
                "sorts its output. Fix the fixture."
            )
    return failures


def check_exclusive_comparison_is_caught():
    """Mutation test: ``<`` instead of ``<=`` drops the boundary sample.

    The comparison lives in ``modules/cutoffs.py``, not in the script, so the
    mutation is applied to a *copied* workflow tree and the copied script is
    run against it. Mutating only the script text would report this guard as
    covered while the operator that decides inclusivity sat one file away
    untested.
    """
    with tempfile.TemporaryDirectory() as tmp:
        meta, sample_files = _write_fixture(tmp)
        good = run_samples_upto(
            meta, sample_files["20A"], Path(tmp) / "good.txt", CUTOFF, "20A"
        )
        scripts = _mutated_workflow(
            tmp,
            {
                "modules/cutoffs.py": (
                    "return [s for s in strains if strain_dates[s] <= bound]",
                    "return [s for s in strains if strain_dates[s] < bound]",
                )
            },
        )
        broken = run_samples_upto(
            meta, sample_files["20A"], Path(tmp) / "bad.txt", CUTOFF, "20A",
            scripts_dir=scripts,
        )
    if "a_boundary" in broken:
        return [
            "an exclusive '<' comparison still kept the boundary sample, so this "
            "test cannot detect an off-by-one at the cutoff"
        ]
    if len(broken) != len(good) - 1:
        return [
            f"the exclusive mutant dropped {len(good) - len(broken)} sample(s), "
            f"expected exactly 1 (the boundary sample)"
        ]
    return []


def check_sorting_the_output_is_caught():
    """Mutation test: a script that sorts its kept list is caught.

    Order is what makes the degenerate arm byte-identical to the full tree, and
    the test data's clade lists are already sorted, so no pipeline run can see
    this.
    """
    with tempfile.TemporaryDirectory() as tmp:
        meta, sample_files = _write_fixture(tmp)
        scripts = _mutated_workflow(
            tmp,
            {
                "modules/cutoffs.py": (
                    "return [s for s in strains if strain_dates[s] <= bound]",
                    "return sorted(s for s in strains if strain_dates[s] <= bound)",
                )
            },
        )
        got = run_samples_upto(
            meta, sample_files["20A"], Path(tmp) / "s.txt", CUTOFF, "20A",
            scripts_dir=scripts,
        )
    if got == EXPECTED_UPTO["20A"]:
        return [
            "sorting the kept list produced the expected order, so this fixture "
            "cannot detect a reordering; the clade sample list must not be sorted"
        ]
    return []


def check_empty_clade_raises_rather_than_writing_nothing():
    """A clade with no pre-cutoff samples must raise, not emit an empty file.

    An empty ``-s`` list makes ``matUtils extract`` fail with an error naming
    neither the clade nor the cutoff, and an empty PCP table downstream would
    read as a true zero. The DAG is supposed to have skipped this clade, so
    reaching the rule at all is a disagreement worth stopping on.
    """
    with tempfile.TemporaryDirectory() as tmp:
        meta, sample_files = _write_fixture(tmp)
        try:
            run_samples_upto(
                meta, sample_files["20B"], Path(tmp) / "e.txt", CUTOFF, "20B"
            )
        except ValueError as exc:
            if "no samples at or before" not in str(exc):
                return [f"raised, but not about the empty selection: {exc}"]
            return []
    return ["a clade with zero pre-cutoff samples wrote a file instead of raising"]


def check_counts_match_the_sample_lists():
    """The checkpoint's counts equal what the per-clade rule actually writes.

    These are two different code paths over the same module function, and the
    roster is only trustworthy if they agree: ``validate_test.py`` checks the
    same equality end-to-end, but only for clades a run happens to include.
    """
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        meta, sample_files = _write_fixture(tmp)
        df = run_counts(
            meta, sample_files, Path(tmp) / "counts.csv", [CUTOFF], min_clade_size=1
        )
        for clade, expected in EXPECTED_UPTO.items():
            row = df[df["clade"] == clade]
            if len(row) != 1:
                failures.append(f"{clade}: expected 1 count row, got {len(row)}")
                continue
            got = int(row["n_samples_upto"].iloc[0])
            if got != len(expected):
                failures.append(
                    f"{clade}: checkpoint counted {got}, the sample list has "
                    f"{len(expected)}"
                )
            if int(row["n_samples_full"].iloc[0]) != len(CLADE_ROWS[clade]):
                failures.append(f"{clade}: n_samples_full disagrees with the fixture")
    return failures


def check_min_clade_size_decides_adequacy():
    """``min_clade_size`` is applied at the cutoff, not on the full tree.

    19A has 3 samples in all but only 2 by the cutoff, so a threshold of 3 must
    mark it inadequate. Testing against the full-tree count instead would keep it
    -- and the clade would then reach ``matUtils`` with fewer samples than the
    pipeline's own definition of a processable clade.
    """
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        meta, sample_files = _write_fixture(tmp)
        df = run_counts(
            meta, sample_files, Path(tmp) / "c3.csv", [CUTOFF], min_clade_size=3
        )
        verdicts = dict(zip(df["clade"], df["adequate_sample_counts"]))
        if verdicts.get("19A") is not False:
            failures.append(
                "19A has 3 samples in total but 2 by the cutoff; at "
                "min_clade_size=3 it must be inadequate (adequacy was read off "
                "the full-tree count)"
            )
        if verdicts.get("20A") is not True:
            failures.append("20A has 3 samples by the cutoff and must be adequate")
        if verdicts.get("20B") is not False:
            failures.append("20B has 0 samples by the cutoff and must be inadequate")
    return failures


def check_adequacy_on_the_full_count_is_caught():
    """Mutation test: scoring adequacy against ``n_samples_full``."""
    original = COUNTS_SCRIPT.read_text()
    mutated = original.replace(
        '"adequate_sample_counts": len(kept) >= min_clade_size,',
        '"adequate_sample_counts": len(strains) >= min_clade_size,',
    )
    assert mutated != original, "could not locate the adequacy test to mutate"
    with tempfile.TemporaryDirectory() as tmp:
        meta, sample_files = _write_fixture(tmp)
        df = run_counts(
            meta,
            sample_files,
            Path(tmp) / "m.csv",
            [CUTOFF],
            min_clade_size=3,
            mutate=lambda t: mutated,
        )
    verdicts = dict(zip(df["clade"], df["adequate_sample_counts"]))
    if verdicts.get("19A") is not True:
        return [
            "the mutant that scores adequacy on the full-tree count did not "
            "change 19A's verdict, so this check cannot detect it"
        ]
    return []


def check_all_inadequate_raises():
    """A cutoff at which nothing is processable is a config error, not a result.

    Without this, every downstream file is simply absent, which is
    indistinguishable from a build that has not run.
    """
    with tempfile.TemporaryDirectory() as tmp:
        meta, sample_files = _write_fixture(tmp)
        try:
            run_counts(
                meta, sample_files, Path(tmp) / "x.csv", [CUTOFF], min_clade_size=99
            )
        except ValueError as exc:
            if "no clade has" not in str(exc):
                return [f"raised, but not about the empty cutoff: {exc}"]
            return []
    return ["a cutoff with no adequate clade produced a table instead of raising"]


def check_counts_file_is_sorted_and_reproducible():
    """The checkpoint's output does not depend on input order.

    Snakemake makes no promise about the order of an expanded input list, and
    ``data-central`` manifests pin a hash per file (pattern 9), so an unsorted
    write would be a spurious diff on every rebuild.
    """
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        meta, sample_files = _write_fixture(tmp)
        first = run_counts(
            meta, sample_files, Path(tmp) / "a.csv", [CUTOFF, "2022-06-01"],
            min_clade_size=1,
        )
        shuffled = dict(reversed(list(sample_files.items())))
        second = run_counts(
            meta, shuffled, Path(tmp) / "b.csv", [CUTOFF, "2022-06-01"],
            min_clade_size=1,
        )
        if (Path(tmp) / "a.csv").read_bytes() != (Path(tmp) / "b.csv").read_bytes():
            failures.append(
                "the counts file changed when the input order changed; it must be "
                "sorted before writing"
            )
        keys = list(zip(first["cutoff"], first["clade"]))
        if keys != sorted(keys):
            failures.append(f"counts rows are not sorted by (cutoff, clade): {keys}")
    return failures


def check_unknown_strain_raises():
    """A sample list naming a strain the metadata does not have must raise.

    Both files derive from the same filtered metadata, so a miss means they do
    not -- most plausibly the raw metadata was wired in, in which case silently
    dropping the strain shrinks every cutoff set by an unreported amount.
    """
    with tempfile.TemporaryDirectory() as tmp:
        meta, sample_files = _write_fixture(tmp)
        rogue = Path(tmp) / "rogue_samples.txt"
        rogue.write_text("d_early\nnot_in_metadata\n")
        try:
            run_samples_upto(meta, str(rogue), Path(tmp) / "r.txt", CUTOFF, "20A")
        except ValueError as exc:
            if "absent from" not in str(exc):
                return [f"raised, but not about the unknown strain: {exc}"]
            return []
    return ["an unknown strain was silently dropped instead of raising"]


def check_config_validation():
    """``validate_cutoffs`` rejects the four ways a cutoff can be wrong.

    Loose ISO, at/before ``min_leaf_date``, past ``max_leaf_date``, and a
    duplicate. Each would otherwise surface as a ``MissingInputException`` naming
    a path nothing asked for, or -- worst -- as a silent full-tree copy under a
    directory name claiming a truncation.
    """
    failures = []
    ok = cutoffs.validate_cutoffs(
        ["2020-06-01", "2022-06-01"], "2019-12-01", "2022-06-01", "t"
    )
    if ok != ["2020-06-01", "2022-06-01"]:
        failures.append(f"valid cutoffs were altered: {ok}")
    cases = {
        "loose ISO": ["2023-1-1"],
        "at or before min_leaf_date": ["2019-12-01"],
        "past max_leaf_date": ["2022-06-02"],
        "duplicate": ["2020-06-01", "2020-06-01"],
    }
    for label, value in cases.items():
        try:
            cutoffs.validate_cutoffs(value, "2019-12-01", "2022-06-01", "t")
        except ValueError:
            continue
        failures.append(f"validate_cutoffs accepted {label}: {value}")
    # A None max_leaf_date means no upper bound, and must not crash.
    try:
        cutoffs.validate_cutoffs(["2030-01-01"], "2019-12-01", None, "t")
    except ValueError as exc:
        failures.append(f"an unbounded tree rejected a late cutoff: {exc}")
    return failures


CHECKS = [
    ("cutoff is inclusive and order-preserving",
     check_cutoff_is_inclusive_and_order_preserving),
    ("mutation test: exclusive '<' is caught", check_exclusive_comparison_is_caught),
    ("mutation test: sorting the output is caught", check_sorting_the_output_is_caught),
    ("empty clade raises", check_empty_clade_raises_rather_than_writing_nothing),
    ("checkpoint counts match the sample lists", check_counts_match_the_sample_lists),
    ("min_clade_size applies at the cutoff", check_min_clade_size_decides_adequacy),
    ("mutation test: adequacy on the full count is caught",
     check_adequacy_on_the_full_count_is_caught),
    ("an all-inadequate cutoff raises", check_all_inadequate_raises),
    ("counts file is sorted and reproducible",
     check_counts_file_is_sorted_and_reproducible),
    ("an unknown strain raises", check_unknown_strain_raises),
    ("config validation rejects four bad forms", check_config_validation),
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


def test_pcp_date_cutoffs():
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
        "\nPASSED: the cutoff is inclusive, order survives, min_clade_size is "
        "applied at the cutoff rather than on the full tree, empty selections "
        "raise, and every guard fires on a deliberately broken input."
    )


if __name__ == "__main__":
    main()
