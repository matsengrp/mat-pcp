#!/usr/bin/env python3
"""Pin the bte floor and the get_haplotype re-mutation contract (#136).

``bte<=0.9.2``'s ``MATree.get_haplotype()`` is wrong at a **twice-mutated
position**: it keeps the *first* mutation and discards later non-reversions,
because the update branch mutates a value copy that is never re-inserted
(``allm.erase`` is only reached on the reversion path). So a site that goes
A->C on one branch and C->G on a descendant reports ``A100C`` instead of
``A100G`` -- a plausible sequence, silently wrong at one base.

Two shipped call sites depend on it:

* ``scripts/pcp_from_tree.py:206`` builds every PCP's *parent* haplotype;
* ``modules/windowing.py:397`` builds every window founder's mutation set.

So the blast radius is the PCP corpus DASM2 trains on and the windowed fitness
trajectory. ``gisaid-and-public_2026-04-22_v9`` is **not** affected -- its resolved
env is ``bte-0.9.3-py310h184ae93_1`` on python 3.10.19 -- but that was luck of the
solve: ``bte.yml`` pinned ``bte>=0.1`` with ``python>=3.9``, which left the solver
free to land on 0.9.0. A downstream consumer did land on it
(``dasm2-experiments-2/experiments/2026-07-21-dasm-tree-score``): 898 of 1,332 nodes
carried at least one wrong nucleotide across 169 positions, 50 of them in the spike
CDS.

This file has two halves, and the split is deliberate rather than a convenience.

**Half A -- the pin (always runs).** A pure YAML check on ``workflow/envs/bte.yml``:
the bte floor is >= 0.9.3, and the python floor is >= 3.10 with a < 3.11 ceiling.
The python floor is not an independent preference -- bioconda's bte recipe is
``skip: True  # [py < 310]`` at 0.9.3, so a py3.9 solve *cannot* satisfy the bte
floor and would resolve backwards. The two pins are therefore checked as a
coupled pair: a future edit that relaxes python to >=3.9 while leaving
``bte>=0.9.3`` in place looks harmless and reopens the hole.

**Half B -- the behaviour (runs where bte exists).** Builds a three-deep MAT with
a twice-mutated site, a reverted site, and a singly-mutated site, and asserts
``get_haplotype`` agrees with an explicit root->node walk.

Half B cannot run in the pixi environment: ``bte`` lives only in the Snakemake
conda env ``workflow/envs/bte.yml`` creates, and adding it to ``pixi.toml`` would
re-solve the shared env for a test. Rather than skip -- which is exactly the
silent-pass this suite exists to fight (``CLAUDE.md`` patterns 8 and 10) -- the
test **locates the bte env and re-executes itself inside it**. On a fresh CI runner
that env does not exist yet at unit-test time, so:

* plain invocation: half A runs, half B reports NOT RUN and says how to run it;
* ``--require-bte``: a missing bte env is a **failure**.

CI wires both -- the plain form with the other unit tests, and ``--require-bte``
after the pipeline run, when ``.snakemake/conda`` is populated. Without the second
step half B would never execute anywhere, which is the vacuity failure mode.

Usage:
    python tests/test_bte_haplotype.py                # half A + half B if available
    python tests/test_bte_haplotype.py --require-bte  # half B is mandatory
"""

import os
import re
import sys
import traceback
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
BTE_ENV_YAML = REPO_ROOT / "workflow" / "envs" / "bte.yml"
SNAKEMAKE_CONDA = REPO_ROOT / ".snakemake" / "conda"

# Guard against an infinite re-exec if the located interpreter still cannot
# import bte (a half-built conda env).
REEXEC_FLAG = "MATPCP_BTE_TEST_REEXECED"

MIN_BTE = (0, 9, 3)
MIN_PYTHON = (3, 10)
MAX_PYTHON_EXCLUSIVE = (3, 11)

# The fixture. A path root -> n1 -> n2 -> {leaf1, leaf2}; bte renames internal
# nodes to node_1/node_2/node_3 on load, which is why the ids below are bte's and
# not the newick's -- discovered by inspection, and asserted before use.
NEWICK = "(((leaf1:1,leaf2:1)n2:1)n1:1)root;"
MUTATIONS = {
    # site 100 is the re-mutation: A->C here, C->G on the child. A pre-0.9.3
    # get_haplotype keeps A100C at node_3 and below.
    # site 200 is reverted on the child, which exercises the `allm.erase` path
    # that pre-0.9.3 handled correctly -- included so a failure distinguishes
    # "re-mutations are broken" from "haplotypes are broken generally".
    # site 300 is an ordinary single mutation.
    "node_2": ["1:A100C", "1:A200T", "1:A300G"],
    "node_3": ["1:C100G", "1:T200A"],
}
EXPECTED_HAPLOTYPES = {
    "node_2": {"A100C", "A200T", "A300G"},
    "node_3": {"A100G", "A300G"},
    "leaf1": {"A100G", "A300G"},
    "leaf2": {"A100G", "A300G"},
}

_MUT = re.compile(r"^([ACGT-])(\d+)([ACGT-])$")


def _parse_floor(spec, package):
    """Extract the ``>=`` floor for ``package`` from a dependency string."""
    match = re.search(rf"{package}\s*>=\s*([0-9.]+)", spec)
    if not match:
        return None
    return tuple(int(p) for p in match.group(1).rstrip(".").split("."))


def check_env_pins():
    """Half A: the bte floor and the python floor that makes it enforceable."""
    failures = []
    text = BTE_ENV_YAML.read_text()

    bte_lines = [
        ln.strip().lstrip("- ").strip()
        for ln in text.splitlines()
        if re.match(r"\s*-\s*bte\b", ln)
    ]
    if len(bte_lines) != 1:
        failures.append(
            f"expected exactly one bte dependency line in {BTE_ENV_YAML.name}, "
            f"found {len(bte_lines)}: {bte_lines}"
        )
    else:
        floor = _parse_floor(bte_lines[0], "bte")
        if floor is None:
            failures.append(
                f"bte dependency {bte_lines[0]!r} declares no >= floor, so the "
                f"solver may resolve bte 0.9.0, whose get_haplotype silently "
                f"drops re-mutations (#136)"
            )
        elif floor < MIN_BTE:
            failures.append(
                f"bte floor is {'.'.join(map(str, floor))}, must be at least "
                f"{'.'.join(map(str, MIN_BTE))}: get_haplotype drops re-mutations "
                f"at or below 0.9.2 (#136)"
            )

    python_lines = [
        ln.strip().lstrip("- ").strip()
        for ln in text.splitlines()
        if re.match(r"\s*-\s*python\b", ln)
    ]
    if len(python_lines) != 1:
        failures.append(
            f"expected exactly one python dependency line, found "
            f"{len(python_lines)}: {python_lines}"
        )
    else:
        spec = python_lines[0]
        py_floor = _parse_floor(spec, "python")
        if py_floor is None or py_floor < MIN_PYTHON:
            shown = ".".join(map(str, py_floor)) if py_floor else "none"
            failures.append(
                f"python floor is {shown}, must be at least "
                f"{'.'.join(map(str, MIN_PYTHON))}. This is not an independent "
                f"preference: bioconda's bte recipe is `skip: True  # [py < 310]` "
                f"at 0.9.3, so a py3.9 solve cannot satisfy bte>={MIN_BTE[0]}."
                f"{MIN_BTE[1]}.{MIN_BTE[2]} and resolves backwards instead. "
                f"Relaxing this while leaving the bte floor in place looks "
                f"harmless and reopens #136."
            )
        ceiling = re.search(r"python[^#]*<\s*([0-9.]+)", spec)
        if not ceiling:
            failures.append(
                f"python dependency {spec!r} declares no upper bound; bte is "
                f"incompatible with python >= 3.11"
            )
        else:
            parts = tuple(int(p) for p in ceiling.group(1).rstrip(".").split("."))
            if parts > MAX_PYTHON_EXCLUSIVE:
                failures.append(
                    f"python ceiling is < {ceiling.group(1)}, but bte is "
                    f"incompatible with python >= "
                    f"{'.'.join(map(str, MAX_PYTHON_EXCLUSIVE))}"
                )

    if not failures:
        print(
            f"  bte.yml pins bte>={'.'.join(map(str, MIN_BTE))} together with "
            f"python>={'.'.join(map(str, MIN_PYTHON))},"
            f"<{'.'.join(map(str, MAX_PYTHON_EXCLUSIVE))}"
        )
    return failures


def _reference_walk(tree, node_id):
    """The haplotype of ``node_id`` by walking root -> node and replaying branches.

    Deliberately independent of ``get_haplotype``: it reads each node's own
    ``mutations`` (which is what matUtils wrote, relative to that node's parent)
    and composes them by hand, tracking the reference base and the running base
    per site. That composition is the thing ``get_haplotype`` is supposed to do.

    Also validates the chain: each mutation's stated "from" base must equal the
    running base at that site, so a fixture that does not actually re-mutate --
    which would make this test pass vacuously -- raises instead.
    """
    path = list(reversed(list(tree.rsearch(node_id, True))))  # root -> node
    ref_base = {}
    cur_base = {}
    for node in path:
        for mutation in node.mutations:
            match = _MUT.match(mutation)
            if not match:
                raise AssertionError(
                    f"cannot parse mutation {mutation!r} on node {node.id}"
                )
            from_base, site, to_base = match.group(1), int(match.group(2)), match.group(3)
            if site not in ref_base:
                ref_base[site] = from_base
                cur_base[site] = from_base
            if cur_base[site] != from_base:
                raise AssertionError(
                    f"mutation {mutation!r} on node {node.id} claims to start from "
                    f"{from_base} but the running base at site {site} is "
                    f"{cur_base[site]}; the fixture or the tree is inconsistent"
                )
            cur_base[site] = to_base
    return {
        f"{ref_base[s]}{s}{cur_base[s]}"
        for s in ref_base
        if cur_base[s] != ref_base[s]
    }


def check_get_haplotype_matches_walk():
    """Half B: get_haplotype must equal an explicit root->node walk."""
    import bte

    failures = []
    tree = bte.MATree(nwk_string=NEWICK)

    node_ids = [n.id for n in tree.depth_first_expansion()]
    expected_ids = ["node_1", "node_2", "node_3", "leaf1", "leaf2"]
    if node_ids != expected_ids:
        return [
            f"bte named the fixture's nodes {node_ids}, expected {expected_ids}; "
            f"the mutation map below targets bte's auto-generated internal node "
            f"ids, so this test would apply mutations to the wrong branches"
        ]

    tree.apply_mutations(MUTATIONS)

    # The fixture must actually contain a re-mutation, or half B proves nothing.
    child_sites = {
        int(_MUT.match(m).group(2)) for m in tree.get_node("node_3").mutations
    }
    parent_sites = {
        int(_MUT.match(m).group(2)) for m in tree.get_node("node_2").mutations
    }
    twice_mutated = child_sites & parent_sites
    if not twice_mutated:
        return [
            "the fixture has no site mutated on both node_2 and node_3, so it "
            "does not exercise the bte<=0.9.2 defect at all"
        ]

    for node_id, expected in EXPECTED_HAPLOTYPES.items():
        actual = set(tree.get_haplotype(node_id))
        walked = _reference_walk(tree, node_id)

        if walked != expected:
            failures.append(
                f"{node_id}: the independent root->node walk produced {sorted(walked)} "
                f"but this test expects {sorted(expected)}; the fixture's hardcoded "
                f"expectation and its own walk disagree, so neither can be trusted"
            )
        if actual != walked:
            failures.append(
                f"{node_id}: get_haplotype returned {sorted(actual)} but the "
                f"root->node walk gives {sorted(walked)}. Missing: "
                f"{sorted(walked - actual)}; spurious: {sorted(actual - walked)}. "
                f"bte {getattr(bte, '__version__', '(unknown version)')} is "
                f"dropping or mis-composing branch mutations -- see #136."
            )

    # State the discriminating case explicitly: it is the one assertion whose
    # failure identifies the bug rather than merely reporting a mismatch.
    site = sorted(twice_mutated)[0]
    node_3 = set(tree.get_haplotype("node_3"))
    at_site = {m for m in node_3 if int(_MUT.match(m).group(2)) == site}
    if at_site != {"A100G"}:
        failures.append(
            f"site {site} is mutated on both node_2 (A->C) and node_3 (C->G), so "
            f"node_3's haplotype must carry exactly {{'A100G'}}; got {sorted(at_site)}. "
            f"{{'A100C'}} means bte kept the first mutation and discarded the "
            f"second -- the bte<=0.9.2 defect (#136)."
        )

    if not failures:
        version = getattr(bte, "__version__", "(unknown version)")
        print(
            f"  get_haplotype agrees with the root->node walk at all "
            f"{len(EXPECTED_HAPLOTYPES)} nodes, including the twice-mutated site "
            f"{site} (bte {version})"
        )
    return failures


def _normalize_env(text):
    """Comparable form of an env YAML: comments and blank lines stripped."""
    return [
        line.split("#", 1)[0].rstrip()
        for line in text.splitlines()
        if line.split("#", 1)[0].strip()
    ]


def _locate_bte_python():
    """A python interpreter that can import bte, or ``None``.

    Snakemake writes ``<hash>_.yaml`` next to the prefix directory ``<hash>_``, and
    a repo that has run the pipeline across more than one version of
    ``workflow/envs/bte.yml`` has **several** bte envs sitting side by side --
    including ones solved from the old ``bte>=0.1`` spec. Picking any of them (say,
    the first in sorted order) would let half B pass against a stale env while the
    shipped spec resolves to something else entirely, which is precisely the
    "tests something adjacent to the real thing" failure mode.

    So an env whose recorded YAML matches the current ``bte.yml`` is required
    first, and a non-matching env is used only as a last resort and announced.
    """
    if not SNAKEMAKE_CONDA.is_dir():
        return None

    wanted = _normalize_env(BTE_ENV_YAML.read_text())
    exact, fallback = [], []
    for yaml_path in sorted(SNAKEMAKE_CONDA.glob("*.yaml")):
        try:
            text = yaml_path.read_text()
        except OSError:
            continue
        if not re.search(r"^\s*-\s*bte\b", text, re.M):
            continue
        candidate = yaml_path.with_suffix("") / "bin" / "python"
        if not candidate.is_file():
            continue
        (exact if _normalize_env(text) == wanted else fallback).append(candidate)

    if exact:
        return exact[0]
    if fallback:
        print(
            f"  WARNING: no conda env matches the current {BTE_ENV_YAML.name}; "
            f"falling back to {fallback[0]}, which was solved from a different "
            f"spec. Re-run the pipeline to build the current env."
        )
        return fallback[0]
    return None


def _bte_importable():
    try:
        import bte  # noqa: F401
    except ImportError:
        return False
    return True


def run_all_checks(require_bte):
    failures = []

    print("Checking bte.yml pins the floor that closes #136...")
    try:
        failures.extend(check_env_pins())
    except Exception:
        failures.append(f"env pin check raised:\n{traceback.format_exc()}")

    print("Checking get_haplotype composes re-mutations correctly...")
    if _bte_importable():
        try:
            failures.extend(check_get_haplotype_matches_walk())
        except Exception:
            failures.append(
                f"get_haplotype check raised:\n{traceback.format_exc()}"
            )
    elif require_bte:
        failures.append(
            "--require-bte was passed but bte is not importable and no "
            f"Snakemake bte env was found under {SNAKEMAKE_CONDA}. Half B is the "
            "only check that exercises the actual defect, so a missing env is a "
            "failure here rather than a skip. Run the test-data pipeline first: "
            "pixi run snakemake -j 4 --use-conda --configfile "
            "config/config_test_data.yaml"
        )
    else:
        print(
            "  NOT RUN: bte is importable neither here nor from any Snakemake\n"
            "  conda env (bte lives only in workflow/envs/bte.yml's env, not in\n"
            "  the pixi env). Half A above still pinned the floor. To run half B:\n"
            "    pixi run snakemake -j 4 --use-conda --configfile config/config_test_data.yaml\n"
            "    pixi run python tests/test_bte_haplotype.py --require-bte"
        )
    return failures


def test_bte_haplotype():
    """Pytest entry point.

    This repo runs tests as standalone scripts (``pixi run python tests/<f>.py``)
    and does not depend on pytest. This wrapper exists so that anyone who does
    reach for ``pytest tests/`` gets a real result rather than a green run that
    collected nothing.
    """
    failures = run_all_checks(require_bte=False)
    assert not failures, "\n".join(failures)


def main():
    require_bte = "--require-bte" in sys.argv[1:]

    # If bte is missing here but a Snakemake bte env exists, re-run inside it so
    # half B actually executes rather than reporting NOT RUN.
    if not _bte_importable() and not os.environ.get(REEXEC_FLAG):
        interpreter = _locate_bte_python()
        if interpreter is not None:
            print(f"bte is not importable here; re-executing under {interpreter}")
            os.environ[REEXEC_FLAG] = "1"
            os.execv(str(interpreter), [str(interpreter), __file__, *sys.argv[1:]])

    failures = run_all_checks(require_bte)
    if failures:
        print(f"\nFAILED ({len(failures)} problem(s)):")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    print(
        "\nPASSED: bte.yml pins a version whose get_haplotype composes "
        "re-mutations, and (where bte was available) it demonstrably does."
    )


if __name__ == "__main__":
    main()
