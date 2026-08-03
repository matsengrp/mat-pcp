#!/usr/bin/env python3
"""Pin the clade-founder guards that stop a silent fallback to Wuhan (#129).

``clade_founder_fasta.py`` used to write the *unmodified reference* as a clade's
founder whenever it could not resolve that clade's defining mutations -- exit 0,
logged as ordinary informational output. In v8 that silently affected 24H, 24I,
25A, 25B and 25C. Nothing downstream noticed, because the wrong founder
propagated consistently: every consumer compared against the same wrong
sequence, so cross-clade analyses were measuring a shift towards Wuhan without
any internal inconsistency to trip on.

Two guards now cover the two routes to that outcome, and this test pins both,
plus the exemptions that must keep working:

* **source missing** -- clade in none of Neher / the Roemer map / Roemer direct.
  This is the v8 bug exactly.
* **result reference-identical** -- a source resolves but yields no *effective*
  mutations (empty substitution list, or every mutation dropped by the
  reference-base mismatch guard). Different route, same wrong founder.
* **reference clade exempt** -- 19A's founder legitimately *is* the reference,
  so it must not raise. Without this the guard would fail every run.
* **config allowlist exempt** -- ``allow_reference_identical_founders`` records a
  deliberate case, so a future legitimate clade needs no code change.

``clade_founder_fasta.py`` is a Snakemake script: it reads a module-level
``snakemake`` object rather than taking arguments. Rather than refactor it just
to be testable, this runs it the way Snakemake does, with a stub ``snakemake``
injected into its globals. That keeps the *shipped* code path under test --
including the params wiring, which a unit test of an extracted helper would miss.

Usage:
    python tests/test_clade_founder_guards.py
"""

import json
import sys
import tempfile
import traceback
import types
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "workflow" / "scripts" / "clade_founder_fasta.py"


def _ensure_biopython():
    """Make ``import Bio.SeqIO`` work, with a stub if biopython is absent.

    The script under test reads its reference with ``Bio.SeqIO.read``, but
    biopython lives in ``workflow/envs/default.yml`` rather than the pixi
    environment -- this repo deliberately keeps science dependencies in the conda
    envs (see the "scipy and statsmodels are in conda envs, not pixi env" note in
    the CI workflow). Rather than widen the pixi env for one test, fall back to a
    minimal single-record FASTA reader so the test runs anywhere.

    Real biopython is preferred when present, so a conda-env run exercises the
    actual parser. FASTA parsing is not what this test pins -- the guards are --
    so the stub costs no coverage that matters here.
    """
    try:
        import Bio.SeqIO  # noqa: F401

        return "biopython"
    except ModuleNotFoundError:
        pass

    class _Record:
        def __init__(self, seq):
            self.seq = seq

    def _read(path, fmt):
        assert fmt == "fasta", f"stub only handles fasta, got {fmt!r}"
        lines = Path(path).read_text().splitlines()
        if not lines or not lines[0].startswith(">"):
            raise ValueError(f"{path} is not a FASTA file")
        body = [ln.strip() for ln in lines[1:] if ln.strip() and not ln.startswith(">")]
        return _Record("".join(body))

    bio = types.ModuleType("Bio")
    seqio = types.ModuleType("Bio.SeqIO")
    seqio.read = _read
    bio.SeqIO = seqio
    sys.modules.setdefault("Bio", bio)
    sys.modules.setdefault("Bio.SeqIO", seqio)
    return "stub"

# A short stand-in for the genome (1-indexed, as mutation strings are):
#
#   pos 1234567890...
#       ATGCGGTACCATGCAATTGGCC...
#
# so position 5 is G, 7 is T, 9 is C, 11 is A. Mutation strings below must state
# the base that is actually there, or apply_mutations skips them -- which is a
# *tested* behaviour (see the 24I case), so getting it wrong by accident turns a
# success case into a spurious raise rather than a silent pass.
REF_SEQ = "ATGCGGTACCATGCAATTGGCCATGCATGCATGCGGTTACCAATTGGCCATG"

# Keyed on Nextstrain clade name, as the real Neher JSON is.
NEHER = {
    "19A": {"nuc": []},          # reference clade: legitimately no mutations
    "19B": {"nuc": ["C9T", "A11G"]},
    "20A": {"nuc": ["G5A"]},
}

# Keyed on Pango lineage, as the real Roemer JSON is.
ROEMER = {
    "XFG": {"nucSubstitutions": ["T7G", "A11C"]},
    "LF.7": {"nucSubstitutions": []},          # resolves, but no substitutions
    "MV.1": {"nucSubstitutions": ["G7C"]},     # position 7 is T, not G -> skipped
}

MAP = {"25C": "XFG", "24H": "LF.7", "24I": "MV.1"}


class _Namespace:
    """Stand-in for Snakemake's input/output/params/wildcards objects."""

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def get(self, key, default=None):
        return self.__dict__.get(key, default)


def _run_script(tmpdir, clade, roemer_map=None, allow_identical=None,
                reference_clade="19A"):
    """Execute clade_founder_fasta.py the way Snakemake does.

    Returns the written founder sequence on success; raises whatever the script
    raises on failure.
    """
    tmpdir = Path(tmpdir)
    neher_path = tmpdir / "neher.json"
    roemer_path = tmpdir / "roemer.json"
    ref_path = tmpdir / "ref.fa"
    out_path = tmpdir / f"{clade}.fa"
    log_path = tmpdir / f"{clade}.log"

    neher_path.write_text(json.dumps(NEHER))
    roemer_path.write_text(json.dumps(ROEMER))
    ref_path.write_text(f">ref\n{REF_SEQ}\n")

    stub = _Namespace(
        input=_Namespace(
            neher_json=str(neher_path),
            roemer_json=str(roemer_path),
            ref_fasta=str(ref_path),
        ),
        output=_Namespace(fasta=str(out_path)),
        params=_Namespace(
            roemer_nextstrain_to_pango=MAP if roemer_map is None else roemer_map,
            reference_clade=reference_clade,
            allow_reference_identical_founders=allow_identical or [],
        ),
        wildcards=_Namespace(clade=clade),
        log=[str(log_path)],
    )

    # exec rather than import: the script runs top-level statements against the
    # `snakemake` global, which is exactly how Snakemake invokes it.
    script_globals = {"__name__": "__main__", "__file__": str(SCRIPT), "snakemake": stub}
    exec(compile(SCRIPT.read_text(), str(SCRIPT), "exec"), script_globals)

    written = out_path.read_text().splitlines()
    return "".join(written[1:])


def _expect_raise(tmpdir, label, must_mention, **kwargs):
    """The script must raise, and the message must be actionable."""
    try:
        _run_script(tmpdir, **kwargs)
    except ValueError as exc:
        message = str(exc)
        missing = [m for m in must_mention if m not in message]
        if missing:
            return [
                f"{label}: raised as required, but the message omits {missing}. "
                f"An operator has to be able to act on it. Got:\n{message}"
            ]
        print(f"  {label}: raised ValueError as required")
        return []
    except Exception:  # noqa: BLE001 - any other failure is itself a failure
        return [
            f"{label}: raised the wrong exception type:\n{traceback.format_exc()}"
        ]
    return [
        f"{label}: did NOT raise. This is the #129 failure mode -- the reference "
        "would be written as this clade's founder and every downstream consumer "
        "would agree with it."
    ]


def _expect_success(tmpdir, label, expect_differs, **kwargs):
    """The script must succeed, and the founder must (or must not) differ."""
    try:
        seq = _run_script(tmpdir, **kwargs)
    except Exception:  # noqa: BLE001
        return [f"{label}: unexpectedly failed:\n{traceback.format_exc()}"]

    differs = seq != REF_SEQ
    if differs != expect_differs:
        return [
            f"{label}: founder {'differs from' if differs else 'is identical to'} "
            f"the reference, expected the opposite"
        ]
    print(
        f"  {label}: succeeded, founder "
        f"{'differs from' if differs else 'matches'} the reference as expected"
    )
    return []


def run_all_checks():
    failures = []

    if not SCRIPT.is_file():
        return [f"missing script {SCRIPT}"]

    print(f"  FASTA reader: {_ensure_biopython()}")

    with tempfile.TemporaryDirectory() as tmpdir:
        print("Checking an unresolvable clade fails loudly (#129)...")
        # 25B is in none of NEHER, MAP, or ROEMER -- the v8 scenario.
        failures += _expect_raise(
            tmpdir,
            "clade absent from every source",
            [
                "25B",
                "roemer_nextstrain_to_pango",
                "#129",
                "own copy",  # each config has its own map; fixing one isn't enough
            ],
            clade="25B",
        )

        print("Checking a clade dropped from the map fails loudly...")
        # 25C resolves normally, but not once the map entry is removed -- this is
        # what "a config falls behind" looks like.
        failures += _expect_raise(
            tmpdir,
            "clade missing from the pango map",
            ["25C", "roemer_nextstrain_to_pango"],
            clade="25C",
            roemer_map={},
        )

        print("Checking a resolved-but-empty source fails loudly...")
        # 24H -> LF.7, which exists in Roemer with an empty substitution list.
        # The source check passes; only the result check catches this.
        failures += _expect_raise(
            tmpdir,
            "source resolves with no substitutions",
            ["24H", "identical to the reference", "#129"],
            clade="24H",
        )

        print("Checking a fully-skipped mutation list fails loudly...")
        # 24I -> MV.1 claims G at position 7, which is A. apply_mutations skips
        # it with a warning, leaving the founder equal to the reference.
        failures += _expect_raise(
            tmpdir,
            "every mutation skipped on reference mismatch",
            ["24I", "identical to the reference", "WARNING"],
            clade="24I",
        )

        print("Checking legitimate founders still build...")
        failures += _expect_success(
            tmpdir, "neher-sourced clade", True, clade="19B"
        )
        failures += _expect_success(
            tmpdir, "roemer-sourced clade via map", True, clade="25C"
        )

        print("Checking the exemptions still work...")
        # Without this the guard would fail every run: 19A's founder *is* the
        # reference by definition.
        failures += _expect_success(
            tmpdir, "reference clade 19A", False, clade="19A"
        )
        # A non-default reference clade must be honoured, not hardcoded.
        failures += _expect_success(
            tmpdir,
            "reference clade honoured from config",
            False,
            clade="24H",
            reference_clade="24H",
        )
        # And the config allowlist records a deliberate case without a code edit.
        failures += _expect_success(
            tmpdir,
            "config allowlist exemption",
            False,
            clade="24H",
            allow_identical=["24H"],
        )

    return failures


def test_clade_founder_guards():
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
    print("\nPASSED: clade-founder guards reject a silent fallback to the reference.")


if __name__ == "__main__":
    main()
