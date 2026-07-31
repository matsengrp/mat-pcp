#!/usr/bin/env python3
"""Pin the PCP gene coordinate override and the codon alignment it feeds (#133).

Two changes, both of which fail silently if wrong:

1. **``pcp_gene_cds_blocks``** lets a PCP-only gene (ORF1b = nsp12-nsp16) be
   given coordinates the shared reference GTF does not carry. A mis-specified
   block does not crash -- it produces a frame-shifted protein that still
   translates and still looks like valid output. So every block is validated,
   and the override is add-only: shadowing a GTF gene would leave the PCP slice
   and every fitness output disagreeing about one gene's coordinates.

2. **``pcp_gene_to_fitness_gene``** aligns a PCP gene with the ``aamut_fitness``
   gene it is validated against. Through v8 the validation filtered on
   ``gene == {gene}``, and ``aamutfit.get_coding`` drops ORF1a, so the ORF1a
   check compared against **zero rows** and reported no discrepancies because it
   had nothing to compare.

The alignment is derived from genomic positions, not from the configured offset,
and check_codon_correspondence_beats_a_bare_offset demonstrates why: under a
bare ``codon_offset: 0``, ORF1a's codons 4402-4405 would map onto ORF1ab codons
4402-4405, which are *ORF1b's* codons 10-13 -- four spurious discrepancies
indistinguishable from real ones. The configured offset is kept as a cross-check
so a stale number fails loudly.

The coordinate identity is asserted on a synthetic genome on purpose: that
ORF1b's blocks reproduce ORF1ab's interior slice is pure positional arithmetic
and holds for any sequence. The genome-specific facts (2704 aa, no internal
stops, the nsp12 N-terminus) are dataset invariants and live in
``validate_test.py``.

Usage:
    python tests/test_pcp_gene_coordinates.py
"""

import random
import sys
import tempfile
import traceback
from importlib import util as importlib_util
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# Real coordinates, as verified in #133 against the shipped ORF1ab founders.
ORF1AB_BLOCKS = [(266, 13468), (13468, 21552)]
ORF1A_BLOCKS = [(266, 13480)]
ORF1B_BLOCKS = [(13442, 13468), (13468, 21552)]
S_BLOCKS = [(21563, 25381)]
GENOME_LENGTH = 29903

# ORF1b is ORF1ab codons 4393-7096, i.e. nt offsets [13176:21288] of the
# reconstructed ORF1ab CDS.
ORF1B_SLICE_START_NT = 13176
ORF1B_SLICE_END_NT = 21288
ORF1B_OFFSET = 4392
ORF1A_ORPHAN_CODONS = [4402, 4403, 4404, 4405]


def _load(module_name):
    """Load a workflow module by path, as the pipeline scripts do."""
    spec = importlib_util.spec_from_file_location(
        module_name, REPO_ROOT / "workflow" / "modules" / f"{module_name}.py"
    )
    module = importlib_util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _synthetic_gtf(path, genes):
    """Write a minimal GTF carrying ``genes`` as {name: [(start, end), ...]}."""
    with open(path, "w") as f:
        for gene, blocks in genes.items():
            for start, end in blocks:
                f.write(
                    f'NC_045512.2\ttest\tCDS\t{start}\t{end}\t.\t+\t0\t'
                    f'gene_id "{gene}"; transcript_id "{gene}";\n'
                )
    return path


def _reference_gtf(tmpdir):
    """A GTF with the real ORF1ab/ORF1a/S blocks and deliberately no ORF1b."""
    return _synthetic_gtf(
        Path(tmpdir) / "ref.gtf",
        {"ORF1ab": ORF1AB_BLOCKS, "ORF1a": ORF1A_BLOCKS, "S": S_BLOCKS},
    )


def check_override_validation():
    """Every malformed override must raise, naming the gene."""
    failures = []
    gtf = _load("gtf")

    bad = {
        "not divisible by 3": [[100, 110]],
        "unsorted blocks": [[500, 502], [100, 102]],
        "end before start": [[500, 400]],
        "start below 1": [[0, 2]],
        "past the reference end": [[29900, 29910]],
        "empty list": [],
        "not a list": "13442-21552",
        "three-element block": [[1, 2, 3]],
        "non-integer coordinate": [["13442", 13468]],
        "boolean coordinate": [[True, 12]],
    }
    for label, blocks in bad.items():
        try:
            gtf.validate_cds_blocks("ORF1b", blocks, genome_length=GENOME_LENGTH)
        except ValueError as exc:
            if "ORF1b" not in str(exc):
                failures.append(
                    f"{label}: raised without naming the gene: {exc}"
                )
        else:
            failures.append(f"{label}: accepted {blocks!r} without raising")

    # The real blocks must pass, and come back as tuples.
    try:
        parsed = gtf.validate_cds_blocks(
            "ORF1b", [[13442, 13468], [13468, 21552]], genome_length=GENOME_LENGTH
        )
        if parsed != ORF1B_BLOCKS:
            failures.append(f"real ORF1b blocks parsed to {parsed!r}")
    except ValueError as exc:
        failures.append(f"real ORF1b blocks rejected: {exc}")

    # Blocks sharing a nucleotide must be accepted: that shared nt is exactly
    # what keeps the frameshift in frame, and it is why the divisibility check
    # sums the blocks rather than measuring the span. ORF1b's span
    # (21552-13442+1 = 8111) is NOT divisible by 3; its summed length (8112) is.
    span = ORF1B_BLOCKS[-1][1] - ORF1B_BLOCKS[0][0] + 1
    total = sum(end - start + 1 for start, end in ORF1B_BLOCKS)
    if span % 3 == 0:
        failures.append(
            f"ORF1b span {span} is divisible by 3, so this check no longer "
            f"distinguishes summing blocks from measuring the span"
        )
    if total != 8112:
        failures.append(f"ORF1b summed length is {total}, expected 8112")
    if not failures:
        print(f"  10 malformed overrides rejected; real blocks parse "
              f"(sum {total} nt = {total // 3} codons, span {span} nt)")
    return failures


def check_override_is_add_only():
    """An override for a gene the GTF already defines must raise, not shadow."""
    failures = []
    gtf = _load("gtf")
    with tempfile.TemporaryDirectory() as tmpdir:
        gtf_path = _reference_gtf(tmpdir)

        try:
            gtf.resolve_cds_blocks(
                gtf_path, "S", {"S": [[21563, 25381]]}, genome_length=GENOME_LENGTH
            )
        except ValueError as exc:
            if "add-only" not in str(exc):
                failures.append(f"raised, but not about add-only: {exc}")
        else:
            failures.append(
                "an override for S (a GTF gene) was silently accepted, leaving "
                "the PCP slice and fitness outputs free to disagree"
            )

        # ORF1b is absent from the GTF, so its override is the legitimate case.
        try:
            blocks = gtf.resolve_cds_blocks(
                gtf_path,
                "ORF1b",
                {"ORF1b": [[13442, 13468], [13468, 21552]]},
                genome_length=GENOME_LENGTH,
            )
            if blocks != ORF1B_BLOCKS:
                failures.append(f"ORF1b resolved to {blocks!r}")
        except ValueError as exc:
            failures.append(f"ORF1b override rejected: {exc}")
    if not failures:
        print("  overriding a GTF gene raises; adding an absent gene works")
    return failures


def check_override_is_inert_when_unused():
    """With no override, resolution must equal the pre-#133 GTF-only behaviour.

    Acceptance criterion 7: every existing output stays byte-identical when
    ORF1b is absent from ``pcp_genes``.
    """
    failures = []
    gtf = _load("gtf")
    with tempfile.TemporaryDirectory() as tmpdir:
        gtf_path = _reference_gtf(tmpdir)
        for gene in ["ORF1ab", "ORF1a", "S"]:
            baseline = gtf.cds_blocks_for_gene(gtf_path, gene)
            for overrides in [None, {}, {"ORF1b": [[13442, 13468], [13468, 21552]]}]:
                resolved = gtf.resolve_cds_blocks(
                    gtf_path, gene, overrides, genome_length=GENOME_LENGTH
                )
                if resolved != baseline:
                    failures.append(
                        f"{gene} with overrides={overrides!r} resolved to "
                        f"{resolved!r}, not the GTF's {baseline!r}"
                    )
        # A gene in neither the GTF nor the override must still raise clearly.
        try:
            gtf.resolve_cds_blocks(gtf_path, "ORF1b", {}, genome_length=GENOME_LENGTH)
        except ValueError as exc:
            if "No CDS blocks found" not in str(exc):
                failures.append(f"unknown gene raised unexpectedly: {exc}")
        else:
            failures.append("a gene in neither GTF nor override did not raise")
    if not failures:
        print("  GTF genes resolve identically with and without an override present")
    return failures


def check_orf1b_reproduces_the_orf1ab_interior_slice():
    """ORF1b's blocks must rebuild exactly ORF1ab's codons 4393-7096.

    Pure positional arithmetic, so it holds for any sequence -- which is why a
    synthetic genome is the right fixture and pins the coordinates rather than
    the genome.
    """
    failures = []
    sequences = _load("sequences")
    rng = random.Random(1729)
    genome = "".join(rng.choice("ACGT") for _ in range(GENOME_LENGTH))

    orf1ab = sequences.slice_cds(genome, ORF1AB_BLOCKS)
    orf1b = sequences.slice_cds(genome, ORF1B_BLOCKS)
    expected = orf1ab[ORF1B_SLICE_START_NT:ORF1B_SLICE_END_NT]

    if len(orf1ab) != 21288:
        failures.append(f"ORF1ab slice is {len(orf1ab)} nt, expected 21288")
    if len(orf1b) != 8112:
        failures.append(f"ORF1b slice is {len(orf1b)} nt, expected 8112")
    if orf1b != expected:
        failures.append(
            f"ORF1b slice != ORF1ab[{ORF1B_SLICE_START_NT}:{ORF1B_SLICE_END_NT}] "
            f"-- the blocks do not describe the intended interior slice"
        )
    if ORF1B_SLICE_START_NT % 3 != 0:
        failures.append(
            f"slice start {ORF1B_SLICE_START_NT} is not codon-aligned in ORF1ab"
        )
    if not failures:
        print(
            f"  ORF1b == ORF1ab[{ORF1B_SLICE_START_NT}:{ORF1B_SLICE_END_NT}] "
            f"({len(orf1b)} nt = {len(orf1b) // 3} codons), codon-aligned"
        )
    return failures


def check_codon_correspondence_beats_a_bare_offset():
    """The derived map must exclude ORF1a's four post-frameshift codons.

    This is the check that justifies deriving the correspondence instead of
    trusting ``codon_offset``. The mutation half shows what a bare offset does.
    """
    failures = []
    genemap = _load("genemap")

    codon_map, unmapped, offset = genemap.pcp_to_fitness_codons(
        ORF1A_BLOCKS, ORF1AB_BLOCKS
    )
    if offset != 0:
        failures.append(f"ORF1a offset derived as {offset}, expected 0")
    if unmapped != ORF1A_ORPHAN_CODONS:
        failures.append(
            f"ORF1a unmapped codons {unmapped!r}, expected {ORF1A_ORPHAN_CODONS!r}"
        )
    if len(codon_map) != 4401:
        failures.append(f"ORF1a mapped {len(codon_map)} codons, expected 4401")
    if codon_map.get(4401) != 4401:
        failures.append("ORF1a codon 4401 should map to ORF1ab codon 4401")

    # Mutation half: a bare offset would map the orphans onto real ORF1ab
    # codons, which belong to ORF1b's frame. Show those targets are genuinely
    # occupied, so the resulting rows would look like ordinary discrepancies.
    orf1b_map, _, orf1b_offset = genemap.pcp_to_fitness_codons(
        ORF1B_BLOCKS, ORF1AB_BLOCKS
    )
    if orf1b_offset != ORF1B_OFFSET:
        failures.append(f"ORF1b offset derived as {orf1b_offset}, expected 4392")
    if len(orf1b_map) != 2704 or orf1b_map.get(1) != 4393:
        failures.append(
            f"ORF1b mapped {len(orf1b_map)} codons, codon 1 -> "
            f"{orf1b_map.get(1)}; expected 2704 and 4393"
        )
    collides = [c for c in ORF1A_ORPHAN_CODONS if c in set(orf1b_map.values())]
    if collides != ORF1A_ORPHAN_CODONS:
        failures.append(
            f"expected ORF1a's orphan codons {ORF1A_ORPHAN_CODONS!r} to collide "
            f"with ORF1b's frame under a bare offset; got {collides!r}"
        )

    # S against itself is the trivial identity, and must stay so.
    s_map, s_unmapped, s_offset = genemap.pcp_to_fitness_codons(S_BLOCKS, S_BLOCKS)
    if s_offset != 0 or s_unmapped or len(s_map) != 1273:
        failures.append(
            f"S vs S should be the identity over 1273 codons; got "
            f"{len(s_map)} mapped, {len(s_unmapped)} unmapped, offset {s_offset}"
        )
    if not failures:
        print(
            f"  ORF1a: 4401 mapped, {ORF1A_ORPHAN_CODONS} orphaned (they land on "
            f"ORF1b's frame under a bare offset); ORF1b: 2704 mapped at +4392"
        )
    return failures


def check_declared_offset_is_cross_checked():
    """A configured offset that disagrees with the coordinates must raise."""
    failures = []
    genemap = _load("genemap")

    # Correct offsets pass.
    for gene, blocks, offset, n_expected in [
        ("ORF1a", ORF1A_BLOCKS, 0, 4401),
        ("ORF1b", ORF1B_BLOCKS, ORF1B_OFFSET, 2704),
    ]:
        try:
            codon_map, _ = genemap.align_to_fitness_gene(
                gene, "ORF1ab", blocks, ORF1AB_BLOCKS, offset
            )
            if len(codon_map) != n_expected:
                failures.append(
                    f"{gene}: aligned {len(codon_map)} codons, expected {n_expected}"
                )
        except ValueError as exc:
            failures.append(f"{gene}: correct offset {offset} rejected: {exc}")

    # Wrong offsets raise, including off-by-one -- the case a human would most
    # plausibly write and least plausibly notice.
    for wrong in [ORF1B_OFFSET - 1, ORF1B_OFFSET + 1, 0, 9999]:
        try:
            genemap.align_to_fitness_gene(
                "ORF1b", "ORF1ab", ORF1B_BLOCKS, ORF1AB_BLOCKS, wrong
            )
        except ValueError as exc:
            if "codon_offset" not in str(exc):
                failures.append(f"offset {wrong}: raised, but not about the offset")
        else:
            failures.append(
                f"declared offset {wrong} accepted though the coordinates imply "
                f"{ORF1B_OFFSET}"
            )

    # Genes that do not correspond at all must raise rather than align emptily.
    try:
        genemap.align_to_fitness_gene("S", "ORF1ab", S_BLOCKS, ORF1AB_BLOCKS, 0)
    except ValueError as exc:
        if "no codon" not in str(exc).lower():
            failures.append(f"disjoint genes raised unexpectedly: {exc}")
    else:
        failures.append("S vs ORF1ab (disjoint) aligned without raising")
    if not failures:
        print("  correct offsets align; wrong offsets and disjoint genes raise")
    return failures


def check_empty_fitness_slice_is_fatal():
    """PCP counts with no fitness rows must fail, not report a clean run."""
    failures = []
    genemap = _load("genemap")

    # The v8 ORF1a situation: PCP counts exist, the fitness slice is empty.
    try:
        genemap.require_nonempty_fitness_slice(
            "ORF1a", "ORF1a", 1234, 0, "all_aamut_fitness.csv"
        )
    except ValueError as exc:
        if "vacuous" not in str(exc):
            failures.append(f"raised, but not about vacuity: {exc}")
    else:
        failures.append(
            "an empty fitness slice with 1234 PCP count rows was accepted -- "
            "this is exactly the v8 ORF1a gap"
        )

    # Both non-empty is the normal case; both empty is a gene with no PCP data
    # at all, which is not this guard's business.
    for n_pcp, n_fit, label in [(10, 10, "both non-empty"), (0, 0, "both empty")]:
        try:
            genemap.require_nonempty_fitness_slice(
                "S", "S", n_pcp, n_fit, "all_aamut_fitness.csv"
            )
        except ValueError as exc:
            failures.append(f"{label} should pass, but raised: {exc}")
    if not failures:
        print("  empty slice with PCP counts is fatal; other combinations pass")
    return failures


CHECKS = [
    ("override validation", check_override_validation),
    ("override is add-only", check_override_is_add_only),
    ("override is inert when unused", check_override_is_inert_when_unused),
    ("ORF1b reproduces the ORF1ab interior slice",
     check_orf1b_reproduces_the_orf1ab_interior_slice),
    ("codon correspondence beats a bare offset",
     check_codon_correspondence_beats_a_bare_offset),
    ("declared offset is cross-checked", check_declared_offset_is_cross_checked),
    ("empty fitness slice is fatal", check_empty_fitness_slice_is_fatal),
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


def test_pcp_gene_coordinates():
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
        "\nPASSED: the PCP coordinate override is validated and add-only, and "
        "the fitness-gene codon alignment is derived, cross-checked, and "
        "cannot be vacuous."
    )


if __name__ == "__main__":
    main()
