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
        # Ascending starts and a whole-codon total do NOT imply non-overlap:
        # [[1, 100], [50, 201]] satisfies both (252 nt) while yielding the codon
        # (100, 51, 52) -- genomic positions running backwards. Caught in the
        # PR #132 review; a frame-shifted slice still translates, so nothing
        # downstream would have objected.
        "blocks overlapping by more than the frameshift nt": [[1, 100], [50, 201]],
        "second block nested inside the first": [[1, 300], [10, 20]],
        # A shared nucleotide is only in frame on a codon boundary. Both of these
        # pass ascending starts, pass the overlap check (they share exactly one
        # nt) and total a whole number of codons -- caught in the second PR #132
        # review pass, which found the first accepted at 21 nt and the second at
        # 12 nt. [[1, 11], [11, 20]] reads nt 11 twice *within* codon (10, 11, 11).
        "shared nt not on a codon boundary": [[1, 11], [11, 20]],
        "two shared junctions, first out of frame": [[1, 4], [4, 7], [7, 10]],
        # Both junctions here are on codon boundaries, so only the "at most one
        # frameshift" rule rejects it. Without that rule this translates cleanly.
        "two in-frame shared junctions": [[1, 3], [3, 5], [5, 7]],
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

    # Abutting and disjoint blocks must stay legal -- the overlap check above
    # must reject only overlaps *larger* than the shared frameshift nucleotide.
    for label, blocks in {
        "sharing exactly the frameshift nt": [[13442, 13468], [13468, 21552]],
        "abutting": [[1, 3], [4, 6]],
        "disjoint": [[1, 3], [100, 102]],
    }.items():
        try:
            gtf.validate_cds_blocks("OK", blocks, genome_length=GENOME_LENGTH)
        except ValueError as exc:
            failures.append(f"{label}: legal blocks rejected: {exc}")
    if not failures:
        print(f"  {len(bad)} malformed overrides rejected, 3 legal shapes "
              f"accepted; real blocks parse (sum {total} nt = {total // 3} "
              f"codons, span {span} nt)")
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


MASK_SCRIPT = REPO_ROOT / "workflow" / "scripts" / "get_problematic_sites.py"


class _Namespace:
    """Stand-in for Snakemake's input/output/params/wildcards objects."""

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def get(self, key, default=None):
        return self.__dict__.get(key, default)


def _run_mask_script(tmpdir, gene, rows, overrides=None, gtf_genes=None,
                     drop_columns=None):
    """Execute get_problematic_sites.py the way Snakemake does.

    ``rows`` is a list of ``(nt_site, exclude, gene_cell, codon_site_cell)``.
    Returns the codon sites actually written to the mask file; raises whatever
    the script raises.

    exec rather than import, because the script runs top-level statements
    against the ``snakemake`` global -- exactly how Snakemake invokes it. This is
    the point of the check: the first version of this test re-implemented the
    script's nt_site -> codon inversion and asserted against its own copy, so
    three separate breaks of the shipped script (always-empty mask, dropping the
    shared frameshift nt, disabling the resolution guard) all left it green.
    """
    tmpdir = Path(tmpdir)
    counts_path = tmpdir / "counts_all.csv"
    out_path = tmpdir / f"{gene}_codon_sites_to_mask.txt"
    log_path = tmpdir / f"{gene}.log"

    drop = set(drop_columns or ())
    columns = [c for c in ("nt_site", "exclude", "gene", "codon_site",
                           "clade", "nt_mutation") if c not in drop]
    with open(counts_path, "w") as f:
        f.write(",".join(columns) + "\n")
        for index, (nt_site, exclude, gene_cell, codon_cell) in enumerate(rows):
            values = {
                "nt_site": nt_site, "exclude": exclude, "gene": gene_cell,
                "codon_site": codon_cell, "clade": "19A",
                "nt_mutation": f"A{nt_site}T_{index}",
            }
            f.write(",".join(str(values[c]) for c in columns) + "\n")

    gtf_path = _synthetic_gtf(
        tmpdir / "ref.gtf",
        gtf_genes if gtf_genes is not None else {"ORF1ab": ORF1AB_BLOCKS},
    )

    stub = _Namespace(
        input=_Namespace(expected_vs_actual=str(counts_path), gtf=str(gtf_path)),
        output=_Namespace(sites_to_mask=str(out_path)),
        params=_Namespace(pcp_gene_cds_blocks=overrides or {}),
        wildcards=_Namespace(gene=gene),
        log=[str(log_path)],
    )
    script_globals = {
        "__name__": "__main__",
        "__file__": str(MASK_SCRIPT),
        "snakemake": stub,
    }
    exec(compile(MASK_SCRIPT.read_text(), str(MASK_SCRIPT), "exec"), script_globals)
    written = out_path.read_text().strip()
    return [int(line) for line in written.splitlines() if line]


def check_pcp_only_gene_gets_a_nonempty_codon_mask():
    """A `pcp_gene_cds_blocks` gene's codon mask must come from its CDS blocks.

    `counts_all.csv`'s `gene` column is derived from the reference GTF, so a gene
    added only through `pcp_gene_cds_blocks` never appears in it. Matching on the
    name selected zero rows and `get_problematic_sites` wrote an **empty** mask at
    exit 0 -- so ORF1b's DASM-ready PCPs kept every problematic codon site that
    every other gene had removed (35 of them on the v9 test tree). Found in the
    PR #132 review; CI's existing check was `test -f`, which an empty file passes.

    Runs the real script against a synthetic table so it holds without a pipeline
    run, and asserts the codon numbers written to the file -- including that the
    shared frameshift nucleotide yields *two* codons, which is the arithmetic a
    naive one-codon-per-nt inversion would get wrong.
    """
    failures = []

    # nt 13442 -> codon 1; nt 13468 is the shared frameshift nt, so codons 9 and
    # 10 (block 1 contributes 27 nt = 9 codons); nt 21552 -> the last codon 2704.
    # nt 13445 is only *partly* excluded, so it must not appear (codon 2), and
    # nt 100 is outside ORF1b entirely.
    rows = [
        (13442, True, "nan", ""),
        (13468, True, "ORF1a;ORF1ab;ORF1ab", "4401;4401;4402"),
        (21552, True, "ORF1ab", "7096"),
        (13445, True, "ORF1ab", "4394"),
        (13445, False, "ORF1ab", "4394"),
        (100, True, "nan", ""),
    ]
    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            got = _run_mask_script(
                tmpdir, "ORF1b", rows, overrides={"ORF1b": ORF1B_BLOCKS}
            )
        except Exception as exc:  # noqa: BLE001 - reported as a failure
            return [f"ORF1b mask derivation raised: {exc}"]

    expected = [1, 9, 10, 2704]
    if got != expected:
        failures.append(
            f"ORF1b mask is {got}, expected {expected} (codon 1 from nt 13442, "
            f"codons 9 and 10 from the shared frameshift nt 13468, codon 2704 "
            f"from nt 21552; nt 13445 is only partly excluded and nt 100 is "
            f"outside the gene)"
        )
    if not got:
        failures.append(
            "ORF1b mask is empty -- the #132 defect, which a `test -f` check "
            "and a re-implemented inversion both pass"
        )
    if not failures:
        print(f"  the shipped script derives ORF1b codons {got} from nt_site "
              f"alone; the shared frameshift nt yields two codons")
    return failures


def check_mask_script_guards_and_gene_pairing():
    """The mask script's three silent-failure routes must all raise, and it must
    not file another gene's codon number under this one.

    At an overlap site `ref_coding_sites` joins every covering CDS feature, so
    `gene` and `codon_site` are `;`-joined and paired **positionally**
    (`"N;ORF9b"` with `"4;1"`). The script flattened both and appended every
    component, filing ORF9b's codon 1 as a codon of N. ORF9b's codon numbers
    (1-98) are all valid N codon numbers, so no range check could see it -- it
    over-masks a real codon of the wrong gene. Latent rather than firing on the
    v9 test tree (297 sites have a divergent compound `codon_site`, none of them
    fully excluded), and newly reachable because v9 adds N to `pcp_genes`.
    """
    failures = []
    overlap_rows = [
        (28300, True, "N;ORF9b", "4;1"),
        (28400, True, "N", "37"),
    ]
    n_blocks = {"N": [(28274, 29533)], "ORF1ab": ORF1AB_BLOCKS}

    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            got = _run_mask_script(tmpdir, "N", overlap_rows, gtf_genes=n_blocks)
        except Exception as exc:  # noqa: BLE001 - reported as a failure
            failures.append(f"N mask derivation raised: {exc}")
            got = None
    if got is not None and got != [4, 37]:
        failures.append(
            f"N's mask is {got}, expected [4, 37]. [1, 4, 37] means ORF9b's "
            f"codon 1 was filed as a codon of N: the components of 'gene' and "
            f"'codon_site' are positionally paired and must be zipped, not "
            f"flattened."
        )

    # The gene name must match a whole component, not a substring: `ORF1a` is a
    # substring of `ORF1ab`, so a substring match files ORF1ab's codon numbers
    # under ORF1a (#114). Without this case the whole-component match is
    # uncovered -- relaxing it to `masking_gene not in gene_component` survives
    # every other check here, and lands in the production tree's largest gene.
    # It cannot be caught by the test config, whose pcp_genes (S, ORF1ab, E,
    # ORF1b) contain no symbol that is a substring of another; `ORF1a` appears
    # only in the production config.
    substring_rows = [
        (500, True, "ORF1a;ORF1ab", "79;79"),
        (13468, True, "ORF1a;ORF1ab;ORF1ab", "4401;4401;4402"),
    ]
    substring_blocks = {"ORF1a": ORF1A_BLOCKS, "ORF1ab": ORF1AB_BLOCKS}
    for gene, expected, note in [
        ("ORF1a", [79, 4401],
         "4402 present means ORF1ab's codon was filed under ORF1a via a "
         "substring match"),
        ("ORF1ab", [79, 4401, 4402],
         "ORF1ab must collect both 4401 and 4402 at the shared frameshift "
         "nucleotide, where it appears twice in one cell"),
    ]:
        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                got_sub = _run_mask_script(
                    tmpdir, gene, substring_rows, gtf_genes=substring_blocks
                )
            except Exception as exc:  # noqa: BLE001 - reported as a failure
                failures.append(f"{gene} mask derivation raised: {exc}")
                continue
        if got_sub != expected:
            failures.append(
                f"{gene}'s mask is {got_sub}, expected {expected} -- {note}"
            )

    # Each of these writes an empty mask at exit 0 if unguarded.
    guarded = {
        # "pcp_gene_cds_blocks" would be too weak a needle: every
        # validate_cds_blocks message opens with `pcp_gene_cds_blocks[...]` and
        # the add-only message names it too, so it cannot tell this error from
        # several others.
        "gene in neither the GTF nor the override": dict(
            gene="ZZZ", rows=[(13442, True, "nan", "")], overrides={},
            must_mention="appears in neither",
        ),
        "override blocks intersect nothing in the grid": dict(
            gene="ORF1b", rows=[(100, True, "nan", ""), (110, True, "nan", "")],
            overrides={"ORF1b": ORF1B_BLOCKS}, must_mention="intersect",
        ),
        "gene and codon_site disagree on component count": dict(
            gene="N", rows=[(28300, True, "N;ORF9b", "4")], overrides={},
            must_mention="components",
        ),
        # A missing codon_site column returned [] and wrote a 0-byte mask at
        # exit 0 for every GTF-carried gene, while its two siblings (missing
        # 'exclude', missing 'gene') raised. build_count_table.py drops a missing
        # column from counts_all.csv silently, so one upstream schema regression
        # would unmask every gene with nothing to attribute it to.
        "counts table has no codon_site column": dict(
            gene="N", rows=[(28300, True, "N", "4")], overrides={},
            must_mention="codon_site", drop_columns=["codon_site"],
        ),
    }
    for label, case in guarded.items():
        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                _run_mask_script(
                    tmpdir, case["gene"], case["rows"],
                    overrides=case["overrides"], gtf_genes=n_blocks,
                    drop_columns=case.get("drop_columns"),
                )
            except (ValueError, KeyError) as exc:
                if case["must_mention"] not in str(exc):
                    failures.append(
                        f"{label}: raised without mentioning "
                        f"{case['must_mention']!r}: {exc}"
                    )
            else:
                failures.append(f"{label}: wrote a mask instead of raising")

    # The dispatch must key off the config, not off whichever genes the table
    # mentions: a gene with an override that the GTF *also* carries must reach
    # resolve_cds_blocks' add-only guard rather than silently taking the GTF path.
    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            _run_mask_script(
                tmpdir, "ORF1ab", [(13442, True, "ORF1ab", "4393")],
                overrides={"ORF1ab": ORF1AB_BLOCKS},
            )
        except ValueError as exc:
            if "add-only" not in str(exc):
                failures.append(
                    f"a gene in both the GTF and the override raised, but not "
                    f"the add-only error: {exc}"
                )
        else:
            failures.append(
                "a gene defined in both the GTF and pcp_gene_cds_blocks silently "
                "took the GTF path, bypassing the add-only guard"
            )

    if not failures:
        print(f"  codon pairing zipped not flattened, and matched on whole "
              f"components (ORF1a != ORF1ab); all {len(guarded)} empty-mask "
              f"routes raise; dispatch keys off the config")
    return failures


CHECKS = [
    ("override validation", check_override_validation),
    ("PCP-only gene gets a non-empty codon mask",
     check_pcp_only_gene_gets_a_nonempty_codon_mask),
    ("mask script guards and cross-gene codon pairing",
     check_mask_script_guards_and_gene_pairing),
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
