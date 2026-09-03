#!/usr/bin/env python3
"""Unit tests for the per-clade DASM censor table (#144).

**The CI dataset produces an empty table, so this file is the only non-vacuous
coverage of the feature.** ``config_test_data.yaml``'s ``max_clades: 3`` gives
19A/19B/20A and the earliest clade in ``data/usher_masked_sites.yaml`` is 20B, so
a green pipeline run asserts nothing about the coordinates, the directionality or
the guards -- the same situation CLAUDE.md patterns 8, 10 and 11 document for
theirs. In particular no yaml entry falls anywhere in ORF1b's leading block
13442-13468, and 13468 is masked in no clade, so **the frameshift junction case
cannot be exercised by any pipeline run** and is synthetic here.

Runs the real ``workflow/scripts/censor_table.py`` against a stub ``snakemake``
global, never a re-implementation of it. CLAUDE.md pattern 10 records why: an
earlier test of the since-retired ``get_problematic_sites`` re-implemented that
script's nt_site -> codon inversion and asserted against its own copy, so three
separate breaks of the shipped script all left it green. The mutation tests below go
further and exec *deliberately corrupted* copies of the shipped script text,
asserting the guards catch each corruption -- a guard that never fires on a
broken input is indistinguishable from no guard.

Tests are hermetic: a deterministic synthetic genome, synthetic founder slices,
and either the real committed yaml or a synthetic one. Nothing reads a pipeline
output.
"""

import csv
import shutil
import sys
import tempfile
import traceback
from importlib import util as importlib_util
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CENSOR_SCRIPT = REPO_ROOT / "workflow" / "scripts" / "censor_table.py"
USHER_YAML = REPO_ROOT / "data" / "usher_masked_sites.yaml"

# Real coordinates, as verified in #133 against the shipped founders.
S_BLOCKS = [(21563, 25381)]
ORF1AB_BLOCKS = [(266, 13468), (13468, 21552)]
ORF1B_BLOCKS = [(13442, 13468), (13468, 21552)]
GENOME_LENGTH = 29903

# 21M is the T1/T4 witness. Its S-range site masks are 22289-22294 (its other
# site, 28254, is in N and must not appear in S's table); its S-range directional
# masks are T22329C, C22629A and T25000C. 22629 is NOT in 21M's `sites` list, so
# the two mask kinds are disjoint for this clade -- which is what makes it a
# clean witness for "a directional mask was not collapsed into a site mask".
WITNESS_CLADE = "21M"
WITNESS_SITES = [22289, 22290, 22291, 22292, 22293, 22294]
WITNESS_MUTATION = ("C", 22629, "A")


def _load(module_name):
    """Load a workflow module by path, as the pipeline scripts do."""
    spec = importlib_util.spec_from_file_location(
        module_name, REPO_ROOT / "workflow" / "modules" / f"{module_name}.py"
    )
    module = importlib_util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


sequences = _load("sequences")
genemap = _load("genemap")
masking = _load("masking")


class _Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def synthetic_genome(length=GENOME_LENGTH):
    """A deterministic pseudo-random genome.

    The round-trip identity these tests pin is about *coordinates*, not bases, so
    a synthetic genome makes the suite hermetic without weakening anything. The
    formula is deliberately trivial rather than seeded ``random``, so the fixture
    cannot drift with a Python release.
    """
    bases = "ACGT"
    return "".join(bases[(i * 7 + 3) % 4] for i in range(length))


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


def _write_fasta(path, name, seq):
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i : i + 60] + "\n")
    return path


def _founder_for(genome, blocks, yaml_path, clade, orientation="from"):
    """A clade's gene founder slice, with each directional mask's base forced.

    ``orientation='from'`` sets the founder base to the mask's ``from_nt``, which
    is the ordinary case (the mutation runs away from the founder).
    ``orientation='to'`` sets it to ``to_nt`` for every mask, which is the
    signature of a from_nt/to_nt swap and must trip the orientation guard.
    """
    _, clade_mutations = masking.load_usher_masks(yaml_path)
    seq = list(sequences.slice_cds(genome, blocks))
    for mutation in clade_mutations.get(clade, ()):
        from_nt, site, to_nt = mutation[0], int(mutation[1:-1]), mutation[-1]
        for index in sequences.genomic_to_pcp_indices(blocks, site):
            seq[index] = from_nt if orientation == "from" else to_nt
    return "".join(seq)


def run_censor_table(
    tmpdir,
    gene,
    blocks,
    clades,
    yaml_path=USHER_YAML,
    genome=None,
    overrides=None,
    gtf_genes=None,
    orientation="from",
    script_text=None,
):
    """Execute censor_table.py the way Snakemake does; return (rows, roster, log).

    exec rather than import, because the script runs top-level statements against
    the ``snakemake`` global -- exactly how Snakemake invokes it. ``script_text``
    lets a caller exec a deliberately corrupted copy for the mutation tests.
    """
    tmpdir = Path(tmpdir)
    genome = genome if genome is not None else synthetic_genome()

    gtf_path = _synthetic_gtf(
        tmpdir / "ref.gtf", gtf_genes if gtf_genes is not None else {gene: blocks}
    )
    ref_path = _write_fasta(tmpdir / "ref.fa", "NC_045512.2", genome)

    founder_dir = tmpdir / "founders"
    founder_dir.mkdir(exist_ok=True)
    founder_fastas = {}
    for clade in clades:
        path = founder_dir / f"{clade}.fa"
        _write_fasta(
            path, clade, _founder_for(genome, blocks, yaml_path, clade, orientation)
        )
        founder_fastas[clade] = str(path)

    out_table = tmpdir / "censored_mutations.csv"
    out_roster = tmpdir / "censored_mutations.clades.txt"
    log_path = tmpdir / f"{gene}.log"

    stub = _Namespace(
        input=_Namespace(
            usher_masked_sites=str(yaml_path),
            gtf=str(gtf_path),
            ref_fasta=str(ref_path),
            founders=list(founder_fastas.values()),
        ),
        output=_Namespace(table=str(out_table), clades=str(out_roster)),
        params=_Namespace(
            clades=list(clades),
            founder_fastas=founder_fastas,
            pcp_gene_cds_blocks=overrides or {},
        ),
        wildcards=_Namespace(gene=gene),
        log=[str(log_path)],
    )
    script_globals = {
        "__name__": "__main__",
        "__file__": str(CENSOR_SCRIPT),
        "snakemake": stub,
    }
    text = script_text if script_text is not None else CENSOR_SCRIPT.read_text()
    exec(compile(text, str(CENSOR_SCRIPT), "exec"), script_globals)

    with open(out_table, newline="") as f:
        rows = list(csv.DictReader(f))
    roster = [line for line in out_roster.read_text().splitlines() if line]
    log_text = log_path.read_text() if log_path.exists() else ""
    return rows, roster, log_text


def _rows_for(rows, clade):
    return [r for r in rows if r["clade"] == clade]


# ---------------------------------------------------------------------------
# T1 -- known answer over a single block
# ---------------------------------------------------------------------------


def check_single_block_known_answer():
    """21M's six S-range site masks land at the arithmetically known indices.

    Also computes the codon number two independent ways -- ``nt_index // 3 + 1``
    and ``genemap.codon_genomic_positions``, the derivation the retired
    ``get_problematic_sites.py`` used for a PCP-only gene -- because a
    1-vs-0-based slip in either direction produces a plausible number on its
    own.
    """
    failures = []
    tmpdir = tempfile.mkdtemp()
    try:
        rows, _, _ = run_censor_table(tmpdir, "S", S_BLOCKS, [WITNESS_CLADE])

        site_rows = [r for r in _rows_for(rows, WITNESS_CLADE) if r["from_nt"] == ""]
        got = [(int(r["nt_site"]), int(r["nt_index"])) for r in site_rows]
        expected = [(site, site - 21563) for site in WITNESS_SITES]
        if got != expected:
            failures.append(
                f"21M's S site masks should map {expected}, got {got}"
            )

        # 22289 is the first nucleotide of its codon (726 % 3 == 0) and 22291 the
        # last (728 % 3 == 2), so this covers both ends of one codon.
        if 726 % 3 != 0 or 728 % 3 != 2:
            failures.append(
                "the witness sites no longer straddle a codon; pick new ones"
            )

        codons = genemap.codon_genomic_positions(S_BLOCKS)
        for site, index in expected:
            by_arithmetic = index // 3 + 1
            by_position = [
                number
                for number, positions in enumerate(codons, start=1)
                if site in positions
            ]
            if by_position != [by_arithmetic]:
                failures.append(
                    f"site {site}: nt_index//3+1 gives codon {by_arithmetic} but "
                    f"codon_genomic_positions gives {by_position}"
                )
        if 726 // 3 + 1 != 243:
            failures.append("site 22289 should be S codon 243")

        # 28254 is 21M's other masked site and lives in N, not S.
        if any(int(r["nt_site"]) == 28254 for r in rows):
            failures.append("N's site 28254 leaked into S's censor table")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return failures


# ---------------------------------------------------------------------------
# T2 -- the ORF1b frameshift junction is one-to-many
# ---------------------------------------------------------------------------


# Carries a directional mask on the shared nucleotide as well as whole-site
# masks across it. The two emit loops in censor_table.py are separate blocks of
# code with their own `for nt_index in indices`, so a yaml with only `sites`
# pins the one-to-many behaviour for whole-site masks and leaves the directional
# loop free to collapse to a single index -- which would silently censor one of
# the two codons the shared nucleotide sits in.
SYNTHETIC_JUNCTION_YAML = """\
junction_block:
  clades:
    - 21M
  sites:
    - 13467
    - 13468
    - 13469
  mutations:
    - A13468G
"""

# (nt_site, nt_index) for the whole-site masks and for the directional mask.
# 13468 is read at BOTH index 26 and index 27 -- ORF1b's two blocks share it so
# the -1 frameshift stays in frame -- and the two indices fall in codons 9 and
# 10, so collapsing either loop to one index loses a codon.
JUNCTION_SITE_ROWS = [(13467, 25), (13468, 26), (13468, 27), (13469, 28)]
JUNCTION_MUT_ROWS = [(13468, 26), (13468, 27)]


def check_frameshift_junction_emits_two_rows():
    """A mask on the shared nucleotide occupies TWO positions of the PCP sequence.

    ORF1b's blocks share nt 13468 so the -1 frameshift stays in frame, and
    ``slice_cds`` reads it twice -- at index 26 and index 27, which fall in
    *different* codons (9 and 10). A flat ``nt_site - cds_start`` offset yields
    three rows, one of them wrong, and is wrong for every position past the
    junction. Synthetic because no yaml entry falls in 13442-13468 and 13468 is
    masked in no clade, so no pipeline run can reach this.
    """
    failures = []
    tmpdir = tempfile.mkdtemp()
    try:
        yaml_path = Path(tmpdir) / "junction.yaml"
        yaml_path.write_text(SYNTHETIC_JUNCTION_YAML)

        rows, _, _ = run_censor_table(
            tmpdir,
            "ORF1b",
            ORF1B_BLOCKS,
            [WITNESS_CLADE],
            yaml_path=yaml_path,
            overrides={"ORF1b": [list(b) for b in ORF1B_BLOCKS]},
            gtf_genes={"ORF1ab": ORF1AB_BLOCKS},
        )

        if len(rows) != len(JUNCTION_SITE_ROWS) + len(JUNCTION_MUT_ROWS):
            failures.append(
                f"three masked sites plus one directional mask spanning the "
                f"frameshift junction should emit "
                f"{len(JUNCTION_SITE_ROWS) + len(JUNCTION_MUT_ROWS)} rows "
                f"(13468 is read twice by each loop), got {len(rows)}"
            )

        # Checked per emit loop, not pooled: censor_table.py walks whole-site
        # masks and directional masks in two separate blocks, each with its own
        # index expansion, so a pooled count can stay right while one loop
        # collapses and the other over-emits.
        site_rows = sorted(
            (int(r["nt_site"]), int(r["nt_index"])) for r in rows if not r["from_nt"]
        )
        if site_rows != JUNCTION_SITE_ROWS:
            failures.append(
                f"whole-site junction indices wrong: {site_rows} != "
                f"{JUNCTION_SITE_ROWS}"
            )

        mut_rows = sorted(
            (int(r["nt_site"]), int(r["nt_index"])) for r in rows if r["from_nt"]
        )
        if mut_rows != JUNCTION_MUT_ROWS:
            failures.append(
                f"directional junction indices wrong: {mut_rows} != "
                f"{JUNCTION_MUT_ROWS}. A directional mask on the shared "
                f"nucleotide must emit both indices, exactly as a whole-site "
                f"mask does."
            )

        bases = {(r["from_nt"], r["to_nt"]) for r in rows if r["from_nt"]}
        if bases != {("A", "G")}:
            failures.append(
                f"the directional junction rows should carry A->G on both "
                f"indices, got {bases}"
            )

        positions = sequences.cds_nt_positions(ORF1B_BLOCKS)
        if not (positions[26] == positions[27] == 13468):
            failures.append(
                f"ORF1b indices 26 and 27 should both read genomic 13468; got "
                f"{positions[26]} and {positions[27]}"
            )
        if 26 // 3 + 1 == 27 // 3 + 1:
            failures.append(
                "the two indices for the shared nucleotide must land in "
                "different codons (9 and 10)"
            )
        if (26 // 3 + 1, 27 // 3 + 1) != (9, 10):
            failures.append(
                f"shared-nt codons should be 9 and 10, got "
                f"{(26 // 3 + 1, 27 // 3 + 1)}"
            )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return failures


def check_flat_offset_is_caught():
    """Mutation test: reverting to a flat offset must go red, not merely differ.

    The point of the row-count reconciliation guard is that it is *independent*
    of the emit loop. Corrupting the emit loop to take only the first index --
    exactly what a flat ``nt_site - cds_start`` offset computes -- must raise
    rather than write a plausible three-row table.
    """
    failures = []
    tmpdir = tempfile.mkdtemp()
    try:
        yaml_path = Path(tmpdir) / "junction.yaml"
        yaml_path.write_text(SYNTHETIC_JUNCTION_YAML)

        text = CENSOR_SCRIPT.read_text()

        # Both emit loops, separately. They are two distinct blocks of code with
        # their own index expansion, so corrupting only the whole-site one leaves
        # the guard unexercised against a directional mask -- and the directional
        # rows are the half a training loss keys on.
        emit_loops = {
            "whole-site": (
                '        for nt_index in indices:\n'
                '            rows.append((clade, nt_site, nt_index, "", ""))'
            ),
            "directional": (
                "        for nt_index in indices:\n"
                "            rows.append((clade, nt_site, nt_index, from_nt, to_nt))"
            ),
        }
        for label, original in emit_loops.items():
            corrupted = text.replace(
                original, original.replace("indices:", "indices[:1]:")
            )
            if corrupted == text:
                failures.append(
                    f"the mutation test could not find the {label} emit loop to "
                    f"corrupt; it is no longer testing anything. Update the patch."
                )
                continue
            try:
                run_censor_table(
                    tmpdir,
                    "ORF1b",
                    ORF1B_BLOCKS,
                    [WITNESS_CLADE],
                    yaml_path=yaml_path,
                    overrides={"ORF1b": [list(b) for b in ORF1B_BLOCKS]},
                    gtf_genes={"ORF1ab": ORF1AB_BLOCKS},
                    script_text=corrupted,
                )
                failures.append(
                    f"dropping the second index at the frameshift junction in "
                    f"the {label} emit loop did NOT raise; the row-count "
                    f"reconciliation guard is vacuous for that loop"
                )
            except ValueError as exc:
                if "account for" not in str(exc):
                    failures.append(
                        f"flat offset in the {label} emit loop raised, but not "
                        f"from the reconciliation guard: {exc}"
                    )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return failures


# ---------------------------------------------------------------------------
# T3 -- round-trip stated as a sequence identity
# ---------------------------------------------------------------------------


def check_round_trip_is_a_sequence_identity():
    """Every emitted row satisfies ``pcp_seq[nt_index] == genome[nt_site - 1]``.

    Stated as a *sequence* identity rather than index arithmetic, so it also
    catches a blocks-vs-sequence mismatch (the wrong gene's blocks) that
    arithmetic alone would pass. This is also the anti-drift check
    ``sequences.slice_cds`` names in its docstring: the slicer keeps its fast
    per-block string slice while ``cds_nt_positions`` builds the index map a
    different way, and nothing but this pins them together.
    """
    failures = []
    tmpdir = tempfile.mkdtemp()
    genome = synthetic_genome()
    try:
        clade_sites, clade_mutations = masking.load_usher_masks(USHER_YAML)
        clades = sorted(set(clade_sites) | set(clade_mutations))

        for gene, blocks, overrides, gtf_genes in [
            ("S", S_BLOCKS, None, None),
            (
                "ORF1b",
                ORF1B_BLOCKS,
                {"ORF1b": [list(b) for b in ORF1B_BLOCKS]},
                {"ORF1ab": ORF1AB_BLOCKS},
            ),
        ]:
            work = Path(tmpdir) / gene
            work.mkdir()
            rows, _, _ = run_censor_table(
                work,
                gene,
                blocks,
                clades,
                genome=genome,
                overrides=overrides,
                gtf_genes=gtf_genes,
            )
            if not rows:
                failures.append(f"{gene}: no rows, so the round-trip is vacuous")

            pcp_seq = sequences.slice_cds(genome, blocks)
            total = sum(end - start + 1 for start, end in blocks)
            if len(pcp_seq) != total:
                failures.append(
                    f"{gene}: sliced length {len(pcp_seq)} != block total {total}"
                )
            if len(pcp_seq) % 3 != 0:
                failures.append(f"{gene}: sliced length {len(pcp_seq)} is not in frame")

            positions = sequences.cds_nt_positions(blocks)
            if len(positions) != len(pcp_seq):
                failures.append(
                    f"{gene}: cds_nt_positions has {len(positions)} entries but "
                    f"slice_cds produced {len(pcp_seq)} characters -- the slicer "
                    f"and the index map have drifted"
                )
            else:
                drifted = [
                    i
                    for i, (char, site) in enumerate(zip(pcp_seq, positions))
                    if char != genome[site - 1]
                ]
                if drifted:
                    failures.append(
                        f"{gene}: slice_cds and cds_nt_positions disagree at "
                        f"{len(drifted)} indices (first: {drifted[0]})"
                    )

            for row in rows:
                index, site = int(row["nt_index"]), int(row["nt_site"])
                if pcp_seq[index] != genome[site - 1]:
                    failures.append(
                        f"{gene}: row {row} does not round-trip -- "
                        f"pcp_seq[{index}]={pcp_seq[index]!r} vs "
                        f"genome[{site - 1}]={genome[site - 1]!r}"
                    )
                    break
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return failures


# ---------------------------------------------------------------------------
# T4 -- directionality survives
# ---------------------------------------------------------------------------


def check_directionality_survives():
    """A directional mask must not be collapsed into a whole-site mask.

    21M's ``C22629A`` must appear as ``(21M, 22629, 1066, C, A)`` and never as
    ``(A, C)``. 22629 is not in 21M's ``sites`` list, so there must also be no
    blank-``from_nt`` row at index 1066 -- an implementation that widened the
    directional mask to the site would over-censor invisibly, since the row count
    would not change.
    """
    failures = []
    tmpdir = tempfile.mkdtemp()
    try:
        rows, _, _ = run_censor_table(tmpdir, "S", S_BLOCKS, [WITNESS_CLADE])
        from_nt, site, to_nt = WITNESS_MUTATION
        index = site - 21563
        if index != 1066:
            failures.append(f"witness index should be 1066, got {index}")

        at_index = [
            r
            for r in _rows_for(rows, WITNESS_CLADE)
            if int(r["nt_index"]) == index
        ]
        if len(at_index) != 1:
            failures.append(
                f"expected exactly one row at index {index} for 21M, got {at_index}"
            )
        else:
            row = at_index[0]
            if (row["from_nt"], row["to_nt"]) != (from_nt, to_nt):
                failures.append(
                    f"directionality lost: expected ({from_nt}, {to_nt}), got "
                    f"({row['from_nt']}, {row['to_nt']})"
                )
            if (row["from_nt"], row["to_nt"]) == (to_nt, from_nt):
                failures.append("from_nt and to_nt are swapped")
            if row["from_nt"] == "":
                failures.append(
                    f"the directional mask at {site} was collapsed to a "
                    f"whole-site mask -- over-censoring with no change in row count"
                )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return failures


def check_swapped_direction_is_caught():
    """A wholesale from_nt/to_nt swap trips the orientation guard.

    Individually, ``from_nt != founder`` is legitimate: matUtils reports each
    branch's mutations relative to its *parent*, not to the clade founder, so a
    mask can sit on a derived background -- including a reversion back toward the
    founder. Measured on the v11 gisaid-and-public tree, of the 1,126 directional
    masks carried by the 32 clades that have any, 1,083 have founder == from_nt,
    36 have founder == to_nt and 7 have neither. Asserting equality (as #144 originally
    specified) would fail on real data, so the guard is on the *ratio*, which a
    swap inverts wholesale.
    """
    failures = []
    tmpdir = tempfile.mkdtemp()
    try:
        try:
            run_censor_table(
                tmpdir, "S", S_BLOCKS, [WITNESS_CLADE], orientation="to"
            )
            failures.append(
                "founders whose base equals every mask's to_nt did NOT raise; "
                "the orientation guard is vacuous"
            )
        except ValueError as exc:
            if "wrong way round" not in str(exc):
                failures.append(f"swap raised the wrong error: {exc}")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    # ...and the ordinary orientation must NOT raise, or the guard is just noise.
    tmpdir = tempfile.mkdtemp()
    try:
        _, _, log_text = run_censor_table(
            tmpdir, "S", S_BLOCKS, [WITNESS_CLADE], orientation="from"
        )
        if "founder == from_nt" not in log_text:
            failures.append(
                "the founder-orientation breakdown is missing from the log; it "
                "is the diagnostic that makes the ratio guard interpretable"
            )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return failures


# ---------------------------------------------------------------------------
# Roster, emptiness, and the CSV contract
# ---------------------------------------------------------------------------


def check_roster_names_every_clade_including_empty_ones():
    """A clade with nothing censored must still appear in the roster.

    Zero rows and a mistyped clade name are the same thing in the CSV, so the
    consumer cannot implement "unknown clade raises, empty clade passes" from the
    rows alone: ``table[clade]`` raises on every legitimately-empty clade and
    ``table.get(clade, [])`` reinstates the silent failure. This is
    dasm2-experiments#1079's C2 and it is a producer-side requirement.
    """
    failures = []
    tmpdir = tempfile.mkdtemp()
    try:
        clades = ["19A", "19B", "20A", WITNESS_CLADE]
        rows, roster, _ = run_censor_table(tmpdir, "S", S_BLOCKS, clades)
        if roster != sorted(clades):
            failures.append(f"roster should be {sorted(clades)}, got {roster}")
        if any(r["clade"] == "19A" for r in rows):
            failures.append("19A has no UShER mask and must contribute no rows")
        if "19A" not in roster:
            failures.append(
                "19A contributes zero rows, so it is only distinguishable from a "
                "typo via the roster -- and it is missing from it"
            )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return failures


def check_empty_table_is_legal_and_explained():
    """The CI clade set yields an empty table, which must not raise.

    ``config_test_data.yaml``'s ``max_clades: 3`` gives 19A/19B/20A and the
    yaml's earliest clade is 20B, so a "nonzero total" assertion would fail every
    correct CI run. Emptiness must instead be *attributable*: the row count
    reconciles, and the log says why it is zero.
    """
    failures = []
    tmpdir = tempfile.mkdtemp()
    try:
        rows, roster, log_text = run_censor_table(
            tmpdir, "S", S_BLOCKS, ["19A", "19B", "20A"]
        )
        if rows:
            failures.append(f"the CI clade set should censor nothing, got {len(rows)}")
        if roster != ["19A", "19B", "20A"]:
            failures.append(f"roster wrong on the empty case: {roster}")
        if "is empty because" not in log_text:
            failures.append(
                "an empty table must log *why* it is empty; a silent empty file "
                "is indistinguishable from a coordinate error"
            )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return failures


def check_csv_contract_and_reproducibility():
    """Blank cells stay blank under ``keep_default_na=False``, and bytes repeat.

    ``load_usher_masks`` returns *sets*, which are unordered, so without the
    explicit sort the file is not byte-reproducible across runs -- and
    data-central manifests pin md5/xxh3 per file (CLAUDE.md pattern 9). The
    ``keep_default_na`` half is pinned here because this is a cross-repo handoff:
    with pandas' defaults the empty cells become NaN and a ``from_nt == ""`` test
    silently matches nothing.
    """
    failures = []
    try:
        import pandas as pd
    except ImportError:
        return ["pandas is required for the CSV contract check"]

    first = tempfile.mkdtemp()
    second = tempfile.mkdtemp()
    try:
        clades = sorted({WITNESS_CLADE, "22E", "23I", "19A"})
        run_censor_table(first, "S", S_BLOCKS, clades)
        run_censor_table(second, "S", S_BLOCKS, clades)
        a = (Path(first) / "censored_mutations.csv").read_bytes()
        b = (Path(second) / "censored_mutations.csv").read_bytes()
        if a != b:
            failures.append("two runs produced different bytes; the sort is not total")

        table = pd.read_csv(
            Path(first) / "censored_mutations.csv", keep_default_na=False
        )
        if list(table.columns) != ["clade", "nt_site", "nt_index", "from_nt", "to_nt"]:
            failures.append(f"unexpected columns: {list(table.columns)}")

        blank_from = table["from_nt"] == ""
        blank_to = table["to_nt"] == ""
        if not blank_from.any():
            failures.append("no whole-site rows, so the blank-cell contract is vacuous")
        if not (blank_from == blank_to).all():
            failures.append(
                "a row must have both from_nt and to_nt populated, or neither"
            )

        naive = pd.read_csv(Path(first) / "censored_mutations.csv")
        if not naive["from_nt"].isna().any():
            failures.append(
                "the keep_default_na trap this contract documents no longer "
                "exists; re-check the contract note in OUTPUT_SCHEMA.md"
            )

        ordered = sorted(
            zip(table["clade"], table["nt_index"], table["to_nt"])
        )
        if list(zip(table["clade"], table["nt_index"], table["to_nt"])) != ordered:
            failures.append("rows are not sorted by (clade, nt_index, to_nt)")
    finally:
        shutil.rmtree(first, ignore_errors=True)
        shutil.rmtree(second, ignore_errors=True)
    return failures


def check_malformed_masks_raise():
    """An unreadable yaml mask must raise, not silently censor nothing.

    An unparseable entry and an out-of-gene one both contribute zero rows, so
    skipping the first is invisible in the emitted table.
    """
    failures = []
    cases = {
        "unparseable mutation": "b:\n  clades: [21M]\n  mutations: ['G265030']\n",
        "non-ACGT base": "b:\n  clades: [21M]\n  mutations: ['N22629A']\n",
        "from_nt == to_nt": "b:\n  clades: [21M]\n  mutations: ['C22629C']\n",
    }
    for label, body in cases.items():
        tmpdir = tempfile.mkdtemp()
        try:
            yaml_path = Path(tmpdir) / "bad.yaml"
            yaml_path.write_text(body)
            try:
                run_censor_table(
                    tmpdir, "S", S_BLOCKS, [WITNESS_CLADE], yaml_path=yaml_path
                )
                failures.append(f"{label}: did not raise")
            except ValueError:
                pass
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)
    return failures


def check_out_of_gene_masks_are_dropped_not_raised():
    """Most yaml masks lie outside any given gene; that must be routine.

    A gene whose CDS contains no masked site at all still gets a well-formed
    empty table, and the drop counts go to the log so "outside this gene" stays
    distinguishable from "the coordinates are wrong".
    """
    failures = []
    tmpdir = tempfile.mkdtemp()
    try:
        rows, _, log_text = run_censor_table(
            tmpdir, "S", S_BLOCKS, [WITNESS_CLADE]
        )
        if not rows:
            failures.append("21M does have S masks; the fixture is wrong")
        if "Dropped as outside S" not in log_text:
            failures.append("out-of-gene drop counts are missing from the log")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return failures


CHECKS = [
    ("T1 single-block known answer", check_single_block_known_answer),
    ("T2 frameshift junction is one-to-many", check_frameshift_junction_emits_two_rows),
    ("T2 mutation test: flat offset is caught", check_flat_offset_is_caught),
    ("T3 round-trip is a sequence identity", check_round_trip_is_a_sequence_identity),
    ("T4 directionality survives", check_directionality_survives),
    ("T4 mutation test: a swap is caught", check_swapped_direction_is_caught),
    ("roster names every clade", check_roster_names_every_clade_including_empty_ones),
    ("empty table is legal and explained", check_empty_table_is_legal_and_explained),
    ("CSV contract and reproducibility", check_csv_contract_and_reproducibility),
    ("malformed masks raise", check_malformed_masks_raise),
    ("out-of-gene masks are dropped", check_out_of_gene_masks_are_dropped_not_raised),
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


def test_censor_table():
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
        "\nPASSED: the censor table's coordinates round-trip, the frameshift "
        "junction emits both indices, directionality survives, and every guard "
        "fires on a deliberately broken input."
    )


if __name__ == "__main__":
    main()
