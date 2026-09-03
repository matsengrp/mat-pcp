#!/usr/bin/env python3
"""Unit tests for the clade-independent nucleotide mask (#153).

The CI pipeline run *does* build this file non-vacuously -- 23 rows for S, 104
for ORF1ab, 39 for ORF1b, 1 for E on the v12 test tree -- so unlike #144's
censor table it is not empty in CI. What a run still cannot see is everything
that makes the coordinates right:

* **The frameshift junction.** No ProblematicSites VCF entry and no config
  ``sites_to_exclude`` entry falls on ORF1b's shared nucleotide 13468, so the
  one place the genomic->index map is one-to-many is unreachable from any
  pipeline run and must be synthetic (the same situation as #144).
* **The ``site_include_range`` reason.** Every shipped config bounds the genome
  at 266-29674, and no PCP gene has a nucleotide outside that, so the third
  reason value is never emitted by any run. It is exercised here on a small
  synthetic gene.
* **Emptiness.** Every configured gene has a non-empty mask on every shipped
  tree, so the header-only case -- legitimate for a short gene, and the
  observable symptom of #133's silent-empty defect one file over -- is only
  reachable on a fixture.
* **Whether the guards fire at all.** All three pass on both v12 production
  trees, which says nothing about what they do when something is wrong.

Runs the real ``workflow/scripts/masked_nt_sites.py`` against a stub
``snakemake`` global, never a re-implementation of it, and the mutation tests
exec deliberately corrupted copies of the shipped script -- or of
``workflow/modules/masking.py``, where the reason vocabulary actually lives, so
the corruption lands in the file an operator would edit rather than in a
paraphrase of it (CLAUDE.md pattern 10).

Hermetic: a deterministic synthetic genome, a synthetic GTF, a synthetic count
table. Nothing reads a pipeline output.
"""

import csv
import shutil
import sys
import tempfile
import traceback
from importlib import util as importlib_util
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT_REL = "scripts/masked_nt_sites.py"
MASKING_REL = "modules/masking.py"

# Real coordinates, as verified in #133 against the shipped founders.
S_BLOCKS = [(21563, 25381)]
ORF1AB_BLOCKS = [(266, 13468), (13468, 21552)]
ORF1B_BLOCKS = [(13442, 13468), (13468, 21552)]
GENOME_LENGTH = 29903
FULL_RANGE = {"start": 266, "end": 29674}

# A short synthetic gene, for the two cases no real gene can reach: the
# site_include_range reason and an empty mask. 21 nt = 7 codons.
TINY_BLOCKS = [(1000, 1020)]


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


def synthetic_genome(length=GENOME_LENGTH):
    """A deterministic pseudo-random genome.

    The round-trip identity these tests pin is about *coordinates*, not bases,
    so a synthetic genome keeps the suite hermetic without weakening anything.
    Deliberately a formula rather than seeded ``random``, so the fixture cannot
    drift with a Python release.
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


def _mutated_workflow(tmp, edits):
    """A throwaway copy of ``workflow/`` with text edits applied.

    Copied rather than edited in place: the scripts resolve
    ``modules/<name>.py`` relative to their own location, so a mutation to a
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


def _write_counts(path, grid_sites, all_excluded, partly_excluded=()):
    """A minimal ``counts_all.csv``: two clades over ``grid_sites``.

    ``all_excluded`` sites carry ``exclude`` in both clades (so
    ``site_all_excluded`` holds); ``partly_excluded`` sites carry it in one clade
    only, which is how a *clade-specific* mask looks and must NOT reach this
    file. The script reads only ``nt_site`` and ``exclude``; ``clade`` is here so
    the fixture has the shape of the real table.
    """
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["clade", "nt_site", "nt_mutation", "exclude"])
        for clade in ("19A", "20A"):
            for site in sorted(grid_sites):
                excluded = site in set(all_excluded) or (
                    site in set(partly_excluded) and clade == "19A"
                )
                writer.writerow([clade, site, f"A{site}G", excluded])
    return path


def run_masked_nt_sites(
    tmpdir,
    gene,
    blocks,
    config_sites=(),
    vcf_sites=(),
    site_include_range=None,
    overrides=None,
    gtf_genes=None,
    grid_sites=None,
    all_excluded=None,
    partly_excluded=(),
    codon_sites=None,
    genome=None,
    workflow_root=None,
):
    """Execute masked_nt_sites.py the way Snakemake does; return (rows, log).

    exec rather than import, because the script runs top-level statements
    against the ``snakemake`` global -- exactly how Snakemake invokes it.
    ``workflow_root`` points the run at a corrupted copy of ``workflow/`` for
    the mutation tests.
    """
    tmpdir = Path(tmpdir)
    tmpdir.mkdir(parents=True, exist_ok=True)
    genome = genome if genome is not None else synthetic_genome()
    site_include_range = site_include_range or dict(FULL_RANGE)
    root = Path(workflow_root) if workflow_root else REPO_ROOT / "workflow"
    script = root / SCRIPT_REL

    gtf_path = _synthetic_gtf(
        tmpdir / "ref.gtf", gtf_genes if gtf_genes is not None else {gene: blocks}
    )
    ref_path = _write_fasta(tmpdir / "ref.fa", "NC_045512.2", genome)

    gene_sites = sorted(set(sequences.cds_nt_positions(blocks)))
    reasons = masking.global_masked_site_reasons(
        gene_sites,
        sites_to_exclude=config_sites,
        vcf_masked_sites=vcf_sites,
        site_include_range=site_include_range,
    )

    # The default grid is the gene's own CDS plus the masked sites, all-excluded
    # exactly where a reason applies -- i.e. a count table that agrees with the
    # config. Every guard-1 check below overrides one of these to make it
    # disagree.
    if grid_sites is None:
        grid_sites = set(gene_sites)
    if all_excluded is None:
        all_excluded = set(reasons)
    counts_path = _write_counts(
        tmpdir / "counts_all.csv", grid_sites, all_excluded, partly_excluded
    )

    with open(tmpdir / "site_mask.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["site"])
        for site in sorted(set(int(s) for s in vcf_sites)):
            writer.writerow([site])

    # The codon file is an *input*, so the default is the answer the shipped
    # get_problematic_sites.py produces for this fixture: widen each masked
    # nucleotide's index to its codon. Derived with the real ``sequences``
    # module, not with a copy of the script's arithmetic.
    if codon_sites is None:
        codon_sites = sorted(
            {
                index // 3 + 1
                for site in reasons
                for index in sequences.genomic_to_pcp_indices(blocks, site)
            }
        )
    codon_path = tmpdir / f"{gene}_codon_sites_to_mask.txt"
    codon_path.write_text("\n".join(str(c) for c in codon_sites))

    out_table = tmpdir / "masked_nt_sites.csv"
    log_path = tmpdir / f"{gene}.log"

    stub = _Namespace(
        input=_Namespace(
            counts=str(counts_path),
            codon_mask=str(codon_path),
            site_mask=str(tmpdir / "site_mask.csv"),
            gtf=str(gtf_path),
            ref_fasta=str(ref_path),
        ),
        output=_Namespace(table=str(out_table)),
        params=_Namespace(
            sites_to_exclude=list(config_sites),
            site_include_range=site_include_range,
            pcp_gene_cds_blocks=overrides or {},
        ),
        wildcards=_Namespace(gene=gene),
        log=[str(log_path)],
    )
    script_globals = {
        "__name__": "__main__",
        "__file__": str(script),
        "snakemake": stub,
    }
    exec(compile(script.read_text(), str(script), "exec"), script_globals)

    with open(out_table, newline="") as f:
        rows = list(csv.DictReader(f))
    log_text = log_path.read_text() if log_path.exists() else ""
    return rows, log_text


def _triples(rows):
    return [(int(r["nt_site"]), int(r["nt_index"]), r["reason"]) for r in rows]


# ---------------------------------------------------------------------------
# T1 -- known answer over a single block
# ---------------------------------------------------------------------------
#
# Three of the five masked sites (22289/22290/22291) sit in ONE codon, which is
# the whole point of the file: at codon granularity they are indistinguishable
# from each other and from the two innocent neighbours of a lone masked
# nucleotide. 23000 carries two reasons, so it emits two rows at one index.
S_CONFIG_SITES = [22289, 22291, 23000]
S_VCF_SITES = [21563, 22290, 23000, 28000]  # 28000 is in N, not S
S_EXPECTED = [
    (21563, 0, "vcf_masked_sites"),
    (22289, 726, "sites_to_exclude"),
    (22290, 727, "vcf_masked_sites"),
    (22291, 728, "sites_to_exclude"),
    (23000, 1437, "sites_to_exclude"),
    (23000, 1437, "vcf_masked_sites"),
]
S_EXPECTED_CODONS = [1, 243, 480]


def check_single_block_known_answer():
    """S's masked sites land at the arithmetically known indices, in order."""
    failures = []
    tmpdir = tempfile.mkdtemp()
    try:
        rows, log_text = run_masked_nt_sites(
            tmpdir,
            "S",
            S_BLOCKS,
            config_sites=S_CONFIG_SITES,
            vcf_sites=S_VCF_SITES,
            codon_sites=S_EXPECTED_CODONS,
        )
        got = _triples(rows)
        if got != S_EXPECTED:
            failures.append(f"expected {S_EXPECTED}, got {got}")

        # The indices are `site - 21563` for a single block; spelled out so a
        # 1-vs-0-based slip cannot pass by agreeing with a helper.
        for site, index, _ in S_EXPECTED:
            if index != site - 21563:
                failures.append(f"site {site} should be index {site - 21563}")

        if any(int(r["nt_site"]) == 28000 for r in rows):
            failures.append("28000 is outside S and must not appear in its mask")

        # Five nucleotides, three codons: the collapse this file exists to undo.
        codons = {int(r["nt_index"]) // 3 + 1 for r in rows}
        if sorted(codons) != S_EXPECTED_CODONS:
            failures.append(
                f"the fixture should widen to {S_EXPECTED_CODONS}, got "
                f"{sorted(codons)} -- 726/727/728 must share codon 243, or the "
                f"nucleotide-vs-codon distinction is not being exercised"
            )
        if len({r["nt_index"] for r in rows}) != 5:
            failures.append("the fixture should mask 5 distinct nt positions")
        if "over-mask" not in log_text:
            failures.append("the over-mask factor is missing from the log")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return failures


def check_multi_reason_site_emits_one_row_per_reason():
    """A site with two reasons emits two rows, and the vocabulary is closed.

    18 of the 121 globally masked sites across the six v12 production genes are
    in both the ProblematicSites VCF and the config list, so a consumer that
    assumes one row per ``nt_index`` over-counts by that much. The alternative
    -- a ``;``-joined ``reason`` cell -- is what CLAUDE.md pattern 10 documents
    the cost of, so the file is long-form and this pins it.
    """
    failures = []
    tmpdir = tempfile.mkdtemp()
    try:
        rows, _ = run_masked_nt_sites(
            tmpdir,
            "S",
            S_BLOCKS,
            config_sites=S_CONFIG_SITES,
            vcf_sites=S_VCF_SITES,
            codon_sites=S_EXPECTED_CODONS,
        )
        at_1437 = [r for r in rows if int(r["nt_index"]) == 1437]
        if len(at_1437) != 2:
            failures.append(
                f"site 23000 is in both sources, so it must emit two rows; got "
                f"{len(at_1437)}"
            )
        if {r["reason"] for r in at_1437} != {
            "sites_to_exclude",
            "vcf_masked_sites",
        }:
            failures.append(f"wrong reasons at 23000: {[r['reason'] for r in at_1437]}")
        if any(";" in r["reason"] for r in rows):
            failures.append("reason must be one value per row, not a joined cell")
        unknown = {r["reason"] for r in rows} - set(masking.GLOBAL_MASK_REASONS)
        if unknown:
            failures.append(f"reasons outside the closed vocabulary: {unknown}")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return failures


# ---------------------------------------------------------------------------
# T2 -- the ORF1b frameshift junction is one-to-many
# ---------------------------------------------------------------------------
#
# 13468 is read at BOTH index 26 and index 27 -- ORF1b's two blocks share it so
# the -1 frameshift stays in frame -- and the two indices fall in codons 9 and
# 10, so collapsing to a single index loses a codon. Synthetic: no VCF or config
# entry falls in 13442-13468, so no pipeline run can reach this.
JUNCTION_VCF_SITES = [13467, 13468, 13469]
JUNCTION_EXPECTED = [
    (13467, 25, "vcf_masked_sites"),
    (13468, 26, "vcf_masked_sites"),
    (13468, 27, "vcf_masked_sites"),
    (13469, 28, "vcf_masked_sites"),
]


def check_frameshift_junction_emits_two_rows():
    """A mask on the shared nucleotide occupies TWO positions of the sequence."""
    failures = []
    tmpdir = tempfile.mkdtemp()
    try:
        rows, _ = run_masked_nt_sites(
            tmpdir,
            "ORF1b",
            ORF1B_BLOCKS,
            vcf_sites=JUNCTION_VCF_SITES,
            overrides={"ORF1b": [list(b) for b in ORF1B_BLOCKS]},
            gtf_genes={"ORF1ab": ORF1AB_BLOCKS},
        )
        got = _triples(rows)
        if got != JUNCTION_EXPECTED:
            failures.append(
                f"three masked sites spanning the frameshift junction should "
                f"emit {JUNCTION_EXPECTED} (13468 is read twice), got {got}"
            )
        codons = sorted({index // 3 + 1 for _, index, _ in got})
        if codons != [9, 10]:
            failures.append(
                f"the junction's two indices must fall in codons 9 and 10, got "
                f"{codons} -- if they shared a codon, dropping one would be "
                f"invisible to the codon cross-check"
            )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return failures


def check_flat_offset_is_caught():
    """A flat ``nt_site - cds_start`` offset must not survive the guards.

    The obvious wrong inversion: it emits one row where the junction needs two,
    and is wrong for every position past the junction. Mutation-tested rather
    than reasoned about, because the correct and incorrect versions agree on
    every single-block gene -- which is five of the six production genes.
    """
    failures = []
    tmp = tempfile.mkdtemp()
    try:
        root = _mutated_workflow(
            tmp,
            {
                SCRIPT_REL: (
                    "for nt_index in sequences.genomic_to_pcp_indices"
                    "(cds_blocks, nt_site):",
                    "for nt_index in [nt_site - cds_blocks[0][0]]:",
                )
            },
        )
        try:
            run_masked_nt_sites(
                Path(tmp) / "run",
                "ORF1b",
                ORF1B_BLOCKS,
                vcf_sites=JUNCTION_VCF_SITES,
                overrides={"ORF1b": [list(b) for b in ORF1B_BLOCKS]},
                gtf_genes={"ORF1ab": ORF1AB_BLOCKS},
                workflow_root=root,
            )
            failures.append(
                "a flat block offset did NOT raise; the junction is unguarded"
            )
        except ValueError as exc:
            # Named, so the check cannot pass on an unrelated failure: the
            # offset emits three rows where four are due, which is the row-count
            # reconciliation's business rather than the codon or round-trip
            # guard's.
            if "frameshift index" not in str(exc):
                failures.append(f"wrong error for a flat offset: {exc}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return failures


def check_round_trip_is_a_sequence_identity():
    """Guard 3 catches a blocks-vs-sequence mismatch index arithmetic passes.

    Mutated so the reference is sliced with the *wrong* coordinates while the
    indices stay right. Every index is still in range and every codon still
    matches the codon file, so guards 1 and 2 pass -- only the string identity
    sees it. That is why the check is written against the sequence rather than
    as ``0 <= nt_index < len(positions)``.
    """
    failures = []
    tmp = tempfile.mkdtemp()
    try:
        root = _mutated_workflow(
            tmp,
            {
                SCRIPT_REL: (
                    "gene_reference = sequences.slice_cds(reference, cds_blocks)",
                    "gene_reference = sequences.slice_cds("
                    "reference, [(1, len(positions))])",
                )
            },
        )
        try:
            run_masked_nt_sites(
                Path(tmp) / "run",
                "S",
                S_BLOCKS,
                config_sites=S_CONFIG_SITES,
                vcf_sites=S_VCF_SITES,
                codon_sites=S_EXPECTED_CODONS,
                workflow_root=root,
            )
            failures.append(
                "slicing the reference at the wrong coordinates did NOT raise; "
                "the round-trip guard is vacuous"
            )
        except ValueError as exc:
            if "round-trip" not in str(exc):
                failures.append(f"wrong error for a blocks mismatch: {exc}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # ...and the correct slicing must NOT raise, or the guard is just noise.
    tmpdir = tempfile.mkdtemp()
    try:
        run_masked_nt_sites(
            tmpdir,
            "S",
            S_BLOCKS,
            config_sites=S_CONFIG_SITES,
            vcf_sites=S_VCF_SITES,
            codon_sites=S_EXPECTED_CODONS,
        )
    except Exception as exc:
        failures.append(f"the correct fixture raised: {exc}")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return failures


# ---------------------------------------------------------------------------
# Guard 2 -- the codon file is exactly this file's codon-widening
# ---------------------------------------------------------------------------


def check_codon_cross_check_fires_both_ways():
    """A codon file that is not this mask's widening raises, either direction.

    This is the check that makes retiring ``{gene}_codon_sites_to_mask.txt``
    later a mechanical step rather than a re-derivation, so it has to fail on a
    stale file. Driven from the *input* rather than by corrupting the script,
    because a stale codon file is the realistic failure: the two rules read the
    same count table but reach codons by different routes.
    """
    failures = []
    cases = {
        "codon file missing a codon": [1, 243],
        "codon file with an extra codon": [1, 243, 480, 999],
        "codon file off by one": [2, 244, 481],
    }
    for label, codon_sites in cases.items():
        tmpdir = tempfile.mkdtemp()
        try:
            try:
                run_masked_nt_sites(
                    tmpdir,
                    "S",
                    S_BLOCKS,
                    config_sites=S_CONFIG_SITES,
                    vcf_sites=S_VCF_SITES,
                    codon_sites=codon_sites,
                )
                failures.append(f"{label}: did not raise")
            except ValueError as exc:
                if "same mask at two granularities" not in str(exc):
                    failures.append(f"{label}: wrong error: {exc}")
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    tmpdir = tempfile.mkdtemp()
    try:
        _, log_text = run_masked_nt_sites(
            tmpdir,
            "S",
            S_BLOCKS,
            config_sites=S_CONFIG_SITES,
            vcf_sites=S_VCF_SITES,
            codon_sites=S_EXPECTED_CODONS,
        )
        if "Codon cross-check" not in log_text:
            failures.append("the codon cross-check is missing from the log")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return failures


# ---------------------------------------------------------------------------
# Guard 1 -- reconciliation against the count table, both directions
# ---------------------------------------------------------------------------


def check_reconciliation_catches_a_count_table_disagreement():
    """The two derivations of "globally masked" must agree, or the run stops.

    Both directions are checked because they mean different things: a site the
    count table excludes everywhere but no global reason explains is masked
    clade-specifically in every clade (which ``reason`` cannot describe, and
    which belongs in ``censored_mutations.csv``), while the reverse means the
    count table was built from a different config or VCF than this rule was
    given. On every shipped tree the two sets are equal -- 628 sites genome-wide
    on both v12 production trees and the test tree -- so neither branch fires in
    any run.
    """
    failures = []

    # Direction 1: the counts table all-excludes a site with no global reason.
    tmpdir = tempfile.mkdtemp()
    try:
        try:
            run_masked_nt_sites(
                tmpdir,
                "S",
                S_BLOCKS,
                config_sites=S_CONFIG_SITES,
                vcf_sites=S_VCF_SITES,
                codon_sites=S_EXPECTED_CODONS,
                all_excluded=set(S_CONFIG_SITES) | {21563, 22290, 23000, 24000},
            )
            failures.append(
                "an all-clade-excluded site matching no global reason did NOT "
                "raise"
            )
        except ValueError as exc:
            if "reason vocabulary cannot express" not in str(exc):
                failures.append(f"wrong error for the missing-reason case: {exc}")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    # Direction 2: a globally masked site the counts table does not all-exclude.
    # Excluded in one clade only, which is exactly how a clade-specific mask
    # looks -- so this also pins that `transform("all")` is not `.any()`.
    tmpdir = tempfile.mkdtemp()
    try:
        try:
            run_masked_nt_sites(
                tmpdir,
                "S",
                S_BLOCKS,
                config_sites=S_CONFIG_SITES,
                vcf_sites=S_VCF_SITES,
                codon_sites=S_EXPECTED_CODONS,
                all_excluded={21563, 22290, 22289, 22291},
                partly_excluded={23000},
            )
            failures.append(
                "a globally masked site that the counts table excludes in only "
                "one clade did NOT raise"
            )
        except ValueError as exc:
            if "different config" not in str(exc):
                failures.append(f"wrong error for the extra-site case: {exc}")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    # Direction 3: a masked site the grid does not mention at all. Separated
    # from the two above because it makes both of them meaningless rather than
    # false, and because for a GTF-carried gene it would also break the codon
    # cross-check (that file is built from the count table's own rows).
    tmpdir = tempfile.mkdtemp()
    try:
        try:
            run_masked_nt_sites(
                tmpdir,
                "S",
                S_BLOCKS,
                config_sites=S_CONFIG_SITES,
                vcf_sites=S_VCF_SITES,
                codon_sites=S_EXPECTED_CODONS,
                grid_sites=set(range(21563, 22000)),
            )
            failures.append("a masked site absent from the count grid did NOT raise")
        except ValueError as exc:
            if "absent from" not in str(exc):
                failures.append(f"wrong error for the off-grid case: {exc}")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return failures


def check_dropping_a_reason_in_the_module_is_caught():
    """Losing one reason in ``masking.py`` is caught by the reconciliation.

    The mutation lands in the module, not in the script: the reason vocabulary
    lives in ``modules/masking.py`` beside ``compute_exclusions``, which is the
    file an operator changing what "globally masked" means would edit. Patching
    the script instead would report the guard as covered while the definition
    sat one file away (CLAUDE.md pattern 15's note on the same shape).
    """
    failures = []
    tmp = tempfile.mkdtemp()
    try:
        root = _mutated_workflow(
            tmp,
            {
                MASKING_REL: (
                    '        if site in config_sites:\n'
                    '            found.append("sites_to_exclude")\n',
                    "",
                )
            },
        )
        try:
            run_masked_nt_sites(
                Path(tmp) / "run",
                "S",
                S_BLOCKS,
                config_sites=S_CONFIG_SITES,
                vcf_sites=S_VCF_SITES,
                codon_sites=S_EXPECTED_CODONS,
                workflow_root=root,
            )
            failures.append(
                "dropping the sites_to_exclude reason did NOT raise; the "
                "reconciliation against counts_all.csv is vacuous"
            )
        except ValueError as exc:
            if "reason vocabulary cannot express" not in str(exc):
                failures.append(f"wrong error for a dropped reason: {exc}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return failures


# ---------------------------------------------------------------------------
# The two cases no shipped config can produce
# ---------------------------------------------------------------------------


def check_site_include_range_reason_is_emitted():
    """The third reason value, which no pipeline run emits.

    Every shipped config bounds the genome at 266-29674 and no PCP gene has a
    nucleotide outside that, so ``site_include_range`` is present in the schema
    and absent from every real output. If the column were assumed to hold only
    two values, a future config narrowing the range would emit rows nobody had
    ever seen.
    """
    failures = []
    tmpdir = tempfile.mkdtemp()
    try:
        rows, _ = run_masked_nt_sites(
            tmpdir,
            "TINY",
            TINY_BLOCKS,
            site_include_range={"start": 1010, "end": 29674},
        )
        got = _triples(rows)
        expected = [(site, site - 1000, "site_include_range") for site in range(1000, 1010)]
        if got != expected:
            failures.append(f"expected {expected}, got {got}")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return failures


def check_empty_mask_is_legal_and_explained():
    """A gene with nothing masked gets a header-only file and a reason why.

    Legitimate for a short gene -- E contributes exactly one site on the
    production trees, so zero is not far off -- but it also means this gene's
    training applies no global mask at all, which is the observable symptom of
    #133's silent-empty defect one file over. Nothing on any shipped tree is
    empty, so a "nonzero" assertion would look fine and this case would be
    untested.
    """
    failures = []
    tmpdir = tempfile.mkdtemp()
    try:
        rows, log_text = run_masked_nt_sites(tmpdir, "TINY", TINY_BLOCKS)
        if rows:
            failures.append(f"an unmasked gene should emit no rows, got {len(rows)}")
        header = (Path(tmpdir) / "masked_nt_sites.csv").read_text().splitlines()
        if header != ["nt_site,nt_index,reason"]:
            failures.append(f"the empty file must still carry its header: {header}")
        if "header-only" not in log_text:
            failures.append(
                "an empty mask must log *why* it is empty; a header-only file is "
                "otherwise indistinguishable from a coordinate error"
            )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return failures


# ---------------------------------------------------------------------------
# The CSV contract
# ---------------------------------------------------------------------------


def check_csv_contract_and_reproducibility():
    """Column order, total sort, and byte-identical repeat runs.

    The sort is required rather than cosmetic: ``global_masked_site_reasons``
    reads *sets*, so without it the file is not byte-reproducible across runs --
    and data-central manifests pin md5/xxh3 per file (CLAUDE.md pattern 9).
    """
    failures = []
    first = tempfile.mkdtemp()
    second = tempfile.mkdtemp()
    try:
        for tmpdir in (first, second):
            run_masked_nt_sites(
                tmpdir,
                "S",
                S_BLOCKS,
                config_sites=S_CONFIG_SITES,
                vcf_sites=S_VCF_SITES,
                codon_sites=S_EXPECTED_CODONS,
            )
        a = (Path(first) / "masked_nt_sites.csv").read_bytes()
        b = (Path(second) / "masked_nt_sites.csv").read_bytes()
        if a != b:
            failures.append("two runs produced different bytes; the sort is not total")

        lines = a.decode().splitlines()
        if lines[0] != "nt_site,nt_index,reason":
            failures.append(f"unexpected header: {lines[0]!r}")
        keys = [
            (int(line.split(",")[1]), line.split(",")[2]) for line in lines[1:]
        ]
        if keys != sorted(keys):
            failures.append(f"rows are not sorted by (nt_index, reason): {keys}")
        if len(set(keys)) != len(keys):
            failures.append(
                "(nt_index, reason) is repeated, so the sort key is not a total "
                "order and the byte-equality above is luck"
            )
    finally:
        shutil.rmtree(first, ignore_errors=True)
        shutil.rmtree(second, ignore_errors=True)
    return failures


def check_global_masked_sites_still_matches_its_old_formula():
    """The delegation in ``masking.py`` did not change the existing set.

    ``global_masked_sites`` feeds the windowed expected-count path and now
    delegates to ``global_masked_site_reasons`` so the three membership tests are
    written once. This pins the refactor against the formula it replaced,
    including the out-of-genome coordinates the reason enumeration cannot reach.
    """
    failures = []
    cases = [
        ([153, 1149], [21563, 29903], {"start": 266, "end": 29674}, 29903),
        ([], [], {"start": 1, "end": 29903}, 29903),
        ([30500], [40000], {"start": 266, "end": 29674}, 29903),
        ([5], [10], {}, 100),
    ]
    for config_sites, vcf_sites, rng, genome_len in cases:
        start = rng.get("start", 266)
        end = rng.get("end", 29674)
        expected = {int(s) for s in config_sites}
        expected.update(int(s) for s in vcf_sites)
        expected.update(s for s in range(1, genome_len + 1) if s < start or s > end)
        got = masking.global_masked_sites(
            sites_to_exclude=config_sites,
            vcf_masked_sites=vcf_sites,
            site_include_range=rng,
            genome_len=genome_len,
        )
        if got != expected:
            failures.append(
                f"global_masked_sites changed for "
                f"({config_sites}, {vcf_sites}, {rng}, {genome_len}): "
                f"{sorted(got ^ expected)[:10]} differ"
            )
    return failures


CHECKS = [
    ("T1 single-block known answer", check_single_block_known_answer),
    ("T1 one row per reason", check_multi_reason_site_emits_one_row_per_reason),
    ("T2 frameshift junction is one-to-many", check_frameshift_junction_emits_two_rows),
    ("T2 mutation test: flat offset is caught", check_flat_offset_is_caught),
    ("T3 round-trip is a sequence identity", check_round_trip_is_a_sequence_identity),
    ("guard 2: codon cross-check fires", check_codon_cross_check_fires_both_ways),
    ("guard 1: count-table reconciliation", check_reconciliation_catches_a_count_table_disagreement),
    ("guard 1 mutation test: a dropped reason is caught", check_dropping_a_reason_in_the_module_is_caught),
    ("site_include_range reason is emitted", check_site_include_range_reason_is_emitted),
    ("empty mask is legal and explained", check_empty_mask_is_legal_and_explained),
    ("CSV contract and reproducibility", check_csv_contract_and_reproducibility),
    ("global_masked_sites is unchanged", check_global_masked_sites_still_matches_its_old_formula),
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


def test_masked_nt_sites():
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
        "\nPASSED: the global nucleotide mask's coordinates round-trip, the "
        "frameshift junction emits both indices, a multi-reason site emits both "
        "reasons, and every guard fires on a deliberately broken input."
    )


if __name__ == "__main__":
    main()
