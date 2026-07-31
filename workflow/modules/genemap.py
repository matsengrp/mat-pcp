"""Map a PCP gene's codons onto the gene it is validated against.

``validate_pcp_vs_aamut.py`` compares PCP-derived amino-acid counts against
``aamut_fitness``. That comparison needs a codon correspondence, because a PCP
gene need not appear in the fitness table under its own name or numbering:

* ``ORF1a`` is **absent** from ``aamut_fitness`` entirely -- ``aamutfit.get_coding``
  drops it (``aamutfit.py:80``, "Remove ORF1a (keep only ORF1ab)"). Comparing
  ``gene == "ORF1a"`` therefore compared against zero rows and reported no
  discrepancies because it had nothing to compare (#133 Part 2).
* ``ORF1b`` (nsp12-nsp16) is a PCP-only gene by design, an interior slice of
  ``ORF1ab`` starting at ORF1ab codon 4393.

The correspondence is derived from **genomic positions**, not from a declared
codon offset. Two codons correspond iff they occupy the same three reference
nucleotides. That is the definition, and deriving it has one concrete payoff a
bare offset does not give:

``ORF1a``'s GTF block is ``(266, 13480)`` -- 4405 codons -- while ``ORF1ab``'s
first block ends at 13468, 4401 codons in. So ORF1a codons 1-4401 sit on the
same nucleotides as ORF1ab codons 1-4401, but **ORF1a codons 4402-4405 have no
ORF1ab counterpart at all**: they run 12 nt past the -1 ribosomal frameshift in
ORF1a's own frame, whereas ORF1ab's codon 4402 starts at nt 13468 and is in the
shifted frame. Under a bare ``codon_offset: 0`` those four codons would be
mapped onto ORF1ab codons 4402-4405 -- which are *ORF1b's* codons 10-13, wholly
different residues -- and would surface as four spurious discrepancies
indistinguishable from real ones. Deriving the map excludes them by
construction and lets the report name them.

The configured ``codon_offset`` is retained as a **cross-check**: this module
returns the offset it derived, and the caller asserts the two agree. A config
number that no longer matches the coordinates is then a loud failure instead of
a silent mis-comparison.
"""


def codon_genomic_positions(blocks):
    """Genomic nucleotide positions of each codon in ``blocks``, as 3-tuples.

    ``blocks`` are 1-based-inclusive ``(start, end)`` pairs in ascending genome
    order, as returned by ``gtf.resolve_cds_blocks``. Blocks that share a
    nucleotide (the -1 ribosomal frameshift in ORF1ab and ORF1b) read that
    nucleotide twice, matching ``sequences.slice_cds``, so the codon tuples stay
    in frame across the join.

    A trailing partial codon is dropped rather than raising; callers that need
    whole-codon blocks enforce that separately (``gtf.validate_cds_blocks``).

    >>> codon_genomic_positions([(1, 6)])
    [(1, 2, 3), (4, 5, 6)]
    >>> codon_genomic_positions([(1, 2), (2, 4)])  # shared nt 2 read twice
    [(1, 2, 2)]
    >>> codon_genomic_positions([(1, 4)])  # trailing partial codon dropped
    [(1, 2, 3)]
    """
    nts = [pos for start, end in blocks for pos in range(start, end + 1)]
    n_codons = len(nts) // 3
    return [tuple(nts[i * 3 : i * 3 + 3]) for i in range(n_codons)]


def pcp_to_fitness_codons(pcp_blocks, fitness_blocks):
    """Correspondence between a PCP gene's codons and a fitness gene's.

    Returns ``(codon_map, unmapped, offset)`` where

    * ``codon_map`` maps 1-based PCP-gene codon number -> 1-based fitness-gene
      codon number, for codons occupying identical reference nucleotides;
    * ``unmapped`` is the sorted list of PCP-gene codon numbers with no
      counterpart (see the module docstring for ORF1a's four);
    * ``offset`` is the single integer ``fitness_codon - pcp_codon`` shared by
      every mapped pair, or ``None`` if the mapped pairs do not share one.

    ``offset is None`` means the two genes interleave rather than one being a
    contiguous in-frame slice of the other. No configured gene does that, so
    callers should treat it as a configuration error rather than trying to
    represent it.
    """
    fitness_index = {
        positions: number + 1
        for number, positions in enumerate(codon_genomic_positions(fitness_blocks))
    }
    codon_map = {}
    unmapped = []
    for number, positions in enumerate(codon_genomic_positions(pcp_blocks), start=1):
        if positions in fitness_index:
            codon_map[number] = fitness_index[positions]
        else:
            unmapped.append(number)

    offsets = {fitness - pcp for pcp, fitness in codon_map.items()}
    offset = offsets.pop() if len(offsets) == 1 else None
    return codon_map, unmapped, offset


def align_to_fitness_gene(
    pcp_gene, fitness_gene, pcp_blocks, fitness_blocks, declared_offset
):
    """Codon correspondence for a configured ``pcp_gene_to_fitness_gene`` entry.

    Returns ``(codon_map, unmapped)``. Raises if the two genes do not correspond
    in a way a single ``codon_offset`` can describe, or if the offset implied by
    the coordinates disagrees with ``declared_offset``.

    Kept here rather than inline in ``validate_pcp_vs_aamut.py`` so the three
    failure modes are unit-testable; the script itself only reads config and
    reports.
    """
    codon_map, unmapped, derived_offset = pcp_to_fitness_codons(
        pcp_blocks, fitness_blocks
    )
    if not codon_map:
        raise ValueError(
            f"pcp_gene_to_fitness_gene maps {pcp_gene!r} to {fitness_gene!r}, "
            f"but no codon of {pcp_gene} occupies the same reference nucleotides "
            f"as any codon of {fitness_gene}. Blocks: {pcp_blocks} vs "
            f"{fitness_blocks}."
        )
    if derived_offset is None:
        raise ValueError(
            f"{pcp_gene!r} codons do not sit at a single constant offset within "
            f"{fitness_gene!r}; the two genes interleave rather than one being a "
            f"contiguous in-frame slice of the other. A single codon_offset "
            f"cannot describe this."
        )
    # The configured offset is a cross-check on the coordinates, not the source
    # of truth: a stale config number becomes a loud failure instead of a silent
    # comparison shifted by a few codons, which would read as real discrepancies.
    if derived_offset != declared_offset:
        raise ValueError(
            f"pcp_gene_to_fitness_gene[{pcp_gene!r}].codon_offset is "
            f"{declared_offset}, but the CDS blocks put {pcp_gene} codon 1 at "
            f"{fitness_gene} codon {derived_offset + 1} (offset "
            f"{derived_offset}). Fix whichever is wrong; do not adjust one to "
            f"silence the other."
        )
    return codon_map, unmapped


def require_nonempty_fitness_slice(
    pcp_gene, fitness_gene, n_pcp_rows, n_fitness_rows, source
):
    """Raise if a gene with PCP counts has no fitness rows to compare against.

    The anti-vacuity guard (#133). Comparing against an empty slice reports "no
    discrepancies" while detecting nothing, which is how the ORF1a check passed
    through v8. This is the single check that would have caught that gap without
    any of the rest of the machinery, so it is deliberately unconditional.
    """
    if n_pcp_rows > 0 and n_fitness_rows == 0:
        raise ValueError(
            f"No aamut_fitness rows for PCP gene {pcp_gene!r} (looked for "
            f"gene == {fitness_gene!r} in {source}), but {n_pcp_rows} "
            f"PCP-derived count rows exist. An empty slice would make this "
            f"validation vacuous: it would report zero discrepancies having "
            f"compared nothing. Either add a pcp_gene_to_fitness_gene entry "
            f"mapping {pcp_gene!r} onto the gene it appears as, or drop it from "
            f"pcp_genes."
        )
