"""Shared sequence-slicing helpers.

Stdlib-only so this module can be loaded (via the by-path
``load_workflow_module`` pattern) from scripts running in any conda env,
including ``bte.yml`` (``pcp_from_tree.py``) and ``default.yml``
(``clade_founder_fasta.py``).

Coordinates are **1-based and inclusive on both ends**, relative to the
Wuhan-Hu-1 reference (NC_045512.2), matching the GTF coordinate conventions
used elsewhere in the pipeline.
"""


def slice_1based_inclusive(seq, start, end):
    """Return the ``[start, end]`` slice of ``seq``, 1-based inclusive on both ends.

    >>> slice_1based_inclusive("ABCDEF", 2, 4)
    'BCD'
    >>> slice_1based_inclusive("ABCDEF", 1, 6)
    'ABCDEF'
    """
    return seq[start - 1 : end]


def _check_cds_blocks(cds_blocks):
    """Shared shape validation for the slicer and the two coordinate helpers.

    Deliberately weaker than ``gtf.validate_cds_blocks``: that one enforces the
    frame/overlap rules a *configured* override must satisfy, while this one only
    rejects blocks that cannot be walked at all. Keeping it here means the three
    functions below cannot disagree about what they accept.
    """
    if not cds_blocks:
        raise ValueError("cds_blocks must be non-empty")
    if any(end < start for start, end in cds_blocks):
        raise ValueError(f"each CDS block must have end >= start: {cds_blocks}")
    starts = [start for start, _ in cds_blocks]
    if starts != sorted(starts):
        raise ValueError(f"cds_blocks must be in ascending genome order: {cds_blocks}")


def slice_cds(seq, cds_blocks):
    """Concatenate 1-based-inclusive ``(start, end)`` blocks, in genome order.

    Single-block genes pass ``[(start, end)]`` and reduce exactly to
    :func:`slice_1based_inclusive`. ORF1ab passes its two GTF CDS blocks, which
    share one nucleotide (the -1 ribosomal frameshift): that nucleotide is read
    twice so the concatenated coding sequence stays in frame.

    Implemented as a per-block string slice rather than a gather over
    :func:`cds_nt_positions`, even though the latter would make the slicer and
    the inverter share one line of code. This runs **twice per PCP** (parent and
    child) over millions of PCPs, and the gather is ~900x slower on ORF1ab
    (1.3 us vs 1166 us; ~500x on S) because it builds a 21,555-element index list
    and concatenates character by character. Caching the position list does not
    rescue it (886 us) -- the cost is the per-character join, not the list. The
    anti-drift guarantee is therefore enforced by test rather than by shared
    code: ``tests/test_censor_table.py`` pins
    ``slice_cds(seq, b)[i] == seq[cds_nt_positions(b)[i] - 1]`` over the real
    gene blocks. Do not "simplify" this to a gather.

    >>> slice_cds("ABCDEF", [(2, 4)])          # single block == plain slice
    'BCD'
    >>> slice_cds("ABCDEF", [(1, 2), (2, 3)])  # shared nt read twice
    'ABBC'
    >>> slice_cds("ABCDEFGH", [(1, 2), (4, 5), (7, 8)])  # three blocks
    'ABDEGH'
    """
    _check_cds_blocks(cds_blocks)
    return "".join(slice_1based_inclusive(seq, start, end) for start, end in cds_blocks)


def cds_nt_positions(cds_blocks):
    """Genomic position of every nucleotide of the concatenated CDS, in order.

    Position ``i`` of the list is the reference coordinate that supplies
    character ``i`` of ``slice_cds(genome, cds_blocks)``, so this is the explicit
    index->site map that :func:`slice_cds` applies implicitly. A block pair that
    shares the -1 frameshift nucleotide lists that coordinate **twice**, matching
    the slicer.

    >>> cds_nt_positions([(2, 4)])
    [2, 3, 4]
    >>> cds_nt_positions([(1, 2), (2, 3)])  # shared nt appears twice
    [1, 2, 2, 3]
    """
    _check_cds_blocks(cds_blocks)
    return [pos for start, end in cds_blocks for pos in range(start, end + 1)]


def genomic_to_pcp_indices(cds_blocks, nt_site):
    """0-based indices into the PCP sequence that read genomic ``nt_site``.

    The inverse of :func:`cds_nt_positions`, and the direction that is *not* a
    function: a site outside the CDS yields ``[]``, an ordinary coding site
    yields one index, and the nucleotide shared by two blocks at a -1 ribosomal
    frameshift yields **two**. ORF1b's blocks
    ``[(13442, 13468), (13468, 21552)]`` read genomic 13468 at both index 26 and
    index 27, which land in different codons (9 and 10) -- so a caller that
    assumes one index per site drops a real position, and one that computes a
    flat ``nt_site - cds_start`` offset is wrong for every position past the
    junction (#133).

    Computed from block offsets rather than by searching
    :func:`cds_nt_positions`, so the round-trip test that pins the two against
    each other compares two independent derivations rather than one function
    against itself.

    >>> genomic_to_pcp_indices([(21563, 25381)], 22289)   # S: single block
    [726]
    >>> genomic_to_pcp_indices([(13442, 13468), (13468, 21552)], 13468)
    [26, 27]
    >>> genomic_to_pcp_indices([(21563, 25381)], 100)     # outside the CDS
    []
    """
    _check_cds_blocks(cds_blocks)
    indices = []
    offset = 0
    for start, end in cds_blocks:
        if start <= nt_site <= end:
            indices.append(offset + nt_site - start)
        offset += end - start + 1
    return indices
