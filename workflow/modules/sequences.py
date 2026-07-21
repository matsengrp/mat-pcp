"""Shared sequence-slicing helpers.

Stdlib-only so this module can be loaded (via the by-path
``load_workflow_module`` pattern) from scripts running in any conda env,
including ``bte.yml`` (``pcp_from_tree.py``) and ``default.yml``
(``clade_founder_fasta.py``).

Coordinates are **1-based and inclusive on both ends**, relative to the
Wuhan-Hu-1 reference (NC_045512.2), matching the ``pcp_gene_region`` and GTF
conventions used elsewhere in the pipeline.
"""


def slice_1based_inclusive(seq, start, end):
    """Return the ``[start, end]`` slice of ``seq``, 1-based inclusive on both ends.

    >>> slice_1based_inclusive("ABCDEF", 2, 4)
    'BCD'
    >>> slice_1based_inclusive("ABCDEF", 1, 6)
    'ABCDEF'
    """
    return seq[start - 1 : end]


def slice_cds(seq, cds_blocks):
    """Concatenate 1-based-inclusive ``(start, end)`` blocks, in genome order.

    Single-block genes pass ``[(start, end)]`` and reduce exactly to
    :func:`slice_1based_inclusive`. ORF1ab passes its two GTF CDS blocks, which
    share one nucleotide (the -1 ribosomal frameshift): that nucleotide is read
    twice so the concatenated coding sequence stays in frame.

    >>> slice_cds("ABCDEF", [(2, 4)])          # single block == plain slice
    'BCD'
    >>> slice_cds("ABCDEF", [(1, 2), (2, 3)])  # shared nt read twice
    'ABBC'
    >>> slice_cds("ABCDEFGH", [(1, 2), (4, 5), (7, 8)])  # three blocks
    'ABDEGH'
    """
    return "".join(slice_1based_inclusive(seq, start, end) for start, end in cds_blocks)
