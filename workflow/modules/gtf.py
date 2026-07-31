"""GTF helpers: derive per-gene CDS block coordinates.

Loaded via the by-path ``load_workflow_module`` pattern from scripts running in
both ``default.yml`` and ``bte.yml``; pandas is available in both.
"""

import pandas as pd


def _read_cds(gtf_path):
    """Parse ``gtf_path``'s CDS features into a frame with a ``gene`` column.

    Mirrors the CDS parsing in ``ref_coding_sites.py`` (filter ``feature ==
    'CDS'``, extract ``gene_id``).
    """
    return (
        pd.read_csv(
            gtf_path,
            sep="\t",
            header=None,
            comment="#",
            names=[
                "seqName",
                "source",
                "feature",
                "start",
                "end",
                "score",
                "strand",
                "frame",
                "attribute",
            ],
        )
        .query("feature == 'CDS'")
        .assign(gene=lambda x: x["attribute"].str.extract(r'^gene_id "(\w+)"'))
    )


def gtf_cds_genes(gtf_path):
    """Gene symbols with at least one CDS feature in ``gtf_path``."""
    return set(_read_cds(gtf_path)["gene"].dropna())


def cds_blocks_for_gene(gtf_path, gene):
    """Return the CDS blocks for ``gene`` as 1-based-inclusive ``(start, end)``
    tuples, sorted by genomic start.

    Single-block genes return one tuple; ORF1ab returns two blocks that share
    nucleotide 13468 (the -1 ribosomal frameshift), e.g.
    ``[(266, 13468), (13468, 21552)]``.
    """
    cds = _read_cds(gtf_path)
    blocks = cds[cds["gene"] == gene].sort_values("start")
    result = [(int(s), int(e)) for s, e in zip(blocks["start"], blocks["end"])]
    if not result:
        raise ValueError(f"No CDS blocks found for gene {gene!r} in {gtf_path}")
    return result


def validate_cds_blocks(gene, blocks, genome_length=None):
    """Validate explicit CDS blocks for ``gene``, returning them as tuples.

    Every check names ``gene``, because a mis-specified block does not fail
    loudly downstream -- it produces a frame-shifted protein that still looks
    like valid output, which is the failure mode this repository keeps hitting
    (#129, #130, #133). Blocks must be 1-based inclusive, ascending, within the
    genome, and total a whole number of codons.

    ``genome_length`` is optional only because not every caller has the
    reference to hand; supply it whenever possible.
    """
    if not isinstance(blocks, (list, tuple)) or not blocks:
        raise ValueError(
            f"pcp_gene_cds_blocks[{gene!r}] must be a non-empty list of "
            f"[start, end] pairs; got {blocks!r}"
        )
    parsed = []
    for block in blocks:
        if not isinstance(block, (list, tuple)) or len(block) != 2:
            raise ValueError(
                f"pcp_gene_cds_blocks[{gene!r}] entries must be [start, end] "
                f"pairs; got {block!r}"
            )
        start, end = block
        if isinstance(start, bool) or isinstance(end, bool):
            raise ValueError(
                f"pcp_gene_cds_blocks[{gene!r}] coordinates must be integers; "
                f"got {block!r}"
            )
        if not isinstance(start, int) or not isinstance(end, int):
            raise ValueError(
                f"pcp_gene_cds_blocks[{gene!r}] coordinates must be integers; "
                f"got {block!r}"
            )
        if start < 1:
            raise ValueError(
                f"pcp_gene_cds_blocks[{gene!r}] is 1-based, so start must be "
                f">= 1; got {block!r}"
            )
        if end < start:
            raise ValueError(
                f"pcp_gene_cds_blocks[{gene!r}] blocks must have end >= start; "
                f"got {block!r}"
            )
        if genome_length is not None and end > genome_length:
            raise ValueError(
                f"pcp_gene_cds_blocks[{gene!r}] block {block!r} ends past the "
                f"reference length {genome_length}"
            )
        parsed.append((start, end))

    starts = [start for start, _ in parsed]
    if starts != sorted(starts):
        raise ValueError(
            f"pcp_gene_cds_blocks[{gene!r}] must be in ascending genome order; "
            f"got {parsed!r}"
        )

    # Blocks may deliberately share a nucleotide (ORF1ab's and ORF1b's -1
    # ribosomal frameshift reads nt 13468 twice), so the total length is the
    # sum over blocks and NOT the span -- that shared nucleotide is exactly why
    # the divisibility check has to be on the sum.
    total = sum(end - start + 1 for start, end in parsed)
    if total % 3 != 0:
        raise ValueError(
            f"pcp_gene_cds_blocks[{gene!r}] totals {total} nt, not a whole "
            f"number of codons (blocks {parsed!r}). A frame-shifted slice "
            f"still translates, so this cannot be caught downstream."
        )
    return parsed


def resolve_cds_blocks(gtf_path, gene, overrides=None, genome_length=None):
    """Return ``gene``'s CDS blocks, preferring an explicit config override.

    ``overrides`` is the ``pcp_gene_cds_blocks`` config mapping. It exists so a
    PCP-only gene can be given coordinates the shared reference GTF does not
    carry -- ORF1b (nsp12-nsp16) is the motivating case. It is deliberately
    *not* fed to ``edit_ref_gtf`` / ``ref_coding_sites``, so fitness outputs are
    untouched: adding ORF1b to the GTF would make every site in an 8 kb stretch
    multi-gene and, absent a matching ``gene_overlaps.retain`` entry, would
    silently drop all of nsp12-nsp16 from ``aamut_fitness`` (see #133).

    The override is **add-only**: a gene that already has CDS features in the
    GTF raises rather than being silently shadowed. Overriding a GTF gene would
    mean the PCP slice and every fitness output disagreed about the same gene's
    coordinates with nothing to flag it, and nothing needs that today.
    """
    overrides = overrides or {}
    if gene not in overrides:
        return cds_blocks_for_gene(gtf_path, gene)
    blocks = validate_cds_blocks(gene, overrides[gene], genome_length)
    if gene in gtf_cds_genes(gtf_path):
        raise ValueError(
            f"pcp_gene_cds_blocks defines {gene!r}, but {gene!r} already has "
            f"CDS features in {gtf_path}. The override is add-only: it exists "
            f"to name genes the GTF lacks, not to shadow ones it has, which "
            f"would leave the PCP slice and the fitness outputs disagreeing "
            f"about the same gene."
        )
    return blocks
