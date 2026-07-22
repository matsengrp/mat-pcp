"""GTF helpers: derive per-gene CDS block coordinates.

Loaded via the by-path ``load_workflow_module`` pattern from scripts running in
both ``default.yml`` and ``bte.yml``; pandas is available in both.
"""

import pandas as pd


def cds_blocks_for_gene(gtf_path, gene):
    """Return the CDS blocks for ``gene`` as 1-based-inclusive ``(start, end)``
    tuples, sorted by genomic start.

    Mirrors the CDS parsing in ``ref_coding_sites.py`` (filter ``feature ==
    'CDS'``, extract ``gene_id``). Single-block genes return one tuple; ORF1ab
    returns two blocks that share nucleotide 13468 (the -1 ribosomal
    frameshift), e.g. ``[(266, 13468), (13468, 21552)]``.
    """
    cds = (
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
    blocks = cds[cds["gene"] == gene].sort_values("start")
    result = [(int(s), int(e)) for s, e in zip(blocks["start"], blocks["end"])]
    if not result:
        raise ValueError(f"No CDS blocks found for gene {gene!r} in {gtf_path}")
    return result
