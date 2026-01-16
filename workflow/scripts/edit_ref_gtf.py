#!/usr/bin/env python3
"""Edit the reference GTF with manual additions for missing genes."""

import csv
import pandas as pd

# Get inputs from snakemake
input_gtf = snakemake.input.gtf
output_gtf = snakemake.output.gtf
edits = snakemake.params.edits
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log(f"Loading GTF from {input_gtf}")

gtf = pd.read_csv(
    input_gtf,
    sep="\t",
    header=None,
    comment="#",
    names=[
        "seqname",
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

log(f"Loaded {len(gtf)} GTF entries")

seqname = gtf["seqname"].unique()
assert len(seqname) == 1, f"Expected single seqname, got: {seqname}"
seqname = seqname[0]

source = gtf["source"].unique()
assert len(source) == 1, f"Expected single source, got: {source}"
source = source[0]

# Verify transcript starts are sorted
transcript_starts = gtf.query("feature == 'transcript'")["start"].tolist()
assert transcript_starts == sorted(transcript_starts), "transcript starts not sorted"

# Add each gene from edits
for gene, coords in edits.items():
    start, end = coords
    log(f"Adding gene {gene} at positions {start}-{end}")

    after_gtf = gtf.query("feature == 'transcript'").query("start > @start")
    if len(after_gtf):
        i = after_gtf.index[0]
    else:
        i = len(gtf)

    records = []
    for feature in ["transcript", "exon", "CDS", "start_codon", "stop_codon"]:
        # Calculate coordinates based on feature type
        if feature == "stop_codon":
            feat_start = end - 2
            feat_end = end
        elif feature == "start_codon":
            feat_start = start
            feat_end = start + 2
        elif feature == "CDS":
            feat_start = start
            feat_end = end - 3  # Exclude stop codon from CDS
        else:  # transcript, exon
            feat_start = start
            feat_end = end

        # Build attribute string
        if feature == "transcript":
            attr = f'gene_id "{gene}"; transcript_id "{gene}";'
        else:
            attr = (
                f'gene_id "{gene}"; transcript_id "{gene}"; '
                f'exon_number "1"; exon_id "{gene}";'
            )

        records.append(
            (
                seqname,
                source,
                feature,
                feat_start,
                feat_end,
                ".",
                "+",
                "." if feature in ["transcript", "exon"] else 0,
                attr,
            )
        )

    add_gtf = pd.DataFrame.from_records(records, columns=gtf.columns)
    gtf = pd.concat([gtf.loc[: i - 1], add_gtf, gtf.loc[i:]], ignore_index=True)
    log(f"  Added {len(records)} features for {gene}")

# Save output
gtf.to_csv(output_gtf, index=False, header=None, sep="\t", quoting=csv.QUOTE_NONE)
log(f"Saved edited GTF with {len(gtf)} entries to {output_gtf}")
