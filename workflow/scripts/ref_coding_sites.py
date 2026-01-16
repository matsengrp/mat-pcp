#!/usr/bin/env python3
"""Get all sites in the GTF that correspond to a coding sequence."""

import Bio.SeqIO
import pandas as pd

# Get inputs from snakemake
gtf_file = snakemake.input.gtf
fasta_file = snakemake.input.fasta
output_csv = snakemake.output.csv
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log(f"Loading GTF from {gtf_file}")

# Get CDS features
cds_features = (
    pd.read_csv(
        gtf_file,
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
    .assign(
        frame=lambda x: x["frame"].astype(int),
        gene=lambda x: x["attribute"].str.extract(r'^gene_id "(\w+)"'),
        genelength=lambda x: x["end"] - x["start"] + 1,
        shift=lambda x: x.groupby("gene")["genelength"].shift(fill_value=0).astype(int),
    )
    .drop(columns=["seqName", "source", "score", "feature", "attribute"])
)

log(f"Found {len(cds_features)} CDS features")
log(f"Genes: {sorted(cds_features['gene'].unique())}")

# Validate CDS features
assert all(cds_features["frame"] == 0), (
    f"Non frame-0 CDSs found:\n{cds_features[cds_features['frame'] != 0]}"
)
assert all(cds_features["strand"] == "+"), (
    f"Non-positive-sense CDSs found:\n{cds_features[cds_features['strand'] != '+']}"
)

# Get sites with codon position assigned
log("Mapping sites to codon positions...")
sites = (
    cds_features.assign(
        site=lambda x: x.apply(lambda r: list(range(r["start"], r["end"] + 1)), axis=1),
    )
    .explode("site", ignore_index=True)
    .assign(
        codon_position=lambda x: x.apply(
            lambda r: (r["site"] - r["start"]) % 3 + 1,
            axis=1,
        ),
        codon_site=lambda x: x.apply(
            lambda r: (r["site"] - r["start"] + r["shift"]) // 3 + 1,
            axis=1,
        ),
    )
    .drop(columns=["strand", "frame", "start", "end"])
    .groupby("site", as_index=False)
    .aggregate(
        codon_position=pd.NamedAgg(
            "codon_position", lambda s: ";".join(map(str, s.values))
        ),
        codon_site=pd.NamedAgg("codon_site", lambda s: ";".join(map(str, s.values))),
        gene=pd.NamedAgg("gene", lambda s: ";".join(s)),
        genelength=pd.NamedAgg("genelength", lambda s: ";".join(map(str, s.values))),
        shift=pd.NamedAgg("shift", lambda s: ";".join(map(str, s.values))),
    )
)

log(f"Mapped {len(sites)} coding sites")

# Add non-coding sites
log(f"Loading reference FASTA from {fasta_file}")
ref_record = Bio.SeqIO.read(fasta_file, "fasta")
length = len(ref_record)
log(f"Reference length: {length} bp")

noncoding_sites = sorted(set(range(1, length + 1)) - set(sites["site"]))
log(f"Non-coding sites: {len(noncoding_sites)}")

sites = pd.concat(
    [
        sites,
        pd.DataFrame(
            {
                "site": noncoding_sites,
                "codon_position": "noncoding",
                "codon_site": "noncoding",
                "gene": "noncoding",
                "genelength": "noncoding",
                "shift": "noncoding",
            }
        ),
    ]
).sort_values("site")

# Select final columns
sites = sites[["site", "codon_position", "codon_site", "gene"]]

sites.to_csv(output_csv, index=False)
log(f"Saved {len(sites)} sites to {output_csv}")
