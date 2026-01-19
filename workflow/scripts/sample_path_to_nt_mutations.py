#!/usr/bin/env python3
"""Convert matUtils sample paths to mutations leading to each node."""

import pandas as pd

# Get inputs from snakemake
input_tsv = snakemake.input.tsv
output_csv = snakemake.output.csv
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log(f"Processing sample paths from {input_tsv}")

# Read sample paths file
# Format: tip_name<tab>path
# Where path is: node1:muts1 node2:muts2 ...
raw_df = pd.read_csv(input_tsv, sep="\t", names=["tip", "path"])
log(f"Read {len(raw_df)} sample paths")

if len(raw_df) == 0:
    log("No sample paths found, creating empty output")
    nodes_and_muts = pd.DataFrame(columns=["node_id", "nt_mutations"])
else:
    # Parse node paths
    nodes_and_muts = (
        raw_df.assign(node_paths=lambda x: x["path"].str.split())
        .drop(columns=["tip", "path"])
        .explode("node_paths")
        .query("node_paths.notnull()")
        .assign(
            node_id=lambda x: x["node_paths"].str.split(":").str[0],
            nt_mutations=lambda x: x["node_paths"].str.split(":").str[1],
        )
        .drop(columns="node_paths")
    )

    log(f"Extracted {len(nodes_and_muts)} node entries")

    # Ensure each node has a unique mutation path
    unique_nodes = nodes_and_muts.drop_duplicates()
    n_unique_nodes = nodes_and_muts["node_id"].nunique()

    if n_unique_nodes != len(unique_nodes):
        log(f"WARNING: {len(unique_nodes) - n_unique_nodes} duplicate node entries")

    # Get just unique entry per node and clean up
    nodes_and_muts = unique_nodes.drop_duplicates(subset=["node_id"]).assign(
        nt_mutations=lambda x: x["nt_mutations"].str.replace(",", ";")
    )

    log(f"Output: {len(nodes_and_muts)} unique nodes")

nodes_and_muts.to_csv(output_csv, index=False)
log(f"Wrote output to {output_csv}")
log("Done!")
