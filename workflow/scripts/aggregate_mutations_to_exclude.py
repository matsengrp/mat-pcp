#!/usr/bin/env python3
"""Compile comprehensive exclusion list from all masking sources."""

import itertools

import Bio.SeqIO
import pandas as pd
import yaml

# Get snakemake variables
ref_fasta = snakemake.input.ref_fasta
site_mask_csv = snakemake.input.site_mask
usher_masked_sites_yaml = snakemake.input.usher_masked_sites
clade_counts_csv = snakemake.input.clade_counts
output_csv = snakemake.output.csv

site_include_range = snakemake.params.site_include_range
sites_to_exclude_param = snakemake.params.sites_to_exclude
clade_column = snakemake.params.clade_column

log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


nts = ["A", "C", "G", "T"]

# Get clades from clade_counts file
clades_df = pd.read_csv(clade_counts_csv)
clades = clades_df[clades_df["adequate_sample_counts"] == True][clade_column].tolist()

log(f"Processing exclusions for {len(clades)} clades")

# Get sequence length
seqlength = len(Bio.SeqIO.read(ref_fasta, "fasta"))
log(f"Reference sequence length: {seqlength}")

# Sites out of range
start_range = site_include_range.get("start", 266)
end_range = site_include_range.get("end", 29674)

out_of_range_sites = list(range(1, start_range)) + list(
    range(end_range + 1, seqlength + 1)
)

sites_out_of_range = pd.DataFrame(
    [
        (clade, site, f"{nt1}{site}{nt2}", False)
        for nt1, site, nt2, clade in itertools.product(
            nts, out_of_range_sites, nts, clades
        )
    ],
    columns=["clade", "site", "mutation", "masked_in_usher"],
)
log(f"Sites out of range: {len(sites_out_of_range)} entries")

# Configured sites to exclude
sites_to_exclude_list = sites_to_exclude_param if sites_to_exclude_param else []
sites_to_exclude = pd.DataFrame(
    [
        (clade, site, f"{nt1}{site}{nt2}", False)
        for nt1, site, nt2, clade in itertools.product(
            nts, sites_to_exclude_list, nts, clades
        )
    ],
    columns=["clade", "site", "mutation", "masked_in_usher"],
)
log(f"Configured sites to exclude: {len(sites_to_exclude)} entries")

# Sites from problematic sites VCF
masked_site_list = pd.read_csv(site_mask_csv)["site"].tolist()
sites_masked = pd.DataFrame(
    [
        (clade, site, f"{nt1}{site}{nt2}", True)
        for nt1, site, nt2, clade in itertools.product(
            nts, masked_site_list, nts, clades
        )
    ],
    columns=["clade", "site", "mutation", "masked_in_usher"],
)
log(f"VCF masked sites: {len(sites_masked)} entries")

to_exclude = pd.concat([sites_out_of_range, sites_to_exclude, sites_masked])

# UShER per-clade masked sites
with open(usher_masked_sites_yaml) as f:
    usher_masked_sites = yaml.safe_load(f)

for mask_name, mask_dict in usher_masked_sites.items():
    # Only include clades that are in our processing list
    mask_clades = [c for c in mask_dict.get("clades", []) if c in clades]
    if not mask_clades:
        continue

    mask_df = pd.DataFrame(
        [
            (clade, site, f"{nt1}{site}{nt2}", True)
            for nt1, site, nt2, clade in itertools.product(
                nts, mask_dict.get("sites", []), nts, mask_clades
            )
        ],
        columns=["clade", "site", "mutation", "masked_in_usher"],
    )
    to_exclude = pd.concat([to_exclude, mask_df])
    log(f"UShER mask '{mask_name}': {len(mask_df)} entries")

# Deduplicate and save
result = (
    to_exclude.groupby(["clade", "site", "mutation"], as_index=False)
    .aggregate(masked_in_usher=pd.NamedAgg("masked_in_usher", "any"))
    .sort_values(["clade", "site", "mutation"])
)

result.to_csv(output_csv, index=False)
log(f"Wrote {len(result)} exclusion entries to {output_csv}")
log("Done!")
