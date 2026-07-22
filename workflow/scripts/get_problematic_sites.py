#!/usr/bin/env python3
"""Identify problematic codon sites in a specific gene for masking."""

import pandas as pd

# Get snakemake variables
expected_vs_actual_csv = snakemake.input.expected_vs_actual
output_file = snakemake.output.sites_to_mask
masking_gene = snakemake.wildcards.gene
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log(f"Identifying problematic codon sites in gene: {masking_gene}")

# Load expected vs actual counts
expected_v_actual = pd.read_csv(expected_vs_actual_csv)
log(f"Loaded {len(expected_v_actual)} count entries")

# Find problematic codon sites to mask in the target gene.
#
# A codon site is masked only if it is excluded across ALL clades -- i.e. every
# row at that ``nt_site`` has ``exclude == True``. This keeps the
# clade-independent masks (out of ``site_include_range``, config
# ``sites_to_exclude``, and ProblematicSites VCF sites -- each of which sets
# ``exclude`` on every row at the site) and drops the clade-specific ones (UShER
# per-clade ``sites``/``mutations`` masks and founder->reference reversions),
# which should not clade-specifically mask positions of the target gene.
if "exclude" not in expected_v_actual.columns:
    raise KeyError("counts table is missing the required 'exclude' column")

excluded = expected_v_actual["exclude"] == True  # noqa: E712
site_all_excluded = excluded.groupby(expected_v_actual["nt_site"]).transform("all")

mask = site_all_excluded
if "gene" in expected_v_actual.columns:
    mask = mask & expected_v_actual["gene"].str.contains(masking_gene, na=False)

problematic_sites = expected_v_actual[mask]

# Get unique codon sites
codon_sites = []
if "codon_site" in problematic_sites.columns:
    raw_sites = problematic_sites["codon_site"].unique()
    # Handle potential semicolon-separated values (e.g., "123;456")
    for site in raw_sites:
        if pd.isna(site):
            continue
        site_str = str(site)
        if ";" in site_str:
            for s in site_str.split(";"):
                s = s.strip()
                if s.replace(".", "").replace("-", "").isdigit():
                    codon_sites.append(int(float(s)))
        elif site_str.replace(".", "").replace("-", "").isdigit():
            codon_sites.append(int(float(site_str)))

codon_sites = sorted(set(codon_sites))

log(f"Found {len(codon_sites)} problematic codon sites")

# Write output
with open(output_file, "w") as f:
    f.write("\n".join([str(c) for c in codon_sites]))

log(f"Wrote {len(codon_sites)} sites to {output_file}")
log("Done!")
