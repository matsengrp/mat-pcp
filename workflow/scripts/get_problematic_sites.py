#!/usr/bin/env python3
"""Identify problematic codon sites in a specific gene for masking."""

import pandas as pd

# Get snakemake variables
expected_vs_actual_csv = snakemake.input.expected_vs_actual
output_file = snakemake.output.sites_to_mask
masking_gene = snakemake.params.get("masking_gene", "S")
additional_sites_to_ignore = snakemake.params.get("additional_sites_to_ignore", [])
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

# Find problematic codon sites:
# - Sites marked as exclude
# - Sites masked in UShER
# - Sites in the additional ignore list

# Build filter conditions
conditions = []

if "gene" in expected_v_actual.columns:
    conditions.append(expected_v_actual["gene"].str.contains(masking_gene, na=False))

# Get sites that should be masked
mask = conditions[0] if conditions else pd.Series([True] * len(expected_v_actual))

if "exclude" in expected_v_actual.columns:
    mask_exclude = mask & (
        (expected_v_actual["exclude"] == True)  # noqa: E712
        | (expected_v_actual.get("masked_in_usher", pd.Series([False] * len(expected_v_actual))) == True)  # noqa
    )
elif "masked_in_usher" in expected_v_actual.columns:
    mask_exclude = mask & (expected_v_actual["masked_in_usher"] == True)  # noqa
else:
    mask_exclude = pd.Series([False] * len(expected_v_actual))

# Also include additional sites to ignore
if additional_sites_to_ignore and "nt_site" in expected_v_actual.columns:
    mask_additional = mask & expected_v_actual["nt_site"].isin(additional_sites_to_ignore)
    mask_exclude = mask_exclude | mask_additional

problematic_sites = expected_v_actual[mask_exclude]

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
