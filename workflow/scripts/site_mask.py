#!/usr/bin/env python3
"""Parse problematic sites VCF to extract site positions."""

import pandas as pd
import urllib.request
import tempfile
import os

# Get params from snakemake
url = snakemake.params.get("url")
output_csv = snakemake.output.csv
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log(f"Downloading problematic sites VCF from {url}")

# Download VCF
with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as tmp:
    tmp_path = tmp.name

try:
    urllib.request.urlretrieve(url, tmp_path)
    log(f"Downloaded to {tmp_path}")

    # Parse VCF
    sites_df = pd.read_csv(
        tmp_path,
        sep="\t",
        comment="#",
        names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"],
    )

    log(f"Found {len(sites_df)} sites in VCF")

    # Extract just the site positions
    result = sites_df.rename(columns={"POS": "site"})[["site"]]
    result.to_csv(output_csv, index=False)

    log(f"Wrote {len(result)} sites to {output_csv}")

finally:
    if os.path.exists(tmp_path):
        os.unlink(tmp_path)

log("Done!")
