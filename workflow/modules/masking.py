"""Shared masking logic: the single source of truth for per-mutation exclusion.

The pipeline formerly computed exclusions in two places (``count_mutations.py``
and ``aggregate_mutations_to_exclude.py``), which duplicated the logic and could
drift. This module centralizes it. It is applied exactly once, in
``build_count_table.py``, against the full mutation grid.

``compute_exclusions`` mirrors the union of the two former implementations:

  * ``sites_to_exclude``    -- global config sites            -> exclude
  * ``vcf_masked_sites``    -- ProblematicSites VCF sites      -> exclude + masked_in_usher
  * ``usher_clade_sites``   -- {clade: set(sites)} from yaml   -> exclude + masked_in_usher
  * ``site_include_range``  -- {start, end}; outside range     -> exclude

``exclude`` is the union of all reasons; ``masked_in_usher`` flags the
artifactual-site masks (VCF + UShER yaml sites) only.
"""

import pandas as pd
import yaml


def load_usher_masked_sites(path):
    """Parse the UShER masked-sites YAML into ``{clade: set(int sites)}``.

    A block applies to a clade iff the clade is listed in the block's
    ``clades``. Only the site-level ``sites`` masks are read here.
    """
    with open(path) as f:
        spec = yaml.safe_load(f) or {}
    clade_sites = {}
    for block in spec.values():
        sites = block.get("sites", []) or []
        for clade in block.get("clades", []) or []:
            clade_sites.setdefault(clade, set()).update(int(s) for s in sites)
    return clade_sites


def compute_exclusions(
    df,
    *,
    sites_to_exclude,
    vcf_masked_sites,
    usher_clade_sites,
    site_include_range,
):
    """Return ``(exclude, masked_in_usher)`` boolean Series aligned to ``df``.

    ``df`` must have columns ``clade`` and ``nt_site`` (int-coercible). See the
    module docstring for the source semantics.
    """
    site = df["nt_site"].astype(int)
    clade = df["clade"]

    config_sites = {int(s) for s in (sites_to_exclude or [])}
    vcf_sites = {int(s) for s in (vcf_masked_sites or [])}

    start = site_include_range.get("start", 266)
    end = site_include_range.get("end", 29674)

    # Per-(clade, site) membership in the UShER yaml masks.
    usher_pairs = [
        (c, int(s)) for c, sites in usher_clade_sites.items() for s in sites
    ]
    in_usher = pd.MultiIndex.from_arrays([clade, site]).isin(usher_pairs)

    in_vcf = site.isin(vcf_sites).to_numpy()
    in_config = site.isin(config_sites).to_numpy()
    out_of_range = ((site < start) | (site > end)).to_numpy()

    masked_in_usher = in_vcf | in_usher
    exclude = out_of_range | in_config | masked_in_usher

    return (
        pd.Series(exclude, index=df.index),
        pd.Series(masked_in_usher, index=df.index),
    )
