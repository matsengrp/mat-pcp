"""Shared masking logic: the single source of truth for per-mutation exclusion.

The pipeline formerly computed exclusions in two places (``count_mutations.py``
and ``aggregate_mutations_to_exclude.py``), which duplicated the logic and could
drift. This module centralizes it. It is applied exactly once, in
``build_count_table.py``, against the full mutation grid.

``compute_exclusions`` applies every masking reason in one place:

  * ``sites_to_exclude``      -- global config sites             -> exclude
  * ``vcf_masked_sites``      -- ProblematicSites VCF sites       -> exclude + masked_in_usher
  * ``usher_clade_sites``     -- {clade: set(sites)} from yaml     -> exclude + masked_in_usher
  * ``usher_clade_mutations`` -- {clade: set(mutations)} from yaml -> exclude + masked_in_usher
  * ``reference_nt`` column   -- founder->reference reversions     -> exclude
  * ``site_include_range``    -- {start, end}; outside range      -> exclude

``exclude`` is the union of all reasons; ``masked_in_usher`` flags the
artifactual masks from the UShER yaml (its ``sites`` and per-clade
``mutations``) and the ProblematicSites VCF only.

The UShER yaml carries two kinds of mask per block: ``sites`` masks the whole
nucleotide site clade-wide, while the optional ``mutations`` list masks only a
specific directional reversion (e.g. ``"G26530A"``, in ``nt_mutation`` format)
without masking the rest of the site.

A **founder->reference reversion** is any grid row whose destination base equals
the reference base at that site (``to_nt == reference_nt``); because the grid
only contains rows with ``to_nt != clade_founder_nt``, this implies the clade
founder differs from the reference there (a clade-defining site). These
reversions are enriched for artifacts and are always excluded. They set
``exclude`` but NOT ``masked_in_usher`` (they are not a UShER/VCF artifact mask).
"""

import pandas as pd
import yaml


def load_usher_masks(path):
    """Parse the UShER masked-sites YAML in one pass.

    Returns ``(clade_sites, clade_mutations)``:

      * ``clade_sites``     -- ``{clade: set(int sites)}`` from each block's
        ``sites`` (whole-site clade-wide masks).
      * ``clade_mutations`` -- ``{clade: set(mutation str)}`` from each block's
        optional ``mutations`` (directional reversions in ``nt_mutation``
        format, e.g. ``"G26530A"``); blocks without a ``mutations`` key
        contribute nothing.

    A block applies to a clade iff the clade is listed in the block's
    ``clades``.
    """
    with open(path) as f:
        spec = yaml.safe_load(f) or {}
    clade_sites = {}
    clade_mutations = {}
    for block in spec.values():
        sites = block.get("sites", []) or []
        mutations = block.get("mutations", []) or []
        for clade in block.get("clades", []) or []:
            clade_sites.setdefault(clade, set()).update(int(s) for s in sites)
            clade_mutations.setdefault(clade, set()).update(str(m) for m in mutations)
    return clade_sites, clade_mutations


def global_masked_sites(
    *,
    sites_to_exclude,
    vcf_masked_sites,
    site_include_range,
    genome_len,
):
    """Return the set of *globally* (clade-independent) masked nucleotide sites.

    This is the subset of :func:`compute_exclusions`'s reasons that does not
    depend on a clade: the config ``sites_to_exclude`` (``in_config``), the
    ProblematicSites VCF sites (``in_vcf``), and everything outside
    ``site_include_range`` (``out_of_range``). The per-clade UShER-yaml site and
    mutation masks and the founder->reference reversions are intentionally NOT
    included -- they are clade/mutation specific and have no meaning on the
    clade-free full tree used by the sliding-window path.

    :func:`compute_exclusions` remains the per-row authority (it must keep
    ``in_vcf`` separable to build ``masked_in_usher``); this helper is the shared
    definition of the global set for callers that need it as a plain site list
    (currently the windowed expected-count path). Sites are 1-based.
    """
    start = site_include_range.get("start", 266)
    end = site_include_range.get("end", 29674)
    masked = {int(s) for s in (sites_to_exclude or [])}
    masked.update(int(s) for s in (vcf_masked_sites or []))
    masked.update(s for s in range(1, genome_len + 1) if s < start or s > end)
    return masked


def compute_exclusions(
    df,
    *,
    sites_to_exclude,
    vcf_masked_sites,
    usher_clade_sites,
    usher_clade_mutations=None,
    site_include_range,
):
    """Return ``(exclude, masked_in_usher)`` boolean Series aligned to ``df``.

    ``df`` must have columns ``clade``, ``nt_site`` (int-coercible),
    ``nt_mutation``, and ``reference_nt``. See the module docstring for the
    source semantics.
    """
    if "reference_nt" not in df.columns:
        raise KeyError(
            "compute_exclusions requires a 'reference_nt' column to flag "
            "founder->reference reversions; none was supplied. The reference "
            "base is always available at the call site (build_count_table)."
        )

    site = df["nt_site"].astype(int)
    clade = df["clade"]

    config_sites = {int(s) for s in (sites_to_exclude or [])}
    vcf_sites = {int(s) for s in (vcf_masked_sites or [])}

    start = site_include_range.get("start", 266)
    end = site_include_range.get("end", 29674)

    # Per-(clade, site) membership in the UShER yaml site masks.
    usher_pairs = [
        (c, int(s)) for c, sites in usher_clade_sites.items() for s in sites
    ]
    in_usher_site = pd.MultiIndex.from_arrays([clade, site]).isin(usher_pairs)

    # Per-(clade, nt_mutation) membership in the optional UShER yaml mutation
    # masks (directional reversions masked at the mutation level, not the site).
    usher_mut_pairs = [
        (c, m)
        for c, muts in (usher_clade_mutations or {}).items()
        for m in muts
    ]
    in_usher_mut = pd.MultiIndex.from_arrays(
        [clade, df["nt_mutation"]]
    ).isin(usher_mut_pairs)

    # Founder->reference reversion: destination base equals the reference base.
    # The grid guarantees to_nt != clade_founder_nt, so this only fires at
    # clade-defining sites (founder != reference). Sets exclude, not masked_in_usher.
    to_nt = df["nt_mutation"].str[-1]
    if to_nt.isna().any():
        raise ValueError(
            "compute_exclusions: 'nt_mutation' has NaN/empty values, so the "
            "reversion destination base cannot be determined; the grid must "
            "supply a well-formed nt_mutation (e.g. 'G26530A') for every row."
        )
    in_reversion = (to_nt == df["reference_nt"]).to_numpy()

    in_vcf = site.isin(vcf_sites).to_numpy()
    in_config = site.isin(config_sites).to_numpy()
    out_of_range = ((site < start) | (site > end)).to_numpy()

    masked_in_usher = in_vcf | in_usher_site | in_usher_mut
    exclude = out_of_range | in_config | masked_in_usher | in_reversion

    return (
        pd.Series(exclude, index=df.index),
        pd.Series(masked_in_usher, index=df.index),
    )
