"""The pre-Omicron / Omicron era partition, defined once (#137).

``SARS2-mut-fitness-v2`` fits **two** neutral rate models — one on pre-Omicron
clades, one on Omicron+ — and dispatches per clade. mat-pcp fits one global model
and applies a per-clade scalar; the ``pre_omicron_or_omicron`` covariate did not
survive the port, while the light-switch boundary covariate from the same
upstream notebook did. This module is the mapping that restores it.

**Two upstream details are deliberately not ported.**

1. *Upstream dispatches on a positional slice.* ``clades_pre_o = clades[0:11]``
   indexes ``count_df.clade.unique()``, which is correct only because the counts
   file happens to arrive in Nextstrain clade order. Nothing asserts that
   ordering, and a row-order change would silently relabel eleven clades. Every
   consumer here dispatches through :func:`era_for_clade` on the label.

2. *Upstream's 11-clade list is not an enumeration of pre-Omicron clades*, it is
   the subset upstream analyzed. ``gisaid-and-public_2026-04-22`` carries 22
   pre-Omicron clades — upstream's 11 plus ``19A``, ``19B``, ``20D``, ``20F``,
   ``21B``, ``21D``, ``21E``, ``21F``, ``21G``, ``21H`` — every one of which
   upstream's ``else 'omicron'`` fallback would label Omicron. :data:`CLADE_ERA`
   is therefore exhaustive and :func:`era_for_clade` **raises** on an
   unrecognized clade rather than defaulting to either era. A default here is the
   #129/#130/#133 silent-drop shape: the wrong era produces a plausible fit, not
   an error.

**Why the mapping is audited rather than derived.** :data:`CLADE_ERA` is 58
hand-maintained entries whose failure mode is a clade in the wrong list, which no
internal cross-check can see. :func:`era_from_label` is a total predicate over
the Nextstrain naming convention (year of emergence + letter) that
:func:`audit_clade_era_mapping` compares against the mapping entry by entry. It
audits and never substitutes: a fallback would silently label a hypothetical
``21N`` and retire the raise that item 2 exists to provide.

Four properties of that predicate are load-bearing, and each is asserted by
``tests/test_era_mapping.py``:

* **It is two-sided.** The ``21K``/``21L``/``21M`` exclusion on the pre-Omicron
  half is what makes it so, and the failure it prevents is worse than a missed
  error: without the exclusion the audit *demands* the wrong answer. Drop it and
  ``era_from_label`` calls every ``21x`` clade pre-Omicron, so auditing the
  shipped mapping raises on ``21K``, ``21L`` and ``21M`` — three failures whose
  obvious "fix" is to move all three into ``pre_omicron``, after which the audit
  passes over all 58 entries with the three Omicron founders in the wrong era.
  (Moving ``21K`` alone does *not* pass: ``21L``/``21M`` then disagree. The
  vacuity needs all three, which is exactly what a maintainer clearing three
  audit failures would do.)
* **It asserts the label format.** The clade label space is *not* guaranteed to be
  ``NN[A-Z]``. ``checkpoint get_clades`` normalizes ``"21K (BA.1)"`` to ``"21K"``
  by accident of awk field splitting (``split($0, a, " ")`` then ``a[2]``,
  ``Snakefile:768-770``), which passes single-token non-clades through intact —
  ``config_test_data.yaml``'s clade column carries ``recombinant`` and an empty
  value. A malformed label raises rather than being bucketed.
* **It parses the year as an integer.** Two-digit years make lexicographic
  comparison agree with numeric comparison by coincidence, not by construction.
* **It holds forward.** Nextstrain names a clade for its year of emergence, so a
  clade minted ``22``+ cannot be pre-Omicron.

**The parenthetical annotation cannot serve as the era key**, which is worth
recording because it looks like a better check. On the test data it nearly does
(``21K (Omicron)``, ``21L (Omicron)``, ``21M (Omicron)``), but on the production
tree ``Nextstrain_clade_usher`` gives WHO names through ``21J`` and Pango
lineages from ``21K`` on — ``21K (BA.1)``, ``21L (BA.2)``, ``22A (BA.4)``,
``23I (BA.2.86)`` — and ``21M (Omicron)`` is the only one of 58 values that names
the era at all. The awk normalization above strips it in any case.
"""

import re

PRE_OMICRON = "pre_omicron"
OMICRON = "omicron"
ERAS = (PRE_OMICRON, OMICRON)

#: Every ``Nextstrain_clade_usher`` value on ``gisaid-and-public_2026-04-22``
#: (52 clades, ``mat/clade_counts.csv``), plus the six ``25D``-``25I`` that the
#: 2026-08-06 trees add. ``21A`` and ``25I`` are mapped but not processed by the
#: pipeline — neither resolves to a clade founder (``21A`` is the one 19A-22F
#: clade absent from the Neher founder JSON; ``25I``'s lineage ``BA.3.2`` is
#: absent from the frozen Roemer JSON), so ``clades_to_exclude`` drops both
#: (#129). They are kept here so the audit covers them and so re-adding a founder
#: does not silently reintroduce an unmapped clade.
CLADE_ERA = {
    # Pre-Omicron: 19A-19B, 20A-20J, 21A-21J (22 clades, 21 processed).
    "19A": PRE_OMICRON,
    "19B": PRE_OMICRON,
    "20A": PRE_OMICRON,
    "20B": PRE_OMICRON,
    "20C": PRE_OMICRON,
    "20D": PRE_OMICRON,
    "20E": PRE_OMICRON,
    "20F": PRE_OMICRON,
    "20G": PRE_OMICRON,
    "20H": PRE_OMICRON,
    "20I": PRE_OMICRON,
    "20J": PRE_OMICRON,
    "21A": PRE_OMICRON,
    "21B": PRE_OMICRON,
    "21C": PRE_OMICRON,
    "21D": PRE_OMICRON,
    "21E": PRE_OMICRON,
    "21F": PRE_OMICRON,
    "21G": PRE_OMICRON,
    "21H": PRE_OMICRON,
    "21I": PRE_OMICRON,
    "21J": PRE_OMICRON,
    # Omicron: 21K (BA.1), 21L (BA.2), 21M (basal B.1.1.529), then 22A+ (30).
    "21K": OMICRON,
    "21L": OMICRON,
    "21M": OMICRON,
    "22A": OMICRON,
    "22B": OMICRON,
    "22C": OMICRON,
    "22D": OMICRON,
    "22E": OMICRON,
    "22F": OMICRON,
    "23A": OMICRON,
    "23B": OMICRON,
    "23C": OMICRON,
    "23D": OMICRON,
    "23E": OMICRON,
    "23F": OMICRON,
    "23G": OMICRON,
    "23H": OMICRON,
    "23I": OMICRON,
    "24A": OMICRON,
    "24B": OMICRON,
    "24C": OMICRON,
    "24D": OMICRON,
    "24E": OMICRON,
    "24F": OMICRON,
    "24G": OMICRON,
    "24H": OMICRON,
    "24I": OMICRON,
    "25A": OMICRON,
    "25B": OMICRON,
    "25C": OMICRON,
    # New in the 2026-08-06 trees. 25I is excluded from processing (see above).
    "25D": OMICRON,
    "25E": OMICRON,
    "25F": OMICRON,
    "25G": OMICRON,
    "25H": OMICRON,
    "25I": OMICRON,
}

#: The three 2021-named clades that are Omicron. ``21K`` is BA.1, ``21L`` is
#: BA.2, ``21M`` is basal B.1.1.529 — all descend from the November 2021
#: emergence, so the year prefix alone cannot separate them from Delta-era
#: clades sharing the ``21`` prefix.
OMICRON_21 = frozenset({"21K", "21L", "21M"})

#: Nextstrain clade labels are a two-digit year of emergence plus one uppercase
#: letter. Anchored with ``fullmatch`` so ``recombinant``, ``""``, and
#: ``"21K (BA.1)"`` all fail rather than partially matching.
_CLADE_LABEL = re.compile(r"(\d{2})([A-Z])")


def era_for_clade(clade):
    """The era of ``clade`` according to :data:`CLADE_ERA`.

    Raises ``KeyError`` on an unmapped clade. This is the *only* function callers
    should use to bucket a clade: it is the explicit mapping, not the year-prefix
    heuristic, and an unmapped clade must stop the run rather than land in a
    default bucket.
    """
    try:
        return CLADE_ERA[clade]
    except KeyError:
        raise KeyError(
            f"clade {clade!r} has no era in CLADE_ERA. Add it to "
            f"workflow/modules/era.py rather than defaulting: an unmapped clade "
            f"fitted into the wrong era produces a plausible rate table, not an "
            f"error. Known clades: {sorted(CLADE_ERA)}"
        ) from None


def era_from_label(clade):
    """The era *implied by* a Nextstrain clade label's year prefix.

    An auditor for :data:`CLADE_ERA`, never a substitute for it — see the module
    docstring. Raises ``ValueError`` on a label that is not ``NN[A-Z]``.
    """
    match = _CLADE_LABEL.fullmatch(clade)
    if match is None:
        raise ValueError(
            f"clade label {clade!r} is not a Nextstrain clade of the form "
            f"NN[A-Z] (two-digit year of emergence plus one uppercase letter). "
            f"get_clades normalizes '21K (BA.1)' to '21K' but passes "
            f"single-token non-clades through unchanged, so this is reachable "
            f"from a config whose clade column carries e.g. 'recombinant'."
        )
    year = int(match.group(1))
    if year <= 20:
        return PRE_OMICRON
    if year == 21:
        return OMICRON if clade in OMICRON_21 else PRE_OMICRON
    return OMICRON


def audit_clade_era_mapping(mapping=None):
    """Check every entry of ``mapping`` against :func:`era_from_label`.

    Raises ``ValueError`` listing every disagreement. Called at fit time so a
    mis-edited mapping fails the run rather than shifting a rate table.
    """
    if mapping is None:
        mapping = CLADE_ERA
    disagreements = []
    for clade, era in sorted(mapping.items()):
        if era not in ERAS:
            disagreements.append(f"{clade}: era {era!r} is not one of {ERAS}")
            continue
        implied = era_from_label(clade)
        if implied != era:
            disagreements.append(
                f"{clade}: mapped {era!r} but the year prefix implies {implied!r}"
            )
    if disagreements:
        raise ValueError(
            "CLADE_ERA disagrees with the Nextstrain year-prefix convention:\n  "
            + "\n  ".join(disagreements)
        )
    return len(mapping)


def partition_clades(clades):
    """Split ``clades`` into ``{era: [clade, ...]}``, preserving input order.

    Only the eras actually present appear as keys, which is what lets the
    Snakefile request exactly the era tables it can fit: a single-era tree (the
    ``max_clades: 3`` test config, and the out-of-scope BA.2/JN.1/XFG trees)
    yields one key, so the absent era's rule is never built (#138 decision 9).
    Raises on an unmapped clade via :func:`era_for_clade` — an empty era is a
    legitimate skip, an unknown clade never is.
    """
    partition = {}
    for clade in clades:
        partition.setdefault(era_for_clade(clade), []).append(clade)
    return partition
