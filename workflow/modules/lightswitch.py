"""The light-switch boundary convention, defined once (#125).

The GLM's ``global_context`` factor splits each mutation type into two regimes at
a configured genome position (``glm_lightswitch_boundaries``). Three places have
to agree on which side the boundary site itself belongs to:

1. **Training** — ``modules/glm.py::_compute_global_context_feature`` computes
   ``global_context = (nt_site > boundary)``, so ``site == boundary`` is fitted
   into the *before* regime.
2. **The rate table's row identity** — ``modules/rates.py`` predicts each grid
   row through a sentinel position (``nt_site = 0 if before_boundary else
   s_max``), which makes ``nt_site_before_boundary=True`` mean exactly the
   ``global_context == 0`` bucket from (1).
3. **The serving label** — ``scripts/annotate_counts.py`` (clade path) and
   ``modules/expected_counts.py`` (windowed trajectory path) each label every
   mutation with ``nt_site_before_boundary``, which is what the rate lookup keys
   on.

Through v8 the serving label used ``site < boundary``. Since ``not (site <
boundary)`` is ``site >= boundary`` rather than ``site > boundary``, the boundary
site was **fitted into the before-bucket and served the after-bucket's rate**.
On ``gisaid-and-public_2026-04-22_v8`` that exposed one real pair (A→T at 21562:
51 rows, 82 counts); the other three (C→T at 13467, C→G and G→C at 21562) had no
observed rows.

``before_boundary`` below is therefore ``site <= boundary``, aligning the label
to the bucket that produced the rate. It is used by *both* label producers so
they cannot drift apart: fixing only one would leave the clade path and the
windowed path disagreeing at the boundary site, which is precisely the
clade/window comparability the trajectory work was built to preserve (the
``evol_time`` normalization note in ``CLAUDE.md``).

**Deliberate divergence from the reference.** ``SARS2-mut-fitness-v2`` carries
this same inconsistency — ``modules/glm.py:59`` trains with ``>`` while
``notebook/counts_by_clade.py.ipynb`` labels with ``<`` — so after this fix
mat-pcp labels site 21562 differently from upstream. That is intentional. Note
that ``tests/test_glm_reference_equivalence.py`` does *not* detect it: its
fixture is a frozen CSV carrying a pre-computed ``nt_site_before_boundary``
column, so it exercises the rate *lookup* over stored labels, not the label
computation. Do not "resolve" the divergence by re-baselining that fixture --
it exists to catch a dropped synonymous filter (#127), and diluting it to
absorb an unrelated convention change would retire the one thing it guards.
"""


def before_boundary(site, boundary):
    """Whether ``site`` falls in the *before-boundary* light-switch bucket.

    Accepts scalars or pandas Series/NumPy arrays for either argument, so the
    per-row clade path and the vectorized windowed path share one definition.

    A mutation type with no configured boundary must come out ``False``
    everywhere: ``glm.py`` adds no ``global_context`` factor for it and
    ``rates.py`` builds only ``nt_site_before_boundary=False`` grid rows, so a
    ``True`` label would have no rate row to match. Callers encode that by
    passing ``boundary = 0`` (what ``scripts/master_table.py`` writes into
    ``nt_site_boundary`` for unbounded types) or by short-circuiting before
    calling; sites are 1-based, so ``site <= 0`` is ``False`` for every real
    site.
    """
    return site <= boundary
