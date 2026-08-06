# Vendored upstream per-era rate tables

`master_table_pre_omicron.csv` and `master_table_omicron.csv`, copied verbatim from
[`neherlab/SARS2-mut-fitness-v2`](https://github.com/neherlab/SARS2-mut-fitness-v2/tree/main/results/master_tables)
at commit `6dc1853c127a` (2025-01-09). 512 rows each — the full
`(mut_type, motif, unpaired, nt_site_before_boundary)` grid — with columns
`mut_type, motif, unpaired, nt_site_boundary, nt_site_before_boundary, rate,
predicted_count, residual`. mat-pcp's `glm/master_table.csv` has the same schema
except it names the second-to-last column `expected_count`.

Vendored rather than fetched at runtime so the comparison is reproducible without
network access, and pinned to a commit so a silent upstream refit cannot change our
numbers underneath us. 38 KB each.

Compare with:

```bash
pixi run python workflow/scripts/compare_rates_to_upstream_eras.py \
    --pooled results/gisaid-and-public_2026-04-22_v9/glm/master_table.csv \
    --outdir analysis/137-era-offsets/upstream
```

## This is not the same check as `tests/test_glm_reference_equivalence.py`

That test runs **upstream's code on mat-pcp's fixture input** and compares at
floating-point-noise tolerance: it asks *is the port faithful*, and any real
difference fails it. `tests/data/glm_reference_rates.csv` is that fixture's output —
pooled, with no `era` column, and on a completely different rate scale (~1e-12)
because it was fit on the tiny test fixture.

The tables here are **upstream's published rates, fit on upstream's own data** — a
different tree, date cutoff and filtering. Agreement can only be approximate, so the
comparison reports statistics and asserts nothing about absolute equality. It does
assert the grids join 1:1 on all 512 cells, because a partial join would compute
every statistic over an unstated subset and still report healthy correlations.

## Results: mat-pcp v9 pooled vs upstream's two eras

| comparison | Pearson (log) | Spearman | median \|log ratio\| | p95 |
|---|---|---|---|---|
| mat-pcp pooled vs upstream pre-Omicron | 0.9972 | 0.9978 | 0.0866 | 0.3040 |
| mat-pcp pooled vs upstream Omicron | 0.9985 | 0.9990 | 0.0703 | 0.2743 |
| **upstream pre-Omicron vs upstream Omicron** | 0.9928 | 0.9947 | **0.1317** | 0.5444 |

Two things follow, and both matter for v10.

**The port is sound.** mat-pcp's pooled rates track upstream's cell-by-cell to within
~7–9% median, across 512 cells, despite different trees and cutoffs. Pooled also lies
*between* the two upstream era rates in **70.1%** of cells, which is where a pooled fit
should sit.

**The era effect is larger than the noise floor.** The two upstream eras differ from
*each other* (median |log ratio| 0.1317) by more than mat-pcp's pooled fit differs from
either of them (0.087, 0.070). So the era signal exceeds the combined
implementation-plus-data difference between two independent codebases — which is the
strongest justification the split has, and does not depend on #137 Task 1's separation
criterion at all.

## Upstream's own era difference per mutation type

`log(pre_Omicron rate / Omicron rate)`, mean over each type's grid cells:

| mut_type | mean | median | cells |
|---|---|---|---|
| **GT** | **+0.6161** | +0.6134 | 32 |
| GC | +0.2707 | +0.2574 | 64 |
| AC | +0.2410 | +0.2496 | 32 |
| AT | +0.1294 | +0.1302 | 64 |
| CG | +0.0972 | +0.0462 | 64 |
| TG | +0.0767 | +0.0539 | 32 |
| AG | +0.0100 | +0.0474 | 32 |
| CT | −0.0142 | −0.0124 | 64 |
| CA | −0.0883 | −0.0768 | 32 |
| TC | −0.1342 | −0.1371 | 32 |
| TA | −0.1550 | −0.1178 | 32 |
| GA | −0.1628 | −0.1550 | 32 |

**G→T ranks first by 2.3×**, independently reproducing Task 1's conclusion from a
different quantity — Task 1 measured synonymous residuals against a pooled model,
this compares two fitted rate tables. Task 1 put G→T's era effect at +0.732; upstream's
rates give +0.616. Same sign, same rank, same order of magnitude.

**One caveat on #140's framing.** Upstream's rates put G→C (+0.271) and A→C (+0.241) at
roughly 40% of G→T — appreciable second-tier era structure. Task 1's weighted-offset
analysis rated both as indistinguishable from zero (A→C at z ≈ 1.2), and
[#140](https://github.com/matsengrp/mat-pcp/issues/140) was closed on that basis. The
two measures are not the same quantity and upstream's figures carry no error bars here,
so this is not a contradiction — but *"eleven of twelve types show no era structure"* is
a stronger claim than upstream's tables support. Re-running this script era-against-era
after #137 Task 2 and the v10 production run settles it directly, which is the intended
post-run sanity check.
