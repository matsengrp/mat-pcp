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

## Results: mat-pcp **v10 per-era** vs upstream's two eras

Run era-against-era after the v10 production run (#137 success criterion 4). This is
the check the pooled comparison below was standing in for.

| comparison | Pearson (log) | Spearman | median \|log ratio\| | p95 |
|---|---|---|---|---|
| **mat-pcp pre-Omicron vs upstream pre-Omicron** | 0.99994 | 0.99992 | **0.0153** | 0.0427 |
| **mat-pcp Omicron vs upstream Omicron** | 0.99989 | 0.99985 | **0.0157** | 0.0654 |
| mat-pcp pre-Omicron vs mat-pcp Omicron | 0.9927 | 0.9946 | 0.1301 | 0.5266 |
| upstream pre-Omicron vs upstream Omicron | 0.9928 | 0.9947 | 0.1317 | 0.5444 |

Two things to read off this.

**The per-era fits agree with upstream ~5x better than the pooled fit did** — median
|log ratio| 0.015 against the pooled 0.087 / 0.070 below. Since the two codebases run
on different trees with different date cutoffs and filtering, 1.5% median disagreement
across 512 cells is close to the floor set by the data alone.

**mat-pcp's own era difference now reproduces upstream's**: 0.1301 against 0.1317, on
tables fitted independently. That is the era effect measured the same way by two
implementations, and it is what the split was for.

### Criterion 4, evaluated

| requirement | result |
|---|---|
| (a) each era table has 512 rows on the same grid | **PASS** — joined 1:1 on all 512 |
| (b) sign of `mean log(rate_pre/rate_om)` matches upstream for all 10 types with \|value\| > 0.05 | **PASS** — 10/10 |
| (c) G→T mean log ratio in [0.4, 1.0] | **PASS** — **+0.6158** (upstream +0.6161) |

Per-type, mat-pcp's fitted era difference against upstream's — every one of the twelve
agrees in sign, and the largest disagreement in magnitude is 0.03:

| mut_type | upstream | mat-pcp v10 |
|---|---|---|
| **GT** | **+0.6161** | **+0.6158** |
| GC | +0.2707 | +0.2889 |
| AC | +0.2410 | +0.2295 |
| AT | +0.1294 | +0.1306 |
| CG | +0.0972 | +0.0646 |
| TG | +0.0767 | +0.0887 |
| AG | +0.0100 | +0.0149 |
| CT | −0.0142 | −0.0154 |
| CA | −0.0883 | −0.0859 |
| TC | −0.1342 | −0.1312 |
| TA | −0.1550 | −0.1634 |
| GA | −0.1628 | −0.1562 |

**This settles the #140 caveat below.** mat-pcp's own fitted tables put G→C at +0.289
and A→C at +0.230 — real second-tier era structure, confirming upstream rather than
Task 1's read. Complete separation of the *synonymous offset* is a stricter test than
"has era structure", so the two measures were never in conflict; the wording *"eleven
of twelve types show no era structure"* was simply too strong, and should be read as
"eleven of twelve do not completely separate".

## Results: mat-pcp v9 pooled vs upstream's two eras (superseded by the above)

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
