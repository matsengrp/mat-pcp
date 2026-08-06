# #137 Task 1 — the all-clade synonymous offset table

Reproduce with:

```bash
pixi run python workflow/scripts/era_offset_table.py \
  --fitness results/gisaid-and-public_2026-04-22_v9/fitness/ntmut_fitness_aggregated.csv \
  --outdir analysis/137-era-offsets
```

Runs in **9.6 s** over the 1.26 GB / 4,559,104-row v9 fitness table. Read-only; it
touches no pipeline output.

`offset(clade, mut_type) = log(Σ actual_count / Σ expected_count)` over synonymous
nucleotide mutations of that type in that clade. Synonymous mutations are excluded
from the amino-acid fitness analyses, so under a well-calibrated neutral model the
offset should be ~0; a model needing an era split shows the offset separating by
era.

## Coverage

| | |
|---|---|
| clades × mutation types | 51 × 12 = **612 cells** |
| cells with `Σ actual_count ≥ 100` (enter the verdict) | **492** |
| cells below threshold (`sub_threshold_cells.csv`) | **118** |
| cells with zero counts, offset `-inf` (`zero_count_cells.csv`) | **2** — `21E`/`GC`, `24I`/`GC` |
| rows read / synonymous non-pooled kept | 4,559,104 / 963,481 |
| rows excluded as `source_clade == "all"` | **88,717** |

492 + 118 + 2 = 612: every cell is accounted for. The 120 excluded cells are listed
rather than dropped silently, and they are concentrated exactly where you would
expect — the smallest clades (`21E` at 588 samples contributes 10, `20F` 9, `24I` 8)
and the rarest transversions (`CG` 32, `GC` 22, `AC` 17).

The `source_clade != "all"` filter matched #137's stated 88,717 rows. The script
asserts no NaN `clade` survives it, so the correctness does not rest on pandas
dropping NaN group keys by default.

## The verdict: G→T passes, and it is the only one

Decision rule, fixed before the table was read: **complete separation**,
`min(offset over pre-Omicron cells) > max(offset over Omicron cells)`.

| mut_type | n pre-O | n O | pre-O range | Omicron range | gap | separates |
|---|---|---|---|---|---|---|
| **GT** | 19 | 29 | **+0.429 … +1.123** | **−0.344 … +0.036** | **+0.393** | **yes** |
| CT | 21 | 30 | +0.002 … +0.132 | −0.107 … +0.104 | −0.103 | no |
| TG | 13 | 25 | −0.005 … +1.134 | −0.111 … +0.121 | −0.126 | no |
| AG | 20 | 30 | −0.178 … +0.125 | −0.173 … +0.026 | −0.205 | no |
| AT | 15 | 28 | +0.036 … +0.410 | −0.049 … +0.323 | −0.288 | no |
| TC | 20 | 30 | −0.273 … −0.020 | −0.090 … +0.186 | −0.459 | no |
| GA | 19 | 29 | −0.305 … −0.079 | −0.124 … +0.223 | −0.527 | no |
| TA | 16 | 28 | −0.014 … +0.635 | +0.058 … +0.551 | −0.565 | no |
| CA | 12 | 28 | −0.286 … +0.109 | −0.076 … +0.343 | −0.629 | no |
| GC | 8 | 19 | +0.503 … +0.885 | +0.494 … +1.144 | −0.641 | no |
| AC | 11 | 23 | −0.133 … +0.264 | −0.362 … +0.550 | −0.683 | no |
| CG | 7 | 12 | +0.484 … +1.300 | +0.213 … +1.788 | −1.304 | no |

**GATE PASSED.** G→T separates with a gap of **+0.393** — reproducing #137's
pre-stated value and ranges exactly. The v10 era split is justified.

## The finding that differs from what #137/#138 predicted

#138's risk section says "**two** of twelve mutation types are expected to fail
complete separation (C→G, G→C)". In fact **eleven of twelve fail**. C→G and G→C do
fail as predicted — they are the two worst, at −1.304 and −0.641 — but so does
everything except G→T.

The prediction came from the four-clade `spread` table in #137, which is a
different quantity: `max − min` across 21J/21L/23A/24A measures *dispersion*, not
*era structure*. A type can have a large spread with the extremes both inside one
era (which is precisely what #137 notes for C→G, whose extreme is 21L, and T→G,
whose extreme is 21F). Complete separation over 19–21 pre-Omicron cells against
28–30 Omicron cells is a far stricter bar than a four-clade spread suggests: one
overlapping cell out of ~48 defeats it.

The count-weighted era means show the same thing without the strictness:

| mut_type | pre-Omicron | Omicron | difference | Σ actual |
|---|---|---|---|---|
| **GT** | **+0.638** | **−0.094** | **+0.732** | 225,699 |
| AC | +0.174 | +0.051 | +0.124 | 37,381 |
| TG | +0.135 | +0.028 | +0.108 | 50,243 |
| GC | +0.773 | +0.686 | +0.087 | 13,793 |
| AT | +0.187 | +0.113 | +0.074 | 70,298 |
| CT | +0.050 | +0.022 | +0.028 | 2,590,660 |
| CG | +1.105 | +1.091 | +0.015 | 14,648 |
| AG | −0.059 | −0.046 | −0.012 | 456,799 |
| CA | −0.011 | +0.100 | −0.112 | 68,884 |
| TC | −0.134 | +0.013 | −0.148 | 949,370 |
| GA | −0.151 | +0.022 | −0.173 | 289,100 |
| TA | +0.109 | +0.351 | −0.242 | 96,725 |

**G→T's era effect is 5.9× the next largest** (+0.732 vs A→C's +0.124), and five
types have an era difference of the *opposite* sign. So the honest description of
this release is that **the era split is a G→T fix**: G→T is the one type where the
pooled model is badly miscalibrated *and* the miscalibration is an era effect.

Two consequences worth carrying into Task 2 rather than treating as noise:

1. The follow-up issue #138 anticipated for "C→G and G→C" is filed as
   [#140](https://github.com/matsengrp/mat-pcp/issues/140) and covers **eleven** types.
   Its framing changes too: the question is not "why do these two resist an
   era split" but "G→T aside, is there any era structure in this model at all?"
2. `glm.py` fits each mutation type independently, so a two-era split halves the
   training data for eleven types that show no era effect, buying variance for no
   bias reduction. That is not an argument against Task 2 — parity with upstream
   and the pooled control arm both argue for the straightforward split — but it is
   an argument for reporting per-type rate uncertainties in the post-run
   validation, and a candidate for the "era as a factor in one pooled fit"
   follow-up #138 already lists as out of scope.

## Files

| file | contents |
|---|---|
| `offset_cells.csv` | all 612 cells: clade, mut_type, summed counts, offset, era |
| `offset_matrix.csv` | 51 × 12 offset matrix in clade order, era column first |
| `separation_verdict.csv` | per-type ranges, gap, pass/fail |
| `era_weighted_means.csv` | count-weighted and unweighted mean offset per (era, type) |
| `sub_threshold_cells.csv` | the 118 cells below `--min-count` |
| `zero_count_cells.csv` | the 2 cells with `-inf` offset |
