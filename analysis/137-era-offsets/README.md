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

> **The gate is genuinely pre-registered, and also genuinely fragile. Both matter.**
>
> *Pre-registered:* `ISSUE-137.md` fixes `Σ actual_count ≥ 100` and states "492 of
> 612 cells clear that threshold" and the +0.393 gap **before** this table was
> produced, so the cutoff was not chosen to make the gate pass.
>
> *Fragile:* the verdict flips at `--min-count 52`. At any threshold ≤ 51,
> **zero of twelve** types separate and the script exits 1. The sole cause is one
> cell — **20F G→T, 51 counts, offset −0.0955** — a pre-Omicron cell sitting inside
> the Omicron range (−0.344 … +0.036). It is listed in `sub_threshold_cells.csv`,
> as the issue required, rather than hidden; but it is a real counterexample that
> the cutoff excludes by one count, not a noise cell (its Poisson SE is ≈0.14, so
> it sits ~5 SE below the pre-Omicron mean of +0.64).
>
> Complete separation is a min-vs-max statistic and therefore maximally sensitive
> to exactly the cells a count cutoff removes. **The release argument should rest on
> the weighted era difference below (G→T at z ≈ 14, nothing else above 1.4σ), not on
> "GATE PASSED".** The substantive conclusion is unaffected and survives every
> confound tested — the era effect is a clean step at the Omicron boundary (flat at
> +0.64 across 2019–2021, flat at −0.13 across 2021–2025), uncorrelated with clade
> size or calendar time within era (|Spearman| ≤ 0.06), and unaffected by per-clade
> UShER masking (8 of 1.48 M G→T counts fall on masked pairs).

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

**G→T is the only era difference distinguishable from zero.** Its weighted
difference is +0.732; adding a between-clade SE (the weights are heavily
concentrated — the largest clade holds 48–61% of pre-Omicron weight, so the Kish
effective sample size is only ~3 per era) puts it at **z ≈ 14**, with no other
positive difference above **1.4σ**. So the honest description of this release is
that **the era split is a G→T fix**: G→T is the one type where the pooled model is
badly miscalibrated *and* the miscalibration is an era effect.

> **Two earlier framings of this table were wrong and are corrected here**, because
> both read as independent support and neither is.
>
> *"G→T's era effect is 5.9× the next largest."* The denominator is noise: A→C's
> +0.124 is z ≈ 1.2. A ratio against a quantity indistinguishable from zero is
> unstable and would change materially on another tree. Use the z column, not the
> ratio.
>
> *"Five types have an era difference of the opposite sign."* This is largely
> **mechanically forced**, not five findings. `rates.py:204` sets
> `T = N_syn / (gen_comp @ rates)` per clade, so total expected synonymous equals
> total observed synonymous *by construction* — the 12 offsets within a clade are a
> near-fixed-sum composition (measured: per-clade Σactual/Σexpected = 1.025 ± 0.011
> against cell offsets spanning ±1.3). G→T's ~+60k pre-Omicron excess therefore
> *must* be paid for by deficits elsewhere, and T→C / G→A / A→G are where it lands.
> They look like era effects because G→T is one. Do not count them as corroboration.

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

---

## v10 outcome: what the split fixed, and what it did not

Measured on `gisaid-and-public_2026-04-22_v10` against the same `_v9` tree, over the
same population (cells with `Σ actual_count ≥ 1000`). Raw per-cell tables:
`gt_offsets_v9.csv`, `gt_offsets_v10.csv`.

### The era structure in G→T is essentially gone — which is what the split was for

| | pre-Omicron mean | Omicron mean | era gap | separation |
|---|---|---|---|---|
| v9 (pooled fit) | +0.6245 | −0.1217 | **+0.3928** | **COMPLETE** |
| v10 (per-era fits) | +0.2762 | +0.1817 | −0.2586 | overlapping |

The complete separation that justified this issue is destroyed: G→T's era difference
falls **+0.7463 → +0.0945, an 87% reduction**, and the two eras' ranges now overlap
heavily. Across all twelve types, mean |era difference| falls **0.1740 → 0.0903**.

### Criterion 1 was not met

> *"Over cells with `Σ actual_count ≥ 1000` (28 G→T cells), `max |offset|` for G→T
> drops below **0.20**, against a current value of 0.796."*

| | max \|offset\|, G→T | verdict |
|---|---|---|
| v9 | 0.7961 | (baseline, reproduces the pre-registered 0.796 exactly) |
| v10 | **0.4482** | **FAIL** — 44% better, but the bar was 0.20 |

**The residual is not era structure and not the site-restriction artifact.** Two
measurements locate it:

- *Not era structure*: 87% of the G→T era difference is gone (above). What remains is
  common to both eras — v10's Omicron cells moved **up** (−0.122 → +0.182) as its own
  fit stopped being dragged by pre-Omicron data, so both eras now sit positive.
- *Not the site restriction* (`rates.py`'s `nt_site.isin(sites)`, which upstream lacks
  and which ISSUE-137 identified as the source of the 0.9–5.5% pooled residue): the
  G→T offset persists **on** the training sites, at +0.417 max over 27 cells, where
  criterion 2's identity holds at **2.8e-15**.

Those two coexist because `T` is a **single scalar per clade** balancing *all*
synonymous types together. G→T can be under-predicted while other types are
over-predicted and the clade total still nets to zero by construction. So the residual
is a **within-era, per-mutation-type** bias in the GLM's relative rates — a different
defect from the one this issue addressed, and one an era split cannot reach.

The likely mechanism, flagged in ISSUE-137's own criterion 6: `glm.py` fits
`log(count + 0.5)` by least squares and exponentiates the prediction, so it recovers a
conditional *geometric* mean where the offset compares *arithmetic* sums. By Jensen
that under-predicts, worst for the most over-dispersed type — which G→T is. Worth a
follow-up issue; it is not a regression, and v10 is strictly better than v9 on this
measure.

### Four types gained era-structured offsets

`AT` moved +0.0037 → −0.1740, and `CG`, `AC`, `AG` moved similarly if less. These all
started near zero, so percentage reductions are meaningless for them. Note mat-pcp's
*fitted* AT era ratio (+0.1306) matches upstream's (+0.1294) almost exactly — so by
upstream's standard these rates are right, and the new offset says the pooled model
had been accidentally well-calibrated for AT. Same follow-up.

### The other criteria

| criterion | result |
|---|---|
| 2 — training-site identity survives | **PASS**, 2.8e-15 (threshold 1e-6) |
| 4 — era difference agrees with upstream | **PASS**, all three parts; see `upstream/README.md` |
| 5 — single-era degeneracy exact | **PASS**, byte-identical on the test tree (stronger than the 1e-10 asked) |
| 6 — runtime ≤ 2× | **PASS**, 3.27 s pooled → 3.28 s + 3.32 s per era |
| 1 — G→T calibrates below 0.20 | **FAIL**, 0.448 |
| 3 — held-out fit improves | not measured |
