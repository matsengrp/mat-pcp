# Eleven of twelve mutation types show no era structure — the v10 split is a G→T fix

Follow-up from #137 Task 1, filed per [#138](https://github.com/matsengrp/mat-pcp/issues/138)'s
instruction to record non-separating mutation types rather than widen #137.

## What Task 1 measured

Task 1 computed the synonymous offset `log(Σ actual / Σ expected)` for all 51 clades ×
12 mutation types on `gisaid-and-public_2026-04-22_v9` and applied the pre-stated
complete-separation criterion, `min(offset over pre-Omicron) > max(offset over Omicron)`.
Full table and reproduction script: `analysis/137-era-offsets/`.

**G→T passes with a gap of +0.393**, reproducing #137's stated value and ranges exactly
(pre-Omicron +0.429…+1.123 over 19 cells, Omicron −0.344…+0.036 over 29). The v10 era
split is justified and #137 Task 2 should proceed.

**Eleven of twelve types fail.** #137 and #138 both predicted that *two* would (C→G and
G→C). Those two are indeed the worst — gaps −1.304 and −0.641 — but every type except
G→T fails.

| | separates | gap |
|---|---|---|
| GT | **yes** | **+0.393** |
| CT, TG, AG, AT | no | −0.103 … −0.288 |
| TC, GA, TA, CA, GC, AC | no | −0.459 … −0.683 |
| CG | no | −1.304 |

## Why the prediction was wrong, which matters for how this is read

The "two of twelve" expectation came from #137's four-clade `spread` table
(`max − min` across 21J/21L/23A/24A). That is dispersion, not era structure. A type can
have a large spread with both extremes inside one era — #137 itself notes exactly this
for C→G (extreme at 21L, Omicron) and T→G (extreme at 21F, pre-Omicron). And complete
separation over ~19 pre-Omicron cells against ~29 Omicron cells is far stricter than a
four-clade spread implies: a single overlapping cell out of ~48 defeats it.

So this is not evidence that Task 1 found something wrong. It is evidence that the
criterion was doing more work than the prediction accounted for.

Count-weighted era means make the shape clearer than the pass/fail does:

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

**G→T's era effect is 5.9× the next largest**, and five types have an era difference of
the opposite sign. The accurate one-line description of the v10 era split is therefore:
*G→T is the one mutation type where the pooled neutral model is badly miscalibrated and
the miscalibration is an era effect.*

## Two questions this raises

**1. C→G and G→C are badly miscalibrated in both eras, and an era split will not touch
them.** Their offsets are large and *positive* everywhere: C→G averages +1.105
pre-Omicron and +1.091 Omicron; G→C +0.773 and +0.686. An offset of +1.1 means ~3× more
synonymous C→G observed than the neutral model expects, in every clade. That is a
first-order calibration failure in a covariate the model *does* have, not a missing one.
Worth noting that C→G and G→C are two of the three types carrying a light-switch boundary
(the third is C→T), and #125 established that mat-pcp's boundary convention diverges from
upstream at site 21562 — where all three shipped `AT`/`CG`/`GC` boundaries sit. Whether
that is related is unknown and is the first thing to check.

**2. Splitting eleven types that show no era effect costs variance for no bias
reduction.** `glm.py` fits each mutation type independently, so a two-era partition
halves the training data per fit. For G→T that is clearly worth it. For A→G (difference
−0.012 over 456,799 counts) it is a pure variance cost. This is *not* an argument against
#137 Task 2 as scoped — parity with `SARS2-mut-fitness-v2` and the pooled control arm
(#138 decision 3) both argue for the straightforward two-fit split, and the pooled table
ships alongside so the comparison is available. But it argues for:

- reporting per-type rate uncertainties in #137's post-run validation, so the variance
  cost is measured rather than assumed negligible;
- taking seriously the "era as a factor in one pooled fit" option #138 lists as out of
  scope. Notes for whoever scopes it: `glm.py` has no formula interface (hand-built numpy
  design matrix in `create_data_matrix`), models are fit independently per mutation type
  so `era × mut_type` is free while `era × motif` needs new interaction columns, and the
  `interacting` constructor argument is dead code (stored at `glm.py:39`, never read).

## Success criteria

1. C→G and G→C's era-independent positive offset explained or attributed — with the
   light-switch boundary interaction at site 21562 checked first and ruled in or out.
2. Per-type rate uncertainty reported for the v10 two-era fits, so the variance cost of
   splitting the eleven no-effect types is a measured number.
3. A decision recorded on whether v11 keeps the split for all twelve types, restricts it
   to types with demonstrated era structure, or moves to era-as-a-factor in one fit.

## Out of scope

- Anything that would delay v10. #137 Task 2 proceeds as scoped; this issue is about what
  the measurement revealed, not a blocker.
- Re-baselining #127's frozen fixture to "match upstream" on the boundary convention. The
  one-site divergence is deliberate (#125); if criterion 1 implicates it, that is a
  finding to write down, not a fixture to change.

## Related

- **#137** — the era split itself. Task 1 is what produced these numbers.
- **#138** — the v10 EPIC. Its risk item "two of twelve mutation types are expected to
  fail complete separation" should be corrected to eleven.
- **#125** — the light-switch boundary convention, relevant to criterion 1.
- **#135** — changes which observations reach the GLM, so it can move every number here.
  Any re-measurement must state which count table it used.
