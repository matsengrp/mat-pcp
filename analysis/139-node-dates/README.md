# #139: internal-node date QC on the v9 production tree

The report and binned artifacts from running the shipped
`workflow/scripts/chronumental_qc.py` against
`results/gisaid-and-public_2026-04-22_v9` (20,961,614 nodes; 425 s, 43.5 GB max RSS).

These are committed because [#141](https://github.com/matsengrp/mat-pcp/issues/141)
picks thresholds from these distributions, and a pre-merge review pointed out that
the numbers otherwise existed **only in a commit message** — unverifiable without
re-running a 7-minute, 43 GB job against a tree that is not in the repo. The
per-tip `node_dates_residuals.csv` is *not* here: it is 1.6 GB.

Reproduce with the real rule (writes into the run directory):

```bash
pixi run snakemake -j 1 --use-conda --configfile <the 2026-04-22 config> \
  --forcerun chronumental_qc
```

| file | contents |
|---|---|
| `node_dates_qc.tsv` | the four gates (all PASS) and the four #139 REPORT rows |
| `node_dates_branch_dt.csv` | `(mut_count_bin, dt_days, count)` — check 5's per-bin day histograms |
| `node_dates_internal_checks.csv` | `(check, days, count)` — checks 6 and 8's day histograms |

## What the tree shows

**The clock converged.** Median Δt rises monotonically with branch mutation count
(16, 23, 36, 47, 58, 68 d for 0–5 mutations): a **median incremental slope of
11 d/mutation on a 16 d zero-mutation intercept**, against a nominal ~15.3. That was
the `--steps 1000` concern and it is answered. Note the high bins *saturate*
(11–20 mut → 80 d, below the 6–10 bin's 84 d) on ~650 branches — plausibly chronic
infections, and `max_nt_mutations` branch QC applies to the counting path, not the
timetree.

**The one real anomaly, and what it actually is.** Check 8 reports
`VIOLATING (>0) = 0.0314 (125031 of 3982424)`: internal nodes dated *after* the
earliest observed sample beneath them. But this is **not independent evidence about
internal placement.** Because `time_monotonicity` forces
`pred(ancestor) ≤ pred(tip)`:

```
0 < pred(a) - obs(t)  <=  pred(t) - obs(t)  =  residual(t)
```

so a positive gap *requires* that tip's own residual to be at least as large. Check
4 now reports its positive tail alongside the median the gate uses, and the two
agree exactly where the inequality says they must:

| | value |
|---|---|
| tips with `predicted > observed` | 208,995 (1.2309%) |
| by more than 30 d | 11,552 |
| **max tip residual** | **+1663 d** |
| **max check-8 gap** | **+1663 d** |

One tip supplies both maxima. So the finding is: **the `tip_residual` gate is a
median over 16,979,190 tips (reading 1.0 d against a 30 d threshold, p95 = 0 d) and
is insensitive by construction to a 1.23% tail reaching 1663 days.** Check 8's value
is that it converts that tail into the quantity that matters downstream — how many
*internal nodes*, hence how many branches, land in the wrong `windows/` bucket —
which is not readable off the residual distribution.

That reframes #141: the unit to attribute is ~209 K **tips**, not 125 K internal
nodes, and the likely fix is a tail statistic on residuals already computed rather
than a new gate on a derived measure.

## The algorithm test: Δt to the nearest descendant (check 9)

Checks 6 and 8 both compare an *inferred* date against an *observed* one, so a
violation is a three-way inconsistency between the inferred date, the topology and a
collection date — and check 8 is bounded by the tip residuals, as above. Neither is a
clean test of the dating algorithm. Check 9 compares **two chronumental outputs**:
given the topology it was handed, how much time does the fit put between an internal
node and its closest descendant?

Because monotonicity makes predicted dates nondecreasing along every path, the
nearest descendant is always an immediate **child**, so this is min-over-children of
the branch Δt — one scatter-min over the edges check 5 already builds, with no
propagation pass and no dependence on node ordering.

Over all 3,982,424 internal nodes:

| | value |
|---|---|
| **Δt = 0** (simultaneous with nearest descendant) | **177,348 (4.45%)** |
| Δt ≤ 1 d | 504,403 (12.67%) |
| p50 | 8 d |
| p95 | 37 d |
| p99 | 67 d |
| p99.9 | 130 d |
| max | 862 d |
| > 365 d | 42 (0.0011%) |

Nothing here is impossible — `frac<0 = 0.0000`, i.e. no internal node is dated after
its nearest descendant — and the shape is what a working clock looks like: most
internal nodes sit a week or two above their closest child. The two features worth a
follow-up judgement are the **4.45% dated simultaneously** with a descendant (which
put a branch's mutations on a window boundary in `modules/windowing.py`) and the
0.0011% with more than a year of unsampled time below them.

## Branch Δt by mutation count (check 5), full table

| mutations | n | p05 | median | p95 | max | frac Δt=0 |
|---|---|---|---|---|---|---|
| 0 | 12,103,637 | 1 | **16** | 71 | 1500 | 2.78% |
| 1 | 5,636,409 | 4 | **23** | 74 | 1031 | 0.63% |
| 2 | 1,819,065 | 8 | 36 | 95 | 1072 | 0.22% |
| 3 | 726,900 | 12 | 47 | 115 | 1098 | 0.12% |
| 4 | 322,823 | 15 | 58 | 135 | 839 | 0.08% |
| 5 | 154,602 | 17 | 68 | 154 | 781 | 0.07% |
| 6–10 | 172,560 | 17 | 84 | 191 | 946 | 0.07% |
| 11–20 | 24,961 | **12** | **80** | 260 | 1321 | 0.02% |
| 21–50 | 603 | 19 | 101 | 489 | 828 | — |
| 51+ | 53 | 27 | 135 | 800 | 1114 | — |
| **all** | **20,961,613** | 2 | **22** | 83 | 1500 | 1.80% |

A clean ~11 d/mutation increment on a 16 d intercept through the 6–10 bin, then it
breaks: the 11–20 bin's median (80 d) falls **below** the 6–10 bin's (84 d) and its
p05 drops to 12 d. High-mutation branches include some very short ones — consistent
with chronic infections or artifacts rather than clock time, and note
`max_nt_mutations` branch QC applies to the counting path, not the timetree.

**All 12,103,637 zero-mutation branches are terminal** — zero are internal, so check
6 covers that whole population rather than a leaf subset. Check 9 reports the split
explicitly so this stops being an inference from two matching row counts.

Children are never dated *earlier* than parents (`frac<0 = 0` on all 20,961,613
branches, `min = 0` in every bin), but **1.80% are exactly simultaneous**. That is
close to forced rather than discovered: chronumental parameterizes non-negative
durations, and `EPS_DAYS = 1` means a −1 day edge would pass `time_monotonicity`
undetected. Every bin's `min` of 0 does settle it, but the gate is not independent
evidence.

## Two report values that must not be quoted as measurements

- **`zero_mut_leaf_predicted_dt`'s `frac<0 = 0.0000`** is forced by
  `time_monotonicity`, a hard gate. It is an invariant, not a confirmed prediction.
  Read `min` instead: `EPS_DAYS = 1` admits a −1 day edge past that gate, and a
  4-decimal fraction over 12,103,637 observations hides up to ~600 negatives.
  (`min = 0d` here, which does settle it.)
- **Checks 6 and 8 use opposite sign conventions** — check 6 is
  `observed − predicted` (negative is bad), check 8 is `predicted − observed`
  (positive is bad). Both rows now print `frac<0`, `frac>0`, and an explicit
  `VIOLATING` side. Comparing the bare `frac<0` values across the two rows gives the
  opposite of the right conclusion: 2.30% bad versus 96.08% good.
