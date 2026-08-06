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
