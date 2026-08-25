# #146: chronumental date-prior sweep on `public_2026-08-06_v11`

Four `--variance_dates` arms at `--steps 1000` against
`results/public_2026-08-06_v11` (10,177,516 branches), plus the traceback from the
two horseshoe arms, which **cannot run in the pinned chronumental 0.0.65** (see
below). Produced by `workflow/scripts/chronumental_sweep.py`; each arm reads the
chronumental inputs already in the output tree and writes to a scratch directory, so
no pipeline rule ran and no production output was touched.

Committed for the same reason as [#139](https://github.com/matsengrp/mat-pcp/issues/139)'s
artifacts: the decision in
[#146](https://github.com/matsengrp/mat-pcp/issues/146) is made from these numbers,
and they would otherwise exist only in a commit message, unverifiable without ~7.5 h
of compute against a tree that is not in the repo. The per-arm `node_dates.tsv`
(351 MB each) and `node_dates_residuals.csv` (284 MB each) are **not** here.

| file | contents |
|---|---|
| `sweep_table.md` / `sweep_metrics.csv` | one row per arm, every metric |
| `<arm>/node_dates_qc.tsv` | the shipped `chronumental_qc` report for that arm |
| `<arm>/node_dates_branch_dt_summary.csv` | per-bin Δt summary for that arm |
| `<arm>/chronumental_steps.tsv` | every 100th step line, plus the final one |
| `<arm>/cmd.txt` | the exact command |
| `B0_horseshoe_traceback.txt` | the upstream crash, verbatim |

Reproduce:

```bash
pixi run python workflow/scripts/chronumental_sweep.py run \
  --tree-dir results/public_2026-08-06_v11 --out-dir <scratch> --steps 1000 --jobs 4
pixi run python workflow/scripts/chronumental_sweep.py qc \
  --tree-dir results/public_2026-08-06_v11 --out-dir <scratch>
pixi run python workflow/scripts/chronumental_sweep.py report --out-dir <scratch>
```

## Results

| arm | RTT/fitted | length_cor | falling | med \|res\| | p95 \|res\| | 0-mut time % | branch-years | >90 d % |
|---|---|---|---|---|---|---|---|---|
| A0 baseline (`--variance_dates 0.3`) | **4.31x** | 0.3964 | **99/100** | 1 d | 1 d | 45.9% | 772,669 | 2.87% |
| A1 `--variance_dates 3` | 4.29x | 0.4377 | 34/100 | 1 d | 1 d | 43.4% | 774,346 | 3.13% |
| A2 `--variance_dates 10` | 3.62x | 0.5571 | 18/100 | 4 d | 7 d | 33.2% | 652,137 | 2.95% |
| A3 `--variance_dates 30` | **1.48x** | 0.8476 | 11/100 | 10 d | 35 d | **1.4%** | 268,546 | 0.51% |

"falling" is the number of the 100 recorded intervals over which `length_cor`
decreased; "0-mut time %" is the share of total branch-time held by zero-mutation
branches. Every arm's tip-residual gate **passes** (the gate is on the median,
threshold 30 d) and every arm's root date lands inside check 7's window.

**The baseline arm reproduces the shipped fit exactly** — `length_cor:0.3964`,
`mutation_rate:9.3575`, identical to
`results/public_2026-08-06_v11/log/chronumental.stdout`. That is what licenses
reading anything off the other three.

**`--variance_dates 30` removes essentially all of the inflation.** The tree's total
elapsed time falls from 772,669 to 268,546 branch-years — a 65% cut — and the
zero-mutation intercept goes from 16 d to **0 d**, taking the zero-mutation share of
branch-time from 45.9% to 1.4%. The aggregate rate moves from 0.23x to 0.68x
chronumental's own root-to-tip estimate, and branches longer than a 90 d window fall
from 2.87% to 0.51%.

**The decay also stops.** At baseline `length_cor` falls at 99 of 100 intervals and
is still falling at step 999 — the fit is actively optimizing away from the
strict-clock solution, so more steps make it worse. At `--variance_dates 30` only 11
of 100 intervals fall. This matters for the five configs that inherit chronumental's
`--steps 20000` default: the baseline's problem compounds with steps, and this arm's
does not appear to.

**The cost is the residual tail, and it is not small.** Median tip residual goes
1 d → 10 d and p95 1 d → 35 d, and the share of tips predicted *later* than their
observed date goes from 1.6% to **49.3%** (234,851 of them by more than 30 d). The
median still passes the 30 d gate with 3x margin, but the distribution has moved
from "essentially exact, with a 1.6% tail" to "centred near the gate". Whether that
is acceptable is a judgement about what the windowed outputs need, not something
these numbers settle on their own.

Note the trade-off curve is steep between arms A2 and A3 (3.62x → 1.48x) and no arm
between 10 and 30 was run, so the knee is uncharacterised.

## The horseshoe arms cannot run in chronumental 0.0.65

Both `--model HorseShoeLike` arms crash after ~10 minutes of setup, at the first
logging step, before any fitting:

```
File ".../chronumental/__main__.py", line 380, in <listcomp>
    f"{name}:{value:.4f}"
ValueError: Unknown format code 'f' for object of type 'str'
```

This is an upstream bug, not a configuration problem, and it is specific to data
like ours. `HorseShoeLike.get_logging_results` sets
`results['argmax_variances_param'] = self.terminal_names[argmax(...)]` — a **strain
name**. The printer at `__main__.py:379-384` formats any value whose `str()`
contains a `.` with `:.4f`:

```python
f"{name}:{value:.4f}" if "." in str(value) else f"{name}:{value}"
```

Every strain name in these trees contains a dot, because they carry a GenBank
accession with a version suffix (`CHN/YN-0306-466/2020|MT396241.1|2020-03-06`) —
5000 of 5000 sampled. So the horseshoe model is unusable on SARS-CoV-2 MAT data in
this version without patching chronumental. It presumably works on flu because those
strain names have no dots.

The fix is one line upstream (format floats as floats and leave strings alone), and
it touches only the logging path, not the model. Testing the horseshoe here
therefore requires running a patched chronumental, which is a decision about the
pinned dependency rather than a sweep parameter.
