# #135 — what the 46,284-row "absent from the founder grid" drop actually is

Reproduce with:

```bash
pixi run python workflow/scripts/audit_off_grid_drop.py \
  --run-dir results/gisaid-and-public_2026-04-22_v9 --outdir analysis/135-off-grid-drop/v9
pixi run python workflow/scripts/audit_off_grid_drop.py \
  --run-dir results/gisaid-and-public_2026-04-22_v8 --outdir analysis/135-off-grid-drop/v8
```

**Measurement only** — no pipeline behaviour changed (EPIC #138 decision 2). The
per-category drop/log/raise recommendation is criterion 3, presented below as a
reviewable table.

The script reconstructs the drop from the same two inputs
`build_count_table.py` uses and **reconciles exactly against the run's own log
line**: 46,284 rows / 210,111 counts.

## The answer: 100% is off-founder background, and 97.8% of it restores the founder

**All 46,284 rows** (100.00% of rows, 100.00% of counts) are observations whose
*from*-base differs from the clade founder base. Within that, the great majority
restore the founder base — a reversion *toward* the founder, which the grid cannot
hold by construction:

| | rows | % of rows | counts | % of counts |
|---|---|---|---|---|
| **off-founder at the nucleotide level** | **46,284** | **100.00** | **210,111** | **100.00** |
| — of which restore the founder base | 43,596 | 94.19 | 205,566 | **97.84** |
| — of which do not | 2,688 | 5.81 | 4,545 | 2.16 |

> **Do not read the 97.84 / 2.16 as a split between two causes.** An earlier
> version of this README did. The two cascade categories are **strictly nested**,
> not alternatives: `reversion_to_clade_founder` is `to_nt == founder`, and
> `off_founder_nt_background` is `from_nt != founder`. Since self-mutations are
> consumed earlier in the cascade, `to_nt == founder` *implies*
> `from_nt != founder` — verified at 100.00% of its rows on both the v8 and v9
> trees. Swapping two `assign` calls in the script reports "100% off-founder, 0%
> reversions" over the identical set of rows. The union is the finding; the
> reversion share is a sub-decomposition of it.

> **"Zero unexplained" is a theorem, not a measurement**, and the earlier version
> of this README quoted it as evidence of completeness. Off-grid is *defined* as
> `founder.isna() | (from != founder) | (to == founder)`. Category 1 consumes the
> first disjunct and categories 5-6 are exactly the other two, so their union is
> the whole remainder for **any** input — `unexplained` is empty by construction
> and the partition assertion is tautological. It is retained only as a tripwire
> on the mask definition in `load_run`. Quoting it as "nothing is unexplained,
> therefore the drop is fully understood" is CLAUDE.md pattern 10's vacuity: a
> check that cannot fail, reported as a check that passed.

What genuinely *is* established, and could have come out otherwise: the four
data-checking categories (missing founder row, non-ACGT founder base, non-ACGT
mutation base, self-mutation) are all **empirically zero** on this tree. Those can
fire. And #135's "open an issue for the unexplained remainder" does not apply —
not because `unexplained` is zero, but because the union above accounts for every
row against a predicate that was independently reconciled with the run's own log
line.

Grid membership is a three-part predicate, which is what makes this tractable.
`(clade, nt_mutation)` is on the grid iff `clade_founder_nts.csv` has a row for
`(clade, nt_site)`, the mutation's *from* base equals that founder base, and its
*to* base differs from it. Two facts collapse most of the possibility space:

- **Every (clade, site) has a founder row.** `clade_founder_nts.csv` is
  1,525,053 rows = 51 clades × 29,903 sites exactly. No site is missing.
- **Every founder base is A/C/G/T** (T 490,251 / A 456,763 / G 298,838 / C 279,201).
  No N, no ambiguity, no gap.

So the grid always offers three rows per (clade, site), and an observation can only
fail to match on the *bases*.

### `reversion_to_clade_founder` (97.8% of counts) — excluded by design, not missing

The observation's **to-base is the clade founder base**. The grid enumerates
mutations *away from* the founder, so a founder→founder row cannot exist. These are
branches that mutated back toward the founder after an ancestor mutated away from it.

Worked example, the single largest dropped observation:

```
clade 20I, site 28095, clade founder base = A

counts_all.csv grid rows at (20I, 28095):
  A28095C  actual_count     7  exclude=True
  A28095G  actual_count     3  exclude=True
  A28095T  actual_count  1118  exclude=True      <- the forward mutation

dropped observation:  T28095A  count 5974        <- the reversion back to founder A
```

`T28095A` is not a mutation away from 20I's founder; it is the return trip. There is
no A→A row for it to land in and there should not be. Dropping it is correct.
(The reversion outnumbering the forward mutation 5.3:1 at a site flagged
`exclude=True` is itself unsurprising — this is a masked, hypermutable site.)

The base-change profile is consistent throughout: the top rows are
T→C where the founder is C (17,272 rows / 86,326 counts), C→T where the founder is
T (6,163 / 25,136), T→G where the founder is G (4,801 / 21,397) — i.e. the
reversion of the commonest forward substitutions.

### `off_founder_nt_background` (2.2% of counts) — #130's case, one level down

Neither the from-base nor the to-base is the founder base: a lineage went
founder→X and then X→Y. This is the *nucleotide-level* analogue of the
*codon-level* filter #130 added, and it should be read next to it, not separately.
Also correctly dropped — it is not an observation of a mutation off the founder
background.

### Three causes named in #135 are structurally impossible

Worth recording so nobody re-investigates them:

- **Outside `site_include_range`** and **masked site**: `compute_exclusions` runs in
  step 3 of `build_count_table.py`, **after** the join, and sets the `exclude` /
  `masked_in_usher` flags on grid rows. It never removes a row, so it cannot move an
  observation off the grid. (The example above shows all three grid rows present
  with `exclude=True`.)
- **Non-coding**: the founder grid covers all 29,903 sites including non-coding ones,
  so being non-coding withholds nothing. 1,431 off-grid rows *are* at non-coding
  sites, but they are off-grid for the base reasons above, not for being non-coding.
- **Clade-founder mismatch between sources**: `aggregated.csv` carries its own
  `clade_founder_nt` column, and it agrees with `clade_founder_nts.csv` on **all
  1,271,502 observations — zero disagreements**.

## Criterion 4: the count is stable, and where it moved is fully attributed

| | v8 | v9 | change |
|---|---|---|---|
| off-grid rows | 46,757 | 46,284 | **−473 (−1.0%)** |
| off-grid counts | 221,539 | 210,111 | **−11,428 (−5.2%)** |
| `reversion_to_clade_founder` rows | 43,810 | 43,596 | −214 |
| `off_founder_nt_background` rows | 2,947 | 2,688 | −259 |
| unexplained | 0 | 0 | — |

**Not growing — shrinking, and 100.0% of the change is in the five clades #129
fixed.** `24H`, `24I`, `25A`, `25B`, `25C` account for −473 of −473 rows and
−11,428 of −11,428 counts. Every one of the other 46 clades is identical between v8
and v9 to the row and to the count.

That is a strong independent corroboration of #129: those five clades had the *bare
reference* as their founder in v8, which made a large number of ordinary
within-clade mutations look like they came off a non-founder base. Correcting the
founders removed exactly those, and touched nothing else. It also means the drop is
deterministic and tracks founder correctness — which is the reassuring answer for a
quantity nobody had checked.

v8's log has no drop line at all, so in v8 this was entirely silent. That is what
#135 was filed about.

## Criterion 3 (the gate): recommended per-category disposition

Not implemented — this is the reviewable table EPIC #138 decision 2 asks for.

| category | today | recommend | why |
|---|---|---|---|
| `reversion_to_clade_founder` | dropped, logged as "absent from the founder grid" | **log under its own name** | Correct to drop, but the current label reads as *something is missing*, when the grid provably cannot hold the row. 97.8% of the counts in that line are this, so the line as written overstates the problem by ~45×. Rename, and report it as a distinct quantity. |
| `off_founder_nt_background` | dropped, same log line | **log next to #130's codon-level line** | Correct to drop. It is the nucleotide-level twin of the filter #130 added; separating them makes the two comparable (2,688/4,545 nt-level vs 7,910/14,113 codon-level). |
| `no_founder_row` | dropped silently | **raise** | Zero today, and cannot legitimately occur — the founder table is a complete 51 × 29,903 grid. If it ever fires, the founder table is truncated. |
| `founder_nt_not_acgt` | dropped silently | **raise** | Zero today. Would mean a founder FASTA carrying N or a gap, i.e. #129's failure mode in a new form. |
| `mutation_nt_not_acgt` | dropped silently | **raise** | Zero today. Would mean indel or ambiguity notation reaching a table that assumes substitutions. |
| `self_mutation` | dropped silently | **raise** | Zero today. A from-base equal to the to-base is meaningless input. |

**The recommendation that matters is the four `raise`s.** The two benign categories
account for 100% of today's drop, and the four impossible ones are currently
*indistinguishable from them* — all six vanish into one log line. Any of the four
becoming non-zero would be a real corruption of the founder table or the counts,
and the present code would absorb it silently. That is the #129/#130/#133 shape,
and it is cheap to close because each is provably zero right now.

None of this changes an output, so if adopted it can land in v10 without moving
#137's baseline.

## The ORF1a question #132's review left open

#132 asked whether ORF1a's 18,596 "only in PCPs" rows
(`validation/ORF1a/pcp_vs_counts_report.txt`) could be attributed between this drop
and #130's filter, since both present identically as `actual_count == 0` on a grid
row. **They can be, partially, and the answer is that this drop dominates.**

Counting distinct `(clade, aa_mutation)` keys over ORF1a-covering rows that have
**no surviving observation** after each filter (both drop masks reproduced exactly —
46,284/210,111 and 7,910/14,113):

| source | keys touched | keys with no surviving observation |
|---|---|---|
| this off-grid drop | 18,109 | **9,792** |
| #130's codon-background filter | 2,499 | **1,865** |
| union | — | **11,634** |

So the two drops account for **11,634 of 18,596 (62.6%)**, and the pre-existing
off-grid drop outweighs #130's filter roughly **5:1**. The `exclude` flag adds at
most 971 more keys (ORF1a-covering keys whose every observed grid row is
`exclude=True`), so ~5,991 (~32%) remains unattributed.

Stated plainly rather than guessed at: the residual is **not** explained from the
count-table side. It needs the PCP-side enumeration — `validate_pcp_vs_aamut.py`
derives AA mutations from actual parent-child pairs through the ORF1a→ORF1ab codon
remapping, with its own branch QC, which is a different enumeration from
`aggregated.csv`. That is a separate investigation, and the numbers above are a
bound on this audit's contribution to it, not a decomposition of the whole
discrepancy.

## Files

Per tree (`v9/`, `v8/`):

| file | contents |
|---|---|
| `by_category.csv` | the decomposition, with rows, counts, percentages, and the reason each category exists |
| `by_clade_category.csv` | category × clade |
| `by_gene_category.csv` | category × gene (grid gene, `;`-joined at overlaps) |
| `by_base_change.csv` | category × (from, to, founder) base triple |
| `by_clade.csv` | totals per clade, used for the v8/v9 comparison |
| `top_200_by_count.csv` | the largest individual dropped observations |
