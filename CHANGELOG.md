# Changelog

This file tracks changes to the mat-pcp pipeline outputs and schema that affect downstream consumers.

## Version Policy

- **Simple integer versioning**: v1, v2, v3, etc.
- **Version bumps**: Any change to output format, filtering logic, or schema that affects downstream consumers should bump the version
- **Breaking changes**: In our use case, essentially any change is a breaking change since models trained on one version's outputs may not work with another's
- **Version in outputs**: The pipeline version is automatically appended to the output directory (e.g., `results/public_2022-06-01_v1`)

## How to Bump Version

1. Update the version number in the `VERSION` file in the repository root
2. Document the changes in this CHANGELOG under a new version heading
3. Update `OUTPUT_SCHEMA.md` if the schema has changed
4. Commit the changes and create a git tag matching the version (e.g., `git tag v2`)

## Versions

### v11

**The reason for the bump is `data/usher_masked_sites.yaml`: the UShER clade mask
is regenerated against the 2026-08-06 tree, so masking flags in `counts_all.csv`
move, and with them the GLM training set, `fitness/` and `aamut_fitness/`.** Models
trained on v10 fitness outputs should not be mixed with v11 ones.

The mask file is a generated artifact (see
[`s2-clade-mask`](https://github.com/matsengrp/s2-clade-mask)) whose masks are
accumulated down the tree topology, so regenerating it against a newer MAT can move
masks for clades that already existed, not only add the new ones. The vendored UCSC
spec `branchSpecificMask.yml` is unchanged (25 lineages, same pinned commit) — every
difference comes from the tree.

This is what forces the bump rather than the two script fixes below, neither of
which moves an existing output on its own: once the mask file differs, the
`results/*_v10` trees are no longer reproducible from `HEAD`, and the version is the
only thing distinguishing them. `gisaid-and-public_2026-04-22` remains a **v10**
artifact, built with the v10 mask file, and is not rebuilt.

- **Clade era mapping**: `workflow/modules/era.py`'s `CLADE_ERA` gains `25D`-`25I`
  (all Omicron). The 2026-08-06 trees carry 58 Nextstrain clades against the
  2026-04-22 tree's 52, and `era_for_clade` raises on an unmapped clade rather than
  defaulting (#137), so the run stopped at `era_fit_summary` until they were added.
  No existing clade's era changes. `tests/test_era_mapping.py` now checks a frozen
  fixture per tree, with the "maps nothing extra" direction taken against the union
  so it stays two-sided.
- **Unparseable dates no longer abort the run**: `filter_metadata.py`'s
  `robust_fromisoformat` guards all four branches, not only the exact-length-10 one.
  The 2026-08-06 metadata carries 231 rows per tree with the Solr-style range
  `[2022 TO 2024]`, which raised `ValueError` and killed the run in `filter_metadata`.
  `determine_resolution` already classifies such a value `MISSING`, which every
  config lists in `invalid_date_resolutions`, so the intended disposition was always
  to drop it. The count and the distinct offending values are now reported, because
  an unparseable date means an unfamiliar metadata format rather than an ordinary
  missing value. Verified over every distinct date value in all three trees: no
  previously-parsed value changes.
  A **config** date is held to the opposite standard: `config_date` parses
  `min_leaf_date` and raises rather than returning `NaT`. Widening the row
  parser's guard had made a malformed `min_leaf_date` silently disable the
  min-date filter -- it is used only as `date < parsed`, and every comparison
  against `NaT` is False, so zero samples were dropped at exit 0 -- while
  `max_leaf_date`, parsed with a bare `datetime.fromisoformat`, still failed
  loudly. The accepted formats are unchanged. Both halves are pinned by the
  new `tests/test_date_parsing.py`, which execs the real script.
- **New trees**: `config/config-2026-08-06-six-genes.yaml` adds
  `public_2026-08-06` and `gisaid-and-public_2026-08-06`. Kept as a separate config
  from the 2026-04-22 one — `roemer_nextstrain_to_pango` is passed whole as a rule
  param to `clade_founder_fasta`, so adding an entry for a new clade invalidates
  every clade founder of every tree sharing the config. `25I` (BA.3.2) is excluded:
  it resolves in neither founder source, and since #129 that is a hard failure
  rather than a silent fallback to the reference.

### v10 (EPIC [#138](https://github.com/matsengrp/mat-pcp/issues/138))

**The reason for the bump is [#137](https://github.com/matsengrp/mat-pcp/issues/137):
the neutral rate model is now fitted once per era instead of once globally, so
`fitness/` and `aamut_fitness/` move on any tree spanning both eras.** Models
trained on v9 fitness outputs should not be mixed with v10 ones.

The other three items in this release change no output on their own — #135 and
#139 are measurement-only, #136 is an environment pin that reproduced v9's md5s
exactly — and ride the coordinated bump rather than carrying one each.

- **[#136](https://github.com/matsengrp/mat-pcp/issues/136)**: `workflow/envs/bte.yml`
  now pins `bte>=0.9.3` with `python>=3.10,<3.11`. `bte<=0.9.2`'s
  `MATree.get_haplotype()` keeps the *first* mutation at a twice-mutated position and
  discards later non-reversions, which would corrupt every PCP parent haplotype
  (`pcp_from_tree.py:206`) and every window founder mutation set
  (`windows/founder_muts.csv`). v9 was unaffected only by luck of the solve.
  **No output change**: the re-solved env moved (bte build `py310h184ae93_1` →
  `py310h1153794_2`, python 3.10.19 → 3.10.20) and a full test-data run still
  reproduced `data/test_data/expected_summary.json` exactly, including every md5.
  So this needs no bump of its own (#138 decision 4). New guard:
  `tests/test_bte_haplotype.py`.
- **[#135](https://github.com/matsengrp/mat-pcp/issues/135)** (measurement only, no
  behaviour change): `workflow/scripts/audit_off_grid_drop.py` +
  `analysis/135-off-grid-drop/` decompose `build_count_table`'s 46,284-row /
  210,111-count "absent from the founder grid" drop. **100%** of it is off-founder
  background at the nucleotide level (the twin of #130's codon-level filter), of
  which **97.8% of counts restore the founder base** -- a reversion the grid cannot
  represent by construction. Those two are strictly *nested*, not competing causes,
  so their ratio is decided by cascade order; and `unexplained` is unreachable by
  construction, so its being zero is a tautology rather than a completeness result.
  Both are stated as such in the README and the script docstring.
  Stable v8 -> v9 (-473 rows), with 100% of the change in the five
  clades #129 fixed. Criterion 3's disposition table is in the README; nothing
  adopted yet.
- **[#139](https://github.com/matsengrp/mat-pcp/issues/139)** (additive outputs,
  report-only): `chronumental_qc` now checks the **internal** node dates, not just
  tips. Five new `REPORT` rows in `mat/node_dates_qc.tsv` plus three new binned
  artifacts: `mat/node_dates_branch_dt.csv`,
  `mat/node_dates_branch_dt_summary.csv` (per-bin median / IQR / p95 for branch Δt
  binned by branch length in mutations, plus a pooled `all` row) and
  `mat/node_dates_internal_checks.csv`. They gate nothing — the existing
  `tip_residual` gate and its threshold are untouched. Internal dates are what
  `modules/windowing.py` uses to bucket every branch into a time window, and until
  now nothing checked them; `tip_residual` passing at `p95|res|=1.0d` over 17 M tips
  is near-circular, since tip collection dates are the priors chronumental is fit
  to. The five rows split into two kinds, which the report now labels: the checks
  comparing an inferred date to an *observed* one measure metadata consistency and
  downstream impact (and `internal_date_vs_min_descendant_observed` is *bounded* by
  the tip residuals, so it is an amplified view of their positive tail rather than
  independent evidence about internal placement), while
  `internal_dt_to_closest_descendant` compares two inferred dates and so tests the
  algorithm. On the v9 tree the latter reads median 8 d, p95 37 d, max 862 d, with
  4.45% of internal nodes dated *simultaneously* with their nearest descendant.
  Production distributions are committed under `analysis/139-node-dates/`, since
  #141 picks thresholds from them. New guard: `tests/test_node_date_qc.py`.
- **[#137](https://github.com/matsengrp/mat-pcp/issues/137) Task 1** (diagnostic only,
  no pipeline rule reads it yet): `workflow/modules/era.py` carries the explicit
  52-clade pre-Omicron/Omicron mapping and its year-prefix audit;
  `analysis/137-era-offsets/` records the 51 × 12 synonymous offset table. G→T
  separates completely by era (gap +0.393), so the era split is justified. Eleven of
  twelve types do *not* separate — see
  [#140](https://github.com/matsengrp/mat-pcp/issues/140).
- **[#137](https://github.com/matsengrp/mat-pcp/issues/137) Tasks 2 and 3 — the
  neutral model is now fitted per era.** This is the change the v10 bump is for:
  it moves `fitness/` and `aamut_fitness/` on any multi-era tree.

  **New outputs** (per tree, per era present): `glm/curated_counts_{era}.csv`,
  `glm/master_table_{era}.csv`, `rates_table_{era}.csv`, and
  `glm/era_fit_summary.csv`. The pooled `glm/master_table.csv`,
  `glm/curated_counts.csv` and `rates_table.csv` are **retained unchanged** — the
  pooled fit is the control arm of the per-era comparison, not legacy baggage.
  Three separate files rather than one with an `era` column: DASM2's
  `TabularRatesNeutralModel` reads rates in a row loop, so an unrecognized column
  would be applied last-row-wins with no error.

  **Changed outputs**: `fitness/expected_counts_by_clade.csv` and everything
  derived from it (`fitness/ntmut_fitness*`, `aamut_fitness/*`) now take each
  clade's `expected_count` and `tau_squared` from its **own era's** fit, dispatched
  on the clade label via `modules/era.py`. Upstream dispatches on the positional
  slice `clades[0:11]`, correct only while the counts file happens to arrive in
  Nextstrain clade order; mat-pcp does not.

  **Unchanged on a single-era tree, byte-for-byte.** `add_expected_count` is now
  `add_expected_count_by_era` with one era, so the pooled path *is* the degenerate
  per-era path. A full test-data run reproduced every v9 md5 —
  `curated_counts.csv`, `master_table.csv`, `rates_table.csv`,
  `expected_counts_by_clade.csv`, `ntmut_fitness_aggregated.csv`,
  `all_aamut_fitness.csv` — and `data/test_data/expected_summary.json` is
  unchanged. That makes the CI run a regression check on the split rather than a
  new path with no baseline.

  **Empty era → the rule is skipped, and the skip is recorded.** The GLM cannot
  fit an empty partition, and single-era trees are normal (the test config's three
  clades are all pre-Omicron; the BA.2/JN.1/XFG trees are single-era by
  construction). `eras_present` resolves the `get_clades` checkpoint and requests
  only non-empty eras, rather than emitting a header-only placeholder that would
  look like a fit. Because a missing file cannot distinguish "legitimately absent"
  from "dropped by a wiring mistake", `glm/era_fit_summary.csv` names every era
  with `fitted` true/false and its clade list.

  **Conserved sites are computed over all clades before the partition**, matching
  upstream. Computing conservation within era would give each era a different
  training site set and break comparability with the pooled fit.

  **`window_expected_counts` deliberately stays on the pooled table** (EPIC #138
  decision 7). The windowed path keys on `founder_id` with no clade to dispatch on,
  and a per-window date label is wrong exactly across the transition. `T_window`
  absorbs any constant rate factor, so only the split's relative structure leaks
  through. Documented at the call site; era dispatch there is a follow-up.

  **`OUTPUT_SCHEMA.md` now documents `glm/master_table.csv`**, which had shipped
  since v1 with no schema entry.

  New guards: `tests/test_era_rate_dispatch.py` (the *only* two-era coverage — no
  shipped config is two-era, so CI's pipeline run cannot exercise the dispatch; it
  execs the real `curate_counts.py` and `expected_counts_by_clade.py`, and is
  mutation-tested against four separate breaks of them) and
  `validate_test.py::validate_era_rate_tables`.

### v9

**Three of these fix outputs that were silently wrong in v8, so v9 is not an optional
upgrade for anyone consuming v8 fitness or count data.** All five clade founders affected by
[#129](https://github.com/matsengrp/mat-pcp/issues/129) and every clade affected by
[#130](https://github.com/matsengrp/mat-pcp/issues/130) produced plausible-looking output
with no error, and the two estimators were corrupted by different mechanisms, so no internal
cross-check disagreed.

**Changes**:

- **An unresolvable clade founder is now a hard failure, not the Wuhan reference
  ([#129](https://github.com/matsengrp/mat-pcp/issues/129)).**
  `clade_founder_fasta.py` builds each clade's founder by applying its defining mutations to
  the reference. When a clade resolved in none of its three sources it applied *zero*
  mutations and wrote the unmodified reference — exit 0, logged as ordinary output. Nothing
  downstream caught it because the wrong founder propagated consistently. It surfaced only
  when an outside consumer noticed implausible cross-clade correlations (24H/24I/25A/25B/25C
  in v8). Two guards now cover the two routes to that output: a missing source, and a
  resolved source yielding no effective mutations. Both raise. `allow_reference_identical_founders`
  is a new (normally empty) config key for clades whose founder legitimately equals the
  reference; the `reference_clade` is exempt automatically.
- **Every config now covers clades 24H–25C
  ([#129](https://github.com/matsengrp/mat-pcp/issues/129), absorbing
  [PR #124](https://github.com/matsengrp/mat-pcp/pull/124)).** `roemer_nextstrain_to_pango`
  gained entries for `24H`–`25C`, and the BA.2/JN.1/XFG tree configs landed. Each config
  carries its own copy of that map, which is *how* #129 happened: adding a clade to one
  config does not add it to the others.
- **Counts include only mutations off the clade founder's *codon*
  ([#130](https://github.com/matsengrp/mat-pcp/issues/130)).** matUtils reports each branch's
  mutations relative to its parent, so the same nucleotide substitution also appears on codons
  that had already mutated — a different amino-acid change sharing one nucleotide event.
  `build_count_table.py` joined observed counts onto the founder grid on
  `(clade, nt_mutation)` alone, dropping `codon_change` (the only column distinguishing the
  two), so off-founder observations were relabelled onto the founder codon and kept as extra
  rows for the same mutation: **6,415 duplicated keys / 13,146 rows in v8, every clade
  affected, concentrated in the highest-count mutations.** `delta_fitness_bayesian`
  inverse-variance-averaged the phantom row in while `delta_fitness_ratio` summed its expected
  count, so reported amino-acid uncertainty got *smaller*. The `from_founder_codon` filter and
  a one-row-per-mutation invariant are restored; a structural mismatch raises rather than
  reading as off-founder.
  **`actual_count` semantics therefore change for affected mutations.**
  **The v8 blast radius is wider than the duplicated keys.** `genome_composition` counts grid
  rows, so v8's 6,731 extra rows (4,581,890 − 4,575,159 = 13,146 − 6,415) inflated the
  denominator of every clade's `evol_time` normalization. Measured v8→v9 on the
  `2026-04-22` tree, per-clade `T` moves from **+5.91% (21J)** down to **−1.26% (21M)**, so
  **every expected count in an affected clade was scaled wrong, not just the 6,415 duplicated
  mutations.** For 21J S:A67V the fitted rate moved only +0.49% while `T` moved +5.91%,
  giving the +6.43% change in `expected_count` (105.629 → 112.416) — i.e. the normalization,
  not the rate table, is the dominant term. Treat v8 fitness values as unusable clade-wide
  rather than mutation-wise.
- **The light-switch boundary site is served the rate it was fitted with
  ([#125](https://github.com/matsengrp/mat-pcp/issues/125)).** `modules/glm.py` trains with
  `global_context = (nt_site > boundary)`, so `site == boundary` is fitted into the
  before-regime, but both label producers used `site < boundary` — and `not (site < b)` is
  `site >= b`, not `site > b`. The boundary site was fitted into one bucket and served the
  other's rate. The convention now has one definition,
  `modules/lightswitch.py::before_boundary` (`site <= boundary`), used by both the clade path
  and the windowed trajectory path. Affects one site per bounded mutation type (A→T at 21562:
  51 rows / 82 counts in v8). mat-pcp deliberately diverges from `SARS2-mut-fitness-v2` at
  site 21562, which carries the same inconsistency between its two codepaths.
- **PCP outputs are gzipped ([#134](https://github.com/matsengrp/mat-pcp/issues/134)).** The
  four PCP families — `{clade}_pcps`, `all_pcps`, `{clade}_pcps_dasm_ready`,
  `all_pcps_dasm_ready` — are now `.csv.gz`. **Every PCP path gained a `.csv.gz` suffix;
  contents and column schemas are unchanged.** `pandas.read_csv` and
  `dasm2.data.registry.load_dataframe` both handle it transparently, the latter by extension.
  Motivation: on `gisaid-and-public_2026-04-22_v8`, `pcps/` was 504 GB apparent (470 GiB), 93% of
  the whole 542 GB output tree; measured gzip level 9 gives ~125× on ORF1a and ~110× on S. On the
  v9 run `pcps/` is **8.3 GB apparent** — across *six* genes where v8 had two, so the ratio
  understates the per-file gain. `mtime=0` is pinned so re-runs are byte-identical, which matters because
  `data-central` manifests pin `md5`/`xxh3` per file.
  `{clade}_filtered_branches.csv` and `filtered_branch_ids.csv` stay plain — **note that this is
  what dominates `pcps/` for a short gene**: the gzipped PCP data for E is 7 MB, while its
  uncompressed rejected-branch ID lists are 592 MB, 98.8% of `pcps/E/`. Because
  `pcp_rows + reject_rows` is the same branch universe for every gene (8,408,804 on the
  `2026-04-22` tree), those lists grow as the gene shrinks, and across six genes they come to
  ~3.0 GB that no known consumer reads. Candidate for a follow-up.
- **ORF1b PCPs, and a validation that was passing vacuously
  ([#133](https://github.com/matsengrp/mat-pcp/issues/133)).** New optional config
  `pcp_gene_cds_blocks` gives explicit CDS blocks for PCP genes the reference GTF does not
  carry, adding `ORF1b` (nsp12–nsp16) as `[[13442, 13468], [13468, 21552]]` — 8112 nt / 2704
  codons, spanning the −1 ribosomal frameshift rather than starting after it, since nsp12's
  first 9 residues are translated in ORF1a's frame. **New additive outputs:** `pcps/ORF1b/`,
  `clade_founders/ORF1b/`, `ORF1b_codon_sites_to_mask.txt`, `validation/ORF1b/`. Fitness
  outputs are untouched — the override is deliberately not fed to `edit_ref_gtf` /
  `ref_coding_sites`.
  Separately, `validate_pcp_vs_aamut.py` filtered `aamut_fitness` on `gene == {gene}`, but
  `aamutfit.get_coding` drops ORF1a, so **through v8 the ORF1a validation compared against
  zero rows and reported no discrepancies because it had nothing to compare.** New
  `pcp_gene_to_fitness_gene` names the fitness gene and codon offset; the correspondence is
  derived from genomic positions (so ORF1a's four post-frameshift codons, 4402–4405, are
  excluded rather than mis-mapped onto ORF1b's frame) and an empty fitness slice is now fatal.
  **`validation/{gene}/pcp_vs_counts_report.txt` contents change**: the report names the
  fitness gene compared against, and for remapped genes the codon offset and corresponding-codon
  count; it also prints the compared row counts on both sides, and reports `VACUOUS` rather
  than `SUCCESS` when zero discrepancies coincide with zero matches.
- **The production config extracts every non-accessory protein.**
  `pcp_genes` for the `2026-04-22` tree is `[S, ORF1a, ORF1b, M, N, E]`: the four structural
  proteins plus the two replicase polyproteins. The accessory ORFs (3a, 6, 7a, 7b, 8, 9b, 10)
  are deliberately out of scope. `M`, `N` and `E` are single-block, codon-aligned GTF features
  (666 / 1257 / 225 nt), so only `ORF1b` needs `pcp_gene_cds_blocks`, and none of the three
  needs a `pcp_gene_to_fitness_gene` entry or a `gene_overlaps` entry.
  **Additive only** — each new protein adds `pcps/{gene}/`, `clade_founders/{gene}/`,
  `{gene}_codon_sites_to_mask.txt` and `validation/{gene}/`, reading no existing schema and
  changing no existing file.
- Guarded the GLM's synonymous-only training set and its equivalence to
  `SARS2-mut-fitness-v2`'s reference rate table
  ([#127](https://github.com/matsengrp/mat-pcp/issues/127)). Tests and assertions only; no
  output change. The neutral-rate restriction lives in `modules/load.py`, not in curation, so
  `glm/curated_counts.csv` deliberately holds every mutation class.
- `validate_test.py` no longer silently validates the previous version's outputs when the
  current `VERSION` directory does not yet exist
  ([#88](https://github.com/matsengrp/mat-pcp/issues/88)). It warns loudly and always prints
  which directory it validated. This triggers exactly at a version bump, which is when it
  matters.

**Not in v9**: `window_expected_counts`'s optimization
([#123](https://github.com/matsengrp/mat-pcp/issues/123)) shipped inside the v8 tag and its
entry stays under v8; v9 carries only its deferred cleanups (a helper extraction and a rename,
with the trajectory output bit-identical).

**Corrections to the released `### v8` entry for
[#123](https://github.com/matsengrp/mat-pcp/issues/123).** Recorded here rather than by
editing a tagged version's section, which stays as published:

- It says the rule was cut "from ~20 h to tens of minutes on the 2026-04-22 tree." The
  measured figure is **1:32:14** (5,534.6 s, 99.4% mean CPU, 9.5 GB max RSS, 1,495,938
  founder-windows). The ~20 h was a *projection*, never a measured baseline, so the speedup is
  ~13× against a projection.
- It calls the output "bit-for-bit identical to the old output" while citing a ~1e-16 max
  difference, which cannot both be true. Delta accumulation reorders floating-point sums, so
  the correct claim is *numerically identical to float precision* (`assert_frame_equal` at
  `rtol=1e-6`, `atol=1e-9`). v9's own cleanup of that code **is** bit-identical, having moved
  statements without reordering arithmetic — hence the care to distinguish the two.
- It names the fast path `ExpectedCounts._founder_rec`; v9 renames it to
  `_recompute_founder_touched_rows`.

### v8

**Always run metadata QC; remove `skip_qc_filtering` (issue
[#121](https://github.com/matsengrp/mat-pcp/issues/121)).** The `skip_qc_filtering`
option — a development convenience that made `filter_metadata` a pass-through to match
legacy SARS2-mut-fitness behavior — is removed. Metadata QC now always runs, so every
tree is filtered by date cutoff (`min_leaf_date` / `max_leaf_date`), date resolution
(MISSING/MONTH/YEAR dropped), host, and genome length before tree filtering. This is a
breaking change: production configs previously ran with QC skipped, so the set of
sequences entering the tree — and therefore all downstream counts, fitness, and windowed
outputs — changes.

**Changes**:

- Removed the `skip_qc_filtering` config key and the corresponding branch in
  `workflow/scripts/filter_metadata.py` / the `filter_metadata` rule params. Leftover
  `skip_qc_filtering` keys in any config are now silently ignored (unknown global keys are
  not validated).
- `max_leaf_date` (per-tree) now actually takes effect; it was a no-op whenever
  `skip_qc_filtering: true`. `subset_frac` likewise now applies on every run.
- Fixes the `chronumental_qc` crash on production configs: the QC path emits the
  `date_resolution_str` column that `chronumental_qc` requires, which the skip path never
  created (`ValueError: Usecols do not match columns` on the missing column).
- Removed the redundant `pcp_max_nt_mutations` / `pcp_max_reversions_to_ref` /
  `pcp_max_reversions_to_clade_founder` config keys. PCP extraction now reuses the shared
  Mutation-Counting keys (`max_nt_mutations` / `max_reversions_to_ref` /
  `max_reversions_to_clade_founder`), so PCP extraction and clade/window mutation counting
  filter branches with the same thresholds. (Config values were already identical, so this
  is a no-op on the shipped configs.)
- Fixed a crash in `window_expected_counts` (`modules/expected_counts.py`) on founders
  whose touched sites are all noncoding: `_touched_sites` now intersects against the coding
  keep-sites, so such a founder short-circuits to the reference table instead of building an
  empty frame whose codon/merge-key columns come back `float64` and raise
  `AttributeError: Can only use .str accessor with string values!`. Behavior-preserving on
  founders that touch any coding site (noncoding rows were always dropped downstream); it
  only surfaced on full-scale production trees, not the test data.
- Performance: `window_expected_counts` now scales to full production trees. The driver
  previously rebuilt and re-reduced the full ~90k-row genome-wide rate table for **every**
  founder-window (~1.5M on a full tree → ~20 h). It now (a) accumulates each window's
  per-ID rate as a *delta* on a once-computed reference base (`N·base − touched-site
  reference rows + founder recomputed rows`), and (b) recomputes each founder's touched
  codons via O(#touched-sites) dict lookups (`ExpectedCounts._founder_rec`) instead of
  per-founder pandas. On the test data this is bit-for-bit identical to the old output
  (max abs diff ~1e-16 in `expected_count` / `delta_fitness_ratio`), and the per-founder
  rate-table build drops ~50× (≈19 ms → ≈0.38 ms), cutting the rule from ~20 h to tens of
  minutes on the 2026-04-22 tree. Output schema is unchanged. (#123)

### v7

**Time-resolved fitness (EPIC [#115](https://github.com/matsengrp/mat-pcp/issues/115)).**
The coordinated version bump for the three phases below. Each phase was purely
additive (new output files, reading no existing schema and read by no current
downstream consumer) and shipped no individual bump per maintainer decision; this
is the single bump that lands the whole EPIC. The new outputs are `mat/node_dates*`
+ `mat/chronumental_timetree.nwk` (Phase 1), `windows/subtree_counts.csv` +
`windows/founder_muts.csv` (Phase 2), and `windows/fitness_trajectory.csv`
(Phase 3). Everything that existed at v6 is byte-identical apart from the version
dir name.

**chronumental node dating** — Phase 1 of EPIC
[#115](https://github.com/matsengrp/mat-pcp/issues/115) (issue
[#116](https://github.com/matsengrp/mat-pcp/issues/116)). Adds time estimates for
every node of the filtered tree; nothing existing changes.

**Changes**:

- New per-tree outputs under `mat/`: `node_dates.tsv` (chronumental-inferred
  date for every tip and internal node), `chronumental_timetree.nwk` (the dated
  time-scaled tree), plus the always-on QC artifacts `node_dates_qc.tsv` and
  `node_dates_residuals.csv`.
- New rules `export_newick`, `chronumental`, `chronumental_qc` and a new conda
  env `workflow/envs/chronumental.yml`.
- New config keys `chronumental_kwargs` and
  `chronumental_max_median_abs_tip_residual_days`.

**sliding-window subtree counts** — Phase 2 of EPIC
[#115](https://github.com/matsengrp/mat-pcp/issues/115) (issue
[#117](https://github.com/matsengrp/mat-pcp/issues/117)). Buckets observed branch
mutations into overlapping time windows over the full dated tree and emits, per
window, each subtree founder plus per-mutation actual counts. No expected counts
or `delta_fitness_ratio` yet (Phase 3).

**Changes**:

- New per-tree outputs under `windows/`: `subtree_counts.csv` (per-window,
  per-mutation actual counts — coding only, with the same multi-gene `;`-collapse
  as `aggregated.csv`) and `founder_muts.csv` (per-window subtree founders +
  reconstructed haplotypes).
- New module `workflow/modules/windowing.py` (`WindowedMAT` + pure helpers) and
  thin wrapper `workflow/scripts/window_subtree_counts.py`; new rule
  `window_subtree_counts` (bte env). Windows are anchored at the end date and
  tiled backwards, so the most recent window is always the full `window_size`.
- Branch QC is identical to the clade path (`count_mutations.py`): it reuses the
  existing `max_nt_mutations` / `max_reversions_to_ref` /
  `max_reversions_to_clade_founder` keys and the unconditional
  drop-branches-with->1-nt-mutation-per-codon filter.
- New config keys `window_size`, `window_step`, `window_start_time`,
  `window_end_time`.
- New correctness harness: `workflow/scripts/validate_windowing.py` (pure +
  output-CSV checks + test-data count regression sentinels, base env) and the
  CI-only `validate_windowing` rule
  (`workflow/scripts/validate_windowing_tree.py`, bte env). The tree-level rule
  is **not** in `results_files`.

**windowed expected counts + `delta_fitness_ratio`** — Phase 3 of EPIC
[#115](https://github.com/matsengrp/mat-pcp/issues/115) (issue
[#118](https://github.com/matsengrp/mat-pcp/issues/118)). Turns the Phase-2
per-window founders + actual counts into a mutation × time-window fitness
trajectory: for each window's subtree founder it enumerates all possible coding
mutations from the founder sequence and looks up mat-pcp's GLM neutral rate for
each mutation's founder-sequence context. Expected counts use mat-pcp's **own
clade-path `evol_time` normalization** (`rates.py`) applied **per window** — the
same neutral-model application as the clade-based `aamut_fitness`, so the two are
directly comparable — rather than s2trajectory's `rate × subtree_mut_count`. Then
`delta_fitness_ratio = log((actual + p) / (expected + p))` per mutation per
window (matching `aamutfit.ratio_fitness`).

**Changes**:

- New per-tree output `windows/fitness_trajectory.csv` (per `(window_start,
  amino-acid ID)`: actual count, expected count, `delta_fitness_ratio`). Coding
  only; ORF1a dropped in favor of ORF1ab (aamutfit convention).
- New module `workflow/modules/expected_counts.py` (a pure-`pandas`/`numpy`
  founder-context enumeration ported from s2trajectory's `PossibleMutations` — no
  `bte`) and thin wrapper `workflow/scripts/window_expected_counts.py`; new rule
  `window_expected_counts` (default env, after the GLM stage).
- Expected counts come **only** from the GLM (via `glm/master_table.csv`); there
  is no synonymous-average fallback. Per window,
  `T = (observed synonymous count) / (rate summed over all possible synonymous
  mutations)` and `expected(ID) = T × Σ_founders rate(context)`, so total expected
  synonymous equals total observed synonymous **by construction** (the calibration
  that makes synonymous `delta_fitness_ratio` center at 0). Because `T` divides out
  the rate's absolute scale, mat-pcp's `rate` is used as-is (no `n_sites`/fraction
  renormalization).
- New shared helper `masking.global_masked_sites` — the single definition of the
  clade-independent global site mask (config `sites_to_exclude` ∪ ProblematicSites
  VCF ∪ outside `site_include_range`), reused as the windowed path's "bad sites".
- **Founder→reference reversions are masked** on the windowed path, mirroring the
  clade path's `to_nt == ref_nt` exclusion (`build_count_table`): a possible
  mutation whose destination is the reference base (arising only at
  founder-mutated sites) is dropped from both the expected enumeration and the
  actual counts, so these artifact-prone reversions are treated identically to
  the clade-based path and windowed fitness stays comparable to `aamut_fitness`.
- New config key `window_min_subtree_size` (default 0 — keep all subtrees).
- `validate_windowing.py` gains a per-window synonymous-calibration proof
  (`sum(actual_syn) == sum(expected_syn)` per window, exact by construction), a
  synonymous-centering sanity check, a formula check, and test-data regression
  sentinels for the trajectory.

### v6

**Multi-protein pipeline** — generalizes PCP extraction, founder slicing, DASM
codon masking, DASM-ready filtering, and PCP-vs-counts validation from a single
hardcoded gene (Spike) to a config-driven list of proteins. Tracks EPIC
[#84](https://github.com/matsengrp/mat-pcp/issues/84).

**Changes**:

- **`pcp_genes` list replaces the scalar `pcp_gene_region` + `masking_gene`
  config keys.** A single Snakemake invocation fans PCP extraction, gene-sliced
  founders, codon masks, DASM-ready filtering, and validation over every protein
  in `pcp_genes` (must be non-empty). Production `config/config.yaml`
  ships `pcp_genes: [S]` to keep v5->v6 Spike outputs byte-identical; the
  test configs exercise the multi-protein path with `[S, ORF1ab, E]`.
  Per-protein coordinates are derived from the reference GTF
  (`cds_blocks_for_gene`) rather than a hardcoded region.

- **Per-gene output layout.** Outputs that were previously Spike-only now live
  under a per-gene subdirectory or filename (gene symbol verbatim):
  - `pcps/{gene}/{clade}_pcps.csv`, `pcps/{gene}/all_pcps.csv`, and the
    `pcps/{gene}/*_pcps_dasm_ready.csv` / `*_dasm_filter_report.txt` /
    `*_dasm_filter_stats.json` / `*_filtered_branches.csv` /
    `filtered_branch_ids.csv` / `pcp_curation_report.txt` siblings.
  - `clade_founders/{gene}/{clade}.fa` (was `clade_founders/spike/{clade}.fa`;
    the `GENE_FOUNDER_DIRNAMES` `S`->`spike` mapping is retired).
  - `{gene}_codon_sites_to_mask.txt` (was `spike_codon_sites_to_mask.txt`).
  - `validation/{gene}/pcp_vs_counts_report.txt` (+ `pcp_aa_counts.csv`,
    `pcp_vs_counts_discrepancies.csv`), one set per protein.

- **ORF1ab is handled as one frameshift-reconstructed CDS**, not split into
  nsps: its two GTF CDS blocks share nt 13468, giving a single in-frame
  21288 nt = 7096 codon founder slice. Gene-sliced founders are now emitted
  natively by the pipeline, superseding the post-hoc `slice_spike_founders.sh`
  step ([#97](https://github.com/matsengrp/mat-pcp/issues/97)).

- **Delta-fitness columns renamed in the amino-acid fitness outputs**
  ([#89](https://github.com/matsengrp/mat-pcp/issues/89)): `naive_delta_fitness`
  -> **`delta_fitness_ratio`** (the primary per-mutation signal) and
  `delta_fitness` -> **`delta_fitness_bayesian`** (inverse-variance-weighted
  estimate). Pure rename (values unchanged), applied consistently across the
  per-clade `{clade}_aamut_fitness.csv` files and the aggregate
  `all_aamut_fitness.csv`.

- **Consolidated production config**
  ([#113](https://github.com/matsengrp/mat-pcp/issues/113)): the four per-tree
  production configs are merged into a single four-tree `config/config.yaml`.
  `min_clade_size` now defaults to `0` for every tree; the `public_2022-06-01`
  tree was previously `10000`, which biased its GLM rates
  ([#81](https://github.com/matsengrp/mat-pcp/issues/81)) — its v6 rate outputs
  therefore differ from v5. Removed the `masking_gene` and `clade_synonyms`
  config keys and the soft `min_clade_size` Snakefile guardrail.

- **`get_problematic_sites` uses exact gene-membership matching**
  ([#114](https://github.com/matsengrp/mat-pcp/issues/114)): compound `gene`
  cells (e.g. `S;ORF3a`) are split on `;` and matched exactly, so a gene symbol
  can no longer bleed into a longer one (e.g. `ORF1a` into `ORF1ab`) once fanned
  over `pcp_genes`. Latent-correctness fix; no shipped config had colliding
  symbols, so the S mask is byte-identical (17 sites).

- **Internal cleanups** (no output-schema impact): consolidated the by-path
  `load_workflow_module` helper into a single `workflow/modules/workflow_loader.py`
  ([#111](https://github.com/matsengrp/mat-pcp/issues/111)); removed the orphan
  `validate_pcp_vs_counts.py` script and a dead `ref_fasta` rule input
  ([#94](https://github.com/matsengrp/mat-pcp/issues/94)); `validate_test.py`'s
  default output dir is VERSION-aware
  ([#88](https://github.com/matsengrp/mat-pcp/issues/88)).

- **Spike outputs are byte-identical to v5** for every Spike-specific artifact
  (PCPs, gene-sliced founders, codon mask, DASM-ready PCPs) apart from the delta
  column renames above and the version dir name. Regression anchors held: 19A
  Spike founder md5 `b948f2a1f20edbce64cfff1dc540c119`, S DASM-ready aggregate
  md5 `5e538bd088e676b485a11a04745221e9`, 17-site S codon mask.

### v5

**Run multiple trees from one config** — tracks issue
[#104](https://github.com/matsengrp/mat-pcp/issues/104).

**Changes**:

- **New `trees:` config schema**: the pipeline now processes a mapping of named
  trees in a single Snakemake invocation. Shared parameters (reference genome,
  GLM, genes, masking, thresholds, fitness, clade settings, and all
  `qc_filter_options`) stay at the top level; a new top-level `trees:` mapping is
  keyed by a short tree **label**, and each entry supplies that tree's `mat` and
  `metadata` (both required). `max_leaf_date` is the only per-tree-overridable
  key; it is optional and defaults to `None` (no upper date cutoff). When a
  global `qc_filter_options.max_leaf_date` is set, a tree without its own value
  inherits it.

- **`output_dir` removed; derived automatically**: output paths are now
  `{output_root}/{label}_{version}` (e.g. `results/public_2022-06-01_v5`), where
  `output_root` is an optional global key defaulting to `results`. The scalar
  `output_dir`, `mat`, and `metadata` top-level keys no longer exist.

- **Parameter renames** (parallel naming for the two leaf-date bounds):
  `qc_filter_options.datetime_cutoff` → `max_leaf_date` and
  `qc_filter_options.min_date` → `min_leaf_date`.

- **Site masking now applied to the fitness outputs**
  ([#106](https://github.com/matsengrp/mat-pcp/issues/106)): masked mutations
  (`exclude == True` in the count table — the union of UShER/ProblematicSites-VCF
  masks, config `sites_to_exclude`, sites outside `site_include_range`, and
  founder→reference reversions) are now dropped in `ntmut_counts_by_clade.py`
  before nucleotide and amino-acid fitness are computed. Previously the fitness
  branch flagged but never removed these rows, so `ntmut_fitness_aggregated.csv`
  and `all_aamut_fitness.csv` included artifactual counts from masked sites; the
  masking had only ever been applied to the GLM training data (`curate_counts`).
  Rows for masked mutations no longer appear in the fitness outputs. The
  filtering is per `(clade, mutation)`, so the `"all"` cross-clade aggregate
  still counts a mutation from clades where it is not masked.

- **`"all"` cross-clade nucleotide aggregation now pools by mutation identity.**
  `ntmut_counts_by_clade` previously grouped the `"all"` cluster by *every*
  non-count column, including clade-founder-derived annotations (`motif`,
  `clade_founder_codon`/`aa`, `mutant_codon`/`aa`, `aa_mutation`, `synonymous`,
  `four_fold_degenerate`, `mut_class`) that can differ between clades for the
  same nucleotide mutation. That silently split a mutation into multiple `"all"`
  rows (one per context) and fragmented its pooled counts. The `"all"` cluster
  now keys strictly on `(nt_site, nt_mutation)`, summing counts across clades
  into one row per mutation; the differing annotation columns keep a
  representative (first contributing clade's) value. Per-clade rows are
  unchanged. In the test data this collapses 36 fragmented mutations (all
  zero-count) in `ntmut_fitness_aggregated.csv`.

- **`exclude` / `masked_in_usher` columns removed from the fitness outputs.**
  After masking is applied, every surviving row has `exclude == False`, so these
  now-constant flags are dropped in `ntmut_counts_by_clade` rather than carried
  as uninformative columns into `fitness/ntmut_counts/*`, `fitness/ntmut_fitness/*`,
  and `fitness/ntmut_fitness_aggregated.csv`. (They remain in `counts_all.csv`.)

- **`tau_squared` removed from the amino-acid fitness outputs**
  (`aamut_fitness/{clade}_aamut_fitness.csv`, `aamut_fitness/all_aamut_fitness.csv`).
  `tau_squared` is a nucleotide-level log-scale variance of the neutral expected
  count, consumed only by `probfit` when estimating nucleotide fitness. It never
  fed the amino-acid `delta_fitness` (an inverse-variance weighted average of the
  nucleotide estimates) or `naive_delta_fitness` (from summed counts); it was
  merely summed over the constituent nucleotide mutations and carried as an
  unused column. It remains in the nucleotide outputs
  (`fitness/ntmut_fitness_aggregated.csv`), where it is used.

- **Per-tree output layout is byte-identical to v4 for everything except the
  fitness outputs**, apart from the version dir name: a single-tree run writes
  the same files under `{output_root}/{label}_v5/...` that a v4 run wrote under
  `{output_dir}_v4/...`, except that the fitness files
  (`fitness/ntmut_counts/*`, `fitness/ntmut_fitness/*`,
  `fitness/ntmut_fitness_aggregated.csv`, `aamut_fitness/*`) now exclude masked
  mutations, drop the `exclude`/`masked_in_usher` columns, and pool the `"all"`
  cluster by nucleotide identity (see above).

- **Multi-tree CI coverage** (test-infra, no schema impact):
  `config/config_test_data_multitree.yaml` and a step in the Test Pipeline job
  build two trees from one config and assert their `get_clades` checkpoints
  resolve independently — guarding against a regression where a per-tree
  `max_leaf_date` override is ignored.

### v4

**Overhaul UShER masking, mask all founder→reference reversions, and make the codon mask
clade-agnostic** — tracks issues [#99](https://github.com/matsengrp/mat-pcp/issues/99) /
[#100](https://github.com/matsengrp/mat-pcp/pull/100) and
[#103](https://github.com/matsengrp/mat-pcp/issues/103).

**Changes**:

- **New UShER mask source**: `data/usher_masked_sites.yaml` is now the generated,
  Nextstrain-clade-keyed file produced by the standalone
  [`s2-clade-mask`](https://github.com/matsengrp/s2-clade-mask) repository (from the UCSC
  branch-specific mask), replacing the hand-maintained file. Each block carries whole-site
  `sites` masks and an optional per-clade `mutations` list of directional reversions;
  `masking.compute_exclusions` consumes both (site masks and mutation masks set
  `exclude` + `masked_in_usher`).

- **Mask all founder→reference reversions**: any count-grid mutation whose destination base
  equals the reference base at that site (`to_nt == reference_nt`) — i.e. a clade-defining
  base reverting to the reference — is now excluded (`exclude = True`, but **not**
  `masked_in_usher`). `build_count_table` gains the reference FASTA as an input to supply the
  per-site reference base. This is unconditional (no config toggle). The GLM is unaffected
  (these sites are clade-defining, hence non-conserved, and were never in GLM training).

- **Clade-agnostic codon mask**: `get_problematic_sites` now emits a codon site in
  `spike_codon_sites_to_mask.txt` only if it is excluded **across all clades** (every row at
  the `nt_site` has `exclude == True`). This keeps clade-independent masks
  (`site_include_range`, `sites_to_exclude`, ProblematicSites VCF) and drops clade-specific
  ones (UShER per-clade masks and founder→reference reversions), which should not
  clade-specifically mask the target gene.

- **Config consolidation**: `additional_sites_to_ignore` is removed; its sites are merged into
  `sites_to_exclude` (de-duplicated), so all globally-excluded sites follow one downstream
  path (`exclude = True`). Original upstream comments from `jbloomlab/SARS2-mut-fitness` are
  restored on `sites_to_exclude`. The merged sites now show `exclude = True` in
  `counts_all.csv` (previously absent from that flag). GLM training set and the codon mask are
  otherwise unchanged.

### v3

**Unify the count table, make the GLM the single neutral expectation, and simplify masking**
— tracks issue [#101](https://github.com/matsengrp/mat-pcp/issues/101).

**Changes**:

- **Single unified count table** (`build_count_table`, `workflow/scripts/build_count_table.py`):
  - New `counts_all.csv` replaces `expected_vs_actual_mut_counts.csv`. It enumerates the full
    grid of all possible mutations per clade (supplying the zero-`actual_count` rows the GLM
    trains on), joins the aggregated actual counts, and computes `exclude` / `masked_in_usher`
    **once**.
  - Removed rules/scripts: `synonymous_mut_rates`, `expected_mut_counts`,
    `merge_expected_and_actual_counts`, `aggregate_mutations_to_exclude`. Removed intermediates:
    `synonymous_mut_rates/rates_by_clade.csv`, `expected_mut_counts/expected_counts.csv`,
    `mutations_to_exclude.csv`, `expected_vs_actual_mut_counts.csv`.

- **Masking computed once** (`workflow/modules/masking.py`):
  - `count_mutations.py` no longer computes masking; `mutation_counts/{clade}.csv` and
    `mutation_counts/aggregated.csv` **lose the `exclude` column** (raw counts only).
  - `exclude` / `masked_in_usher` are computed in `build_count_table` via
    `masking.compute_exclusions`, from `sites_to_exclude`, the ProblematicSites VCF, the UShER
    yaml `sites`, and `site_include_range`.
  - The `exclude_ref_to_founder_muts` config key is removed. The count grid enumerates only
    mutations *away from the clade founder* (`nt_mutation` = founder_nt + site + to_nt), so
    clade-defining ref→founder mutations (`ref_nt` + site + founder_nt) can never appear in it
    and need no explicit exclusion — including in the GLM training set.
  - The now-unused `clade_founders/{clade}_ref_to_founder_muts.csv` output is removed;
    `clade_founder_fasta` still applies those mutations in-memory to build the founder FASTA.

- **The GLM is the single neutral-expectation estimator**:
  - The synonymous 4-fold `expected_count` is gone. The GLM output column `predicted_count` is
    **renamed `expected_count`** everywhere, and is used as the neutral baseline for **both** the
    pseudocount (`naive_delta_fitness`) and the Bayesian (`delta_fitness`) estimators.
  - Fitness outputs (`fitness/*`, `aamut_fitness/*`) carry `expected_count` in place of
    `predicted_count`; the separate synonymous `expected_count` column no longer exists.
  - The `predicted_counts_by_clade` rule and its output `fitness/predicted_counts_by_clade.csv`
    are renamed to `expected_counts_by_clade` / `fitness/expected_counts_by_clade.csv`.
  - The `synonymous_spectra_min_counts` config key is removed.

- **Behavior:** `naive_delta_fitness` now compares actual counts against the GLM expectation
  (previously the synonymous-rate expectation); `delta_fitness` shifts slightly because
  ref→founder mutations are now excluded from GLM training.

### v2

**Bundled v2 cut** — tracks issues [#75](https://github.com/matsengrp/mat-pcp/issues/75), [#76](https://github.com/matsengrp/mat-pcp/issues/76), and [#80](https://github.com/matsengrp/mat-pcp/issues/80). [#81](https://github.com/matsengrp/mat-pcp/issues/81) ships the first three production runs against this version.

**Changes**:

- **Fused PCP extraction and curation** ([#75](https://github.com/matsengrp/mat-pcp/issues/75)):
  - The `curate_pcps_by_clade` rule, `workflow/scripts/curate_pcps_by_clade.py`, and `workflow/modules/pcp_validation.py` are removed.
  - `extract_pcps_by_clade` now writes the final curated `pcps/{clade}_pcps.csv` directly (no `mat/clades/{clade}_pcps_raw.csv` intermediate) and emits `pcps/{clade}_filtered_branches.csv` alongside it.
  - The identical-parent/child filter is applied inside `pcp_from_tree.py` after the gene-region slice, replacing the separate curation step.
  - Output schema is unchanged: `pcps/{clade}_pcps.csv`, `pcps/all_pcps.csv`, `pcps/{clade}_filtered_branches.csv`, `pcps/filtered_branch_ids.csv`, and `pcps/pcp_curation_report.txt` all keep their v1 columns and content.
  - The unit-test file `tests/test_pcp_validation.py` is removed along with the module it covered.

- **`min_clade_size` guardrail** ([#75](https://github.com/matsengrp/mat-pcp/issues/75)):
  - The Snakefile prints a `WARNING` to stderr at config-load time if any config matching `config-public-*.yaml` sets `min_clade_size > 0`. v1 production configs shipped with `min_clade_size: 10000` and biased the GLM rate estimation; the guardrail makes a recurrence of that mistake noisy.

- **DASM-ready PCP outputs** ([#80](https://github.com/matsengrp/mat-pcp/issues/80)):
  - New `pcps/{clade}_pcps_dasm_ready.csv` and `pcps/all_pcps_dasm_ready.csv` outputs land alongside the existing `*_pcps.csv` files. These absorb the post-processing that `data-central/scripts/filter_viral_pcps.sh` does today (invalid-PCP filtering, masked-only-mutation filtering, 90/10 train/val split).
  - Schema: all columns from `*_pcps.csv` are preserved verbatim; one new column `in_train` (boolean) is appended. `True` means the row is in the train split, `False` means val.
  - Seed convention: per-clade `seed=1729`, aggregate `seed=42` (matches the guardrail in `filter_viral_pcps.sh`).
  - Byte-equivalence to `dasm2 preprocess` is asserted on the test pipeline: the aggregate has md5 `5e538bd088e676b485a11a04745221e9` (sentinel in `validate_test.py`).
  - The filter is a NumPy-only reimplementation of `dasm2.preprocess_command.preprocess`, vendored into `workflow/scripts/curate_pcps_for_dasm.py` so the new rule stays in `envs/default.yml` (no torch/dasm2 dep).
  - Companion files emitted alongside each output: `{clade}_dasm_filter_report.txt` (human-readable per-filter drop counts) and `{clade}_dasm_filter_stats.json` (machine-readable).
  - Once v2 ships, these outputs make `data-central/scripts/filter_viral_pcps.sh` redundant.

- **Gene-sliced clade founder FASTAs** ([#76](https://github.com/matsengrp/mat-pcp/issues/76)):
  - `rule clade_founder_fasta` now also emits a gene-sliced founder FASTA at `clade_founders/{gene_dir}/{clade}.fa` whenever `pcp_gene_region` is configured. For the canonical Spike slice this is `clade_founders/spike/{clade}.fa` (3819 nt, positions 21563-25381, 1-based inclusive, excludes the trailing stop codon).
  - The slice uses the same `seq[start-1:end]` convention as PCP extraction, so the founder slice byte-aligns with the parent/child sequences in `pcps/{clade}_pcps.csv`.
  - Cross-version integrity sentinel: the 19A spike slice has md5 `b948f2a1f20edbce64cfff1dc540c119`; this is asserted in `validate_test.py`.
  - The `pcp_gene_region` block in all production configs and the slicing helper in `workflow/scripts/pcp_from_tree.py` are now documented as 1-based inclusive on both ends.
  - Once v2 ships, this replaces the post-hoc `data-central/scripts/slice_spike_founders.sh` step.

### v1

**Initial release** - Pipeline combining mat-pcp, SARS2-mut-fitness, and SARS2-mut-fitness-v2

**Output Schema**: See [OUTPUT_SCHEMA.md](OUTPUT_SCHEMA.md) for complete schema documentation.

**Key Features**:
- Parent-child pair (PCP) extraction from mutation-annotated trees
- Expected vs actual mutation count computation
- Neutral mutation rate estimation via GLM
- Nucleotide and amino acid fitness estimation
- Per-clade and aggregated outputs

**Main Output Files**:
- `pcps/all_pcps.csv` - Curated parent-child pairs across all clades
- `aamut_fitness/all_aamut_fitness.csv` - Amino acid fitness estimates
- `fitness/ntmut_fitness_aggregated.csv` - Nucleotide fitness estimates
- `rates_table.csv` - GLM-estimated neutral mutation rates
- `validation/pcp_vs_counts_report.txt` - Internal validation report

**Configuration Parameters**: Pipeline behavior is controlled by configuration files (see `config/config.yaml` for full parameter documentation).
