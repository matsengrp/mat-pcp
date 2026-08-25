# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`mat-pcp` is a unified Snakemake pipeline for analyzing SARS-CoV-2 genomic evolution from mutation-annotated trees (MATs). The pipeline combines three previously separate analysis workflows into a single coherent pipeline:

1. **mat-pcp** - Parent-child pair extraction from MATs
2. **SARS2-mut-fitness** - Expected vs actual mutation count computation
3. **SARS2-mut-fitness-v2** - Neutral mutation rate estimation via GLM

### Pipeline Stages

| Stage | Description | Key Rules |
|-------|-------------|-----------|
| 1-2 | Input Processing & Reference | `filter_metadata`, `filter_mat_tree`, `export_newick`, `chronumental`, `chronumental_qc` (tips **and** internal nodes, #139), `chronumental_rtt_plot` (#146), `window_subtree_counts`, `get_ref_fasta`, `ref_coding_sites` |
| 3 | Clade Processing & PCPs | `get_clades` (checkpoint), `mat_clade_subset`, `extract_pcps_by_clade`, `censor_table` (per-clade DASM censor table, #144) |
| 4 | Mutation Counting | `count_mutations_clade`, `aggregate_mutation_counts` |
| 5 | Unified Count Table | `build_count_table` (full mutation grid + actual counts + masking flags, via `modules/masking.py`) |
| 6 | GLM Rate Estimation | `annotate_counts`, `curate_counts`, `master_table`, `get_rates_table`, plus their per-era counterparts `curate_counts_era`, `master_table_era`, `get_rates_table_era` and `era_fit_summary` (#137), `window_expected_counts` (windowed expected counts + `delta_fitness_ratio`; stays on the **pooled** table) |
| 7 | Fitness Estimation | `ntmut_fitness_by_clade`, `aamut_fitness_by_clade`, `aggregate_*_fitness` |
| 8 | Validation | `extract_mutations_from_pcps`, `validate_pcp_vs_counts` |

## Versioning

The pipeline uses simple integer versioning (v1, v2, v3, etc.) for reproducibility and output tracking.

### Version Management

- **Version File**: `VERSION` in the repository root contains the current pipeline version
- **Automatic Versioning**: The Snakefile derives each tree's output directory as `{output_root}/{label}_{version}` (e.g., `results/public_2022-06-01_v11`)
- **Git Tags**: Each version should be tagged in git (e.g., `git tag v1`) when the branch is ready to merge
- **Changelog**: All output schema changes must be documented in `CHANGELOG.md`

### When to Bump Version

Any change that affects downstream consumers requires a version bump:
- Changes to output file formats or schemas (see `OUTPUT_SCHEMA.md`)
- Changes to filtering logic or curation criteria
- Changes to fitness estimation algorithms
- Changes to GLM rate estimation

In general, any change to the pipeline is a breaking change since models trained on one version's outputs may not work with another's.

**Additive-output carve-out.** A change that only *adds* new output files — reading no existing schema and read by no current downstream consumer — may ship without a version bump; the coordinated bump then rides with the next schema-breaking change. (Example: the additive outputs of EPIC [#115](https://github.com/matsengrp/mat-pcp/issues/115) — `mat/node_dates*` (Phase 1), `windows/subtree_counts.csv` + `founder_muts.csv` (Phase 2), `windows/fitness_trajectory.csv` (Phase 3) — each shipped no-bump; per maintainer decision the single coordinated bump for the whole EPIC landed as **v7** with Phase 4 (#119).)

### How to Bump Version

1. Update the version number in the `VERSION` file (e.g., from `v1` to `v2`)
2. Document all changes in `CHANGELOG.md` under a new version heading
3. Update `OUTPUT_SCHEMA.md` if the schema has changed
4. Update any affected documentation (README.md, this file)
5. Commit the changes
6. Create a git tag matching the version: `git tag v2`
7. Push the tag: `git push origin v2`

> **Also audit output-path defaults.** Any script with its own default output directory (currently `workflow/scripts/validate_test.py`) must track `VERSION` the same way the Snakefile does, or it will silently read a stale `results/<label>_v{N-1}` directory after a bump (see #88).

### Development with Versioning

When making changes that will affect outputs:
1. Work on a feature branch as usual
2. Update VERSION, CHANGELOG.md, and OUTPUT_SCHEMA.md as part of your changes
3. When the PR is ready to merge, the new version tag should be created

## Development Commands

### Running the Pipeline

Every run requires exactly one explicit `--configfile` (there is no default
`configfile:` directive; see the Snakefile header for why).

`config/config.yaml` lists **all four production trees**, one of which
(`gisaid-and-public_2023-01-01`) needs restricted GISAID-derived inputs. A
bare `rule all` run builds every tree; to build a single tree, target its
per-tree output path (see "Run specific target" below).

```bash
# Run with pixi (recommended)
pixi run snakemake -j 4 --use-conda --configfile config/config.yaml

# Dry run
pixi run snakemake -n --configfile config/config.yaml

# Run specific target (per-tree output path)
pixi run snakemake -j 4 --use-conda --configfile config/config.yaml \
  results/public_2022-06-01_v11/aamut_fitness/all_aamut_fitness.csv

# Clean all outputs
pixi run snakemake --delete-all-output --configfile config/config.yaml
```

### Testing
```bash
# Unit tests over frozen/synthetic fixtures (no pipeline run needed; CI runs
# all twelve before the pipeline steps)
pixi run python tests/test_glm_reference_equivalence.py
pixi run python tests/test_clade_founder_guards.py
pixi run python tests/test_founder_codon_background.py
pixi run python tests/test_lightswitch_boundary.py
pixi run python tests/test_pcp_compression.py
pixi run python tests/test_pcp_gene_coordinates.py
pixi run python tests/test_era_mapping.py
pixi run python tests/test_era_rate_dispatch.py
pixi run python tests/test_bte_haplotype.py   # half A (env pin) only without bte
pixi run python tests/test_node_date_qc.py
pixi run python tests/test_date_parsing.py
pixi run python tests/test_censor_table.py

# Half B of the bte test needs the Snakemake-built bte env, so it runs after a
# pipeline run; --require-bte makes a missing env a failure rather than a skip.
pixi run python tests/test_bte_haplotype.py --require-bte

# Run on test data (fast, uses max_clades=3)
pixi run snakemake -j 2 --use-conda \
  --configfile config/config_test_data.yaml

# Validate results
pixi run python workflow/scripts/validate_test.py
```

> **`CLAUDE.md`, `CHANGELOG.md` and `workflow/Snakefile` use CRLF line endings.** Edit them in place. A
> Python `read_text()` / `write_text()` round-trip normalizes them to LF, which
> rewrites every line and buries a three-line content edit in a ~2,200-line
> whitespace diff. Nothing errors and the rendered text is identical, so it is only
> visible in `git diff --stat`. Check with `file CLAUDE.md` after scripted edits.

```bash

# Multi-tree smoke test: two trees in one config that resolve their
# get_clades checkpoints independently (the earlier max_leaf_date yields
# fewer clades). CI builds only each tree's clade_counts.csv and asserts
# the counts differ.
pixi run snakemake -j 2 --use-conda \
  --configfile config/config_test_data_multitree.yaml
```

## Architecture

### Directory Structure
```
mat-pcp/
├── config/
│   ├── config.yaml                  # Default config (all public/GISAID trees)
│   ├── config_test_data.yaml        # Test data configuration
│   └── config_test_data_multitree.yaml # Multi-tree smoke-test configuration
├── workflow/
│   ├── Snakefile            # Main workflow definition
│   ├── envs/                # Conda environments
│   │   ├── default.yml      # Python/pandas/biopython
│   │   ├── usher.yml        # UShER tools (matUtils)
│   │   └── bte.yml          # BTE library for tree traversal
│   ├── modules/             # Python modules
│   │   ├── glm.py           # General Linear Model for rate estimation
│   │   ├── rates.py         # Rate calculations
│   │   ├── probfit.py       # Probabilistic fitness estimation
│   │   ├── aamutfit.py      # Amino acid fitness aggregation
│   │   └── load.py          # Data loading utilities
│   └── scripts/             # Snakemake scripts
├── data/                    # Static data files
│   ├── usher_masked_sites.yaml
│   ├── lan_2022/           # RNA structure data
│   └── test_data/          # Test datasets and expected results
├── results/                 # Active pipeline outputs (gitignored)
│   ├── test_data/          # Test runs
│   └── {dataset}/          # Production runs
└── archives/                # Archived results (gitignored)
    └── README.md           # Archive documentation (tracked)
```

### Key Technologies
- **Snakemake**: Workflow management with checkpoints for dynamic clade processing
- **UShER/matUtils**: MAT tree manipulation and filtering
- **BTE**: Efficient tree traversal for PCP extraction
- **Pandas/NumPy**: Data processing and analysis
- **BioPython**: Sequence handling
- **SciPy/Statsmodels**: GLM rate estimation (in conda envs)

### Important Patterns

1. **Checkpoint pattern**: The `get_clades` rule is a checkpoint that determines the list of clades to process at runtime. All downstream per-clade rules depend on this checkpoint.

2. **Three conda environments**:
   - `default.yml` - Python, pandas, numpy, biopython
   - `usher.yml` - matUtils for tree manipulation
   - `bte.yml` - BTE library for tree traversal

3. **max_clades parameter**: Set `max_clades` in config to limit the number of clades processed. Useful for testing (CI uses `max_clades=3` for ~4 minute runs).

4. **Validation stage**: Cross-validates two independent mutation counting approaches (PCP-derived vs matUtils-derived) to ensure pipeline consistency.

5. **The GLM's synonymous-only filter lives in the loader, not in curation.**
   `glm/curated_counts.csv` deliberately holds *all* mutation classes (~22%
   synonymous on the v7/v8 test tree; the share is tree-specific) because
   downstream rules need the full grid.
   The neutral-rate restriction is applied when that file is *read*:
   `master_table.py` and `expected_counts_by_clade.py` both go through
   `modules/load.py::load_curated_muts`, an alias for `load_synonymous_muts`
   (`synonymous == True`), mirroring `SARS2-mut-fitness-v2`'s
   `master_tables.py.ipynb`. Inspecting `curated_counts.csv` therefore makes it
   *look* as though the GLM is fit on everything; check
   `log/master_table.log`'s "Loaded N curated mutations for training" instead.
   Replacing either loader call with a bare `pd.read_csv` would silently train
   the neutral model on mutations under purifying selection. Guarded by
   `tests/test_glm_reference_equivalence.py` (unit, frozen fixtures) and
   `validate_test.py::validate_glm_synonymous_training` (end-to-end). See
   [#126](https://github.com/matsengrp/mat-pcp/issues/126).

6. **An unresolvable clade founder is a hard failure, not a fallback.**
   `clade_founder_fasta.py` builds each clade's founder by applying that clade's
   defining mutations to the Wuhan reference. If the clade resolves in none of its
   three sources (Neher JSON → Roemer via `roemer_nextstrain_to_pango` → Roemer
   keyed directly), it used to apply *zero* mutations and write the unmodified
   reference — exit 0, logged as ordinary output. Nothing downstream caught it
   because the wrong founder propagated *consistently*: every consumer compared
   against the same wrong sequence, so no internal cross-check disagreed. It
   surfaced only when an outside consumer noticed implausible cross-clade
   correlations (24H/24I/25A/25B/25C in v8).
   Two guards now cover the two routes to that output — a missing *source*, and a
   resolved source that yields no *effective* mutations (empty substitution list,
   or every mutation dropped by the reference-base mismatch check). Both raise.
   The second check is deliberately **genome-wide, not per-gene**: at gene level
   reference-identity is legitimate (19B's founder equals 19A's for S, since
   C8782T and T28144C both fall outside spike), so a per-gene assertion would need
   an allowlist growing with every clade/gene pair — the same shape as the bug it
   guards. Guarded by `tests/test_clade_founder_guards.py`. See
   [#129](https://github.com/matsengrp/mat-pcp/issues/129).

7. **A clade's counts include only mutations off its founder's *codon*, not just
   its founder's nucleotide.** matUtils reports each branch's mutations relative
   to its parent, so the same nucleotide substitution also appears on codons that
   had already mutated — a different amino-acid change sharing one nucleotide
   event. `build_count_table.py` joined observed counts onto the founder grid on
   `(clade, nt_mutation)` alone, dropping `codon_change` (the only column that
   distinguishes the two) and using an unvalidated `how="left"` merge. Off-founder
   observations were therefore relabelled onto the founder codon and kept as
   *extra rows for the same mutation*: 6,415 duplicated keys / 13,146 rows in v8,
   every clade affected, concentrated in the highest-count mutations. Both
   estimators were corrupted by *different* mechanisms — `delta_fitness_bayesian`
   inverse-variance-averages the phantom row in, `delta_fitness_ratio` sums its
   expected count — so they never diverged enough to trip a cross-check, and the
   reported amino-acid uncertainty got *smaller*.
   The filter and the one-row-per-mutation invariant are ports of upstream's
   `merge_expected_and_actual_counts.py.ipynb`, lost when three rules were fused.
   Codon strings are `;`-joined across genes at overlap sites and compared
   whole, which requires the same gene order on both sides — true, and verified
   against the out-of-frame N/ORF9b overlap. A *structural* disagreement (one
   side coding, or a different component count) raises rather than reading as
   off-founder, because that would silently zero a legitimate count and the
   row-count invariant cannot see it. `aamutfit.aa_fitness` also rejects a
   repeated `nt_mutation` within one amino-acid group as independent defence.
   Guarded by `tests/test_founder_codon_background.py` (unit; the test dataset
   has no off-founder backgrounds, so it is the *only* coverage) and
   `validate_test.py::validate_founder_codon_background` (dataset invariant).
   See [#130](https://github.com/matsengrp/mat-pcp/issues/130).

8. **The light-switch boundary site belongs to the *before* bucket, and that
   convention has exactly one definition.** `modules/glm.py` trains with
   `global_context = (nt_site > boundary)`, so `site == boundary` is fitted into
   the before-regime, and `modules/rates.py`'s sentinel (`nt_site = 0 if
   before_boundary else s_max`) makes `nt_site_before_boundary=True` mean exactly
   that bucket. Through v8 both *label producers* used `site < boundary` — and
   `not (site < b)` is `site >= b`, not `site > b` — so the boundary site was
   fitted into one bucket and **served the other's rate**. The convention now
   lives in `modules/lightswitch.py::before_boundary` (`site <= boundary`), called
   by both producers: `scripts/annotate_counts.py` (clade path) and
   `modules/expected_counts.py` (windowed trajectory path). Fixing only one would
   leave the clade and windowed fitness paths disagreeing at the boundary site,
   which is the comparability the `evol_time` normalization exists to preserve.
   Two things make this class of bug invisible, which is why it needs a unit test
   rather than a pipeline assertion: `tests/test_glm_reference_equivalence.py`'s
   fixture is a frozen CSV carrying a *pre-computed* `nt_site_before_boundary`
   column, so it exercises the rate lookup over stored labels and not the label
   computation; and every shipped config sets `AT`/`CG`/`GC` to the same 21562, so
   the separate bug where `CG` and `GC` both read the `AT` key (honored per-type
   in training, ignored at serving) is unobservable from any config — that same
   equivalence test's config check *pins* those values, hiding it rather than
   catching it. mat-pcp therefore diverges from `SARS2-mut-fitness-v2` at site
   21562 deliberately; upstream carries the same inconsistency between its
   `modules/glm.py:59` (`>`) and `notebook/counts_by_clade.py.ipynb` (`<`). Do not
   resolve it by re-baselining the #127 fixture. Guarded by
   `tests/test_lightswitch_boundary.py` (unit, and the *only* coverage: the test
   dataset has just 3 rows at a boundary site — AT @ 21562, one per clade — all
   with `actual_count` 0, so a test-data run has nothing asserting which rate
   they are served).
   **Both producers must be executed, not just named.** That test's
   two-producer check originally called `lightswitch.before_boundary` directly
   for the windowed side, so it compared `annotate_counts.py` against the shared
   module and never ran `expected_counts.py` at all — reverting that call site to
   `site < boundary` left the whole suite green. It now drives the real
   `ExpectedCounts._add_context_features` (built with `object.__new__` and the two
   attributes the method reads, since `__init__` demands a coding-sites table
   irrelevant to the label). A test named for two producers that exercises one is
   the same vacuity as pattern 10's `gene == {gene}` filter. See
   [#125](https://github.com/matsengrp/mat-pcp/issues/125).

9. **PCP outputs are gzipped, and the DASM md5 is a *payload* hash.** The four
   PCP families (`{clade}_pcps`, `all_pcps`, `{clade}_pcps_dasm_ready`,
   `all_pcps_dasm_ready`) are `.csv.gz` as of v9; `{clade}_filtered_branches.csv`
   and `filtered_branch_ids.csv` deliberately stay plain. On the v8 production
   tree `pcps/` was **504 GB apparent (470 GiB), 93% of the whole 542 GB output
   tree** — quote apparent or allocated consistently, never a mix (allocated is
   403/415 GiB, i.e. 97%) — and
   measured gzip level 9 gives ~125× on ORF1a / ~110× on S. The policy lives in
   `modules/pcpio.py::PCP_CSV_COMPRESSION` rather than at each `to_csv`, because a
   writer that omitted `mtime=0` would break reproducibility for that one file
   with nothing to attribute it. Level 9 is *not* pinned explicitly: it is
   pandas' default, and it measured both smaller and faster than level 6 (~125×
   at ~97 MiB/s vs ~17× at 29 MiB/s — at level 9 gzip's full match search finds
   the long parent/child and adjacent-row repeats, so it emits far less).
   Three consequences fail silently and so are unit-tested rather than left to
   the pipeline run: (a) `validate_test.py`'s `DASM_READY_ALL_SENTINEL_MD5` is a
   *cross-implementation* check against `dasm2 preprocess`, which writes plain
   CSV, so `test_summary.py` hashes the **decompressed** stream and the sentinel
   value is unchanged by compression — re-baselining it to a compressed-file md5
   would silently retire the check; (b) pandas' gzip output embeds `mtime=now`
   *and* the output filename (the FNAME flag), so `mtime=0` is required for a
   re-run to be byte-identical, which matters because `data-central` manifests pin
   `md5`/`xxh3` per file; (c) `aggregate_clade_pcps.py` derives the `clade` column
   from each input's filename, and its old `replace("_pcps.csv", "")` would have
   yielded `"19A.gz"` — every row of every aggregate mislabelled identically, with
   the row count unchanged and no cross-check able to see it. Guarded by
   `tests/test_pcp_compression.py`. See
   [#134](https://github.com/matsengrp/mat-pcp/issues/134).

10. **nsp12 begins *before* the ribosomal frameshift, and a PCP gene is not
    necessarily an `aamut_fitness` gene.** Two facts, each of which makes an
    incorrect implementation look right.

    *ORF1b is nsp12–nsp16, not the post-slippage product.* `orf1ab_to_nsps` puts
    nsp12 (RdRp) at ORF1ab codons 4393–5324, but the −1 frameshift is at codon
    4402 — so **nsp12's first 9 residues are translated in ORF1a's frame**.
    Defining ORF1b as the post-slippage block alone (`(13468, 21552)`, one block,
    2695 codons) is the obvious reading and silently clips 9 residues off the
    N-terminus of the protein the work is chiefly for. ORF1b is therefore
    `[[13442, 13468], [13468, 21552]]` — 8112 nt / 2704 codons, sharing nt 13468
    exactly as ORF1ab does. Note the summed length (8112) is a whole number of
    codons while the *span* (8111) is not; that shared nucleotide is why
    `validate_cds_blocks` sums blocks rather than measuring the span. It is added
    via `pcp_gene_cds_blocks` rather than the GTF because `edit_ref_gtf` feeds
    `ref_coding_sites` → `counts_all.csv` → GLM → fitness: a GTF entry would make
    every site in an 8 kb stretch multi-gene, and absent a matching
    `gene_overlaps.retain` entry, `aamutfit.get_coding`'s `exclude_overlap`
    default would silently drop all of nsp12–nsp16 from `aamut_fitness` — the
    #129/#130 silent-drop shape again.

    *The codon correspondence is derived, not declared.* `aamutfit.get_coding`
    drops ORF1a (`aamutfit.py:80`), so `validate_pcp_vs_aamut.py`'s
    `gene == {gene}` filter matched **zero rows** for ORF1a: through v8 that
    validation reported no discrepancies because it had nothing to compare.
    `modules/genemap.py` now maps codons by **genomic position** — two codons
    correspond iff they occupy the same three reference nucleotides — which
    matters because ORF1a's GTF block ends at 13480 while ORF1ab's first block
    ends at 13468, so **ORF1a codons 4402–4405 have no ORF1ab counterpart at
    all.** Under a bare `codon_offset: 0` those four would map onto ORF1ab codons
    4402–4405, which are *ORF1b's* codons 10–13, and surface as four spurious
    discrepancies indistinguishable from real ones. The configured
    `codon_offset` is kept as a cross-check and disagreement raises.
    `require_nonempty_fitness_slice` makes an empty slice fatal — applied to the
    *per-clade* table as well as the aggregate, since the aggregate being
    non-empty does not imply the per-clade slice is and per-clade carries almost
    all the rows. When both sides are empty nothing raises, so the report says
    **VACUOUS** rather than "SUCCESS: All amino acid mutation counts match
    exactly!", which is true of an empty comparison and is the wording that
    concealed the v8 gap. Every report now also prints the compared row counts on
    both sides, so "agreed" is distinguishable from "never compared".

    *A PCP gene's codon mask has its own two coordinate systems.* A gene reaching
    `pcp_genes` only through `pcp_gene_cds_blocks` never appears in
    `counts_all.csv`'s GTF-derived `gene` column, so `get_problematic_sites.py`
    matched zero rows and wrote an **empty** mask at exit 0 — ORF1b's DASM-ready
    PCPs kept every problematic codon site that every other gene had removed. Its
    codons are now derived from `nt_site` via its own CDS blocks, and which
    derivation applies is read from the *config* (does the gene have an override),
    never from whether the count table happens to mention the gene: branching on
    the data lets a gene defined in both places silently ignore its override and
    bypass `resolve_cds_blocks`' add-only guard.
    Two further traps here. `gene` and `codon_site` are **positionally paired**
    `;`-joined cells (`"N;ORF9b"` with `"4;1"`), so they must be zipped and
    filtered, never flattened — flattening files ORF9b's codon 1 as a codon of N,
    and since ORF9b's numbers (1–98) are all valid N numbers no range check can
    see it. And `validate_cds_blocks` must check that a shared nucleotide falls on
    a **codon boundary**: `[[1, 11], [11, 20]]` shares exactly one nt and totals
    21, so it passes the overlap and divisibility checks while reading nt 11 twice
    *inside one codon*. At most one shared junction is allowed, since
    `[[1, 3], [3, 5], [5, 7]]` is in frame at both and translates cleanly.

    Guarded by `tests/test_pcp_gene_coordinates.py` (unit; every failure mode
    mutation-tested, and it **execs the real `get_problematic_sites.py`** against a
    stub `snakemake` — an earlier version re-implemented the script's
    nt_site→codon inversion and asserted against its own copy, so three separate
    breaks of the shipped script left it green) and
    `validate_test.py::validate_orf1b_founder_slice`, which
    checks every clade's ORF1b founder byte-equals `ORF1ab[13176:21288]` — two
    coordinate sets, one slicer, so a frame error moves one and not the other.
    The test config carries ORF1b so the override runs end-to-end in CI; **ORF1a
    is not in it**, so no pipeline run asserts the vacuity fix itself.
    See [#133](https://github.com/matsengrp/mat-pcp/issues/133).

11. **The neutral model is fitted once per era, and the pooled fit is the
    degenerate one-era case rather than a second implementation.**
    `SARS2-mut-fitness-v2` fits two rate models — pre-Omicron and Omicron+ — and
    dispatches per clade; the `pre_omicron_or_omicron` covariate did not survive
    mat-pcp's port, while the light-switch covariate from the same upstream
    notebook did. v10 restores it. Three things about the shape are load-bearing.

    *`add_expected_count` is `add_expected_count_by_era` with one era.* That makes
    "the split did not move the pre-existing fit" true by construction instead of
    by measurement, and it is why a single-era tree reproduces the v9 outputs
    **byte-for-byte** — verified on the test data across `curated_counts.csv`,
    `master_table.csv`, `rates_table.csv`, `expected_counts_by_clade.csv`,
    `ntmut_fitness_aggregated.csv` and `all_aamut_fitness.csv`. Re-implementing
    the pooled path separately would retire that guarantee silently, since both
    paths would still produce plausible tables.

    *Conserved sites are computed over every clade **before** the era partition*
    (`curate_counts.py`), matching upstream, which computes `conserved_sites` from
    all clades and only then queries by era. Computing conservation within era
    gives each era a different training site set — so the two era fits and the
    pooled fit would no longer be comparable, and the degeneracy guarantee above
    would break — while still writing three well-formed files. The era filter
    applies to *rows*, never to the site-set derivation.

    *An empty era is skipped at DAG-construction time, and the skip is recorded.*
    `populate_rates` raises on an empty partition, and single-era trees are
    normal, not hypothetical: `config_test_data.yaml`'s `max_clades: 3` gives
    19A/19B/20A, all pre-Omicron. `Snakefile: eras_present` resolves the
    `get_clades` checkpoint and requests only non-empty eras, rather than having
    the rule emit a placeholder — a header-only file or a `touch()` would ship a
    rate table that looks fitted and contains nothing. But a skipped era leaves
    *no file*, which is indistinguishable from an era dropped by a wiring
    mistake, so `glm/era_fit_summary.csv` names **every** era with `fitted`
    true/false. The two failure modes stay distinct: *no clades in the era* →
    skip and record; *clades present but routed to an unfitted era* → raise. One
    `try/except` around both turns a broken fit into a silently absent table.

    **Three tables, not one with an `era` column.** DASM2's
    `TabularRatesNeutralModel` assigns `wt_to_mut_rates[wt][mut]` in a row loop,
    so an unrecognized extra column is read *last-row-wins, silently*; and
    `rates.residual_variance` zips `groupby(...).apply(...)` output
    **positionally** against `self.rates.loc[idx]`, guarded only by a
    `len(...) != len(...)` early `return` — correct today only because the grid's
    construction order happens to equal sorted order. The pooled
    `master_table.csv` is retained with its v9 meaning as the **control arm** of
    the per-era comparison, not as backwards compatibility.
    `window_expected_counts` deliberately still reads it: the windowed path keys
    on `founder_id` with no clade to dispatch on, and a per-window date label is
    wrong exactly across the transition.

    **No shipped config is two-era, so CI cannot check the dispatch.** Every
    config resolves to one era, on which the whole path reduces to v9's pooled
    fit — a pipeline run cannot tell whether the `curated_{era}` inputs are read,
    matched to the right era, or dispatched through `era_for_clade` at all.
    `tests/test_era_rate_dispatch.py` is therefore the only two-era coverage, and
    it **execs the real `curate_counts.py` and the real
    `expected_counts_by_clade.py`** against two constructed eras whose G→T rates
    differ by a known factor: one table serving everyone yields a ratio of 1.0
    and a swapped pair inverts it, so both are distinguished from the right
    answer. Checking `add_expected_count_by_era` alone and calling that two-era
    coverage is pattern 8's vacuity — a test named for the dispatch that never
    runs the dispatching code. Mutation-tested against four separate breaks of
    the shipped scripts. `validate_test.py::validate_era_rate_tables` is the
    end-to-end counterpart, reconciling the summary against the output tree in
    both directions.
    See [#137](https://github.com/matsengrp/mat-pcp/issues/137).

12. **UShER masking censors, and a censored mutation is not relative to the
    clade founder.** Two facts about `pcps/{gene}/censored_mutations.csv`, each
    of which makes a plausible guard wrong.

    *Masked sites are censored, not ambiguous.* UShER does not write `N`; it
    **deletes** the mutations from the MAT, so the node inherits its parent's
    state. The state is known, the *event* was destroyed. That rules out
    representing the mask in the sequence — an `N` tokenizes its whole codon to
    `AA_AMBIG_IDX` and would throw away a state we have — so the correction
    belongs in the observation model, and mat-pcp's job is to emit coordinates.
    `{gene}_codon_sites_to_mask.txt` cannot serve: it never reaches the loss
    (only the PCP drop filter and the dashboards read it), it is codon-granular,
    and `get_problematic_sites.py`'s `site_all_excluded` keeps by construction
    only the reasons holding in *every* clade, so the entire clade-specific half
    of the yaml is invisible to training. The two files carry disjoint *reasons*
    (two of `compute_exclusions`' five here, three there) but **not** disjoint
    sites — a site can be both clade-masked in the yaml and in the
    ProblematicSites VCF, 32 such `(clade, nt_site)` pairs on the v11
    gisaid-and-public tree — so a training run applies both and neither file may
    assert disjointness against the other.

    *matUtils reports each branch relative to its parent, not to the founder.*
    A directional mask's `from_nt` therefore need **not** equal the clade
    founder's base: the mask can legitimately sit on a derived background,
    including a reversion back toward the founder, which is exactly the artifact
    class UShER masking targets. Of the 1,126 directional masks carried by the
    32 clades that have any on the v11 gisaid-and-public tree, 1,083 have
    `founder == from_nt`, 36 have `founder == to_nt` and 7 have neither (e.g.
    25A `A22896T`, founder `G`). An
    equality assertion — the obvious guard, and the one #144 originally
    specified — fails on real data. The sensitive property is the *ratio*: a
    `from_nt`/`to_nt` swap inverts the distribution wholesale, so that is what
    raises, with the three-way counts logged as the diagnostic. This is the same
    observation as pattern 7, one level up: a clade's observed mutations are not
    all off its founder.

    Two further traps. The genomic→`nt_index` map is **one-to-many** at a −1
    frameshift (ORF1b reads genomic 13468 at indices 26 *and* 27, in different
    codons), so a mask there emits two rows and a flat `nt_site - cds_start`
    offset is wrong for every position past the junction — pattern 10's trap
    again. And the clade **roster sidecar is part of the deliverable**: a clade
    with nothing censored contributes zero rows, so in the CSV alone it is
    indistinguishable from a mistyped clade name, and the consumer cannot
    implement "unknown clade raises, empty clade passes" without it.

    `sequences.slice_cds` deliberately keeps its per-block string slice rather
    than gathering over `cds_nt_positions`, even though sharing the code would
    be tidier: it runs twice per PCP over millions of PCPs and the gather
    measured ~900× slower on ORF1ab (1.3 µs vs 1166 µs), which caching the
    position list does not rescue. The anti-drift guarantee is a test, not
    shared code.

    Guarded by `tests/test_censor_table.py`, which is the **only** non-vacuous
    coverage: the CI clades are 19A/19B/20A and the yaml's earliest clade is
    20B, so the table is empty for every gene in every CI run, and no yaml entry
    falls in ORF1b's leading block 13442–13468, so the junction case is
    unreachable from any pipeline run. It execs the real script — including
    deliberately corrupted copies, so a guard that never fires is caught —
    because a test that re-implements the inversion and asserts against its own
    copy is pattern 10's vacuity. `validate_test.py::validate_censor_table` is
    the end-to-end counterpart and checks the roster, which stays non-vacuous
    even when the table is empty.
    See [#144](https://github.com/matsengrp/mat-pcp/issues/144).

13. **An incremental slope and an aggregate rate are different statistics, and only
    one of them can see duration inflation.** chronumental's durations are free
    parameters under `dist.Uniform(0, 365 * 10000)` with a `Delta` guide -- nothing
    penalises a long branch -- while every day-resolved leaf date is asserted to
    **sd 0.3 days** (`--variance_dates` defaults to 0.3 and multiplies a per-leaf
    error of 1 day; the flag name is a misnomer, numpyro's second positional is a
    scale). So the fit can buy leaf-date accuracy by stretching internal branches,
    and on both v11 production trees it did: the fitted rate landed **4.3-4.4x below
    chronumental's own pre-fit root-to-tip estimate**, `length_cor` fell
    monotonically from 0.995 to ~0.38 and was still falling at the last step, and
    57% of branches carry zero mutations while holding **46% of the tree's total
    elapsed time**. mat-pcp's learned clock does not protect against this; the rate
    simply follows the durations down (38.2 -> 8.87).

    *Three of the four gates cannot see it, for structural reasons rather than
    threshold reasons.* Tip residuals are the fitting target, so passing them is
    close to circular (`median|res| = 1.0 d` against a 30 d gate); monotonicity holds
    by construction, since durations are drawn non-negative and cumsummed; and
    coverage is unrelated. The one check that does measure branch Delta-t against
    mutation count, check 5, reported only the **incremental** slope between adjacent
    single-count bins -- deliberately intercept-free, because `median(Delta-t)/n`
    carries the zero-mutation intercept and decays 9x across the bins for reasons
    that are arithmetic rather than biological. But the inflation lives *entirely* in
    that intercept, so differencing adjacent medians cancels exactly the term that
    carries it: the slope read 11 d/mutation against a ~15.3 d nominal and looked
    healthy. Check 5 now also reports `Sigma mutations / Sigma branch-years`, which is
    the statistic the slope cannot see. Both are kept -- the slope is still the right
    way to ask whether the clock is linear, and the aggregate is the right way to ask
    whether the total elapsed time is real.

    Three traps in the aggregate itself. It must include the **zero-mutation
    branches** in both sums, since dropping them reproduces the slope's blindness
    (the unit test's fixture is built so the two statistics *must* disagree, and a
    version computing either one the other's way fails it). It is computed over the
    **per-branch arrays**, not reconstructed from `node_dates_branch_dt_summary.csv`,
    whose ranged bins (`6-10`, `11-20`, `21-50`, `51+`) need invented representative
    counts that move the answer ~+/-7%. And the zero-mutation share is reported as a
    share of branch-**time**, not of branches: most branches carrying no mutations is
    normal under dense sampling, whereas those branches holding half the elapsed time
    is duration with no mutational support.

    *The reference matters more than the statistic.* Comparing the aggregate to *the
    rate chronumental reported* is near-identity under a learned clock -- measured
    **1.027** -- so a "ratio near 1 means fine" rule passes on a 4.3x-inflated tree.
    The reference has to be independent: chronumental's own root-to-tip regression
    line from `log/chronumental.stdout`, or `DAYS_PER_MUTATION_NOMINAL`. That log is
    a Snakemake `log:` entry rather than a declared output, so it is passed as a
    **param** (an input would drag the multi-hour dating rule back into the DAG when
    a log is rotated away) and its absence is **reported with the path it looked
    at**, never silently dropped -- a detail line missing the comparison reads
    exactly like one where the comparison passed.

    `scripts/chronumental_sweep.py` measures prior settings by executing the **real**
    `chronumental_qc.py` per arm, including the tip-residual gate verdict, rather
    than re-implementing the metrics: a sweep whose own "median absolute tip
    residual" disagreed with the gate the pipeline enforces would be measuring the
    wrong quantity with total confidence. The QC script writes its report *before*
    raising, so a gate-breaching arm still yields a complete row -- breaching the
    gate is a result, not an error. Note also that adopting a prior setting is **not**
    a free flag change: it moves `mat/node_dates.tsv` -> `windows/subtree_counts.csv`
    -> `windows/fitness_trajectory.csv`, so it needs a `VERSION` bump.
    *The reference is now a picture, not just a number.* Because everything above
    turns on chronumental's pre-fit root-to-tip slope, that slope had better describe
    a line -- and chronumental computes the cloud, prints the slope and discards the
    cloud, so nothing had ever looked. `chronumental_qc` therefore also emits
    `mat/node_dates_rtt.csv` and a `root_to_tip_regression` REPORT row, and
    `chronumental_rtt_plot` renders `mat/node_dates_rtt_scatter.png`. Four things
    about the shape are load-bearing. The histogram is **exact, not a downsample**:
    both axes are integers (a day-resolution date, a mutation count), so it *is* the
    joint distribution and the figure's regression is bit-identical to the 17 M-point
    cloud's -- which is asserted against an expanded point array, because a
    key-encoding collision would move the slope only slightly and no eye could
    falsify it. The root's own mutation list is **excluded** from divergence: on a MAT
    it is the root's difference from the reference, so counting it offsets every tip
    by one constant -- the cloud shifts up as a block and the slope is unchanged,
    which is precisely the error a scatter cannot show you. The two lines run over
    **different tip sets** (day-resolution tips here, every dated tip in
    chronumental's own regression), so the tripwire is a wide 25% ratio rather than an
    equality; a tight tolerance would cry wolf on a correct run. And plotting gets its
    **own conda env**: `chronumental_qc` runs on `bte.yml` and widening a shared env
    changes its hash, handing every rule using it a `software environment has changed`
    trigger -- a full rebuild of a released tree for a figure.

    See [#145](https://github.com/matsengrp/mat-pcp/issues/145) and
    [#146](https://github.com/matsengrp/mat-pcp/issues/146).

### Configuration Sections

The configuration (`config/config.yaml`) is organized into:

1. **Input/Output**: `output_root` (optional, default `results`) and `trees` (a
   mapping of tree **label** → per-tree config). Each tree requires `mat` and
   `metadata`; `max_leaf_date` is the only per-tree-overridable key. Output dirs
   are derived as `{output_root}/{label}_{version}` — there is no `output_dir`
   key. A single run builds every tree in `trees`.
2. **Reference Genome**: `ref_fasta`, `ref_gtf`, `add_to_ref_gtf`
3. **Metadata QC**: `qc_filter_options` (date filtering via `min_leaf_date` /
   `max_leaf_date`, host validation, genome length). Global except that a tree
   may override `max_leaf_date`. `max_leaf_date` is optional; when unset it
   defaults to `None`, meaning no upper date cutoff.
4. **Tree Processing**: `optimize_tree`, `mat_subset_kwargs`
4b. **Node Dating (chronumental)**: `chronumental_kwargs` (passed verbatim to
   `chronumental`, e.g. `--steps`), `chronumental_max_median_abs_tip_residual_days`
   (the `chronumental_qc` gate — the run fails if the median absolute tip-date
   residual exceeds it). The QC rule runs on every tree.
4c. **Window Fitness (sliding windows)** (EPIC #115 Phases 2-3): `window_size`,
   `window_step`, `window_start_time`, `window_end_time`, `window_min_subtree_size`.
   The `window_subtree_counts` rule (Phase 2) buckets branch mutations into
   overlapping time windows (anchored at the end date, tiled backwards) and emits
   the `windows/` actual-count outputs. Its branch QC reuses the Mutation-Counting
   keys (`max_nt_mutations` / `max_reversions_to_ref` /
   `max_reversions_to_clade_founder`), so the clade-based and window-based paths
   filter branches identically. The `window_expected_counts` rule (Phase 3, #118,
   `modules/expected_counts.py` — a pure-pandas port of s2trajectory's
   `PossibleMutations`, no bte) then computes per-window expected counts from the
   GLM neutral model (`glm/master_table.csv`) and emits
   `windows/fitness_trajectory.csv` with
   `delta_fitness_ratio = log((actual + p) / (expected + p))`. Expected counts use
   the GLM only (no synonymous-average fallback) and **mat-pcp's own clade-path
   `evol_time` normalization applied per window** (`T = observed_syn / rate_syn`;
   `expected = T × rate`), NOT s2trajectory's `rate × subtree_mut_count` — so total
   expected synonymous == total observed synonymous per window by construction, and
   windowed fitness is directly comparable to the clade-based `aamut_fitness`.
5. **Clade Configuration**: `clade_column`, `min_clade_size`, `clades_to_exclude`, `max_clades`
6. **PCP Extraction**: `pcp_genes` (PCP branch QC reuses the shared
   Mutation-Counting keys `max_nt_mutations` / `max_reversions_to_ref` /
   `max_reversions_to_clade_founder`), plus two optional maps that exist only
   because a PCP gene need not be a GTF feature *or* an `aamut_fitness` gene
   (#133):
   - `pcp_gene_cds_blocks` — explicit 1-based-inclusive CDS blocks for genes the
     reference GTF does not carry (currently `ORF1b` only). Consulted by
     `gtf.resolve_cds_blocks` ahead of the GTF and deliberately **not** fed to
     `edit_ref_gtf` / `ref_coding_sites`, so fitness outputs are untouched.
     Add-only: overriding a gene the GTF defines raises. Every key must also be
     in `pcp_genes`.
   - `pcp_gene_to_fitness_gene` — PCP gene → `{gene, codon_offset}` naming the
     `aamut_fitness` gene it is validated against. Needed for `ORF1a` (absent
     from `aamut_fitness`) and `ORF1b` (a slice of `ORF1ab`). The offset is
     re-derived from the CDS blocks and cross-checked against the configured
     value, so a stale number fails loudly.
7. **Clade Founders**: `clade_founders` (Neher and Roemer JSON sources),
   `roemer_nextstrain_to_pango` (Nextstrain clade → Pango lineage, for clades not
   in the Neher JSON), `allow_reference_identical_founders` (optional; clades
   whose founder may legitimately equal the reference genome — the
   `reference_clade` is exempt automatically, so this is normally empty).
   **Every config carries its own copy of `roemer_nextstrain_to_pango`**, so
   adding a clade to one config does not add it to the others; that is how #129
   happened. `clade_founder_fasta.py` now raises rather than silently writing the
   reference as an unresolvable clade's founder — see the note under "Important
   Patterns" below.
8. **Mutation Counting**: `max_nt_mutations`, `max_reversions_*`, `site_include_range`
9. **GLM**: `rna_structure_data`, `glm_lightswitch_boundaries`
11. **Fitness**: `fitness_pseudocount`, `fitness_N_f`, `fitness_N_n_exp`
12. **Gene Configuration**: `orf1ab_to_nsps`, `gene_overlaps`, `genes`
13. **Site Masking**: `sites_to_exclude`, `site_include_range`
14. **Threads**: Per-rule thread allocation

### Output Organization

Results are stored per tree in `{output_root}/{label}_{version}/`:
- `mat/`: Filtered trees, clade subtrees
- `windows/`: Sliding-window subtree counts + founders (Phase 2) + per-window expected counts / `delta_fitness_ratio` trajectory (Phase 3)
- `pcps/`: Curated parent-child pairs, plus each gene's `censored_mutations.csv` DASM training censor table and its `.clades.txt` roster (#144)
- `ref/`: Reference genome, coding sites
- `mutation_counts/`: Per-clade and aggregated counts
- `clade_founders/`: Founder sequences and mutations
- `glm/`: Annotated counts, curated counts, master table
- `fitness/`: Predicted counts, nucleotide fitness
- `aamut_fitness/`: Amino acid fitness
- `validation/`: PCP vs counts validation report
- `log/`: Processing logs

### Data Files
- `data/usher_masked_sites.yaml`: UShER masking information for problematic sites by clade
- `data/lan_2022/41467_2022_28603_MOESM11_ESM.txt`: RNA secondary structure data

## Pre-PR Quality Checklist

Before any pull request:

1. **Issue Compliance**: Verify completion of all requirements from the linked issue
2. **Versioning**: If changes affect outputs or downstream consumers:
   - Bump the version number in the `VERSION` file
   - Document changes in `CHANGELOG.md` under the new version
   - Update `OUTPUT_SCHEMA.md` if schema changed
   - Plan to tag the version in git when PR is ready to merge
3. **Documentation Sync**: If changes affect outputs or behavior:
   - Update `OUTPUT_SCHEMA.md` for any changes to output file formats
   - Update `README.md` for new outputs, changed behavior, or configuration options
   - Update this file (`CLAUDE.md`) for architectural changes or new patterns
4. **Test on Test Data**: Run the pipeline on test data to verify changes work:
   ```bash
   pixi run snakemake -j 2 --use-conda --configfile config/config_test_data.yaml
   ```
5. **Validate Results**: Run validation to check for regressions:
   ```bash
   pixi run python workflow/scripts/validate_test.py
   ```
6. **PR Description**: Ensure PR description accurately reflects the scope of changes
