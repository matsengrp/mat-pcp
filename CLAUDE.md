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
| 3b | Date-Cutoff PCPs (#149) | `get_clades_by_cutoff` (checkpoint), `cutoff_pcp_summary`, `samples_upto_cutoff`, `mat_clade_subset_upto`, `extract_pcps_by_clade_upto`, `aggregate_clade_pcps_upto`, `aggregate_filtered_branches_upto`, `curate_pcps_for_dasm_upto`, `curate_all_pcps_for_dasm_upto`, `dasm_ready_pcps_all_clades_upto` |
| 3c | Subsampling Arms (#151) | `subsample_clade_samples`, `subsample_correlations`; everything else an arm rebuilds is an existing rule reached through the optional `{arm_dir}` prefix |
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
# all fourteen, plus a `marimo check` on the analysis notebook, before the pipeline steps)
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
pixi run python tests/test_pcp_date_cutoffs.py
pixi run python tests/test_clade_subsamples.py

# Half B of the bte test needs the Snakemake-built bte env, so it runs after a
# pipeline run; --require-bte makes a missing env a failure rather than a skip.
pixi run python tests/test_bte_haplotype.py --require-bte

# Run on test data (fast, uses max_clades=3)
pixi run snakemake -j 2 --use-conda \
  --configfile config/config_test_data.yaml

# Validate results
pixi run python workflow/scripts/validate_test.py
```

> **`CLAUDE.md`, `CHANGELOG.md`, `workflow/Snakefile`, `config/config.yaml`,
> `config/config_test_data.yaml` and
> `config/config-gisaid-and-public-BA2-JN1-XFG.yaml` use CRLF line endings** (the
> other configs are LF, so a scripted edit across `config/*.yaml` must read each
> file's ending rather than assume one). Edit them in place. A
> Python `read_text()` / `write_text()` round-trip normalizes them to LF, which
> rewrites every line and buries a three-line content edit in a ~2,200-line
> whitespace diff. Nothing errors and the rendered text is identical, so it is only
> visible in `git diff --stat`. Check with `file CLAUDE.md` after scripted edits.
>
> **The `.py` files are mixed too, so never assume.** 29 of the tracked Python
> files are CRLF, including `workflow/scripts/count_mutations.py`,
> `annotate_counts.py`, `expected_counts_by_clade.py`,
> `get_problematic_sites.py` and `workflow/modules/{glm,rates,load,aamutfit,
> probfit}.py`; the rest — including `pcp_from_tree.py`, `validate_test.py` and
> everything under `tests/` — are LF. There is no rule that predicts which, so a
> scripted edit must read each file's ending:
> `raw = path.read_bytes(); crlf = b"\r\n" in raw` → edit the decoded,
> LF-normalized text → write back with the original ending. Enumerate the
> current set with
> `for f in $(git ls-files '*.py'); do file -b "$f" | grep -q CRLF && echo "$f"; done`.

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
    *The fitted clock is drawn as a prediction, not a third slope.* The figure
    also carries chronumental's **fitted** rate, and its anchoring is different
    from the pre-fit line's on purpose. The pre-fit line is slope-only at the
    centroid because chronumental logs no intercept; the fitted line has a real
    one, since a node's date is the root date plus the cumsum of its path's
    durations, so those durations telescope and the per-branch Poisson
    expectation sums to `rate * (tip - root)` along **any** root-to-tip path.
    Drawn from zero divergence on the root's fitted date it is therefore the
    model's own prediction, and the vertical gap to the cloud is #145 measured in
    mutations: 9.36 subs/yr predicts 62 at the last collection date, where the
    cloud's own regression reads 247. Note what did *not* move -- the root's fitted date spans one day
    across the whole sweep (2019-12-10 for the shipped setting, 2019-12-11 for
    `--variance_dates 30`), so the inflation is interior, and
    "the root got pushed back" is the wrong story. Three traps. The rate must come
    from the **last** logged step: `step:0` is the pre-fit estimate after one
    update (39.16 against 40.36), so a first-match parse returns a clock that
    agrees with its own starting point and is 4.3x wrong, with nothing downstream
    able to tell. The anchor is a **sidecar** (`mat/node_dates_rtt_anchor.csv`)
    because finding the root needs the tree, which the plot rule deliberately does
    not load -- and it cannot be decoded from the log's `root_date` field, an
    internal parameterization: the two arms whose fields read -17.5946 and
    -17.5023 wrote root dates 0.81 days apart, where reading the field as days
    predicts 2.2 hours and as years predicts 34 days. And the report states the
    fitted rate against the **pre-fit**
    slope only, never against the tree's implied aggregate: that ratio is
    near-identity by construction, so a "near 1 is fine" rule passes an inflated
    tree. Guarded by two checks in `tests/test_node_date_qc.py`, both
    mutation-tested, and by the anchor reconciliation in
    `validate_test.py::validate_root_to_tip_cloud`.

    *The sweep's verdict, adopted in v12.* `--variance_dates 30` is now in every
    shipped config. On `public_2026-08-06` it moves root-to-tip/fitted from 4.31x to
    **1.48x**, `length_cor` from 0.40 to **0.85**, and the zero-mutation share of
    total branch-time from 45.9% to **1.4%**; branches longer than `window_size`
    fall from 2.87% to 0.51%. What it spends is tip-residual margin -- median
    |residual| 1 d -> **10 d** against the 30 d gate, p95 1 d -> **35 d**, which is
    *above* the bound the median clears. The gate value is deliberately left at 30:
    retuning it in the same version that changes the prior would leave nothing
    measuring whether the change was affordable. The two `HorseShoeLike` arms
    produced no numbers at all, blocked by
    [#147](https://github.com/matsengrp/mat-pcp/issues/147) -- chronumental formats
    a strain name as a float -- so the horseshoe remains untested here rather than
    rejected.

    **Beware the reference when quoting the size of the effect.** The pre-fit
    root-to-tip slope is itself **1.69x the nominal** 23.9 subs/yr (40.38 measured
    on the public tree, r=0.94 over 8.37 M tips), and nothing here explains that
    gap -- plausibly saltational stems plus the tip-weighting of heavily sampled
    post-jump lineages, but that is a hypothesis, not a measurement. So it is not a
    neutral yardstick. Against the nominal instead, v11's aggregate implied rate is
    9.47 subs/yr = **0.40x** and v12's is 27.25 = **1.14x**. Both references pick
    the same arm, so the adoption does not turn on the choice -- but "4.3x too
    slow" overstates it, and the nominal is the fairer reference for magnitude.

    *What the reference choice does not affect: the inflation is where the
    mutations are.* On the v11 public tree, 91.1% of all mutations sit on branches
    carrying 1-5 of them, and the cumulative implied rate over exactly those
    branches is **8.88 subs/yr** -- so the deficit is not a tail artifact that a
    different reference could argue away. The high-count bins run the other way
    (57.9 and 72.0 subs/yr above 20 mutations), i.e. the fit *compresses* the
    saltational branches while parking 45.9% of the tree's elapsed time on
    zero-mutation ones. That is why the punctuated-evolution reading does not
    rescue the default: the Poisson MLE for a shared rate is `sum(n)/sum(t)`, which
    has no robustness at all -- an outlier raises the numerator and pushes the rate
    *up*. The one comparison that is *never* informative is fitted-vs-implied:
    it is near-identity by construction under a learned clock (measured 1.012 at
    the default, 1.014 at `--variance_dates 30`), so a "ratio near 1 is fine" rule
    passes an inflated tree.

    See [#145](https://github.com/matsengrp/mat-pcp/issues/145) and
    [#146](https://github.com/matsengrp/mat-pcp/issues/146).

14. **A date-truncated PCP set must be re-extracted, and its clade roster is part
    of the deliverable.** For each entry in `pcp_date_cutoffs` the whole per-gene
    PCP family is rebuilt under `pcps/{gene}/upto_{cutoff}/` from the samples
    collected **at or before** that date, so DASM can be trained on what was
    knowable then (#149). Four things about the shape are load-bearing.

    *Row-filtering `all_pcps.csv.gz` is not an approximation of this -- it is
    wrong.* `matUtils extract -s` returns the subtree *induced* on the retained
    samples, so internal nodes that lose all but one child are collapsed and their
    mutations merged onto the surviving branch. A cutoff table is therefore not a
    subset of the full one, and this is measurable rather than theoretical: on the
    test data's `2020-06-01` arm, **2 of 7** ORF1ab cutoff PCPs carry a
    (parent, child) pair occurring nowhere in the full table -- sample
    `MEX/SEARCH-3636-TIJ/2020` has parent `node_13` at branch length 1 in the full
    tree and parent `node_1` at branch length **4** at the cutoff. A row filter
    would keep the PCP whose parent is not reconstructible from pre-cutoff samples
    *and* drop three mutations outright, leaving a set of pairs that is the branch
    set of no tree.

    *The cutoff filters the clade's own sample list, never the metadata's clade
    column.* `rule samples_by_clade` matches a clade with the awk regex
    `^clade( |$)` against values like `21L (Omicron)`; re-deriving that in Python
    would be a second implementation agreeing on `19A` and disagreeing on exactly
    the values the regex exists for, with both results looking like plausible clade
    rosters. `modules/cutoffs.py::select_upto` therefore takes the already-written
    `mat/clades/{clade}_samples.txt` and filters it by date alone, **preserving
    order**. That is what makes the degenerate arm exact: a cutoff equal to
    `max_leaf_date` selects everything, so its outputs are byte-identical to the
    full tree's -- verified across all four test genes for both `all_pcps.csv.gz`
    and `all_pcps_dasm_ready.csv.gz`. This is pattern 11's construction, one
    subsystem over: the ordinary path *is* the degenerate case.

    *An inadequate (cutoff, clade) pair is skipped at DAG-construction time and the
    skip is recorded.* The `get_clades_by_cutoff` checkpoint scores every pair and
    `pcps/cutoff_summary.csv` names all of them with `included` true/false, because
    a clade that does not exist yet contributes zero files and is otherwise
    indistinguishable from a wiring mistake -- #137's empty era and #144's empty
    censor table restated. `included` is *not* a copy of the checkpoint's verdict:
    it comes from the Snakefile function that drove the fan-out and is reconciled
    against the checkpoint **in both directions**, raising either way. Requested-
    but-inadequate publishes an empty subtree; adequate-but-not-requested is a
    silently declined clade, and recording it as `included=False` would hide a
    broken fan-out in the artifact meant to expose it.

    *"Adequate" needs a floor of its own, because ``min_clade_size`` is 0 in every
    production config.* On the full-tree path that is harmless -- a clade only
    reaches the roster because metadata rows carry it, so it always has samples and
    `adequate_sample_counts` is vacuously true. At a cutoff it is not: `0 >= 0`
    marked all 26 zero-sample (cutoff, clade) pairs on the v12 2026-08-06 trees
    adequate, the checkpoint committed the fan-out, and `samples_upto_cutoff`'s
    empty-list guard then raised on each -- so the skip-and-record mechanism above
    was **inert on exactly the configs it was written for**, and the failure took
    every downstream aggregate with it. `cutoffs.MIN_SAMPLES_AT_CUTOFF` floors the
    configured value, and the floor is **2, not 1**, for a second and independent
    reason: a one-sample cutoff subtree writes a well-formed `.pb` that `matUtils`
    produces without complaint, and `bte.MATree` then dies with **SIGSEGV** loading
    it (measured on 25C at 2023-01-01, six times, once per gene; subtrees cut to 2
    and 3 leaves from the 5-leaf 23E tree both load in milliseconds, so the boundary
    is exactly one leaf). A signal 11 in a C++ extension leaves no traceback and no
    partial output, so nothing above it can attribute the failure to clade size. The
    floor is structural rather than a threshold -- `min_clade_size` says how big a
    clade the user wants, this says what can be built at all -- and only the sub-2
    cases are floored: a clade with a *handful* of samples is kept even when it
    contributes zero rows for some gene (23E has 5 samples and 0 S rows), because
    that is what the full-tree path already ships under the same setting (21E has 47
    samples and 6 S rows in the released v11 tree). The same latent crash exists on
    the full-tree path, where `get_clades` applies `min_clade_size` unfloored; no
    shipped tree has a one-sample clade (the smallest on v11 public is 25I at 20),
    so it is recorded rather than changed there. **No CI run can see any of this**:
    `config_test_data.yaml` sets `min_clade_size: 2`, which floors the comparison by
    accident. `tests/test_pcp_date_cutoffs.py` therefore drives the real checkpoint
    at `min_clade_size=0` against a fixture carrying both a 0-sample and a 1-sample
    clade, and mutation-tests the floor by editing **`modules/cutoffs.py`**, not the
    script -- patching script text would report the floor as covered while the
    operator implementing it sat one file away. `validate_test.py::validate_cutoff_pcps`
    also now treats a surviving `mat/cutoffs/{clade}_samples.txt` for an
    `included=False` row as a stale intermediate, since checking only the PCP files
    left the excluded clades' subtrees invisible.

    *The cutoff cannot be answered from the metadata date alone, because some of
    those dates are wrong by years.* On the v12 public tree's 2023-01-01 set, 16
    leaves are Nextstrain 24A (JN.1) or 24E (KP.3.1.1) -- lineages that did not
    exist in 2023 -- reported between 2020-01-10 and 2022-12-09 and predicted
    **337 to 570 days past the cutoff**; their pango calls (JN.1.4, KP.3.1.1,
    MC.24) agree with the tree and not with the date. They arrive in batches
    sharing one placeholder value (20-odd `Spain/AR-IACS_LaSARD-*` all reporting
    exactly 2021-07-16 or 2022-10-18, nine `USA/ND-NDDH-260xx` all 2021-09-02),
    which is a submission date written into the collection field.
    `pcp_max_predicted_days_past_cutoff` drops a sample whose **predicted** date
    lands more than that many days after the cutoff. Four things about the shape
    are load-bearing.

    The test is against the **cutoff**, not against the reported date. Dropping
    when the two dates merely disagree by more than N is the obvious alternative
    and it is a strictly broader filter: it also removes records whose dates
    disagree while *both* sit comfortably before the cutoff, which are
    untrustworthy but not leaks. What a temporal holdout can be wrong about is
    admitting a post-cutoff genome, so that is what is tested -- 826 / 2,695
    leaves on the two production trees against 2,879 / 13,839 under the
    discordance reading. It also makes the predicate cutoff-dependent: the same
    strain is legitimately kept at a later cutoff and dropped at an earlier one.

    The tolerance is **non-zero because the predicted date is fitted**. Median
    |predicted - reported| is 10 d on both trees, p95 35/37 d, p99 56/62 d under
    the `--variance_dates 30` prior, so a sample genuinely collected days before
    the cutoff predicts past it on ordinary noise. 30 days removes 0.011% /
    0.018% of leaves while still catching every one of the 36 / 272 that sit
    more than 90 days past. *What it removes is whole clades, not a thin
    margin.* Every clade the filter touched went to **exactly zero**
    pre-cutoff samples rather than merely fewer, so the 2023-01-01 rosters
    fell from 38 to 33 included clades on the public tree and 43 to 34 on
    gisaid-and-public -- 24A (JN.1) had 14 leaves by reported date and none
    survive. That is the shape to expect: a lineage that did not exist yet
    has no legitimate pre-cutoff sample, so a clade surviving with a
    *handful* of leaves would be the surprising outcome, not this. The bound
    is a config value rather than a constant because it is a property of the
    date prior, which pattern 13 has already moved once.

    The comparison is **day-granular on both sides**. `node_dates.tsv` carries a
    time (`2023-12-10 16:57:02.752815`) while the cutoff is a date, so comparing
    raw timestamps would make membership depend on the time of day chronumental
    happened to emit -- a prediction of exactly cutoff+X at 12:00 excluded and
    the same prediction at 00:00 kept.

    Two knock-on effects that are easy to miss. **The cutoff sets now depend on
    the node dating**, so changing `chronumental_kwargs` changes which PCPs
    exist; the filter is confined to the cutoff path precisely because the
    full-tree outputs never read the reported date (clade rosters come from the
    clade column, the windowed path keys on predicted dates), so filtering there
    would move every output while fixing nothing. And **the degenerate arm's
    byte-equality becomes conditional**: classifying an arm by
    `n_samples_upto == n_samples_full` would reclassify it as truncating the
    moment one sample was dropped, and the one check that says the cutoff
    machinery left the ordinary path alone would silently stop running -- pattern
    10's vacuity. `validate_cutoff_pcps` classifies on
    `n_samples_upto + n_predicted_late_dropped` instead, and then asserts
    byte-equality when nothing was dropped **or byte-inequality when something
    was**, the latter being the only end-to-end evidence that the filter reaches
    the extraction rather than merely being counted in the roster.

    A strain chronumental did not date is **kept and counted** (154 on the public
    tree): dropping it would remove data on a missing prediction rather than a
    contradicting one, and keeping it silently would hide a half-missing
    chronumental output. Reported as `n_undated_kept`.

    No CI run can see any of this either: at X=30 the filter removes **nothing**
    on the test data (0 leaves past cutoff+30 in either arm; the degenerate arm
    has none even past cutoff+0), so a run is green whether the filter works, is
    inverted, or never executes. `tests/test_pcp_date_cutoffs.py` therefore
    drives both real scripts against a fixture whose predictions carry a 12:00
    time as production's do, with a row that is 91 days discordant yet 61 days
    *before* the cutoff -- so a mutant implementing the discordance reading gives
    a different answer rather than the same one -- and mutation-tests the filter
    by editing **`modules/cutoffs.py`**, not the script text.

    *The masks are shared across cutoffs, which is a real leak, bounded and
    stated.* `{gene}_codon_sites_to_mask.txt` and `censored_mutations.csv` come from
    the full tree, so a cutoff's DASM-ready files are filtered with a mask that saw
    post-cutoff data. Recomputing them per cutoff means re-running the whole
    counting and GLM path per cutoff. The leak is confined to which sites are
    flagged as *sequencing artifacts*, not to any fitness signal, and the censor
    table is keyed by clade so extra post-cutoff clades in it are inert. Clade
    founders are shared for a stronger reason: a founder is a property of the
    clade's defining mutations, not of which samples have been collected.

    Two things no pipeline run can check, which is why
    `tests/test_pcp_date_cutoffs.py` exists. Every test-data date is *interior* to
    both CI arms, so a comparison written `<` for `<=` moves no sample and the
    counts, the roster and the byte-equality all still agree -- inclusivity is
    unobservable end-to-end. And the byte-equality rests on order preservation,
    but the test data's clade lists are already sorted, so a script that sorted its
    output would pass CI. Both are checked on a fixture built so each failure moves
    a sample, against deliberately corrupted copies of the shipped scripts -- and
    the inclusivity mutant lives in `modules/cutoffs.py`, not in the script, so the
    test copies the whole `workflow/` tree rather than patching script text, which
    would report the guard as covered while the operator sat one file away.
    `validate_test.py::validate_cutoff_pcps` is the end-to-end counterpart: it
    reconciles the roster against the file tree in both directions, byte-compares
    the degenerate arm, and asserts no cutoff PCP has a leaf child collected after
    its cutoff.

    See [#149](https://github.com/matsengrp/mat-pcp/issues/149).

15. **A subsampled arm caps its large clades and reuses everything else, and the
    degenerate arm is the ordinary path rather than a second implementation.**
    Measuring when a pooled model beats per-clade counts needs the *same* clade
    re-estimated at several sampling depths (#151, for dasm2-experiments#1130).
    An **arm** is one per-clade cap: `max_samples: X` thins every clade with
    more than X samples down to exactly X and leaves smaller clades untouched.
    It intervenes on the tree's **composition**, which is what the
    pooled-vs-per-clade question actually asks about. Seven things about the
    shape are load-bearing.

    *The hook is one file, which is what makes the counts arm and the PCP arm
    move together.* `mat/clades/{clade}_samples.txt` is the only input
    determining a clade's subtree: `translate_mat`, `mat_sample_paths` and
    `extract_pcps_by_clade` all read `mat_clade_subset` and nothing below it
    re-reads the full tree. So `subsample_clade_samples` is the *only* rule an
    arm adds — everything else is an existing rule reached through an optional
    `{arm_dir}` prefix, and the two halves cannot drift apart because there is
    no second wiring. This diverges from #149, which duplicated nine `_upto`
    rules; that was affordable for a PCP-only feature and is not for one
    spanning ~21 rules, and it is why #149's `_upto` family and #151's
    `{arm_dir}` coexist as two mechanisms for one shape.

    *The optional prefix is what buys byte-identity by construction.*
    Snakemake has no optional path components, so `arm_dir` is declared as
    `|cap_[A-Za-z0-9._+-]+/` — the empty alternative *is* the full-tree path.
    A cap at or above every clade's size therefore reproduces the full tree
    byte-for-byte through the same rule bodies, rather than through a second
    rule set that has to be kept in sync. Note where that identity does **not**
    come from: an untouched clade is never drawn at all. `resolve_unit` admits
    a clade only when `n > cap`, so the fan-out reads the full tree's file
    directly and `subsample_clade_samples` **raises** if one arrives anyway.
    `select_cap`'s early return for `n <= cap` is unreachable from the
    pipeline and kept only to state the invariant for other callers. Two
    details still make the identity real: the draw sorts its chosen indices
    back into the input file's order, so an arm's list is a *subsequence* of
    the full tree's (a sampler emitting picks in draw order is still a correct
    uniform sample, and no pipeline check would see the difference — CI asserts
    set containment and the retained count); and gzip's FNAME field carries
    only the *basename*, so a same-named file in a different directory
    compresses to the same bytes — verified, and the reason byte-equality is
    available at all given `pcpio.PCP_CSV_COMPRESSION`'s `mtime=0` (pattern 9).

    *One arm spans the depth curve, and its controls are free.* At
    `max_samples: 50000` on the v12 gisaid tree, 30 of 56 clades are capped and
    their effective fractions run 1.15% (21J) to 100%; the 26 untouched clades
    are byte-identical controls **by construction**. That is the whole reason
    this is a cap rather than a fraction. A `frac` kind, which applied one
    sampling fraction to an explicitly *named* unit of clades, was implemented
    first and **removed** in #151. It intervened on **size** — each clade kept
    the same share, so the tree's sampling imbalance survived intact and one
    arm was one point on a depth curve — which answers a different question;
    its controls had to be *chosen*, a clade someone remembered to leave
    unthinned; and naming a unit honestly is harder than it looks, because a
    clade left at full depth sets an evidence floor through founder proximity
    (measured: 21I's spike founder differs from 21J's by one amino acid while
    carrying ~10% of its samples, and 21L alone leaves 2,737,183 samples within
    20 nt of 21J's founder — 139% of that unit — so an honest 21J ladder needs
    a five-clade unit). A cap leaves no large clade at full depth at all, which
    dissolves the problem instead of documenting it. `validate_entries`
    rejects `sampling_fractions` and `unit` **by name** with that explanation,
    rather than as unknown keys: an operator reading "unknown key" is as likely
    to conclude the key is misspelled as that the mechanism is gone.

    *There are two fan-outs, not one, and conflating them thins a clade the
    cap never reached.* Anything sliced out of the arm's own re-aggregated grid
    (`fitness/{arm}/...` downward) takes the arm prefix for **every** clade in
    the roster — `per_clade_inputs`. Anything derived from a clade's **subtree**
    (`mutation_counts`, `mutations/`, `pcps/`) takes it only for *unit* clades
    and reads the full tree's file otherwise — `per_unit_clade_inputs`. Using
    the grid-derived fan-out for the subtree-derived paths requests
    `mat/clades/{arm}/19A_samples.txt` for a clade below the cap, and the
    sampler would have to dutifully thin it: an arm whose "unchanged" clades
    are quietly subsampled, with every count still well-formed. The sampler
    now raises instead, deciding membership from the sample list it has already
    read rather than from a resolved unit param — resolving a checkpoint inside
    `params` would tie the rule's scheduling to it. The split is also most of
    the cost saving: an arm re-*aggregates* all 56 clades because
    `build_count_table` needs the whole grid, but re-*counts* only the unit.

    *`resolve_unit` is the single definition of membership, and it must be.*
    An arm's unit depends on the tree's clade sizes, which only the
    `get_clades` checkpoint knows, so it cannot be written in the config. Every
    fan-out helper, the roster param and the summary read that one function; a
    second derivation would let the fan-out and the sampler disagree about one
    clade. The unit is also resolved against the *processed* roster, not
    `clade_counts.csv`, which names every clade the metadata mentions including
    ones dropped by `clades_to_exclude` — 21A on every production tree — for
    which requesting `mat/clades/{arm}/21A_samples.txt` would ask for a subtree
    that is never extracted. Membership is **strictly** above the cap: a clade
    of exactly `cap` samples is already at the target, so thinning it would
    write an identical file into the arm directory and turn one of the arm's
    own controls into an arm row. `cap_rng` keys on `(cap, seed, clade)` and
    deliberately **omits the unit** — keying on a derived unit would make each
    clade's draw depend on how many *other* clades happened to exceed the cap,
    so a tree gaining one large clade would silently redraw every arm with
    nothing able to attribute it. That key is a string literal, so a refactor
    could rewrite it without any assertion noticing; removing the frac kind is
    checked against it, and all three production arms' sample lists reproduce
    byte-for-byte through the cap-only module.

    *`glm/` straddles the reuse line and the cut has to be deliberate.*
    `glm/{arm}/annotated_counts.csv` is count-bearing and is rebuilt;
    `glm/curated_counts.csv` is the GLM's *training set* and is not. But the
    chain `count_mutations_clade → ... → annotate_counts → curate_counts →
    master_table → rates_table.csv` makes refitting the **default**, so
    `curate_counts`, `curate_counts_era` and `get_problematic_sites` are pinned
    to literal full-tree paths rather than `rules.annotate_counts.output` /
    `rules.build_count_table.output` — those references now carry `{arm_dir}`
    and would drag the whole GLM into every arm. A refit arm would still produce
    well-formed tables whose expected counts were normalized against rates
    fitted to the thinned data, so `validate_subsample_arms` asserts no
    `glm/{arm}/{curated_counts,master_table,rates_table}*` exists. Sharing the
    fit is conservative for the hypothesis under test: it can only make the
    counts estimator look better. The per-clade `evol_time` normalization *does*
    track the subsample, which is the intended effect.

    *A nominal cut is not an evidence cut, and the two thresholds differ in
    kind.* A branch survives whenever **any** sample below it is retained, so
    `E[retained](p) = Σ k_node (1 - (1-p)^n_node)` and retention is concave: on
    the v11 gisaid tree 21J's retained/p runs 1.79x at p=0.20, 2.22x at 0.10,
    2.72x at 0.05 and 4.10x at 0.01, making a 1% draw a ~24x evidence cut
    rather than a 100x one. Measured on the production caps, the tree keeps
    23.2% of its mutation counts from 11.5% of its sequences at X=50000. Every
    downstream number must key to `correlation_summary.csv`'s measured
    `count_fraction`, never to the configured cap — which is not even a
    fraction. Separately: a unit clade thinned below `min_clade_size` is
    **excluded** (`get_clades` applies that threshold to the full metadata and
    so cannot see a clade that falls below it only after thinning; reachable
    exactly when `min_clade_size > cap`), while a draw of **zero** is a **hard failure** —
    `matUtils extract -s` on an empty sample list is undefined. No accepted
    config can produce one, since `max_samples >= 1` is validated and
    membership requires `n > cap`, so that guard is defence in depth against a
    crash inside `matUtils` that nothing above it could attribute to clade size
    (#149's SIGSEGV, one axis over). Collapsing the two into one branch is
    #149's `MIN_SAMPLES_AT_CUTOFF` story restated: there, `min_clade_size: 0`
    made the skip-and-record mechanism inert on exactly the configs it was
    written for.

    *The arms ship no roster file, and that is a deliberate exception to
    #137/#149.* A `subsample_summary.csv` recording every `(arm, clade)` pair
    with `included` true/false, plus `branch_qc/{clade}.csv` sidecars on
    `count_mutations_clade` and `extract_pcps_by_clade` as its only inputs, was
    built and then **removed**. Three facts decided it. Nothing read it —
    `subsample_correlations` declared it not an input by design,
    `notebooks/cap_arm_scatters.py` discovers arms by listing directories, and
    the production `correlation_summary.csv` was produced on the v12 gisaid tree
    without the file ever existing. Its one irreplaceable column, the excluded
    clade, is **unreachable from every production config**: all eight set
    `min_clade_size: 0`, and exclusion needs `min_clade_size > cap`, so the
    auditable-absence mechanism was inert on exactly the trees it was written
    for — `MIN_SAMPLES_AT_CUTOFF`'s failure mode (#149) recurring in the
    machinery built to answer it. And requesting it on a released tree
    cascaded: its retention denominators are the **full tree's** `branch_qc/`
    sidecars, co-declared outputs of two rules that are absent on the v12 trees
    (0 of 56 clades), so building it meant a GLM refit and a 470 GiB PCP
    rewrite for a file with no reader.
    What that costs, stated rather than absorbed: an excluded clade now
    contributes no row and no files, so it is indistinguishable from a clade a
    wiring mistake dropped, and `tests/test_clade_subsamples.py` is its only
    coverage of any kind. What it does not cost: `failure_reason` still
    attributes each dropped branch to one limit in one order, since the
    per-limit tally survives as a **log** line in both producers — but nothing
    keys on the attribution now, so the two orders no longer have to agree and
    neither script constrains the other.

    **No shipped production config sets `clade_subsamples`, so no production run
    can check any of this**, and `config_test_data.yaml`'s clade sizes (19A at 2
    samples, 19B at 5, 20A at 46) are far below the shipped caps of
    50000/20000/5000, so those would thin nothing there. The test config
    therefore carries `max_samples: [1000, 20]`: `cap_1000_0` thins nothing and
    is the **degenerate arm** — the one shape where the fan-out helpers, the
    roster and the summary all have nothing to iterate — while `cap_20_0` thins
    20A and leaves 19A and 19B as controls, so one CI run exercises a thinned
    clade and two controls at once. Two things stay unreachable even there: the
    membership boundary (`n > cap` against `n >= cap`) needs a clade of exactly
    `cap` samples, and the `min_clade_size` exclusion needs
    `min_clade_size > cap`.
    `tests/test_clade_subsamples.py` therefore drives the real
    `subsample_clade_samples.py` and
    `subsample_correlations.py` against fixtures built so each failure moves a
    sample (seven deliberately unsorted samples at a cap of 4, so every draw is
    partial and the ordering is observable), and mutation-tests every guard by
    editing **`workflow/modules/subsample.py`**, not the script text — patching
    the script would report the sampler as covered while the operator
    implementing it sat one file away. `validate_test.py::validate_subsample_arms`
    is the end-to-end counterpart. It reconciles `correlation_summary.csv`
    against the file tree in both directions — that file's `thinned` column
    comes from `_arm_roster_param`, the function that drove the fan-out, so it
    is not a second evaluation of the sizing arithmetic — asserts the arm
    directory holds *only* clades strictly above the cap, byte-compares the
    degenerate arm, and reports **VACUOUS** rather than success when no arm
    thinned a surviving clade.

    *Two invariants make an arm comparable to full depth at all, and both must be
    asserted rather than assumed.* `build_count_table` enumerates the grid from
    `clade_founder_nts.csv` and fills 0, so an arm's `counts_all.csv` must differ
    from the full tree's in the **count columns and nothing else** -- compared
    row-for-row over the clades the arm kept, which is also what shows the shared
    codon mask is a no-op rather than a leak, since `exclude` and
    `masked_in_usher` are among the compared columns. And because
    `genome_composition` counts *grid* rows, the per-clade normalization
    `K = sum(expected) / sum(actual)` over synonymous training rows is a property
    of the grid and the shared rates alone: it must not move with the cap, to
    1e-10 relative. A K that moved would shift every
    `delta_fitness_ratio` by a factor no correlation against full depth could
    attribute -- the arm would look like a degraded estimate rather than a
    renormalized one. Asserted **per clade**, since pooling lets one clade's K
    rise while another's falls. Both were mutation-tested against corrupted
    output (one flipped `exclude` value; `expected_count` scaled by 1.001) rather
    than trusted because they passed. The third leg, that the GLM is not refit, is
    checked two ways: no `glm/{arm}/{curated_counts,master_table,rates_table}*`
    file exists, and a dry run on an arm target schedules none of those rules
    while `expected_counts_by_clade` reads the **full-tree**
    `glm/curated_counts_{era}.csv` -- verified, with the three shared files'
    md5s unchanged across an arm rebuild.

    *The analysis is part of the deliverable, and its population is not the
    grid.* `correlation_summary.csv` reports rho per
    `(arm, clade, estimator, filter)` over the **intersection of scored
    mutations** -- `actual_count > 0` on both sides. Every grid row carries a
    `delta_fitness_ratio` even where nothing was observed (121,348 of 121,427
    rows for 20A on the test tree), and those rows are invariant to the cap by
    construction, so correlating over them measures
    `build_count_table`'s enumeration and returns a confident number describing
    nothing: 0.619 for the ratio estimator where the 29 genuinely scored
    mutations give 1.000. What the thinning cost is reported as
    `count_fraction` -- `sum(actual_count_arm) / sum(actual_count_full)` over
    the whole grid -- and that is the **only** retention fraction reported. A
    fraction of *distinct scored mutations* is the obvious alternative and is
    deliberately absent: counting mutations rather than observations treats one
    falling from 900 counts to 3 as fully retained, so it diverges by up to 30x
    (21J at cap 5,000, 0.0033 against 0.103) and it predicts the fitted slope
    *worse* -- Spearman 0.985 for `count_fraction` against 0.940 for it, over
    the 120 completed cap fits. Two retention fractions differing threefold
    with no stated reason to prefer one is a standing invitation to quote
    whichever is kinder, so only one ships. `n_scored_full` / `n_scored_arm`
    remain as the raw deletion counts, because that channel is real and
    separable -- a mutation resting on one observation is *deleted* rather than
    degraded when it goes, which on the test tree is the more alarming number
    (50% of samples keeps 61.6% of nucleotide *events* but only **37%** of
    20A's scored amino-acid mutations).
    Three further points. `n_mutations` is on **every** row including the empty
    ones, and an unmeasurable population is kept with a blank coefficient rather
    than dropped, since a missing row reads as "not applicable" where `n=0`
    reads as "asked, and nothing to answer with" -- on the test tree the
    `>=10` filter matches zero rows, so this is the ordinary case. Non-unit
    clades are reported as **controls** and must reproduce full depth exactly,
    which is the end-to-end assertion that the two fan-outs stayed separate.
    And `max_actual_count_full` is reported because when it is 1 the ratio
    estimator is collinear with the expected count by construction, so its
    r = 1.000 is arithmetic -- the test tree's whole regime, which is why both
    `validate_test.py` and `notebooks/cap_arm_scatters.py` say
    **NOT A NOISE MEASUREMENT** rather than printing a reassuring coefficient.
    Guarded by six checks in `tests/test_clade_subsamples.py`, whose fixture is
    built so the grid correlation and the ground-truth correlation must
    **disagree in sign** (80 sub-threshold rows agreeing perfectly, 6
    ground-truth rows anti-correlated), so a script correlating the grid is
    caught rather than merely unproven.

    *The population is fixed on the full-depth row, and the obvious alternative
    is survivorship bias.* A mutation is in the ground-truth set when
    `max(actual_count, expected_count) >= 10` **at full depth**, so the set is
    identical at every sampling depth. Reading the threshold off the *arm's*
    row instead — which this script did first — drops exactly the mutations the
    thinning destroyed, so the denominator shrinks with the cut and a degrading
    estimator reports a flat correlation over an ever-easier population. Nothing
    end-to-end can see it: every coefficient is well-formed under the wrong
    reading and the population it selects looks entirely plausible, it is simply
    a different set per arm. The fixture therefore drops the arm's counts below
    the threshold while leaving full depth's above it, so the correct
    implementation compares 6 mutations where the survivorship mutant compares
    0. And **r is never reported alone**, because it is invariant to any linear
    transform of the arm's estimates and so cannot see either a shift or a
    compression: at `cap_5000_0` 21J's regression slope is 0.117 while r still
    reads 0.431, and its estimates' sd falls from 1.77 to 0.48 — an 8.5x
    compression that r prices at a halving.
    `rmse`, `bias`, `regression_slope`, `sd_full` and `sd_arm` sit beside every
    coefficient for that reason.

    *An arm that thins **nothing** is legitimate.* A cap at or above every
    clade's size is the degenerate arm: its unit is empty, every clade is a
    control, and every output must byte-equal the full tree's. It still appears
    in `correlation_summary.csv`, since that file has a row per clade of the
    arm's whole grid and a degenerate arm's grid is all controls — so it is the
    one arm shape whose rows must *all* read `thinned=False`, and
    `validate_subsample_arms` treats the absence of a `mat/clades/{arm}/`
    directory as its defining property.

    See [#151](https://github.com/matsengrp/mat-pcp/issues/151).


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
   residual exceeds it). The QC rule runs on every tree. **Every shipped config
   sets `--variance_dates 30` as of v12** (#146); see pattern 13 for what that
   buys and what it spends.
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
6. **PCP Extraction**: `pcp_date_cutoffs` and
   `pcp_max_predicted_days_past_cutoff` (optional, per-tree overridable; see
   pattern 14), plus `pcp_genes` (PCP branch QC reuses the shared
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
6b. **Clade Subsampling** (#151): `clade_subsamples` — a list of
   `{max_samples, seed}` entries, each expanding to one **arm** per cap. An arm
   thins every clade above the cap down to it and leaves smaller clades
   untouched, so its membership is derived from the tree's clade sizes rather
   than named. Empty in every production config; arms are launched by naming
   their paths. `sampling_fractions` / `unit` belonged to the removed fraction
   arms and are rejected by name. See pattern 15.
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
- `pcps/`: Curated parent-child pairs, plus each gene's `censored_mutations.csv` DASM training censor table and its `.clades.txt` roster (#144), each gene's `upto_{cutoff}/` date-truncated PCP families and the `cutoff_summary.csv` roster (#149), each gene's `{arm}/` reduced-depth PCP families (#151)
- `correlation_summary.csv` (output root): arm-vs-full-depth fitness correlations, keyed to `n_mutations` and `count_fraction` (#151)
- `notebooks/cap_arm_scatters.py`: the arm marimo notebook -- sequences per clade coloured by emergence date and the tree-level sequence and mutation totals, then per clade the count fraction, the count eCDFs and the ground-truth scatter (#151)
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
