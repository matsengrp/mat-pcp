# Generating `data/usher_masked_sites.yaml` from the UCSC branch-specific mask

This document describes how `data/usher_masked_sites.yaml` is produced by
`workflow/scripts/generate_usher_masked_sites.py`, how the generator works
stage by stage, and how the regenerated masks are verified to recover the
old, hand-maintained masking scheme. It is intended as a reference for future
maintenance and for fine-tuning the generation logic.

## The problem the generator solves

Two masking systems with incompatible keys:

- **Old (mat-pcp):** `data/usher_masked_sites.yaml`, keyed by **Nextstrain
  clade** (`21A`, `22B`). Inheritance was *faked* by manually repeating the
  same sites across every descendant clade's block. Masking is by exact string
  match — no tree involved at runtime.
- **New (UCSC `branchSpecificMask.yml`):** keyed by **PANGO lineage** (`BA.4`,
  `BQ.1.1`), each anchored to a **representative sample**. Upstream,
  `matUtils mask` applies each mask to that sample's node *and its entire
  subtree*, so inheritance is implicit in tree topology. This spec replaced the
  old `maskDelta.sh` approach; see
  [yatisht/usher#324](https://github.com/yatisht/usher/issues/324).

`generate_usher_masked_sites.py` is an **offline** bridge: it reads the PANGO
spec, uses a full public MAT to compute real tree ancestry, coarsens
PANGO → Nextstrain, and emits the flat Nextstrain-keyed file the pipeline
already knows how to consume. The pipeline's runtime masking logic is
essentially unchanged.

## How the script works, stage by stage

### 1. Parse the UCSC spec (`parse_branch_spec`)

For each of the 25 PANGO clades it collects the `representative`
(+ `representativeBacktrack`), expands every `ranges: [[lo, hi]]` into individual
integer positions and merges with `sites` → a flat **site set**; validates each
`reversions` string against `^[ACGT][0-9]+[ACGT]$` → a **mutation list**; keeps
`exclusions` aside. It also asserts the spec has a `version` key (fails loudly
if the upstream format changes).

### 2. Get node IDs from the tree

Two `matUtils extract` passes on the full `.pb.gz`:

- `--clade-paths` → `{clade: founder_node_id}` — the founder (MRCA) node of
  *every* annotated clade, both PANGO (`BA.4`) and Nextstrain (`"22A (BA.4)"`).
  This is the **primary** node source (`parse_clade_paths`).
- `--sample-paths` → each sample's root → tip node chain, used to read ancestor
  chains and as a **fallback** representative-sample resolver.

### 3. Resolve each PANGO clade to a node (`resolve_clade_nodes`)

Prefer the annotated founder node from `--clade-paths`. This changed
mid-implementation: 8 of 25 representative samples are GISAID-only names absent
from the *public* tree, so representative-sample resolution alone failed.
Founder nodes are annotated for all 25, so we use those and keep the
representative-sample path (`resolve_representative_nodes`, mirroring upstream's
backtrack logic) only as a fallback.

### 4. Enumerate Nextstrain clades (`load_nextstrain_clade_values`)

From the metadata `Nextstrain_clade` column. Keeps full values (`"22A (BA.4)"`)
for founder lookup; the **output key** is the bare token (`"22A"`) via
`bare_clade`, exactly matching the pipeline's `get_clades` awk split.

### 5. Assign PANGO masks to Nextstrain clades (`assign_masks_by_founder`)

The core logic. Working from **founder nodes** (not individual samples) is
deliberate: it is noise-free, whereas a per-sample membership rule lets
recombinants and mislabeled samples bleed masks across unrelated clades. A
Nextstrain clade `C` inherits PANGO clade `P`'s masks iff:

- **(a)** `founder(P)` lies on the root → `founder(C)` chain — `P` is an ancestor
  of `C`, so all of `C` descends from `P`; the mask applies **exactly**; or
- **(b)** `C` is the *deepest* Nextstrain clade enclosing `founder(P)` (last
  Nextstrain founder on the root → `founder(P)` chain) — `P` is finer than any
  Nextstrain clade, so its mask is **widened** to all of `C` (over-mask
  coarsening).

Rule (b) is what fixed an earlier bug where basal/trunk clades absorbed every
descendant's mask. It attributes a fine PANGO clade to the clade its samples
*actually carry* — e.g. `KP.3.1.1` → `24E`, not its ancestor `24A`; `21L` (BA.2)
gets only BA.2's own masks. Clades with no founder node (non-monophyletic
`recombinant`) are logged and left unmasked.

### 6. Accumulate & emit (`build_clade_masks` → `write_output`)

Union the site sets and mutation lists of all PANGO clades assigned to each
Nextstrain clade; write one `mask_<bare_clade>` block with `clades` / `sites` /
`mutations`, a `# accumulated from PANGO clades: [...]` provenance comment per
block, and a file header (source commit, tree, metadata, dropped-exclusion
note for BA.2.86 site 23599). Result: **43 clade blocks, 6,937 site-entries,
1,126 mutation-entries.**

Mapping of source fields:

- `ranges` + `sites` → **site masks**.
- `reversions` → **mutation masks** (directional `refPOSalt`, matching
  `nt_mutation`, *not* expanded to whole sites).
- `exclusions` → dropped + logged (cannot express per-sample un-masking at
  clade granularity).

## How the new masks match the old ones

The generator enforces an explicit acceptance criterion (`check_recovers_old`):
the new file must be a **superset** — for every `(clade, site)` in the old file,
that site must still be masked for that clade. An old site counts as recovered
if it is either:

1. still masked as a **site**, or
2. now masked at the **mutation level** — i.e. ≥ 1 new mutation sits at that
   position for that clade.

Case 2 is an intentional, *more precise* "downgrade": the old file blunt-masked
whole sites, while the UCSC spec masks specific reversion mutations at some of
those positions. Only sites recovered by **neither** mechanism are true
failures. Running with `--old-yaml` reported:

- **Exit 0 — full recovery.** Every old `(clade, site)` mask is preserved.
- **141 site → mutation downgrades** (in `22B` / `22D` / `22E`): the old
  whole-site masks there are now the finer directional mutation masks; zero of
  them truly missing.
- Everything else preserved as-is, plus large new coverage for post-2022 clades
  the hand-made file never had.

`--strict` makes any *true* shortfall a hard error; `--union-old` would
force-include unrecovered old masks — neither was needed since recovery is
clean.

## Consumer-side changes (backward-compatible)

Both consumer scripts were extended to read the new optional `mutations:` field.
Blocks without `mutations` behave exactly as before.

- `workflow/scripts/count_mutations.py` unions `mutations` into
  `muts_to_exclude` (filtered via `nt_mutation.isin(...)`).
- `workflow/scripts/aggregate_mutations_to_exclude.py` emits one row per
  mutation (position parsed from the middle of the string) rather than expanding
  over all 4 nucleotides.

## Regeneration procedure

Run offline against the **full** public MAT (all clades must be present /
annotated — do not run against a filtered per-run tree):

```bash
python workflow/scripts/generate_usher_masked_sites.py \
    --branch-spec data/branchSpecificMask.yml \
    --tree public-YYYY-MM-DD.all.masked.pb.gz \
    --metadata public-YYYY-MM-DD.metadata.tsv.gz \
    --old-yaml data/usher_masked_sites.yaml \
    --output data/usher_masked_sites.yaml \
    --strict
```

`matUtils extract --sample-paths` on an ~8M-sample tree produces a large paths
file (tens of GB); the generator reads it only until every founder node's
ancestor chain has been seen (early exit), so memory stays modest. Point
`--workdir` at a disk with room, or pass precomputed `--sample-paths` and
`--clade-paths` files to skip `matUtils`.

To update the vendored spec: re-download `branchSpecificMask.yml` from the UCSC
kent repo at a newer commit, update the provenance block in
`data/branchSpecificMask.yml`, then regenerate as above.

## Caveats / accepted approximations

- **Coarsening loss (accepted).** PANGO clades finer than any Nextstrain clade
  are widened to the enclosing Nextstrain clade (over-mask). Logged per clade.
- **Reversions.** Handled at the mutation level, not blanket-masked over the
  site.
- **Exclusions.** Per-sample un-masking of a site cannot be expressed at
  Nextstrain granularity → skipped + logged (currently only BA.2.86 site 23599).
- **Recombinants.** The non-monophyletic `recombinant` clade has no single
  founder node and is left unmasked.
- **Drift.** The UCSC spec grows as new lineages are designated; the vendored
  copy is pinned to a specific commit and should be refreshed periodically.
