# GLM reference-equivalence fixtures (#126)

These two files pin mat-pcp's GLM neutral-rate estimation to the reference
implementation in
[`SARS2-mut-fitness-v2`](https://github.com/matsengrp/SARS2-mut-fitness-v2)
(`modules/load.py` + `modules/rates.py` + `modules/glm.py`, driven by
`notebook/master_tables.py.ipynb`).

| file | what it is |
|---|---|
| `glm_curated_counts_fixture.csv` | A stratified subsample of a real `glm/curated_counts.csv` (6,047 rows; 3,161 synonymous). **Deliberately contains all four `mut_class` values**, so a test that trains on it can tell synonymous-only from pooled training. Not a *proportional* sample: non-synonymous rows are thinned by a further 4x, since they are never trained on and only need to be numerous enough to show that pooling them in would change the fit. |
| `glm_reference_rates.csv` | The 512-row rate table that the *reference* implementation produces from that fixture. |

`tests/test_glm_reference_equivalence.py` feeds the fixture through mat-pcp's
own `modules/load.py` + `modules/rates.py` and asserts the resulting rates
match `glm_reference_rates.csv` exactly. It also fits the grid a second time on
the *unfiltered* fixture and requires that most conditions differ (currently
433 of 512), so the first two checks cannot pass vacuously. A fourth check
asserts the shipped configs still carry the reference light-switch boundaries,
which the pipeline passes explicitly rather than taking from the module default.

Both the fixture fit and the frozen reference table contain negative rates, and
running the test prints a `WARNING: N conditions have negative rates` block even
on success. That is expected: the subsample's thinnest strata (45 synonymous rows
against 8-9 GLM parameters) extrapolate below zero. These fixtures pin
*agreement with the reference*, not the plausibility of the fit.

## Why a frozen fixture rather than pipeline output

The reference pipeline consumes Bloom-lab count tables for a fixed clade set,
while mat-pcp counts mutations off whatever MAT it is pointed at. Absolute
mutation counts therefore differ between the two pipelines and between mat-pcp
trees, and rates are *not* invariant to that: `Rates.populate_rates` rescales by
`n_training_rows / total_mutations` (`df_train.shape[0]`, i.e. site x mut-type
entries, not distinct sites), so tree size changes every rate.

Feeding **both implementations the same curated counts** removes that variable
entirely. The comparison is then exact rather than approximate, and it stays
valid no matter which tree the pipeline is run on.

## Provenance

The committed fixture was subsampled from a **v7** test run
(`results/test_data_v7/glm/curated_counts.csv`, 87,777 rows: 64,183
nonsynonymous / 18,934 synonymous / 4,129 nonsense / 531 noncoding). A v8 run
reproduces that file byte-for-byte (md5 `89fbce7c…`), since the test config
already ran metadata QC before #122 made it unconditional.

Do not rely on that holding in general: the composition is tree-specific, and
production trees — which did run with `skip_qc_filtering: true` before #122 —
will differ. Nothing in the tests depends on these numbers; the fixture is frozen
precisely so the comparison does not move with the pipeline.

## Regenerating

Only needed if the reference implementation itself changes. Requires a local
clone of `SARS2-mut-fitness-v2` and any completed run's curated counts:

```bash
python tests/generate_glm_reference_fixtures.py \
    --reference-repo ../SARS2-mut-fitness-v2 \
    --curated-counts results/<label>_<version>/glm/curated_counts.csv
```

The generator hard-fails if the resulting fixture cannot support the tests: if it
holds fewer than two `mut_class` values, or if its *synonymous* rows do not reach
past the largest light-switch boundary (21562). The second matters because
`populate_rates` takes `s_max` from the training set, so a synonymous subset that
stopped short of the boundary would collapse the `True`/`False` halves of the
AT/CG/GC/CT grid into duplicate rates.

Regenerating is a deliberate act: it re-baselines what "agrees with the
reference" means, so it belongs in its own commit with a stated reason.
