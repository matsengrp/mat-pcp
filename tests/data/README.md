# GLM reference-equivalence fixtures (#126)

These two files pin mat-pcp's GLM neutral-rate estimation to the reference
implementation in
[`SARS2-mut-fitness-v2`](https://github.com/matsengrp/SARS2-mut-fitness-v2)
(`modules/load.py` + `modules/rates.py` + `modules/glm.py`, driven by
`notebook/master_tables.py.ipynb`).

| file | what it is |
|---|---|
| `glm_curated_counts_fixture.csv` | A stratified subsample of a real `glm/curated_counts.csv`. **Deliberately contains all four `mut_class` values**, so a test that trains on it can tell synonymous-only from pooled training. |
| `glm_reference_rates.csv` | The 512-row rate table that the *reference* implementation produces from that fixture. |

`tests/test_glm_reference_equivalence.py` feeds the fixture through mat-pcp's
own `modules/load.py` + `modules/rates.py` and asserts the resulting rates
match `glm_reference_rates.csv` exactly.

## Why a frozen fixture rather than pipeline output

The reference pipeline consumes Bloom-lab count tables for a fixed clade set,
while mat-pcp counts mutations off whatever MAT it is pointed at. Absolute
mutation counts therefore differ between the two pipelines and between mat-pcp
trees, and rates are *not* invariant to that: `Rates.populate_rates` rescales by
`n_synonymous_sites / total_mutations`, so tree size changes every rate.

Feeding **both implementations the same curated counts** removes that variable
entirely. The comparison is then exact rather than approximate, and it stays
valid no matter which tree the pipeline is run on.

## Regenerating

Only needed if the reference implementation itself changes. Requires a local
clone of `SARS2-mut-fitness-v2`:

```bash
python tests/generate_glm_reference_fixtures.py \
    --reference-repo ../SARS2-mut-fitness-v2 \
    --curated-counts results/test_data_v7/glm/curated_counts.csv
```

Regenerating is a deliberate act: it re-baselines what "agrees with the
reference" means, so it belongs in its own commit with a stated reason.
