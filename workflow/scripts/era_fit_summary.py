#!/usr/bin/env python3
"""Record which eras got a neutral-rate fit and which were skipped (#137).

An era with no clades on this tree is skipped: no ``curated_counts_{era}.csv``,
no ``master_table_{era}.csv``, no ``rates_table_{era}.csv``. That is the correct
behaviour -- ``populate_rates`` cannot fit an empty partition, and every shipped
single-era tree (the ``max_clades: 3`` test config, and the BA.2/JN.1/XFG trees)
hits it -- but *an unlogged skip is the #129/#130/#133 silent-drop shape*: the
output tree simply lacks a file, which reads identically to a wiring mistake that
dropped an era that was in fact present.

So the skip is recorded in a tracked output rather than only in a log. This file
names **every** era in ``era.ERAS`` with ``fitted`` true/false, so "Omicron was
absent from this tree" and "Omicron went missing" are distinguishable after the
fact, from the output tree alone.

It also cross-checks the DAG against the era module: the eras whose tables were
requested as inputs must be exactly the non-empty eras of the clade list. Those
two are derived independently -- the Snakefile's ``eras_present`` resolves the
``get_clades`` checkpoint, this re-partitions the clade list handed to it -- so a
disagreement means one of them is wrong and neither can be trusted.
"""

import os
import pandas as pd
from importlib import util as importlib_util

output_csv = snakemake.output.summary
clades = list(snakemake.params.clades)
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


# The by-path module loader is defined once in
# workflow/modules/workflow_loader.py. Because that loader is itself the
# mechanism for importing shared workflow code, it cannot be imported via
# itself -- it is wired up here by absolute path (#111).
_loader_spec = importlib_util.spec_from_file_location(
    "workflow_loader",
    os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "modules",
        "workflow_loader.py",
    ),
)
workflow_loader = importlib_util.module_from_spec(_loader_spec)
_loader_spec.loader.exec_module(workflow_loader)
load_workflow_module = workflow_loader.load_workflow_module

script_dir = os.path.dirname(os.path.abspath(__file__))
workflow_dir = os.path.dirname(script_dir)
era_mod = load_workflow_module("era", workflow_dir)


log(f"Summarizing era neutral-rate fits over {len(clades)} clades")

partition = era_mod.partition_clades(clades)
glm_dir = os.path.dirname(os.path.abspath(output_csv))

# Which eras' tables the DAG actually built, read off the declared inputs rather
# than assumed -- see the module docstring on why these are cross-checked.
built = set()
for path in snakemake.input:
    name = os.path.basename(path)
    for era_name in era_mod.ERAS:
        if name == f"master_table_{era_name}.csv":
            built.add(era_name)

expected = set(partition)
if built != expected:
    raise ValueError(
        f"the DAG built master tables for eras {sorted(built)} but the clade "
        f"list partitions into {sorted(expected)}. The Snakefile's eras_present "
        f"and era.partition_clades must agree; they are derived independently, "
        f"so a disagreement means one is wrong.\n  clades: {clades}"
    )

rows = []
for era_name in era_mod.ERAS:
    era_clades = partition.get(era_name, [])
    fitted = era_name in expected
    rows.append(
        {
            "era": era_name,
            "fitted": fitted,
            "n_clades": len(era_clades),
            "clades": ";".join(era_clades),
            "master_table": (
                f"glm/master_table_{era_name}.csv" if fitted else ""
            ),
            "rates_table": (f"rates_table_{era_name}.csv" if fitted else ""),
        }
    )
    if fitted:
        log(f"  FITTED  {era_name}: {len(era_clades)} clades -> {era_clades}")
    else:
        log(
            f"  SKIPPED {era_name}: no clades on this tree, so no rate table is "
            f"emitted for it"
        )

summary = pd.DataFrame(rows)
summary.to_csv(output_csv, index=False)

# Every clade must land in exactly one era. partition_clades cannot produce an
# overlap (it appends to one bucket per clade), but the total is asserted so a
# future refactor to two list literals cannot silently drop or double-count one.
n_partitioned = int(summary["n_clades"].sum())
if n_partitioned != len(clades):
    raise ValueError(
        f"the era partition covers {n_partitioned} clades but {len(clades)} were "
        f"supplied. The eras must be a partition of the clade list: disjoint and "
        f"covering."
    )

log(f"Wrote era fit summary to {output_csv}")
log(
    f"{int(summary['fitted'].sum())} of {len(era_mod.ERAS)} eras fitted, "
    f"covering all {len(clades)} clades"
)
log("Done!")
