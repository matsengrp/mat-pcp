#!/usr/bin/env python3
"""Count each clade's samples at each date cutoff, and mark the adequate ones.

This is the ``get_clades_by_cutoff`` checkpoint's body. It decides, at
DAG-construction time, which (cutoff, clade) pairs get a truncated PCP set --
the shape #137 used for empty eras. A clade that does not yet exist at a cutoff,
or that sits below ``min_clade_size`` there, is **skipped rather than emitted
empty**: a header-only PCP file would look like a real one, and every count
downstream would read as a true zero.

Its inputs are the per-clade sample lists ``rule samples_by_clade`` already
wrote, so the clade vocabulary here is by construction the same one the
full-tree path uses -- there is no second clade-matching implementation to drift
(see ``modules/cutoffs.py``). Restricting to those files also means
``clades_to_exclude`` and ``max_clades`` are already applied: this file can only
ever be a subset of the full roster, so no clade can appear at a cutoff that the
full tree does not process.

One long-format file rather than one file per cutoff, so the whole decision is a
single checkpoint: N wildcard checkpoint outputs would each have to be resolved
separately, and a cutoff absent from the config would then be indistinguishable
from one whose file had not been built.
"""

import os
from importlib import util as importlib_util

import pandas as pd

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

cutoffs_mod = load_workflow_module("cutoffs")

metadata = snakemake.input.metadata
sample_files = list(snakemake.input.samples)
counts_out = snakemake.output.counts
cutoffs = list(snakemake.params.cutoffs)
min_clade_size = snakemake.params.min_clade_size
clade_column = snakemake.params.clade_column
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None

COLUMNS = [
    "cutoff",
    "clade",
    "n_samples_full",
    "n_samples_upto",
    "adequate_sample_counts",
]


def log(msg):
    if log_file:
        with open(log_file, "a") as handle:
            handle.write(msg + "\n")
    print(msg)


def clade_from_samples_path(path):
    """Derive the clade from a ``{clade}_samples.txt`` filename.

    Anchored on the suffix and raising otherwise, for the reason
    ``aggregate_clade_pcps.clade_from_pcp_path`` documents: a silently wrong
    clade label leaves the row count unchanged and mislabels every row
    identically, so nothing downstream can see it (pattern 9).
    """
    basename = os.path.basename(path)
    suffix = "_samples.txt"
    if not basename.endswith(suffix):
        raise ValueError(
            f"Sample-list input {basename!r} does not end in {suffix!r}; "
            f"cannot derive its clade."
        )
    return basename[: -len(suffix)]


def main():
    strain_dates = cutoffs_mod.load_strain_dates(metadata)
    log(f"Loaded {len(strain_dates)} strain dates from {metadata}")
    log(f"Cutoffs: {cutoffs}; min_clade_size={min_clade_size}")

    rows = []
    for cutoff in cutoffs:
        for path in sorted(sample_files):
            clade = clade_from_samples_path(path)
            strains = cutoffs_mod.read_sample_list(path)
            kept = cutoffs_mod.select_upto(strains, strain_dates, cutoff)
            rows.append(
                {
                    "cutoff": cutoff,
                    "clade": clade,
                    "n_samples_full": len(strains),
                    "n_samples_upto": len(kept),
                    "adequate_sample_counts": len(kept) >= min_clade_size,
                }
            )

    df = pd.DataFrame(rows, columns=COLUMNS)
    # Sorted so the file is byte-reproducible regardless of the order Snakemake
    # happened to hand over the per-clade inputs (pattern 9).
    df = df.sort_values(["cutoff", "clade"]).reset_index(drop=True)
    df.to_csv(counts_out, index=False)

    for cutoff in cutoffs:
        arm = df[df["cutoff"] == cutoff]
        adequate = arm[arm["adequate_sample_counts"]]
        skipped = sorted(arm.loc[~arm["adequate_sample_counts"], "clade"])
        log(
            f"cutoff {cutoff}: {len(adequate)} of {len(arm)} clades adequate "
            f"(>= {min_clade_size} samples); skipped {skipped or 'none'}"
        )
        if arm.empty or adequate.empty:
            raise ValueError(
                f"cutoff {cutoff}: no clade has {min_clade_size} samples at or "
                f"before it, so the cutoff would produce nothing at all. Either "
                f"the cutoff predates the data (check "
                f"qc_filter_options.min_leaf_date) or min_clade_size is too "
                f"high for it. An empty cutoff is a config error, not a result: "
                f"every downstream file would be absent and look like a build "
                f"that had not run."
            )

    log(f"Wrote {len(df)} (cutoff, clade) rows to {counts_out}")
    log(f"Clade column (for reference): {clade_column}")


if __name__ == "__main__":
    main()
