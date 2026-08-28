#!/usr/bin/env python3
"""Restrict one clade's sample list to samples collected at or before a cutoff.

The input is the clade's *already-written* ``mat/clades/{clade}_samples.txt``,
not the metadata's clade column. That is the whole point: ``rule
samples_by_clade`` selects a clade with the awk regex ``^clade( |$)`` against
values like ``21L (Omicron)``, and re-deriving that match here would be a second
implementation of it -- one that agrees on ``19A`` and disagrees on exactly the
values the regex exists for, with both results looking like plausible clade
rosters. Filtering the existing list by date alone also makes the degenerate arm
exact: a cutoff at or after every sample's date returns the file unchanged, so
its PCPs are byte-identical to the full tree's (#149).

Samples whose *predicted* date lands more than
``pcp_max_predicted_days_past_cutoff`` days after the cutoff are dropped as well
-- a metadata date that is wrong by years otherwise admits a post-cutoff genome
into a pre-cutoff set. Note that this makes the degenerate arm's byte-equality
conditional on the filter removing nothing, which ``validate_cutoff_pcps``
checks rather than assumes.
"""

import os
from importlib import util as importlib_util

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

samples_in = snakemake.input.samples
metadata = snakemake.input.metadata
node_dates = snakemake.input.node_dates
samples_out = snakemake.output.samples
max_days_past_cutoff = snakemake.params.max_days_past_cutoff
cutoff = snakemake.wildcards.cutoff
clade = snakemake.wildcards.clade
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    if log_file:
        with open(log_file, "a") as handle:
            handle.write(msg + "\n")
    print(msg)


def main():
    strains = cutoffs_mod.read_sample_list(samples_in)
    strain_dates = cutoffs_mod.load_strain_dates(metadata)
    predicted = (
        None
        if max_days_past_cutoff is None
        else cutoffs_mod.load_predicted_dates(node_dates)
    )
    selection = cutoffs_mod.select_for_cutoff(
        strains, strain_dates, cutoff, predicted, max_days_past_cutoff
    )
    kept = selection.kept

    # One name per line, in the input file's order -- ``matUtils extract -s``
    # reads this, and reordering it would change the output for no reason
    # (pattern 9: data-central manifests pin a hash per file).
    with open(samples_out, "w") as handle:
        for strain in kept:
            handle.write(f"{strain}\n")

    log(
        f"clade {clade} at cutoff {cutoff}: kept {len(kept)} of {len(strains)} "
        f"samples ({selection.n_after_cutoff} collected after the cutoff, "
        f"{selection.n_predicted_late} predicted more than "
        f"{max_days_past_cutoff} days past it, {selection.n_undated} kept with "
        f"no predicted date)"
    )
    if not kept:
        # Not raised here: whether an empty clade is an error is a DAG-level
        # question, and the DAG has already answered it -- the checkpoint marks
        # a clade below ``min_clade_size`` at this cutoff as inadequate and the
        # per-cutoff rules are never requested for it. Reaching this line means
        # a clade was requested that the checkpoint said to skip, so say so
        # loudly rather than handing matUtils an empty -s file, which fails with
        # an error naming neither the clade nor the cutoff.
        raise ValueError(
            f"Clade {clade!r} has no samples at or before {cutoff}, yet its "
            f"cutoff subtree was requested. The checkpoint "
            f"(mat/clade_counts_by_cutoff.csv) should have marked it "
            f"inadequate; this means the DAG and the checkpoint disagree."
        )


if __name__ == "__main__":
    main()
