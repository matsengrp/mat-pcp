#!/usr/bin/env python3
"""Bucket branch mutations into sliding time windows over the full dated tree.

Thin Snakemake wrapper around ``workflow/modules/windowing.py`` (EPIC #115
Phase 2, #117). Reads ``snakemake.*``, builds a ``WindowedMAT``, runs one
window/step layout, and writes the two per-window CSVs:

  * ``windows/founder_muts.csv``   -- per-window subtree founders + haplotype.
  * ``windows/subtree_counts.csv`` -- per-window per-mutation actual counts,
    coding only, with mat-pcp multi-gene column collapse.
"""

import os
from datetime import datetime, timedelta
from importlib import util as importlib_util

import pandas as pd


# The by-path module loader is defined once in
# workflow/modules/workflow_loader.py; wire it up by absolute path (#111).
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
windowing = workflow_loader.load_workflow_module("windowing")


# --- snakemake I/O -----------------------------------------------------------
tree_path = snakemake.input.tree
node_dates_path = snakemake.input.node_dates
gtf_path = snakemake.input.gtf
fasta_path = snakemake.input.fasta
coding_sites_csv = snakemake.input.coding_sites

founder_muts_out = snakemake.output.founder_muts
subtree_counts_out = snakemake.output.subtree_counts

window_size_days = int(snakemake.params.window_size)
window_step_days = int(snakemake.params.window_step)
window_start_time = snakemake.params.window_start_time
window_end_time = snakemake.params.window_end_time
max_nt_mutations = snakemake.params.max_nt_mutations
max_reversions_to_ref = snakemake.params.max_reversions_to_ref
max_reversions_to_founder = snakemake.params.max_reversions_to_clade_founder

log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    if log_file:
        with open(log_file, "a") as f:
            f.write(str(msg) + "\n")
    print(msg)


def parse_time(value):
    """Parse an ISO date/datetime config value, or None to defer to the tree."""
    if value is None or value == "":
        return None
    return datetime.fromisoformat(str(value))


FOUNDER_COLUMNS = [
    "founder_id",
    "window_start",
    "subtree_size",
    "branches_with_mutations",
    "founder_nt_mutations",
]


def main():
    window_size = timedelta(days=window_size_days)
    step_size = timedelta(days=window_step_days)

    mat = windowing.WindowedMAT(
        tree_path,
        node_dates_path,
        gtf_path,
        fasta_path,
        coding_sites_csv,
        max_nt_mutations=max_nt_mutations,
        max_reversions_to_ref=max_reversions_to_ref,
        max_reversions_to_founder=max_reversions_to_founder,
        start_time=parse_time(window_start_time),
        end_time=parse_time(window_end_time),
        log=log,
    )

    subtree_data = mat.aggregate_data(window_size, step_size)

    founder_rows = []
    count_records = []
    for window_start, subtrees in subtree_data:
        window_str = window_start.strftime("%Y-%m-%d")
        for (
            founder_id,
            founder_muts,
            subtree_size,
            branches_with_mutations,
        ), records in subtrees:
            founder_rows.append(
                (
                    founder_id,
                    window_str,
                    subtree_size,
                    branches_with_mutations,
                    " ".join(
                        sorted(founder_muts, key=lambda m: int(m[1:-1]))
                    ),
                )
            )
            for (
                gene,
                nt_index,
                aa_index,
                aa,
                nuc,
                original_codon,
                alternative_codon,
                is_synonymous,
                count,
            ) in records:
                count_records.append(
                    {
                        "window_start": window_str,
                        "founder_id": founder_id,
                        "protein": gene,
                        "nt_site": nt_index,
                        "nt_mutation": nuc,
                        "aa_mutation": aa,
                        "codon_change": f"{original_codon}>{alternative_codon}",
                        "synonymous": is_synonymous,
                        "count": count,
                    }
                )

    founder_df = pd.DataFrame(founder_rows, columns=FOUNDER_COLUMNS)
    founder_df.to_csv(founder_muts_out, index=False)
    log(f"Wrote {len(founder_df)} founder rows to {founder_muts_out}")

    # Collapse per-gene records to one row per (window_start, founder_id,
    # nt_mutation), matching mat-pcp aggregated.csv (;-joined genes).
    counts_df = windowing.collapse_multigene(count_records, log=log)
    counts_df.to_csv(subtree_counts_out, index=False)
    log(f"Wrote {len(counts_df)} mutation-count rows to {subtree_counts_out}")


if __name__ == "__main__":
    main()
