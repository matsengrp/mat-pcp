#!/usr/bin/env python3
"""Filter PCPs into a DASM-training-ready CSV with a train/val split.

Reimplements the filters that ``dasm2 preprocess`` applies (and that
``data-central/scripts/filter_viral_pcps.sh`` invokes today) using only
NumPy + Pandas, so the rule stays in ``envs/default.yml`` rather than
pulling in torch + dasm2 + netam.

Filters applied in order:

1. Empty parent OR empty child sequence
2. Length mismatch parent vs child
3. Length not divisible by 3 (non-codon-aligned)
4. Ambiguous nucleotides (``N``) in either sequence
5. Stop codons (``TAA``, ``TAG``, ``TGA``) anywhere in parent or child
6. Every nucleotide-level mutation falls at a masked nucleotide position

After filtering, an ``in_train`` boolean column is added using a NumPy
RNG (``np.random.default_rng(seed).random(N) < val_ratio``, then
``in_train = ~val_mask``). The call sequence matches dasm2's
``preprocess_command.preprocess`` exactly.

**The equivalence to ``filter_viral_pcps.sh`` covers filters 1-5 and the
train/val RNG, not filter 6.** It did cover all six until #158: dasm2's own
masking is still codon-granular, so where a masked codon contains an
unmasked mutated nucleotide the two implementations now deliberately
disagree -- this one keeps the PCP and dasm2 drops it. Filters 1-5 plus the
split are unchanged and still byte-identical given the same inputs and seed.
Since #134 this output is gzipped, so that equivalence is on the
**decompressed** stream -- the reference writes plain CSV, and gzip headers
additionally embed the filename and (by default) the mtime. It is checked
that way by ``validate_test.py``'s ``DASM_READY_ALL_SENTINEL_MD5``, which
has never exercised filter 6 in either implementation; see
``modules/pcpio.py``.

Filter 6 reads ``pcps/{gene}/masked_nt_sites.csv`` (#153), whose
``nt_index`` column is a **0-based** index into this gene's PCP nucleotide
sequence and is used as-is. It replaced ``{gene}_codon_sites_to_mask.txt``,
which was the same clade-independent mask widened to codons, so one suspect
nucleotide took its two innocent neighbours with it -- a 2.85x over-mask
across the six v12 production genes (121 masked sites became 345 masked
nucleotides). The widening mattered because the *loss* never saw the codon
file: DASM censors at nucleotide granularity from ``masked_nt_sites.csv``
and ``censored_mutations.csv``, so a PCP was discarded here as carrying no
signal outside masked codons while the loss it would have fed would have
seen signal, at an unmasked nucleotide inside one of those codons. One mask
now serves both (#158).

The script also emits a small ``*_dasm_filter_report.txt`` log alongside
the output CSV summarizing per-filter drop counts. This is the
auditability hook from the issue body (#80).
"""

import json
import os
from importlib import util as importlib_util
from pathlib import Path

import numpy as np
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


masking = load_workflow_module("masking")
pcpio = load_workflow_module("pcpio")

# ---------------------------------------------------------------------------
# Snakemake glue
# ---------------------------------------------------------------------------
input_pcps_path = snakemake.input.pcps
mask_sites_path = snakemake.input.mask_sites
output_dasm_ready = snakemake.output.dasm_ready
output_report = snakemake.output.report
output_stats = snakemake.output.stats

seed = int(snakemake.params.seed)
val_fraction = float(snakemake.params.val_fraction)
parent_col = snakemake.params.get("parent_col", "parent")
child_col = snakemake.params.get("child_col", "child")

log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


# ---------------------------------------------------------------------------
# Filters (vendored from dasm2.data.validation + dasm2.sequence.masking)
# ---------------------------------------------------------------------------
STOP_CODONS = frozenset(("TAA", "TAG", "TGA"))


def _contains_stop_codon(nt_seq):
    """Return True iff ``nt_seq`` contains a stop codon in frame."""
    for i in range(0, len(nt_seq), 3):
        if nt_seq[i : i + 3] in STOP_CODONS:
            return True
    return False


def _is_pcp_valid(parent, child):
    """Match dasm2.data.validation.is_pcp_valid semantics on a single row."""
    if not parent or not child:
        return False
    if len(parent) != len(child):
        return False
    if len(parent) % 3 != 0:
        return False
    if "N" in parent.upper() or "N" in child.upper():
        return False
    if _contains_stop_codon(parent) or _contains_stop_codon(child):
        return False
    return True


def _has_unmasked_nt_mutation(parent, child, masked_indices):
    """Return True iff at least one nucleotide differs at an UNMASKED position.

    ``masked_indices`` is a sorted ``np.ndarray`` of 0-based positions into the
    gene's PCP nucleotide sequence -- the distinct ``nt_index`` values of
    ``masked_nt_sites.csv``, used without adjustment. A PCP whose every
    difference sits at a masked position carries no signal once those positions
    are blacked out, so it is dropped.

    Vectorized rather than looped because this is the one rule in the pipeline
    whose cost is unmeasured and whose input is largest: 4.3 M PCPs x 13,203 nt
    on gisaid ORF1a. Comparing the two encoded buffers finds the (typically
    single-digit) differing positions in one pass, so this is *faster* than the
    codon loop it replaced despite testing three times as many positions. The
    plain-Python reference implementation lives in
    ``tests/test_dasm_ready_nt_filter.py``.
    """
    parent_bytes = np.frombuffer(parent.encode(), dtype=np.uint8)
    child_bytes = np.frombuffer(child.encode(), dtype=np.uint8)
    diff_idx = np.nonzero(parent_bytes != child_bytes)[0]
    if diff_idx.size == 0:
        return False
    return bool(np.isin(diff_idx, masked_indices, invert=True).any())


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
log("=" * 60)
log(f"DASM-ready PCP curation: seed={seed}, val_fraction={val_fraction}")
log(f"  input: {input_pcps_path}")
log(f"  mask:  {mask_sites_path}")
log(f"  out:   {output_dasm_ready}")
log("=" * 60)

pcp_df = pd.read_csv(input_pcps_path)
initial_count = len(pcp_df)
log(f"Loaded {initial_count} PCPs")

if initial_count == 0:
    # Empty-input edge case (issue body explicitly calls this out as a valid
    # outcome for clades where every PCP gets filtered). Write the schema
    # header with the in_train column appended and stop.
    out_cols = list(pcp_df.columns) + ["in_train"]
    empty_df = pd.DataFrame(columns=out_cols)
    Path(output_dasm_ready).parent.mkdir(parents=True, exist_ok=True)
    empty_df.to_csv(
        output_dasm_ready, index=False, compression=pcpio.PCP_CSV_COMPRESSION
    )
    with open(output_report, "w") as f:
        f.write("Input was empty; no PCPs to filter.\n")
    with open(output_stats, "w") as f:
        json.dump(
            {
                "initial_count": 0,
                "after_valid_filter": 0,
                "after_mask_filter": 0,
                "final_count": 0,
                "removed_invalid": 0,
                "removed_masked_only": 0,
                "num_masked_sites": 0,
                "seed": seed,
                "val_fraction": val_fraction,
            },
            f,
            indent=2,
        )
    log("Empty input - wrote header-only output.")
    raise SystemExit(0)

# --- Step 1: validity filter ------------------------------------------------
valid_mask = pcp_df.apply(
    lambda row: _is_pcp_valid(str(row[parent_col]), str(row[child_col])),
    axis=1,
)
pcp_df = pcp_df.loc[valid_mask].reset_index(drop=True)
after_valid = len(pcp_df)
log(f"After validity filter: {after_valid} (removed {initial_count - after_valid})")

# --- Step 2: masked-only-mutation filter -----------------------------------
# The reader lives in modules/masking.py so its missing-column guard can be
# mutation-tested there rather than in this script's text (#158). It returns
# DISTINCT indices: masked_nt_sites.csv is long-form, one row per
# (nt_index, reason), and 18 of the 121 masked sites on the six v12 production
# genes carry two reasons.
problematic_sites = masking.read_masked_nt_indices(mask_sites_path)
num_masked_sites = len(problematic_sites)
log(
    f"Loaded {num_masked_sites} distinct masked nucleotide positions from "
    f"{mask_sites_path}"
)

# The `num_masked_sites > 0` arm is deliberate, not a micro-optimization: a
# header-only mask is legitimate (masked_nt_sites.py declares it legal for a
# gene with no globally masked site), and skipping the filter then keeps every
# valid PCP, which is the right answer. Note the two arms disagree about a
# zero-mutation PCP -- the skip keeps it, the predicate drops it -- so do not
# collapse them. That is unreachable today only because pcp_from_tree.py drops
# `parent == child` rows before writing.
if after_valid > 0 and num_masked_sites > 0:
    masked_indices = np.asarray(problematic_sites, dtype=np.int64)
    keep_mask = pcp_df.apply(
        lambda row: _has_unmasked_nt_mutation(
            str(row[parent_col]), str(row[child_col]), masked_indices
        ),
        axis=1,
    )
    pcp_df = pcp_df.loc[keep_mask].reset_index(drop=True)
after_mask = len(pcp_df)
log(f"After masked-only filter: {after_mask} (removed {after_valid - after_mask})")

# --- Step 3: train/val split -----------------------------------------------
# Must match dasm2.preprocess_command.preprocess byte-for-byte:
#   rng = np.random.default_rng(seed=seed)
#   val_mask = rng.random(len(df)) < val_ratio
#   in_train = ~val_mask
rng = np.random.default_rng(seed=seed)
val_mask_arr = rng.random(len(pcp_df)) < val_fraction
pcp_df["in_train"] = ~val_mask_arr
n_train = int(pcp_df["in_train"].sum())
n_val = int((~pcp_df["in_train"]).sum())
log(f"Train/val split: {n_train} train, {n_val} val ({val_fraction:.0%} val target)")

# --- Write outputs ----------------------------------------------------------
Path(output_dasm_ready).parent.mkdir(parents=True, exist_ok=True)
pcp_df.to_csv(
    output_dasm_ready, index=False, compression=pcpio.PCP_CSV_COMPRESSION
)
log(f"Wrote {len(pcp_df)} DASM-ready PCPs to {output_dasm_ready}")

stats = {
    "initial_count": initial_count,
    "after_valid_filter": after_valid,
    "after_mask_filter": after_mask,
    "final_count": len(pcp_df),
    "removed_invalid": initial_count - after_valid,
    "removed_masked_only": after_valid - after_mask,
    "num_masked_sites": num_masked_sites,
    "n_train": n_train,
    "n_val": n_val,
    "seed": seed,
    "val_fraction": val_fraction,
    "mask_file": str(mask_sites_path),
}
with open(output_stats, "w") as f:
    json.dump(stats, f, indent=2)

with open(output_report, "w") as f:
    f.write("DASM-ready PCP filter report\n")
    f.write("=" * 60 + "\n")
    f.write(f"Input:        {input_pcps_path}\n")
    f.write(
        f"Mask:         {mask_sites_path} "
        f"({num_masked_sites} masked nt positions)\n"
    )
    f.write(f"Seed:         {seed}\n")
    f.write(f"val_fraction: {val_fraction}\n")
    f.write("-" * 60 + "\n")
    f.write(f"Initial PCPs:                {initial_count}\n")
    f.write(f"Removed (invalid):           {initial_count - after_valid}\n")
    f.write(f"Removed (only masked muts):  {after_valid - after_mask}\n")
    f.write(f"Final PCPs:                  {len(pcp_df)}\n")
    f.write(f"  train: {n_train}\n")
    f.write(f"  val:   {n_val}\n")
log(f"Wrote report to {output_report}")
log("Done.")
