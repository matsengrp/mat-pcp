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
6. Every codon-level mutation falls in a masked codon position

After filtering, an ``in_train`` boolean column is added using a NumPy
RNG (``np.random.default_rng(seed).random(N) < val_ratio``, then
``in_train = ~val_mask``). The call sequence matches dasm2's
``preprocess_command.preprocess`` exactly, so the new file is
byte-identical to ``filter_viral_pcps.sh`` output when given the same
inputs/seed.

The script also emits a small ``*_dasm_filter_report.txt`` log alongside
the output CSV summarizing per-filter drop counts. This is the
auditability hook from the issue body (#80).
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

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
# Mask file convention: 1-indexed codon positions (matches dasm2 default).
mask_zero_indexed = bool(snakemake.params.get("mask_zero_indexed", False))

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


def _parse_problematic_sites(filepath, zero_indexed=False):
    """Parse mask file (one site per line) into a sorted list of 0-indexed sites.

    Matches dasm2.sequence.masking.parse_problematic_sites_file: # comments and
    blank lines are skipped; the default 1-indexed convention is converted to
    0-indexed.
    """
    sites = []
    with open(filepath) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            site = int(stripped)
            sites.append(site if zero_indexed else site - 1)
    return sorted(set(sites))


def _has_unmasked_codon_mutation(parent, child, masked_codon_set):
    """Return True iff at least one codon-level mutation is at an UNMASKED site.

    ``masked_codon_set`` is the set of 0-indexed codon positions that are
    masked. A PCP whose every mutation falls in masked codons carries no
    signal once those positions are blacked out, so it's dropped.
    """
    n = len(parent)
    for codon_idx, start in enumerate(range(0, n, 3)):
        if codon_idx in masked_codon_set:
            continue
        if parent[start : start + 3] != child[start : start + 3]:
            return True
    return False


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
    empty_df.to_csv(output_dasm_ready, index=False)
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
problematic_sites = _parse_problematic_sites(
    mask_sites_path, zero_indexed=mask_zero_indexed
)
num_masked_sites = len(problematic_sites)
log(f"Loaded {num_masked_sites} masked codon sites from {mask_sites_path}")

if after_valid > 0 and num_masked_sites > 0:
    masked_codon_set = set(problematic_sites)
    keep_mask = pcp_df.apply(
        lambda row: _has_unmasked_codon_mutation(
            str(row[parent_col]), str(row[child_col]), masked_codon_set
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
pcp_df.to_csv(output_dasm_ready, index=False)
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
    f.write(f"Mask:         {mask_sites_path} ({num_masked_sites} sites)\n")
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
