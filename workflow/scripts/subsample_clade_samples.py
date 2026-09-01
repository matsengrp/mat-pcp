#!/usr/bin/env python3
"""Draw one arm's subsample of one clade's sample list (#151).

Takes the already-written ``mat/clades/{clade}_samples.txt`` and thins it,
rather than re-deriving the roster from the metadata's clade column.
``rule samples_by_clade`` matches a clade with the awk regex ``^clade( |$)``
against values like ``21L (Omicron)``; re-implementing that here would be a
second implementation agreeing on ``19A`` and disagreeing on exactly the values
the regex exists for, with both results looking like plausible clade rosters
(#149's observation, one axis over).

The draw preserves the input file's order, so ``sampling_fractions: [1.0]``
reproduces the input byte-for-byte and every output below the arm is byte-equal
to the full tree's.
"""

import os
from importlib import util as importlib_util

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
subsample = workflow_loader.load_workflow_module("subsample")

full_samples_path = snakemake.input.samples
output_path = snakemake.output.samples
clade = snakemake.wildcards.clade
arm = snakemake.wildcards.arm
unit = list(snakemake.params.unit)
seed = int(snakemake.params.seed)
frac = float(snakemake.params.frac)
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


with open(full_samples_path) as f:
    # Keep the lines verbatim, including any trailing whitespace the upstream
    # awk emitted: the degenerate arm's promise is byte-equality, and stripping
    # then re-joining would quietly normalize a file it claims to reproduce.
    samples = f.read().split("\n")

trailing_newline = samples and samples[-1] == ""
if trailing_newline:
    samples = samples[:-1]

if clade not in unit:
    # Defensive: only unit clades are thinned, and a non-unit clade reaching
    # here means a fan-out routed an untouched clade through the sampler --
    # which would silently produce an arm whose "unchanged" clades are not.
    raise ValueError(
        f"Arm {arm!r} was asked to subsample clade {clade!r}, which is not in "
        f"its unit {unit}. Only unit clades are thinned; every other clade "
        f"reads the full tree's files unchanged."
    )

log(f"Arm {arm}: clade {clade}, unit {unit}, seed {seed}, fraction {frac}")
log(f"Read {len(samples)} samples from {full_samples_path}")

if not samples:
    raise ValueError(
        f"{full_samples_path} is empty, so clade {clade!r} has nothing to "
        f"subsample for arm {arm!r}. A clade reaches the roster only because "
        f"metadata rows carry it, so an empty list here means the clade "
        f"extraction disagrees with the checkpoint."
    )

drawn = subsample.select_subsample(samples, frac, unit, seed, clade)
log(f"Retained {len(drawn)} samples (round({frac} * {len(samples)}))")

with open(output_path, "w") as f:
    f.write("\n".join(drawn))
    if trailing_newline:
        f.write("\n")

log(f"Wrote {output_path}")
