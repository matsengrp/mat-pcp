#!/usr/bin/env python3
"""Draw one arm's subsample of one clade's sample list (#151).

Takes the already-written ``mat/clades/{clade}_samples.txt`` and thins it,
rather than re-deriving the roster from the metadata's clade column.
``rule samples_by_clade`` matches a clade with the awk regex ``^clade( |$)``
against values like ``21L (Omicron)``; re-implementing that here would be a
second implementation agreeing on ``19A`` and disagreeing on exactly the values
the regex exists for, with both results looking like plausible clade rosters
(#149's observation, one axis over).

The draw preserves the input file's order, so an arm's list is a
**subsequence** of the full tree's rather than a permutation of a subset. Note
what that does *not* buy: a clade at or below the cap never reaches this script
at all -- the membership check below raises -- because an untouched clade is a
control by reading the full tree's file, not by having a copy of it drawn here.

**Membership is decided from the sample list, not from a param.** An arm's unit
is derived (``n_clade > max_samples``), and the number this script would need to
re-derive it -- the clade's full sample count -- is exactly what it has just
read. Applying the predicate to the file rather than to a resolved list keeps
the sampler independent of the ``get_clades`` checkpoint (resolving one inside
``params`` would tie this rule's scheduling to it) while still failing loudly if
a fan-out routes an untouched clade here: that would silently produce an arm
whose "unchanged" clades are not.
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
spec = dict(snakemake.params.spec)
seed = int(spec["seed"])
cap = int(spec["cap"])
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


with open(full_samples_path) as f:
    # Keep the lines verbatim, including any trailing whitespace the upstream
    # awk emitted. Stripping and re-joining would quietly normalize the sample
    # names, which are the join key every downstream comparison against full
    # depth uses.
    samples = f.read().split("\n")

trailing_newline = samples and samples[-1] == ""
if trailing_newline:
    samples = samples[:-1]

log(f"Arm {arm} ({subsample.depth_label(spec)}): clade {clade}, seed {seed}")
log(f"Read {len(samples)} samples from {full_samples_path}")

# Emptiness is checked *before* membership, and the order is load-bearing. Every
# cap is >= 1, so an empty list is also at or below the cap and would otherwise
# be reported as "this clade is not one of the arm's" -- true but not the
# problem. An empty list means the clade extraction disagrees with the
# checkpoint, which is a different bug with a different owner.
if not samples:
    raise ValueError(
        f"{full_samples_path} is empty, so clade {clade!r} has nothing to "
        f"subsample for arm {arm!r}. A clade reaches the roster only because "
        f"metadata rows carry it, so an empty list here means the clade "
        f"extraction disagrees with the checkpoint."
    )

# Defensive: a clade the arm does not thin reaching here means a fan-out routed
# an untouched clade through the sampler.
if len(samples) <= cap:
    raise ValueError(
        f"Arm {arm!r} was asked to subsample clade {clade!r}, which has "
        f"{len(samples)} samples -- at or below the cap of {cap}, so it is "
        f"not one of this arm's clades. A clade at or below the cap is "
        f"untouched and reads the full tree's files unchanged; writing a "
        f"copy here would make the arm directory claim a clade it did not "
        f"thin."
    )

drawn = subsample.select_cap(samples, cap, seed, clade)
log(
    f"Retained {len(drawn)} of {len(samples)} samples "
    f"({subsample.depth_label(spec)})"
)

with open(output_path, "w") as f:
    f.write("\n".join(drawn))
    if trailing_newline:
        f.write("\n")

log(f"Wrote {output_path}")
