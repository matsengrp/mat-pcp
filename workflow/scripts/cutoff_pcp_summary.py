#!/usr/bin/env python3
"""Record which (cutoff, clade) pairs got a truncated PCP set, and which did not.

A clade below ``min_clade_size`` at a cutoff -- or one that does not exist yet --
is skipped at DAG-construction time, so it contributes no files. In the output
tree alone that is indistinguishable from a clade dropped by a wiring mistake,
which is #137's and #144's problem restated: the absence has to be *recorded*
somewhere, or "no files" and "broken" read the same. This is that record.

It is not a copy of the checkpoint. The ``included`` column comes from the
Snakefile function that actually requested the per-cutoff files, and is
reconciled against the checkpoint's adequacy verdict **in both directions**:

  * requested but not adequate -> raise. This is #137's "clades present but
    routed to an unfitted era" case: the DAG built something the checkpoint said
    to skip, so an empty or malformed subtree is about to be published.
  * adequate but not requested -> raise. A clade the checkpoint approved and the
    DAG silently declined is the failure this file exists to make visible, and
    reporting it as ``included=False`` would hide it in the very artifact meant
    to expose it.

A single ``try/except`` spanning both, or a summary that merely mirrors the
checkpoint's own column, would turn a broken fan-out into a plausible-looking
table -- the shape #137 calls out explicitly.
"""

import pandas as pd

counts_in = snakemake.input.counts
summary_out = snakemake.output.summary
requested = {str(k): sorted(v) for k, v in dict(snakemake.params.requested).items()}
genes = list(snakemake.params.genes)
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None

COLUMNS = [
    "cutoff",
    "clade",
    "output_subdir",
    "n_samples_full",
    "n_samples_upto",
    "adequate_sample_counts",
    "included",
]


def log(msg):
    if log_file:
        with open(log_file, "a") as handle:
            handle.write(msg + "\n")
    print(msg)


def main():
    df = pd.read_csv(counts_in)
    df["included"] = [
        clade in requested.get(str(cutoff), [])
        for cutoff, clade in zip(df["cutoff"], df["clade"])
    ]
    df["output_subdir"] = [f"upto_{cutoff}" for cutoff in df["cutoff"]]

    adequate = set(
        zip(df.loc[df["adequate_sample_counts"], "cutoff"].astype(str),
            df.loc[df["adequate_sample_counts"], "clade"])
    )
    asked = {(cutoff, clade) for cutoff, clades in requested.items() for clade in clades}

    routed_wrong = sorted(asked - adequate)
    if routed_wrong:
        raise ValueError(
            f"{len(routed_wrong)} (cutoff, clade) pair(s) had truncated PCPs "
            f"requested despite failing the checkpoint's adequacy test, e.g. "
            f"{routed_wrong[:5]}. The DAG and mat/clade_counts_by_cutoff.csv "
            f"disagree; the published subtree would be empty or below "
            f"min_clade_size."
        )
    dropped_silently = sorted(adequate - asked)
    if dropped_silently:
        raise ValueError(
            f"{len(dropped_silently)} (cutoff, clade) pair(s) are adequate in "
            f"mat/clade_counts_by_cutoff.csv but were never requested, e.g. "
            f"{dropped_silently[:5]}. Recording these as included=False would "
            f"hide a broken fan-out in the file whose job is to expose it."
        )

    # Cutoffs the config asked for that the checkpoint never scored at all. The
    # two lists come from the same config, so a mismatch means the checkpoint
    # was built against a different one -- a stale output tree.
    scored = set(df["cutoff"].astype(str))
    missing = sorted(set(requested) - scored)
    if missing:
        raise ValueError(
            f"Cutoff(s) {missing} were requested but do not appear in "
            f"{counts_in}. The checkpoint predates the current config; rebuild "
            f"mat/clade_counts_by_cutoff.csv."
        )

    df[COLUMNS].to_csv(summary_out, index=False)

    for cutoff in sorted(scored):
        arm = df[df["cutoff"].astype(str) == cutoff]
        n_in = int(arm["included"].sum())
        skipped = sorted(arm.loc[~arm["included"], "clade"])
        log(
            f"cutoff {cutoff} -> upto_{cutoff}/: {n_in} of {len(arm)} clades "
            f"included across {len(genes)} gene(s) "
            f"({int(arm.loc[arm['included'], 'n_samples_upto'].sum())} samples); "
            f"skipped {skipped or 'none'}"
        )
    log(f"Wrote {len(df)} (cutoff, clade) rows to {summary_out}")


if __name__ == "__main__":
    main()
