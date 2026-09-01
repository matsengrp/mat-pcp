#!/usr/bin/env python3
"""Record what each subsampling arm actually retained (#151).

A nominal sampling fraction is **not** an evidence cut. A branch survives
whenever *any* sample below it is retained, so retention of mutations runs well
above the fraction of samples drawn -- on the v11 gisaid tree, 1% of 21J's
samples keeps 4.10% of its nucleotide mutations and 10% keeps 22.2%, a curve
that is concave rather than proportional. Every correlation computed downstream
therefore has to be keyed to *measured* evidence, which is what this table is
for.

It is also the auditable-absence half of the deliverable. A unit clade thinned
below ``min_clade_size`` contributes zero files, and in the output tree alone
that is indistinguishable from a wiring mistake -- the problem #137's era
summary and #149's cutoff roster each had to solve. So every (arm, clade) pair
is named with ``included`` true or false.

Two reconciliations, both raising rather than reporting:

*Roster against the file tree, in both directions.* ``included`` comes from the
Snakefile function that drove the fan-out, not from a second evaluation of the
sizing arithmetic -- a summary that merely re-derived the verdict could not
detect a fan-out which declined a clade the arithmetic approved, which is the
failure this file exists to expose. Requested-but-absent and
present-but-not-requested both raise.

*The draw against the file it was drawn from.* ``n_samples_arm ==
round(frac * n_samples_full)`` is asserted against the **sample files**, not
against the clade-counts checkpoint the DAG used. The two can disagree: the
checkpoint keys on the metadata clade column's first whitespace-delimited token
while ``samples_by_clade`` matches with ``^clade( |$)``, so a divergence would
silently produce an arm at a slightly different depth than the one recorded.
"""

import os
from importlib import util as importlib_util

import pandas as pd

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

roster = snakemake.params.roster
genes = list(snakemake.params.genes)
min_clade_size = int(snakemake.params.min_clade_size)
outdir = snakemake.params.outdir
output_path = snakemake.output.summary
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None

# Drop reasons, in ``pcp_from_tree.failure_reason``' early-return order. The
# counts arm scores the whole genome and so has only the first and the two
# reversion limits; the PCP arm slices to one gene and adds the rest.
COUNTS_REASONS = [
    "max_nt_mutations",
    "max_reversions_to_ref",
    "max_reversions_to_clade_founder",
]
PCP_REASONS = COUNTS_REASONS + [
    "no_mutations",
    "duplicate_codon_position",
    "identical_slice",
]

COLUMNS = [
    "arm",
    "unit",
    "seed",
    "sampling_fraction",
    "clade",
    "gene",
    "included",
    "n_samples_full",
    "n_samples_arm",
    "sample_fraction_realized",
    "n_nt_mutations_full",
    "n_nt_mutations_arm",
    "nt_mutation_retention",
    "n_branches_scored_full",
    "n_branches_scored_arm",
    "n_branches_kept_full",
    "n_branches_kept_arm",
    "predicted_nt_mutations_arm",
    "predicted_nt_mutation_retention",
    "nt_mutation_deviation_sd",
] + [f"branches_dropped_{r}" for r in PCP_REASONS]


def log(msg):
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


def n_lines(path):
    with open(path) as f:
        return sum(1 for line in f if line.strip())


def read_branch_qc(path):
    """``(outcome_counts, n_scored)`` from one branch-QC sidecar."""
    qc = pd.read_csv(path)
    outcomes = qc[qc["category"] == "outcome"]
    counts = dict(zip(outcomes["key"], outcomes["n_branches"].astype(int)))
    scored = qc.loc[qc["category"] == "total", "n_branches"]
    return counts, int(scored.iloc[0]) if len(scored) else 0


def node_table(outdir, clade):
    """``(n_node, k_node)`` arrays for one clade's branches, at full depth.

    ``n_node`` is the clade samples below the node and ``k_node`` the nucleotide
    mutations on its own branch -- the two inputs the closed form for expected
    retention needs. Both come from files ``count_mutations_clade`` already
    reads: ``leaves_sharing_mutations`` in ``{clade}_mutations.tsv`` *is*
    ``n_node``, and ``{clade}_nt_mutations.csv`` carries the full mutation list
    including noncoding sites (the tsv's own ``nt_mutations`` column is the
    coding subset, with a site repeated once per overlapping gene).

    Deliberately **not** derived from ``{clade}_sample_paths.tsv``, which is one
    row per *tip* rather than per node -- multiple GB on 21J -- where these two
    are one row per node and are already in the DAG.
    """
    tsv = pd.read_csv(
        os.path.join(outdir, "mutations", f"{clade}_mutations.tsv"),
        sep="\t",
        usecols=["node_id", "leaves_sharing_mutations"],
    )
    csv = pd.read_csv(
        os.path.join(outdir, "mutations", f"{clade}_nt_mutations.csv"),
        dtype={"nt_mutations": str},
    )
    merged = tsv.merge(csv, on="node_id", how="inner", validate="one_to_one")
    # ``validate="one_to_one"`` checks that each key is unique on both sides, not
    # that the two key *sets* agree, so an inner join can silently shrink the
    # node list -- and the closed form would then predict retention over fewer
    # branches than ``count_nt_mutations`` measured it over, making the
    # deviation an artifact of the join rather than of the draw.
    if len(merged) != len(csv) or len(merged) != len(tsv):
        raise ValueError(
            f"clade {clade!r}: {len(tsv)} nodes in {clade}_mutations.tsv and "
            f"{len(csv)} in {clade}_nt_mutations.csv join to {len(merged)}. The "
            f"two are written from the same node list, so a disagreement means "
            f"one of them is stale."
        )
    muts = merged["nt_mutations"].fillna("")
    k = (muts.str.count(";").where(muts != "", -1) + 1).astype(int)
    return merged["leaves_sharing_mutations"].astype(int).tolist(), k.tolist()


def count_nt_mutations(path):
    """Total nucleotide mutations across a clade's branches, pre-QC.

    ``{clade}_nt_mutations.csv`` holds one row per unique node with its own
    branch's mutations ``;``-joined, so the branch count is the row count and
    the mutation count is the sum of list lengths. Counted here rather than
    taken from the branch-QC sidecar because the closed form for expected
    retention is over **every** branch, while the sidecar's histogram covers
    only branches that survived parsing into the counting step.
    """
    df = pd.read_csv(path, dtype={"nt_mutations": str})
    if df.empty:
        return 0, 0
    muts = df["nt_mutations"].fillna("")
    per_branch = muts.str.count(";").where(muts != "", -1) + 1
    return int(per_branch.sum()), len(df)


def ratio(num, den):
    return round(num / den, 6) if den else ""


rows = []
problems = []

for arm in sorted(roster):
    spec = roster[arm]
    unit = list(spec["unit"])
    frac = float(spec["frac"])
    seed = int(spec["seed"])
    included = set(spec["included"])
    excluded = set(spec["excluded"])
    unit_label = "-".join(unit)

    if included | excluded != set(unit):
        problems.append(
            f"arm {arm}: roster's included+excluded {sorted(included | excluded)} "
            f"does not partition its unit {sorted(unit)}"
        )

    for clade in unit:
        arm_samples_path = os.path.join(
            outdir, "mat", "clades", arm, f"{clade}_samples.txt"
        )
        full_samples_path = os.path.join(
            outdir, "mat", "clades", f"{clade}_samples.txt"
        )
        n_full = n_lines(full_samples_path)
        is_included = clade in included
        exists = os.path.exists(arm_samples_path)

        # Reconciliation, both directions. An excluded clade whose sample list
        # nonetheless exists is a stale intermediate from an earlier config, and
        # leaving it would let a later run read a subtree the roster says was
        # never built.
        if is_included and not exists:
            problems.append(
                f"arm {arm}, clade {clade}: roster says included, but "
                f"{arm_samples_path} does not exist"
            )
            continue
        if not is_included and exists:
            problems.append(
                f"arm {arm}, clade {clade}: roster says excluded, but "
                f"{arm_samples_path} exists. Stale intermediate from an earlier "
                f"'clade_subsamples'; delete mat/clades/{arm}/ and re-run."
            )

        base = {
            "arm": arm,
            "unit": unit_label,
            "seed": seed,
            "sampling_fraction": frac,
            "clade": clade,
            "included": is_included,
            "n_samples_full": n_full,
        }

        if not is_included:
            expected = subsample.retained_count(n_full, frac)
            log(
                f"{arm}/{clade}: excluded -- round({frac} * {n_full}) = "
                f"{expected}, against min_clade_size {min_clade_size}"
            )
            if expected >= min_clade_size:
                problems.append(
                    f"arm {arm}, clade {clade}: recorded as excluded, but "
                    f"round({frac} * {n_full}) = {expected} >= min_clade_size "
                    f"{min_clade_size}. A clade the arithmetic approves must not "
                    f"be silently declined by the fan-out."
                )
            rows.append({**base, "gene": "", "n_samples_arm": expected})
            continue

        n_arm = n_lines(arm_samples_path)
        expected = subsample.retained_count(n_full, frac)
        if n_arm != expected:
            problems.append(
                f"arm {arm}, clade {clade}: drew {n_arm} samples, expected "
                f"round({frac} * {n_full}) = {expected}"
            )
        if n_arm < min_clade_size:
            problems.append(
                f"arm {arm}, clade {clade}: included with {n_arm} samples, below "
                f"min_clade_size {min_clade_size}"
            )

        muts_arm, branches_arm = count_nt_mutations(
            os.path.join(outdir, "mutations", arm, f"{clade}_nt_mutations.csv")
        )
        muts_full, branches_full = count_nt_mutations(
            os.path.join(outdir, "mutations", f"{clade}_nt_mutations.csv")
        )

        # Check 6: realized retention against the closed form. Reported rather
        # than gated -- the sd is a *scale* for judging a deviation, not a
        # p-value, because sibling survival is negatively correlated under
        # fixed-size sampling and the fixed-size draw makes the mean itself a
        # slight over-prediction. What the number supports is "a same-signed
        # deviation of many sd at two independent seeds indicates a bug".
        n_nodes, k_nodes = node_table(outdir, clade)
        # ``p`` is the *realized* fraction, not the configured one. The draw is
        # fixed-size ``round(frac * n)``, so the effective per-sample retention
        # probability is ``n_arm / n_full`` -- 4/7 rather than 0.5 for a
        # seven-sample clade at 0.5. The two agree to many digits at production
        # clade sizes and differ by 14% at CI sizes, so using the configured
        # value would make the comparison look biased on exactly the data where
        # it is cheap to check.
        predicted, predicted_sd = subsample.expected_mutation_retention(
            n_nodes, k_nodes, n_arm / n_full
        )
        deviation = (
            round((muts_arm - predicted) / predicted_sd, 4)
            if predicted_sd > 0
            else ""
        )

        qc_arm, scored_arm = read_branch_qc(
            os.path.join(outdir, "mutation_counts", arm, "branch_qc", f"{clade}.csv")
        )
        qc_full, scored_full = read_branch_qc(
            os.path.join(outdir, "mutation_counts", "branch_qc", f"{clade}.csv")
        )

        counts_row = {
            **base,
            # Empty on counts-arm rows: ``count_mutations.py`` applies the QC
            # limits once over the whole genome while ``pcp_from_tree.py``
            # applies them to the *gene-sliced* mutation count, so one branch
            # can be dropped for S and kept for ORF1a. No single drop count
            # describes both, which is why this column exists at all.
            "gene": "",
            "n_samples_arm": n_arm,
            "sample_fraction_realized": ratio(n_arm, n_full),
            "n_nt_mutations_full": muts_full,
            "n_nt_mutations_arm": muts_arm,
            "nt_mutation_retention": ratio(muts_arm, muts_full),
            "n_branches_scored_full": scored_full,
            "n_branches_scored_arm": scored_arm,
            "n_branches_kept_full": qc_full.get("kept", ""),
            "n_branches_kept_arm": qc_arm.get("kept", ""),
            "predicted_nt_mutations_arm": round(predicted, 2),
            "predicted_nt_mutation_retention": ratio(predicted, muts_full),
            "nt_mutation_deviation_sd": deviation,
        }
        for reason in COUNTS_REASONS:
            counts_row[f"branches_dropped_{reason}"] = qc_arm.get(reason, "")
        rows.append(counts_row)
        log(
            f"{arm}/{clade}: {n_arm}/{n_full} samples "
            f"({ratio(n_arm, n_full)}), {muts_arm}/{muts_full} nt mutations "
            f"({ratio(muts_arm, muts_full)}); closed form predicts "
            f"{predicted:.1f} +/- {predicted_sd:.1f} "
            f"(deviation {deviation} sd)"
        )

        for gene in genes:
            g_arm, g_scored_arm = read_branch_qc(
                os.path.join(outdir, "pcps", gene, arm, "branch_qc", f"{clade}.csv")
            )
            g_full, g_scored_full = read_branch_qc(
                os.path.join(outdir, "pcps", gene, "branch_qc", f"{clade}.csv")
            )
            gene_row = {
                **base,
                "gene": gene,
                "n_samples_arm": n_arm,
                "sample_fraction_realized": ratio(n_arm, n_full),
                "n_branches_scored_full": g_scored_full,
                "n_branches_scored_arm": g_scored_arm,
                "n_branches_kept_full": g_full.get("kept", ""),
                "n_branches_kept_arm": g_arm.get("kept", ""),
            }
            for reason in PCP_REASONS:
                gene_row[f"branches_dropped_{reason}"] = g_arm.get(reason, "")
            rows.append(gene_row)

if problems:
    raise ValueError(
        "subsample_summary found "
        + str(len(problems))
        + " inconsistency(ies) between the arm roster and the output tree:\n  "
        + "\n  ".join(problems)
    )

summary = pd.DataFrame(rows, columns=COLUMNS)
summary.to_csv(output_path, index=False)

n_arms = summary["arm"].nunique()
n_included = int(summary.loc[summary["gene"] == "", "included"].sum())
n_excluded = int((~summary.loc[summary["gene"] == "", "included"]).sum())
log(
    f"Wrote {len(summary)} rows for {n_arms} arm(s): {n_included} included and "
    f"{n_excluded} excluded (arm, clade) pair(s) -> {output_path}"
)
