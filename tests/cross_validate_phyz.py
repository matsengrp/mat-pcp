#!/usr/bin/env python3
"""Cross-validate phyz pcp against mat-pcp PCP extraction.

Runs both pipelines on the same MAT subtrees and checks that every edge
produced by mat-pcp is also produced by phyz with identical parent/child
sequences. phyz is expected to be a superset because it does not apply
two mat-pcp filters:
  1. Empty branches: mat-pcp drops branches with zero mutations.
  2. Whole-genome coding count: mat-pcp's max_nt_mutations applies to
     coding mutations anywhere in the genome; phyz's --max-mutations
     applies only within the configured gene region.
  3. Residual (ambiguous codon or reversions-to-clade-founder, not
     cross-checked directly here).

Usage:
    PHYZ=/path/to/phyz MATUTILS=/path/to/matUtils \\
        pixi run python tests/cross_validate_phyz.py

Requires that the test pipeline has been run:
    pixi run snakemake -j 4 --use-conda \\
        --configfile config/config_test_data.yaml
"""

import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path
from collections import Counter


def run_phyz(phyz_bin, tree_pb, refseq, region, max_mut, max_rev, out_tsv):
    cmd = [
        phyz_bin, "pcp",
        "--tree-pb", str(tree_pb),
        "--refseq", str(refseq),
        "--gene-region", region,
        "--max-mutations", str(max_mut),
        "--max-reversions-to-ref", str(max_rev),
        "--output", str(out_tsv),
        "--quiet",
    ]
    subprocess.run(cmd, check=True)


def run_matutils_all_paths(matutils_bin, tree_pb, out_dir, out_name):
    """Dump {node_name: [mutation_strings]} for every node in the tree.

    Caveat: phyz renames unlabeled internal nodes with its own preorder
    indexing, so phyz's "node_N" is not guaranteed to match BTE's "node_N".
    This dict indexes by the BTE names that matUtils emits.
    """
    subprocess.run(
        [matutils_bin, "extract", "-i", str(tree_pb),
         "-d", str(out_dir), "-A", out_name],
        check=True, capture_output=True,
    )
    branch_muts = {}
    with open(Path(out_dir) / out_name) as fh:
        for line in fh:
            parts = line.rstrip("\n").split(": ", 1)
            if len(parts) != 2:
                continue
            name, muts = parts
            branch_muts[name] = [m for m in muts.split(",") if m]
    return branch_muts


def load_coding_sites(csv_path):
    sites = set()
    with open(csv_path) as fh:
        for r in csv.DictReader(fh):
            if r["gene"] != "noncoding":
                sites.add(int(r["site"]))
    return sites


def load_matpcp(csv_path):
    """Return list of (parent_seq, child_seq, parent_name, child_name, branch_length)."""
    rows = []
    with open(csv_path) as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            rows.append((r["parent"], r["child"], r["parent_name"],
                         r["child_name"], int(r["branch_length"])))
    return rows


def load_phyz(tsv_path):
    """Return list of (parent_seq, child_seq, parent_name, child_name, n_mutations, is_terminal)."""
    rows = []
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            rows.append((r["parent"], r["child"], r["parent_name"],
                         r["child_name"], int(r["n_mutations"]),
                         r["child_is_leaf"] == "true"))
    return rows


def compare_clade(clade, matpcp_csv, phyz_tsv, branch_muts=None,
                  coding_sites=None, max_coding=4):
    """Compare using multiset semantics on (parent_seq, child_seq) pairs.

    Multiple rows with the same (parent_seq, child_seq) are common because
    phyz expands condensed nodes inline (each collapsed leaf gets its own
    row with identical reconstructed sequences). We compare multisets so a
    phyz row without a matching mat-pcp row shows up as phyz-only.
    """
    mp = load_matpcp(matpcp_csv)
    pz = load_phyz(phyz_tsv)

    mp_multi = Counter((p, c) for (p, c, *_) in mp)
    pz_multi = Counter((p, c) for (p, c, *_) in pz)

    all_keys = set(mp_multi) | set(pz_multi)
    n_agreed_rows = 0
    mp_only_rows = 0
    pz_only_rows = 0
    mp_only_examples = []
    for k in all_keys:
        m, p = mp_multi.get(k, 0), pz_multi.get(k, 0)
        agreed = min(m, p)
        n_agreed_rows += agreed
        if m > p:
            mp_only_rows += m - p
            mp_only_examples.append((k, m, p))
        elif p > m:
            pz_only_rows += p - m

    # Classify the phyz rows that have no matching mat-pcp row.
    phyz_only_zero_mut = 0
    phyz_only_exceeds_whole_genome_coding = 0
    phyz_only_other = []
    mp_remaining = Counter((p, c) for (p, c, *_) in mp)
    for row in pz:
        key = (row[0], row[1])
        if mp_remaining[key] > 0:
            mp_remaining[key] -= 1
            continue
        n_mut = row[4]
        cname = row[3]
        if n_mut == 0:
            phyz_only_zero_mut += 1
            continue
        if branch_muts is not None and coding_sites is not None:
            whole_genome_muts = branch_muts.get(cname, [])
            coding_count = sum(1 for m in whole_genome_muts
                               if int(m[1:-1]) in coding_sites)
            if coding_count > max_coding:
                phyz_only_exceeds_whole_genome_coding += 1
                continue
        phyz_only_other.append((key, n_mut, row[2], row[3]))

    return {
        "clade": clade,
        "n_matpcp": len(mp),
        "n_phyz": len(pz),
        "n_agreed": n_agreed_rows,
        "n_matpcp_only": mp_only_rows,
        "n_phyz_only": pz_only_rows,
        "phyz_only_zero_mut": phyz_only_zero_mut,
        "phyz_only_exceeds_whole_genome_coding":
            phyz_only_exceeds_whole_genome_coding,
        "phyz_only_other": phyz_only_other,
        "matpcp_only_examples": mp_only_examples[:5],
    }


def default_results_dir():
    """Versioned test output dir, tracking ``VERSION`` (#88).

    A literal ``results/test_data_v1`` here went eight bumps stale without
    anything noticing, because a missing directory reads as "the pipeline has
    not been run" rather than as a wrong path.
    """
    try:
        with open("VERSION") as vf:
            return f"results/test_data_{vf.read().strip()}"
    except OSError:
        return "results/test_data"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phyz-bin", default=os.environ.get("PHYZ", "./phyz"))
    ap.add_argument("--results-dir", default=default_results_dir())
    ap.add_argument("--clades", nargs="+", default=["19A", "19B", "20A"])
    ap.add_argument("--region", default="21563-25381")
    ap.add_argument("--max-mutations", type=int, default=4)
    ap.add_argument("--max-reversions", type=int, default=1)
    ap.add_argument("--phyz-out-dir", default="/tmp/phyz_cross_validate")
    ap.add_argument("--matutils-bin", default=os.environ.get("MATUTILS"),
                    help="Optional matUtils path. If set, the script dumps "
                         "branch mutations and classifies phyz-only rows by "
                         "whole-genome coding-mutation count.")
    args = ap.parse_args()

    results_dir = Path(args.results_dir)
    refseq = results_dir / "ref" / "ref.fa"
    coding_sites_csv = results_dir / "ref" / "coding_sites.csv"
    phyz_out_dir = Path(args.phyz_out_dir)
    phyz_out_dir.mkdir(exist_ok=True)

    coding_sites = None
    if args.matutils_bin and coding_sites_csv.exists():
        coding_sites = load_coding_sites(coding_sites_csv)

    summaries = []
    for clade in args.clades:
        tree_pb = results_dir / "mat" / "clades" / f"{clade}_tree.pb"
        # NOTE: this path is stale independently of the extension -- since #92
        # PCPs live under a per-gene subdirectory (``pcps/{gene}/``), so this
        # helper does not currently resolve against a real run. Extension
        # updated for #134; the missing gene component is pre-existing.
        matpcp_csv = results_dir / "pcps" / f"{clade}_pcps.csv.gz"
        phyz_tsv = phyz_out_dir / f"{clade}_phyz_pcps.tsv"

        if not tree_pb.exists():
            print(f"SKIP {clade}: {tree_pb} missing", file=sys.stderr)
            continue

        print(f"Running phyz on clade {clade}...")
        run_phyz(args.phyz_bin, tree_pb, refseq, args.region,
                 args.max_mutations, args.max_reversions, phyz_tsv)

        branch_muts = None
        if args.matutils_bin:
            branch_muts = run_matutils_all_paths(
                args.matutils_bin, tree_pb, phyz_out_dir,
                f"{clade}_all_paths.txt",
            )

        summary = compare_clade(
            clade, matpcp_csv, phyz_tsv,
            branch_muts=branch_muts, coding_sites=coding_sites,
            max_coding=args.max_mutations,
        )
        summaries.append(summary)

    # Report
    print()
    header = (
        f"{'Clade':<6} {'mat-pcp':>8} {'phyz':>6} {'agreed':>7} "
        f"{'mp-only':>8} {'pz-only':>8} {'0mut':>5} {'>4cod':>6} {'other':>6}"
    )
    print("=" * len(header))
    print(header)
    print("-" * len(header))
    for s in summaries:
        print(
            f"{s['clade']:<6} {s['n_matpcp']:>8} {s['n_phyz']:>6} "
            f"{s['n_agreed']:>7} {s['n_matpcp_only']:>8} "
            f"{s['n_phyz_only']:>8} {s['phyz_only_zero_mut']:>5} "
            f"{s['phyz_only_exceeds_whole_genome_coding']:>6} "
            f"{len(s['phyz_only_other']):>6}"
        )
    print("=" * len(header))

    total_mp_only = sum(s["n_matpcp_only"] for s in summaries)
    total_pz_other = sum(len(s["phyz_only_other"]) for s in summaries)

    print(f"\nmat-pcp-only rows (phyz missing): {total_mp_only}")
    print(f"phyz-only rows (0 mutations — mat-pcp rejects as empty branch): "
          f"{sum(s['phyz_only_zero_mut'] for s in summaries)}")
    print(f"phyz-only rows (>4 whole-genome coding muts — mat-pcp's coding "
          f"count filter acts on the whole genome, phyz's on the gene "
          f"region): {sum(s['phyz_only_exceeds_whole_genome_coding'] for s in summaries)}")
    print(f"phyz-only rows (other, need investigation): {total_pz_other}")

    if total_pz_other > 0:
        print("\nphyz-only 'other' rows (possibly ambiguous codon or "
              "reversions-to-founder in mat-pcp):")
        for s in summaries:
            for (key, n_mut, pname, cname) in s["phyz_only_other"][:10]:
                print(f"  [{s['clade']}] {pname} -> {cname} ({n_mut} muts)")

    if total_mp_only > 0:
        print("\n*** FAILURE: mat-pcp produced rows phyz did not. Examples:")
        for s in summaries:
            for (key, mp_count, pz_count) in s["matpcp_only_examples"]:
                print(f"  [{s['clade']}] mp_count={mp_count} pz_count={pz_count} "
                      f"parent[:30]={key[0][:30]}... child[:30]={key[1][:30]}...")
        sys.exit(1)

    print("\nSUCCESS: every mat-pcp edge is present in phyz output with "
          "identical sequences.")


if __name__ == "__main__":
    main()
