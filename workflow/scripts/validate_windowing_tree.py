#!/usr/bin/env python3
"""Tree-level algorithmic checks for the windowing module (EPIC #115 #117).

Runs in the **bte env** (unlike the pure/CSV checks in ``validate_windowing.py``).
Builds a real ``WindowedMAT`` on the test ``.pb`` and, for controlled-topology
cases, overwrites ``mat.times`` with a **synthetic** dict over the real node ids
(the algorithm takes ``times`` as a dict, so no hand-built ``.pb`` is needed).

Wired as a CI-only Snakemake rule (``validate_windowing``) on the test config; it
is a correctness harness, not a per-tree production gate, so it is NOT in
``results_files``. Writes a report sentinel and fails the run on any violation.

Checks:
  1. Window grid: backward-anchored, all windows ``window_size`` wide, newest
     ends exactly at ``end_time``, consecutive starts ``step_size`` apart.
  2. Partition (``step == size``): kept subtrees across windows are node-disjoint.
  3. Overlap multiplicity (``step == size // 2``): >=1 node lies in two windows.
  4. Step-vs-size asymmetry (synthetic times): founder->A->B with
     ``date(A) in (t+step, t+size]`` and ``date(B) > t+size`` -> A is in the
     founder's subtree, B is masked out, A is not independently a founder.
"""

import os
from datetime import datetime, timedelta
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
windowing = workflow_loader.load_workflow_module("windowing")

tree_path = snakemake.input.tree
node_dates_path = snakemake.input.node_dates
gtf_path = snakemake.input.gtf
fasta_path = snakemake.input.fasta
coding_sites_csv = snakemake.input.coding_sites
report_path = snakemake.output.report
window_size_days = int(snakemake.params.window_size)
window_step_days = int(snakemake.params.window_step)
max_nt_mutations = snakemake.params.max_nt_mutations
max_reversions_to_ref = snakemake.params.max_reversions_to_ref
max_reversions_to_founder = snakemake.params.max_reversions_to_clade_founder

results = []
failures = []


def check(name, ok, detail=""):
    status = "PASS" if ok else "FAIL"
    line = f"[{status}] {name}" + (f": {detail}" if detail else "")
    print(line)
    results.append(line)
    if not ok:
        failures.append(f"{name}: {detail}")


def build_mat():
    return windowing.WindowedMAT(
        tree_path,
        node_dates_path,
        gtf_path,
        fasta_path,
        coding_sites_csv,
        max_nt_mutations=max_nt_mutations,
        max_reversions_to_ref=max_reversions_to_ref,
        max_reversions_to_founder=max_reversions_to_founder,
    )


def subtree_node_ids(mat, founder_id, window_start, window_size):
    """Re-derive a founder's kept-subtree node-id set (applies the size mask)."""
    founder_node = mat.ids_to_nodes[founder_id]
    return {
        n.id
        for n in windowing.traverse_mat(
            founder_node,
            mask=lambda n: mat.times[n.id] > window_start + window_size,
        )
    }


def founder_ids_per_window(subtree_data):
    return [
        (window_start, [rec[0][0] for rec in subtrees])
        for window_start, subtrees in subtree_data
    ]


def check_grid(mat, window_size, step_size):
    starts = mat._window_starts(window_size, step_size)
    check("grid non-empty", len(starts) >= 1, f"{len(starts)} windows")
    if not starts:
        return
    check(
        "grid newest window ends at end_time",
        starts[-1] + window_size == mat.end_time,
        f"{starts[-1] + window_size} vs {mat.end_time}",
    )
    diffs = {(b - a) for a, b in zip(starts[:-1], starts[1:])}
    check(
        "grid consecutive starts are step_size apart",
        diffs == {step_size} or len(starts) == 1,
        f"diffs={diffs}",
    )
    check(
        "grid sorted ascending",
        starts == sorted(starts),
    )


def check_partition(mat, window_size):
    # step == size: kept subtrees across windows must be node-disjoint.
    subtree_data = mat.aggregate_data(window_size, window_size)
    seen = {}
    overlaps = 0
    for window_start, founder_ids in founder_ids_per_window(subtree_data):
        window_nodes = set()
        for fid in founder_ids:
            window_nodes |= subtree_node_ids(mat, fid, window_start, window_size)
        for nid in window_nodes:
            if nid in seen and seen[nid] != window_start:
                overlaps += 1
            seen[nid] = window_start
    check(
        "partition: no node in two windows (step == size)",
        overlaps == 0,
        f"{overlaps} cross-window node reuses",
    )


def check_overlap(mat, window_size, step_size):
    subtree_data = mat.aggregate_data(window_size, step_size)
    windows_by_node = {}
    for window_start, founder_ids in founder_ids_per_window(subtree_data):
        for fid in founder_ids:
            for nid in subtree_node_ids(mat, fid, window_start, window_size):
                windows_by_node.setdefault(nid, set()).add(window_start)
    multi = sum(1 for w in windows_by_node.values() if len(w) >= 2)
    check(
        "overlap: >=1 node in two windows (step == size // 2)",
        multi >= 1,
        f"{multi} nodes span >=2 windows",
    )


def check_asymmetry(mat):
    """Synthetic-times case isolating the step/size masking asymmetry."""
    # Find a real chain founder F -> A -> B (A child of F, B child of A).
    chain = None
    for f in mat.nodes:
        if not f.children:
            continue
        for a in f.children:
            if a.children:
                b = a.children[0]
                chain = (f.id, a.id, b.id)
                break
        if chain:
            break
    if chain is None:
        check("asymmetry: found a founder->A->B chain", False, "no depth-2 chain")
        return
    f_id, a_id, b_id = chain

    t = datetime(2021, 1, 1)
    size = timedelta(days=30)
    step = timedelta(days=10)
    # Everything defaults to well before t; set the chain explicitly.
    synthetic = {nid: t - timedelta(days=100) for nid in mat.ids_to_nodes}
    synthetic[f_id] = t                       # founder: in window, parent < t
    parent_f = mat.parent_of.get(f_id)
    if parent_f is not None:
        synthetic[parent_f] = t - timedelta(days=10)
    synthetic[a_id] = t + step + timedelta(days=1)   # in (t+step, t+size]
    synthetic[b_id] = t + size + timedelta(days=1)   # > t+size
    mat.times = synthetic

    def leaf_func(node):
        parent_id = mat.parent_of.get(node.id)
        return mat.times[node.id] >= t and (
            parent_id is None or mat.times[parent_id] < t
        )

    def mask_func(node):
        return mat.times[node.id] > t + step

    founders = {n.id for n in windowing.iter_leaves(mat.root_node, leaf_func, mask_func)}
    subtree = subtree_node_ids(mat, f_id, t, size)

    check("asymmetry: F is a founder of window t", f_id in founders)
    check("asymmetry: A is NOT an independent founder (founder-band mask)",
          a_id not in founders)
    check("asymmetry: A is in F's subtree (extent uses window_size)",
          a_id in subtree)
    check("asymmetry: B is masked out of F's subtree", b_id not in subtree)


def main():
    window_size = timedelta(days=window_size_days)
    step_size = timedelta(days=window_step_days)

    print("=== window grid ===")
    mat = build_mat()
    check_grid(mat, window_size, step_size)

    print("\n=== partition (step == size) ===")
    check_partition(mat, window_size)

    print("\n=== overlap multiplicity (step == size // 2) ===")
    half = timedelta(days=max(window_size_days // 2, 1))
    check_overlap(mat, window_size, half)

    print("\n=== step-vs-size asymmetry (synthetic times) ===")
    # Rebuild so the synthetic-times overwrite doesn't affect other checks.
    mat2 = build_mat()
    check_asymmetry(mat2)

    with open(report_path, "w") as fh:
        fh.write("\n".join(results) + "\n")

    if failures:
        raise SystemExit(
            "validate_windowing failed:\n  " + "\n  ".join(failures)
        )
    print("\nAll tree-level windowing checks passed.")


if __name__ == "__main__":
    main()
