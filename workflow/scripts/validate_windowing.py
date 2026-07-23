#!/usr/bin/env python3
"""Validate the windowing module + its test outputs (EPIC #115 Phase 2, #117).

Runs in the **base env** (no ``bte``): imports only the pure helpers from
``workflow/modules/windowing.py`` (which keep ``import bte`` lazy), asserts their
known-answer / boundary behavior, then -- when the test outputs exist -- runs
structural checks on ``windows/subtree_counts.csv`` + ``windows/founder_muts.csv``.

Invoked like ``validate_test.py`` (``pixi run python workflow/scripts/validate_windowing.py``);
the tree-level algorithmic checks live in the separate ``validate_windowing`` rule
(bte env). Exits non-zero on any failure.
"""

import glob
import os
import sys
from datetime import timedelta
from importlib import util as importlib_util

import pandas as pd
import yaml

# --- by-path module loader ---------------------------------------------------
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

TEST_CONFIG = "config/config_test_data.yaml"
SUBTREE_COUNTS_COLUMNS = [
    "window_start",
    "founder_id",
    "protein",
    "nt_site",
    "nt_mutation",
    "aa_mutation",
    "codon_change",
    "synonymous",
    "count",
]

failures = []


def check(name, ok, detail=""):
    status = "PASS" if ok else "FAIL"
    print(f"[{status}] {name}" + (f": {detail}" if detail else ""))
    if not ok:
        failures.append(f"{name}: {detail}")


# =============================================================================
# (A) Pure-function known-answer / boundary checks
# =============================================================================

def check_apply_muts():
    check(
        "apply_muts known-answer",
        windowing.apply_muts("ACGT", ["A1G", "T4A"]) == "GCGA",
        windowing.apply_muts("ACGT", ["A1G", "T4A"]),
    )
    # A from-base mismatch warns but still applies (C1G on an 'A').
    check(
        "apply_muts mismatch still applies",
        windowing.apply_muts("ACGT", ["C1G"]) == "GCGT",
    )


def check_get_codon():
    # site 5 (1-based) of AAATTTGGG; the codon changed TTT->TCT (offset 1), so
    # the codon around site 5 is positions 4-6 = "TTT".
    codon = windowing.get_codon("AAATTTGGG", "TTT", "TCT", 5)
    check("get_codon known-answer", codon == "TTT", codon)


# Enough codons that the criterion-1 test (which places max_nt+1 mutations at
# sites 1, 4, 7, ... = 3*max_nt + 1) stays within the fixture even if the test
# config's max_nt_mutations is raised well above its current value of 4.
_N_CODONS = 40


def _qc_map():
    """coding_sites + site_to_codon_keys for the qc tests.

    Sites 1..(3*_N_CODONS) = gene G codons 1.._N_CODONS (3 nt each); sites 50,51
    = an overlap of genes G/Gb both in one codon; sites 1000,1001 (well beyond
    the coding range) = noncoding (absent from the map).
    """
    n_sites = 3 * _N_CODONS
    coding_sites = set(range(1, n_sites + 1)) | {50, 51}
    site_to_codon_keys = {
        s: [("G", str((s - 1) // 3 + 1))] for s in range(1, n_sites + 1)
    }
    site_to_codon_keys[50] = [("G", "1"), ("Gb", "10")]
    site_to_codon_keys[51] = [("G", "1"), ("Gb", "10")]
    return coding_sites, site_to_codon_keys


def _qc(nt_mutations, ref, founder, thresholds):
    coding_sites, site_to_codon_keys = _qc_map()
    return windowing.qc_filter(
        nt_mutations,
        founder,
        ref,
        coding_sites,
        site_to_codon_keys,
        *thresholds,
    )


def check_qc_filter(thresholds):
    max_nt, max_ref, max_founder = thresholds
    N = 200
    ref = "C" * N       # reversion-to-ref  == mutation ending in C
    founder = "A" * N   # reversion-to-founder == mutation ending in A

    # Criterion 1: n coding mutations (keep-if <=). Muts end in G (not C, not A)
    # so they are neither ref- nor founder-reversions; distinct codons.
    distinct = [f"T{1 + 3 * k}G" for k in range(max_nt + 1)]  # max_nt+1 muts
    check(
        "qc criterion 1 reject (> max_nt_mutations)",
        _qc(distinct, ref, founder, thresholds) is False,
    )
    check(
        "qc criterion 1 keep (== max_nt_mutations)",
        _qc(distinct[:max_nt], ref, founder, thresholds) is True,
    )

    # Criterion 2: reversions to reference (end in C).
    ref_rev = [f"T{1 + 3 * k}C" for k in range(max_ref + 1)]
    check(
        "qc criterion 2 reject (> max_reversions_to_ref)",
        _qc(ref_rev, ref, founder, thresholds) is False,
    )
    check(
        "qc criterion 2 keep (== max_reversions_to_ref)",
        _qc(ref_rev[:max_ref], ref, founder, thresholds) is True,
    )

    # Criterion 3: reversions to founder (end in A).
    founder_rev = [f"T{1 + 3 * k}A" for k in range(max_founder + 1)]
    check(
        "qc criterion 3 reject (> max_reversions_to_clade_founder)",
        _qc(founder_rev, ref, founder, thresholds) is False,
    )
    check(
        "qc criterion 3 keep (== max_reversions_to_clade_founder)",
        _qc(founder_rev[:max_founder], ref, founder, thresholds) is True,
    )

    # Noncoding excluded from the count: 3 coding + 2 noncoding muts, max_nt>=3.
    mixed = ["T1G", "T4G", "T7G", "T1000G", "T1001G"]
    check(
        "qc noncoding excluded from count",
        _qc(mixed, ref, founder, thresholds) is True,
    )

    # Criterion 4 (unconditional): 2 coding muts in the same codon -> reject;
    # in different codons -> keep.
    check(
        "qc criterion 4 reject (2 muts, same codon)",
        _qc(["T1G", "T2G"], ref, founder, thresholds) is False,
    )
    check(
        "qc criterion 4 keep (2 muts, different codons)",
        _qc(["T1G", "T4G"], ref, founder, thresholds) is True,
    )
    # Overlap-gene site: both genes' codon groups flag together.
    check(
        "qc criterion 4 reject (overlap-gene same codon)",
        _qc(["T50G", "T51G"], ref, founder, thresholds) is False,
    )


def check_collapse_multigene():
    records = [
        {
            "window_start": "2021-01-01",
            "founder_id": "n1",
            "protein": "ORF1a",
            "nt_site": 6080,
            "nt_mutation": "G6080T",
            "aa_mutation": "D1939Y",
            "codon_change": "GAT>TAT",
            "synonymous": False,
            "count": 3,
        },
        {
            "window_start": "2021-01-01",
            "founder_id": "n1",
            "protein": "ORF1ab",
            "nt_site": 6080,
            "nt_mutation": "G6080T",
            "aa_mutation": "D1939Y",
            "codon_change": "GAT>TAT",
            "synonymous": False,
            "count": 3,
        },
    ]
    out = windowing.collapse_multigene(records)
    check("collapse yields one row", len(out) == 1, f"{len(out)} rows")
    if len(out) == 1:
        row = out.iloc[0]
        check("collapse protein joined in order", row["protein"] == "ORF1a;ORF1ab", row["protein"])
        check("collapse aa_mutation joined", row["aa_mutation"] == "D1939Y;D1939Y")
        check("collapse codon_change joined", row["codon_change"] == "GAT>TAT;GAT>TAT")
        check("collapse count is shared (not summed)", int(row["count"]) == 3, row["count"])
        check("collapse columns", list(out.columns) == SUBTREE_COUNTS_COLUMNS)

    # synonymous = .all(): a True/False pair collapses to False.
    mixed = [dict(records[0], synonymous=True), dict(records[1], synonymous=False)]
    out2 = windowing.collapse_multigene(mixed)
    check("collapse synonymous is all()", bool(out2.iloc[0]["synonymous"]) is False)

    # Empty input -> empty frame with the right columns.
    empty = windowing.collapse_multigene([])
    check("collapse empty schema", list(empty.columns) == SUBTREE_COUNTS_COLUMNS)


# =============================================================================
# (B) Output-CSV structural checks (only when the test outputs exist)
# =============================================================================

def default_output_dir():
    """Versioned test output dir (tracks VERSION, like validate_test.py)."""
    try:
        with open("VERSION") as vf:
            versioned = f"results/test_data_{vf.read().strip()}"
        if os.path.isdir(versioned):
            return versioned
    except OSError:
        pass
    candidates = sorted(glob.glob("results/test_data*"))
    return candidates[-1] if candidates else "results/test_data"


def window_grid_dates(node_dates_path, window_size_days, window_step_days):
    """Formatted (YYYY-MM-DD) window-start grid, anchored at max(node_dates).

    Mirrors WindowedMAT._window_starts backward-anchoring without needing the
    tree: the floor uses min(node_dates) (<= the root's start_time), so the grid
    is a superset of any real window start.
    """
    times = windowing.get_time_dict(node_dates_path)
    end_time = max(times.values())
    floor = min(times.values()) - timedelta(days=window_size_days)
    window_size = timedelta(days=window_size_days)
    step = timedelta(days=window_step_days)
    grid = set()
    i = 0
    while True:
        w = end_time - window_size - i * step
        if w + window_size <= floor:
            break
        grid.add(w.strftime("%Y-%m-%d"))
        i += 1
        if i > 100000:  # safety
            break
    return grid


def check_outputs(output_dir, cfg):
    counts_path = os.path.join(output_dir, "windows/subtree_counts.csv")
    founder_path = os.path.join(output_dir, "windows/founder_muts.csv")
    node_dates_path = os.path.join(output_dir, "mat/node_dates.tsv")
    ref_path = os.path.join(output_dir, "ref/ref.fa")

    if not (os.path.exists(counts_path) and os.path.exists(founder_path)):
        print(f"[SKIP] output CSVs not found under {output_dir} (build them first)")
        return

    counts = pd.read_csv(counts_path)
    founders = pd.read_csv(founder_path)

    check("subtree_counts columns", list(counts.columns) == SUBTREE_COUNTS_COLUMNS,
          list(counts.columns))
    check("subtree_counts non-empty", len(counts) >= 1, f"{len(counts)} rows")
    check("founder_muts non-empty", len(founders) >= 1, f"{len(founders)} rows")

    if len(counts):
        check("count >= 1", int(counts["count"].min()) >= 1, str(counts["count"].min()))
        check("nt_site >= 1", int(counts["nt_site"].min()) >= 1, str(counts["nt_site"].min()))
        if os.path.exists(ref_path):
            ref_len = len(windowing.load_reference_from_fasta(ref_path))
            check("nt_site <= reference length",
                  int(counts["nt_site"].max()) <= ref_len,
                  f"max {counts['nt_site'].max()} vs len {ref_len}")

        # Every founder_id in the counts appears in founder_muts.
        missing = set(counts["founder_id"]) - set(founders["founder_id"])
        check("every count founder_id in founder_muts", not missing,
              f"{len(missing)} missing")

        # Every window_start is on the computed backward-anchored grid.
        if os.path.exists(node_dates_path):
            grid = window_grid_dates(
                node_dates_path, cfg["window_size"], cfg["window_step"]
            )
            off_grid = set(counts["window_start"].astype(str)) - grid
            check("every window_start on the window grid", not off_grid,
                  f"off-grid: {sorted(off_grid)[:5]}")


def main():
    with open(TEST_CONFIG) as fh:
        cfg = yaml.safe_load(fh)
    thresholds = (
        cfg["max_nt_mutations"],
        cfg["max_reversions_to_ref"],
        cfg["max_reversions_to_clade_founder"],
    )

    print("=== (A) pure-function checks ===")
    check_apply_muts()
    check_get_codon()
    check_qc_filter(thresholds)
    check_collapse_multigene()

    print("\n=== (B) output-CSV structural checks ===")
    output_dir = sys.argv[1] if len(sys.argv) > 1 else default_output_dir()
    check_outputs(output_dir, cfg)

    print()
    if failures:
        print(f"FAILED ({len(failures)}):")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    print("All windowing validation checks passed.")


if __name__ == "__main__":
    main()
