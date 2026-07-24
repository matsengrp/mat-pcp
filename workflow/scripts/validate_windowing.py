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

# Regression sentinels for the config_test_data windowing outputs. These lock
# Phase 2's counts down while later EPIC #115 phases churn windowing.py, so an
# unintended counting regression surfaces immediately. chronumental is
# deterministic run-to-run (PRNGKey(0); verified byte-identical node dates across
# reruns), so these actual (integer) counts are exact -- matching how
# validate_test.py treats clade `actual_count` (rel_tol=1e-6). Phase 3 will add
# windowed expected-count / fitness sums here; those are floats and get the
# looser rel_tol=0.01, like the clade `delta_fitness_*_sum`.
#
# UPDATE these deliberately when the output intentionally changes (a
# windowing-logic change, OR a chronumental/jax version bump that shifts a node
# across a window boundary). A mismatch with no such intended change is a bug.
WINDOW_COUNT_SENTINELS = {
    "n_windows": 19,
    "n_founder_rows": 370,
    "n_count_rows": 4041,
    "total_count": 4150,
}
COUNT_REL_TOL = 1e-6

# ---- Phase 3 (#118): windowed expected counts + delta_fitness_ratio ---------
FITNESS_TRAJECTORY_COLUMNS = [
    "window_start",
    "ID",
    "gene",
    "aa_site",
    "aa_mutation",
    "wt_aa",
    "mut_aa",
    "synonymous",
    "actual_count",
    "expected_count",
    "delta_fitness_ratio",
]

# The headline correctness proof: per-window synonymous calibration. Expected
# counts are scaled by the evol_time factor T so that, per window, total expected
# synonymous == total observed synonymous BY CONSTRUCTION (rates.py evol_time
# strategy, same as the clade path). So sum(actual_syn)/sum(expected_syn) must be
# ~1 (exactly, up to float). A departure means the T normalization or the
# synonymous bookkeeping is wrong.
SYN_CALIBRATION_REL_TOL = 1e-4
# Secondary sanity check: synonymous delta_fitness_ratio should also center near 0.
SYN_CENTER_ABS_TOL_MEDIAN = 0.1
SYN_CENTER_ABS_TOL_MEAN = 0.3

# Regression sentinels for the windowed fitness trajectory on config_test_data.
# Integer-valued counts are exact (COUNT_REL_TOL); expected_count is a float sum
# and gets the looser FITNESS_FLOAT_REL_TOL, matching how validate_test.py treats
# the clade delta_fitness_*_sum. Update deliberately on an intended change.
# n_windows is 18 (not 19): window 2020-05-23 has no observed synonymous mutations,
# so its evol_time T is undefined and it is skipped (logged), unlike the actual-only
# subtree_counts.csv which keeps all 19 windows.
WINDOW_FITNESS_SENTINELS = {
    "n_traj_rows": 1278390,
    "n_windows": 18,
    "actual_sum": 4165.0,
    "expected_sum": 5631.6126,
}
FITNESS_FLOAT_REL_TOL = 0.01

failures = []


def is_close(actual, expected, rel_tol=COUNT_REL_TOL, abs_tol=1e-9):
    """Float-tolerant equality, mirroring validate_test.py."""
    if expected == 0:
        return abs(actual) <= abs_tol
    return abs(actual - expected) <= max(rel_tol * abs(expected), abs_tol)


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

        # Regression sentinels -- test data only (the values are specific to
        # config_test_data; a production tree would legitimately differ).
        if "test_data" in os.path.basename(os.path.normpath(output_dir)):
            observed = {
                "n_windows": int(counts["window_start"].nunique()),
                "n_founder_rows": int(len(founders)),
                "n_count_rows": int(len(counts)),
                "total_count": int(counts["count"].sum()),
            }
            for key, expected in WINDOW_COUNT_SENTINELS.items():
                check(
                    f"count sentinel: {key}",
                    is_close(observed[key], expected),
                    f"got {observed[key]}, expected {expected} "
                    f"(update WINDOW_COUNT_SENTINELS if intended)",
                )
        else:
            print(f"[SKIP] count sentinels (non-test output dir {output_dir})")


# =============================================================================
# (C) Windowed expected counts + delta_fitness_ratio (#118)
# =============================================================================

def check_expected_outputs(output_dir, cfg):
    import numpy as np

    traj_path = os.path.join(output_dir, "windows/fitness_trajectory.csv")
    if not os.path.exists(traj_path):
        print(f"[SKIP] fitness_trajectory.csv not found under {output_dir}")
        return

    traj = pd.read_csv(traj_path)
    check(
        "fitness_trajectory columns",
        list(traj.columns) == FITNESS_TRAJECTORY_COLUMNS,
        list(traj.columns),
    )
    check("fitness_trajectory non-empty", len(traj) >= 1, f"{len(traj)} rows")
    if not len(traj):
        return

    # ORF1a is dropped in favor of ORF1ab (aamutfit convention).
    check(
        "fitness_trajectory has no ORF1a rows",
        int((traj["gene"] == "ORF1a").sum()) == 0,
        f"{int((traj['gene'] == 'ORF1a').sum())} ORF1a rows",
    )

    syn = traj[traj["synonymous"] == True]  # noqa: E712
    check("synonymous rows present", len(syn) >= 1, f"{len(syn)} rows")

    # Headline correctness proof: per-window synonymous calibration is exact by
    # construction (T = observed_syn / rate_syn), so sum(actual_syn) ==
    # sum(expected_syn) per window (and overall). This is a real, sensitive check
    # (unlike the median, it is NOT dominated by unobserved mutations).
    if len(syn):
        for w, g in syn.groupby("window_start"):
            exp_sum = g["expected_count"].sum()
            act_sum = g["actual_count"].sum()
            check(
                f"synonymous calibration sum(actual)==sum(expected) [{w}]",
                is_close(act_sum, exp_sum, rel_tol=SYN_CALIBRATION_REL_TOL),
                f"actual={act_sum:.4f}, expected={exp_sum:.4f}",
            )
        # Secondary sanity: the log-ratios also center near 0.
        dfr = syn["delta_fitness_ratio"]
        check(
            "synonymous delta_fitness_ratio centered (|median| small)",
            abs(dfr.median()) <= SYN_CENTER_ABS_TOL_MEDIAN,
            f"median={dfr.median():.4f} (tol {SYN_CENTER_ABS_TOL_MEDIAN})",
        )

    # delta_fitness_ratio consistency with the pseudocount log-ratio formula.
    p = cfg.get("fitness_pseudocount", 0.5)
    recomputed = (
        (traj["actual_count"] + p) / (traj["expected_count"] + p)
    ).apply(lambda x: float("nan") if x <= 0 else x)
    ok = np.allclose(
        traj["delta_fitness_ratio"].to_numpy(),
        np.log(recomputed.to_numpy()),
        equal_nan=True,
        atol=1e-9,
    )
    check("delta_fitness_ratio == log((actual+p)/(expected+p))", ok)

    # Regression sentinels -- test data only. Sentinels left as None until a green
    # run fills them in (avoids locking in a stale pre-redesign number).
    if "test_data" in os.path.basename(os.path.normpath(output_dir)):
        observed = {
            "n_traj_rows": (int(len(traj)), COUNT_REL_TOL),
            "n_windows": (int(traj["window_start"].nunique()), COUNT_REL_TOL),
            "actual_sum": (float(traj["actual_count"].sum()), COUNT_REL_TOL),
            "expected_sum": (float(traj["expected_count"].sum()), FITNESS_FLOAT_REL_TOL),
        }
        for key, expected in WINDOW_FITNESS_SENTINELS.items():
            got, tol = observed[key]
            if expected is None:
                print(f"[INFO] fitness sentinel {key} observed = {got} (not yet pinned)")
                continue
            check(
                f"fitness sentinel: {key}",
                is_close(got, expected, rel_tol=tol),
                f"got {got}, expected {expected} "
                f"(update WINDOW_FITNESS_SENTINELS if intended)",
            )
    else:
        print(f"[SKIP] fitness sentinels (non-test output dir {output_dir})")


def check_founder_rate_table_invariant(output_dir, cfg):
    """The per-founder diff-splice must equal a full-genome rebuild.

    ``ExpectedCounts.founder_rate_table`` recomputes GLM rates only for the sites
    a founder's mutations touch and splices them into the shared reference
    enumeration. That optimization is correct **iff** the recomputed site set
    covers every site whose rate can change -- both the mutated codons (for
    translation) and the ±1 motif halo (the 3-nt rate context). A regression here
    (previously: the motif halo was omitted, silently staling the rate of coding
    sites adjacent to a founder mutation across a codon boundary) is invisible to
    the aggregate synonymous calibration, so guard it directly: rebuild the whole
    genome from the founder sequence and require a byte-for-byte-equal rate table.
    """
    import numpy as np

    needed = {
        "master": os.path.join(output_dir, "glm/master_table.csv"),
        "coding": os.path.join(output_dir, "ref/coding_sites.csv"),
        "ref": os.path.join(output_dir, "ref/ref.fa"),
        "site_mask": os.path.join(output_dir, "ref/site_mask.csv"),
        "founders": os.path.join(output_dir, "windows/founder_muts.csv"),
    }
    rna = cfg.get("rna_structure_data")
    if not all(os.path.exists(p) for p in needed.values()) or not (
        rna and os.path.exists(rna)
    ):
        print(f"[SKIP] founder_rate_table invariant (inputs missing under {output_dir})")
        return

    EC = workflow_loader.load_workflow_module("expected_counts")
    masking = workflow_loader.load_workflow_module("masking")

    ref_seq = EC.load_reference_from_fasta(needed["ref"])
    rates_df = EC.prepare_rates_df(needed["master"])
    vcf_sites = pd.read_csv(needed["site_mask"])["site"].tolist()
    bad_sites = masking.global_masked_sites(
        sites_to_exclude=cfg.get("sites_to_exclude", []),
        vcf_masked_sites=vcf_sites,
        site_include_range=cfg.get("site_include_range", {}),
        genome_len=len(ref_seq),
    )
    basepair_df = EC.build_site_unpaired_map(rna, len(ref_seq))
    ec = EC.ExpectedCounts(
        ref_seq, needed["coding"], bad_sites, basepair_df, rates_df
    )

    founders = pd.read_csv(needed["founders"]).fillna({"founder_nt_mutations": ""})
    founders["_n"] = founders["founder_nt_mutations"].str.split().apply(len)
    # The most-mutated founder maximizes codon-boundary coverage (the failure mode).
    top = founders.sort_values("_n", ascending=False).head(1)
    if top.empty or int(top["_n"].iloc[0]) == 0:
        print("[SKIP] founder_rate_table invariant (no mutated founders)")
        return
    fmuts = EC.founder_muts_to_dict(top["founder_nt_mutations"].iloc[0], bad_sites)

    key = ["site", "mut_nt", "gene"]
    opt = ec.founder_rate_table(fmuts).set_index(key)["rate"].sort_index()
    founder_seq = EC.apply_muts(ref_seq, fmuts)
    full_map = ec.reference_site_map.assign(
        wt_nt=[fmuts.get(s, ref_seq[s - 1]) for s in ec.reference_site_map.site]
    )
    truth = ec._build_rate_table(full_map, founder_seq).set_index(key)["rate"].sort_index()
    joined = opt.to_frame("opt").join(truth.to_frame("truth"), how="outer")
    n_mismatch = int(
        (~np.isclose(joined["opt"], joined["truth"], equal_nan=True)).sum()
    )
    check(
        "founder_rate_table == full-genome rebuild (motif halo)",
        n_mismatch == 0,
        f"{n_mismatch} rate mismatches for founder "
        f"{top['founder_id'].iloc[0]} ({int(top['_n'].iloc[0])} muts)",
    )


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

    print("\n=== (C) windowed expected counts + delta_fitness_ratio (#118) ===")
    check_expected_outputs(output_dir, cfg)
    check_founder_rate_table_invariant(output_dir, cfg)

    print()
    if failures:
        print(f"FAILED ({len(failures)}):")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    print("All windowing validation checks passed.")


if __name__ == "__main__":
    main()
