#!/usr/bin/env python3
"""Validate test pipeline outputs against expected results."""

import json
import math
import re
import sys
from pathlib import Path

# Cross-version integrity sentinel for the gene-sliced 19A spike founder.
# 19A spike is byte-stable across SARS-CoV-2 mat-pcp cutoff dates because the
# Wuhan-Hu-1 reference founder for clade 19A doesn't change between runs.
# A mismatch here means either the slice math drifted, the reference FASTA
# changed, or the FASTA writer's header/wrap convention changed.
SPIKE_19A_SENTINEL_MD5 = "b948f2a1f20edbce64cfff1dc540c119"

# Cross-implementation byte-equivalence sentinels for DASM-ready PCP outputs
# on the test dataset (#80). These md5s come from running the reference
# implementation in ``dasm2.preprocess_command.preprocess`` (the same
# function that ``data-central/scripts/filter_viral_pcps.sh`` invokes) on
# the test pipeline's ``all_pcps.csv`` + ``spike_codon_sites_to_mask.txt``
# with the canonical seeds. A mismatch here means either the filter logic
# drifted from dasm2 or the train/val RNG call sequence changed.
DASM_READY_ALL_SENTINEL_MD5 = "5e538bd088e676b485a11a04745221e9"


def is_close(actual, expected, rel_tol=1e-6, abs_tol=1e-9):
    """Check if two floats are close within tolerance."""
    if expected == 0:
        return abs(actual) <= abs_tol
    return abs(actual - expected) <= max(rel_tol * abs(expected), abs_tol)


def load_json(filepath):
    """Load JSON from file."""
    with open(filepath, 'r') as f:
        return json.load(f)


def compare_summaries(actual, expected):
    """Compare actual and expected summaries, returning differences."""
    differences = []

    # Check total PCPs
    if actual['pcp_summary']['total_pcps'] != expected['pcp_summary']['total_pcps']:
        differences.append(
            f"PCP count mismatch: expected {expected['pcp_summary']['total_pcps']}, "
            f"got {actual['pcp_summary']['total_pcps']}"
        )

    # Check columns
    actual_cols = set(actual['pcp_summary']['columns'])
    expected_cols = set(expected['pcp_summary']['columns'])
    if actual_cols != expected_cols:
        differences.append(
            f"Column mismatch: missing {expected_cols - actual_cols}, "
            f"extra {actual_cols - expected_cols}"
        )

    # Check clade counts
    if actual['clade_counts_summary']['total_clades'] != expected['clade_counts_summary']['total_clades']:
        differences.append(
            f"Clade count mismatch: expected {expected['clade_counts_summary']['total_clades']}, "
            f"got {actual['clade_counts_summary']['total_clades']}"
        )

    # Check total sequences
    if actual['clade_counts_summary']['total_sequences'] != expected['clade_counts_summary']['total_sequences']:
        differences.append(
            f"Total sequence count mismatch: expected {expected['clade_counts_summary']['total_sequences']}, "
            f"got {actual['clade_counts_summary']['total_sequences']}"
        )

    # Check clade counts checksum
    if actual['clade_counts_summary']['checksum'] != expected['clade_counts_summary']['checksum']:
        differences.append("Clade counts file has changed (checksum mismatch)")

    # Check for empty clade files
    if len(actual.get('empty_clade_files', [])) > 0:
        differences.append(f"Found empty clade files: {actual['empty_clade_files']}")

    # Check file existence (normalize paths by removing directory prefix)
    def normalize_path(path):
        """Strip the results directory prefix to compare relative paths."""
        # Remove results/test_data or results/test_data_v1 etc prefix
        return re.sub(r'^results/test_data[^/]*/', '', path)

    expected_normalized = {normalize_path(k): v for k, v in expected['files'].items()}
    actual_normalized = {normalize_path(k): v for k, v in actual['files'].items()}

    for filepath_norm, stats in expected_normalized.items():
        if filepath_norm not in actual_normalized:
            differences.append(f"Missing expected file: {filepath_norm}")
        elif actual_normalized[filepath_norm] is None:
            differences.append(f"Expected file does not exist: {filepath_norm}")

    # Check rates_table.csv
    if expected.get('rates_table_summary'):
        if not actual.get('rates_table_summary'):
            differences.append("rates_table.csv summary missing")
        else:
            exp_rates = expected['rates_table_summary']
            act_rates = actual['rates_table_summary']

            if act_rates['num_rows'] != exp_rates['num_rows']:
                differences.append(
                    f"rates_table row count mismatch: expected {exp_rates['num_rows']}, "
                    f"got {act_rates['num_rows']}"
                )
            if set(act_rates.get('mut_types', [])) != set(exp_rates.get('mut_types', [])):
                differences.append(
                    f"rates_table mutation types mismatch: expected {exp_rates['mut_types']}, "
                    f"got {act_rates['mut_types']}"
                )
            if act_rates.get('num_motifs') != exp_rates.get('num_motifs'):
                differences.append(
                    f"rates_table motif count mismatch: expected {exp_rates['num_motifs']}, "
                    f"got {act_rates['num_motifs']}"
                )
            # Check rate_sum with floating-point tolerance
            if 'rate_sum' in exp_rates and 'rate_sum' in act_rates:
                if not is_close(act_rates['rate_sum'], exp_rates['rate_sum'], rel_tol=1e-6):
                    differences.append(
                        f"rates_table rate_sum mismatch: expected {exp_rates['rate_sum']}, "
                        f"got {act_rates['rate_sum']}"
                    )

    # Check clade_founder_nts.csv
    if expected.get('clade_founder_nts_summary'):
        if not actual.get('clade_founder_nts_summary'):
            differences.append("clade_founder_nts.csv summary missing")
        else:
            exp_nts = expected['clade_founder_nts_summary']
            act_nts = actual['clade_founder_nts_summary']

            if act_nts['num_rows'] != exp_nts['num_rows']:
                differences.append(
                    f"clade_founder_nts row count mismatch: expected {exp_nts['num_rows']}, "
                    f"got {act_nts['num_rows']}"
                )
            if set(act_nts.get('clades', [])) != set(exp_nts.get('clades', [])):
                differences.append(
                    f"clade_founder_nts clades mismatch: expected {exp_nts['clades']}, "
                    f"got {act_nts['clades']}"
                )
            if act_nts['checksum'] != exp_nts['checksum']:
                differences.append("clade_founder_nts.csv has changed (checksum mismatch)")

    # Check gene-sliced (spike) founder FASTAs, including the cross-version
    # 19A integrity sentinel (see module-level ``SPIKE_19A_SENTINEL_MD5``).
    if expected.get('spike_founders_summary'):
        if not actual.get('spike_founders_summary'):
            differences.append("spike_founders_summary missing")
        else:
            exp_spike = expected['spike_founders_summary']
            act_spike = actual['spike_founders_summary']

            if act_spike.get('num_files') != exp_spike.get('num_files'):
                differences.append(
                    f"spike_founders file count mismatch: expected "
                    f"{exp_spike.get('num_files')}, got {act_spike.get('num_files')}"
                )
            if set(act_spike.get('files', [])) != set(exp_spike.get('files', [])):
                differences.append(
                    f"spike_founders files mismatch: expected "
                    f"{exp_spike.get('files')}, got {act_spike.get('files')}"
                )
            for clade, exp_md5 in exp_spike.get('checksums', {}).items():
                act_md5 = act_spike.get('checksums', {}).get(clade)
                if act_md5 != exp_md5:
                    differences.append(
                        f"spike_founders {clade} md5 mismatch: expected "
                        f"{exp_md5}, got {act_md5}"
                    )

            # Hard-coded cross-version sentinel: 19A spike slice must match
            # the canonical Wuhan-Hu-1 reference checksum.
            act_19a = act_spike.get('checksums', {}).get('19A')
            if act_19a and act_19a != SPIKE_19A_SENTINEL_MD5:
                differences.append(
                    f"19A spike sentinel md5 mismatch: expected "
                    f"{SPIKE_19A_SENTINEL_MD5}, got {act_19a}"
                )

    # Check DASM-ready PCP outputs (#80), including byte-equivalence sentinel
    # against the dasm2 reference implementation.
    if expected.get('dasm_ready_pcps_summary'):
        if not actual.get('dasm_ready_pcps_summary'):
            differences.append("dasm_ready_pcps_summary missing")
        else:
            exp_dasm = expected['dasm_ready_pcps_summary']
            act_dasm = actual['dasm_ready_pcps_summary']

            exp_agg = exp_dasm.get('aggregate', {})
            act_agg = act_dasm.get('aggregate', {})

            if act_agg.get('num_rows') != exp_agg.get('num_rows'):
                differences.append(
                    f"dasm_ready aggregate row count mismatch: expected "
                    f"{exp_agg.get('num_rows')}, got {act_agg.get('num_rows')}"
                )
            if not act_agg.get('has_in_train'):
                differences.append("dasm_ready aggregate missing in_train column")
            if act_agg.get('n_train') != exp_agg.get('n_train'):
                differences.append(
                    f"dasm_ready aggregate n_train mismatch: expected "
                    f"{exp_agg.get('n_train')}, got {act_agg.get('n_train')}"
                )
            # n_val asserted separately so a split bug that flips
            # membership without changing total count is caught even if
            # row count and n_train both happen to match.
            if act_agg.get('n_val') != exp_agg.get('n_val'):
                differences.append(
                    f"dasm_ready aggregate n_val mismatch: expected "
                    f"{exp_agg.get('n_val')}, got {act_agg.get('n_val')}"
                )
            if act_agg.get('checksum') != exp_agg.get('checksum'):
                differences.append(
                    f"dasm_ready aggregate md5 mismatch: expected "
                    f"{exp_agg.get('checksum')}, got {act_agg.get('checksum')}"
                )

            # Hard-coded cross-implementation sentinel: the aggregate md5
            # must equal the dasm2 reference output (proves byte-equivalence
            # to filter_viral_pcps.sh).
            if act_agg.get('checksum') and act_agg.get('checksum') != DASM_READY_ALL_SENTINEL_MD5:
                differences.append(
                    f"dasm_ready aggregate sentinel md5 mismatch: expected "
                    f"{DASM_READY_ALL_SENTINEL_MD5}, got {act_agg.get('checksum')}"
                )

            # Per-clade checksum diffs
            for clade, exp_md5 in exp_dasm.get('per_clade_checksums', {}).items():
                act_md5 = act_dasm.get('per_clade_checksums', {}).get(clade)
                if act_md5 != exp_md5:
                    differences.append(
                        f"dasm_ready per-clade {clade} md5 mismatch: expected "
                        f"{exp_md5}, got {act_md5}"
                    )

    # Check aamut_fitness outputs
    if expected.get('aamut_fitness_summary'):
        if not actual.get('aamut_fitness_summary'):
            differences.append("aamut_fitness summary missing")
        else:
            exp_aa = expected['aamut_fitness_summary']
            act_aa = actual['aamut_fitness_summary']

            if act_aa['num_rows'] != exp_aa['num_rows']:
                differences.append(
                    f"aamut_fitness row count mismatch: expected {exp_aa['num_rows']}, "
                    f"got {act_aa['num_rows']}"
                )
            if act_aa['num_clade_files'] != exp_aa['num_clade_files']:
                differences.append(
                    f"aamut_fitness clade file count mismatch: expected {exp_aa['num_clade_files']}, "
                    f"got {act_aa['num_clade_files']}"
                )
            if set(act_aa.get('clade_files', [])) != set(exp_aa.get('clade_files', [])):
                differences.append(
                    f"aamut_fitness clade files mismatch: expected {exp_aa['clade_files']}, "
                    f"got {act_aa['clade_files']}"
                )

            # Check numerical sums with appropriate tolerances
            # Integer-like sums (expected_count, actual_count) use tight tolerance
            for key in ['expected_count_sum', 'actual_count_sum']:
                if key in exp_aa and key in act_aa:
                    if not is_close(act_aa[key], exp_aa[key], rel_tol=1e-6):
                        differences.append(
                            f"aamut_fitness {key} mismatch: expected {exp_aa[key]}, "
                            f"got {act_aa[key]}"
                        )
            # Floating-point sums can vary more across platforms, use looser tolerance
            for key in ['predicted_count_sum', 'delta_fitness_sum']:
                if key in exp_aa and key in act_aa:
                    if not is_close(act_aa[key], exp_aa[key], rel_tol=0.01):
                        differences.append(
                            f"aamut_fitness {key} mismatch: expected {exp_aa[key]}, "
                            f"got {act_aa[key]} (rel_diff={abs(act_aa[key] - exp_aa[key]) / abs(exp_aa[key]):.4f})"
                        )

    return differences


def main():
    """Run validation."""
    import os
    import subprocess

    # Determine output directory from args, env var, or default
    output_dir = sys.argv[1] if len(sys.argv) > 1 else os.environ.get('OUTPUT_DIR', 'results/test_data')

    # Generate current summary
    result = subprocess.run(
        ['python', 'workflow/scripts/test_summary.py', output_dir],
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        print(f"Failed to generate test summary: {result.stderr}")
        sys.exit(1)
    
    actual = json.loads(result.stdout)
    
    # Load expected summary
    expected_path = Path('data/test_data/expected_summary.json')
    if not expected_path.exists():
        print(f"Expected summary not found at {expected_path}")
        print("Current summary:")
        print(json.dumps(actual, indent=2))
        sys.exit(1)
    
    expected = load_json(expected_path)
    
    # Compare
    differences = compare_summaries(actual, expected)
    
    if differences:
        print("VALIDATION FAILED!")
        print("\nDifferences found:")
        for diff in differences:
            print(f"  - {diff}")
        print("\nActual summary:")
        print(json.dumps(actual, indent=2))
        sys.exit(1)
    else:
        print("✓ All validation checks passed!")
        print(f"  - Total PCPs: {actual['pcp_summary']['total_pcps']}")
        print(f"  - Total clades: {actual['clade_counts_summary']['total_clades']}")
        print(f"  - Total sequences: {actual['clade_counts_summary']['total_sequences']}")
        if actual.get('rates_table_summary'):
            print(f"  - Rates table rows: {actual['rates_table_summary']['num_rows']}")
        if actual.get('clade_founder_nts_summary'):
            print(f"  - Clade founder nucleotides: {actual['clade_founder_nts_summary']['num_rows']}")
        if actual.get('aamut_fitness_summary'):
            print(f"  - AA mutation fitness rows: {actual['aamut_fitness_summary']['num_rows']}")
            print(f"  - AA mutation fitness clade files: {actual['aamut_fitness_summary']['num_clade_files']}")
        sys.exit(0)


if __name__ == "__main__":
    main()