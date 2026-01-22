#!/usr/bin/env python3
"""Validate test pipeline outputs against expected results."""

import json
import sys
from pathlib import Path


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

    # Check file existence
    for filepath, stats in expected['files'].items():
        if filepath not in actual['files']:
            differences.append(f"Missing expected file: {filepath}")
        elif actual['files'][filepath] is None:
            differences.append(f"Expected file does not exist: {filepath}")

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
            if act_rates['checksum'] != exp_rates['checksum']:
                differences.append("rates_table.csv has changed (checksum mismatch)")

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

            # Check numerical sums (with tolerance for floating point)
            for key in ['expected_count_sum', 'actual_count_sum', 'predicted_count_sum', 'delta_fitness_sum']:
                if key in exp_aa and key in act_aa:
                    if abs(act_aa[key] - exp_aa[key]) > 1e-6:
                        differences.append(
                            f"aamut_fitness {key} mismatch: expected {exp_aa[key]}, "
                            f"got {act_aa[key]}"
                        )

            if act_aa['checksum'] != exp_aa['checksum']:
                differences.append("all_aamut_fitness.csv has changed (checksum mismatch)")

    return differences


def main():
    """Run validation."""
    # Generate current summary
    import subprocess
    result = subprocess.run(
        ['python', 'workflow/scripts/test_summary.py', 'results_test_data'],
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        print(f"Failed to generate test summary: {result.stderr}")
        sys.exit(1)
    
    actual = json.loads(result.stdout)
    
    # Load expected summary
    expected_path = Path('test_data/expected_summary.json')
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