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
        sys.exit(0)


if __name__ == "__main__":
    main()