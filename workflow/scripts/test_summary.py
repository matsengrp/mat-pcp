#!/usr/bin/env python3
"""Generate a summary of test pipeline outputs for validation."""

import sys
import os
import pandas as pd
import hashlib
from pathlib import Path


def get_file_stats(filepath):
    """Get basic statistics about a file."""
    if not os.path.exists(filepath):
        return None
    
    stats = {
        'exists': True,
        'size': os.path.getsize(filepath),
        'lines': None,
        'checksum': None
    }
    
    # Get line count for text files
    if filepath.endswith('.csv') or filepath.endswith('.txt'):
        with open(filepath, 'r') as f:
            stats['lines'] = sum(1 for _ in f)
    
    # Get checksum for small files
    if stats['size'] < 1024 * 1024:  # Less than 1MB
        with open(filepath, 'rb') as f:
            stats['checksum'] = hashlib.md5(f.read()).hexdigest()
    
    return stats


def summarize_pcps(pcp_file):
    """Get summary statistics from the PCP file."""
    df = pd.read_csv(pcp_file)
    
    summary = {
        'total_pcps': len(df),
        'columns': list(df.columns),
        'num_columns': len(df.columns),
        'unique_child_clades': df['child_clade'].nunique() if 'child_clade' in df else None,
        'unique_parent_clades': df['parent_clade'].nunique() if 'parent_clade' in df else None,
    }
    
    # Sample first few rows (without full data)
    if len(df) > 0:
        summary['first_5_child_nodes'] = df['child_node'].head().tolist() if 'child_node' in df else None
    
    return summary


def summarize_clade_counts(clade_counts_file):
    """Get summary of clade counts."""
    df = pd.read_csv(clade_counts_file)
    
    return {
        'total_clades': len(df),
        'total_sequences': int(df['count'].sum()),
        'clades_with_adequate_samples': int(df['adequate_sample_counts'].sum()) if 'adequate_sample_counts' in df else None,
        'checksum': hashlib.md5(open(clade_counts_file, 'rb').read()).hexdigest()
    }


def main(results_dir):
    """Generate test summary for the given results directory."""
    results_path = Path(results_dir)
    mat_path = results_path / "mat"
    
    summary = {
        'results_directory': str(results_path),
        'files': {},
        'pcp_summary': None,
        'clade_counts_summary': None
    }
    
    # Check key output files
    key_files = [
        mat_path / "all_pcps.csv",
        mat_path / "clade_counts.csv",
        mat_path / "filtered_metadata.tsv.gz",
        mat_path / "filtered_tree_raw.pb"
    ]
    
    for filepath in key_files:
        summary['files'][str(filepath)] = get_file_stats(str(filepath))
    
    # Get PCP summary
    pcp_file = mat_path / "all_pcps.csv"
    if pcp_file.exists():
        summary['pcp_summary'] = summarize_pcps(str(pcp_file))
    
    # Get clade counts summary  
    clade_counts_file = mat_path / "clade_counts.csv"
    if clade_counts_file.exists():
        summary['clade_counts_summary'] = summarize_clade_counts(str(clade_counts_file))
    
    # Count clade-specific files
    clades_dir = mat_path / "clades"
    if clades_dir.exists():
        clade_files = list(clades_dir.glob("*_samples.txt"))
        summary['num_clade_sample_files'] = len(clade_files)
        
        # Check if all clades have non-empty sample files
        empty_clade_files = []
        for clade_file in clade_files:
            if os.path.getsize(clade_file) == 0:
                empty_clade_files.append(clade_file.name)
        summary['empty_clade_files'] = empty_clade_files
    
    # Print summary in readable format
    import json
    print(json.dumps(summary, indent=2))
    
    return summary


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python test_summary.py <results_directory>")
        sys.exit(1)
    
    main(sys.argv[1])