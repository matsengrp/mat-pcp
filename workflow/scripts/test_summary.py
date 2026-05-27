#!/usr/bin/env python3
"""Generate a summary of test pipeline outputs for validation."""

import hashlib
import json
import os
import sys
from pathlib import Path

import pandas as pd

# Must match ``GENE_FOUNDER_DIRNAMES["S"]`` in ``workflow/Snakefile``. Pulled
# into a named constant so a v3 multi-gene refactor that changes the dirname
# can be found by grep rather than a bare-string sweep across the codebase.
SPIKE_FOUNDER_DIRNAME = "spike"


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


def summarize_rates_table(rates_file):
    """Get summary of rates_table.csv."""
    if not os.path.exists(rates_file):
        return None

    df = pd.read_csv(rates_file)

    return {
        'num_rows': len(df),
        'columns': list(df.columns),
        'mut_types': sorted(df['mut_type'].unique().tolist()) if 'mut_type' in df else None,
        'num_motifs': df['motif'].nunique() if 'motif' in df else None,
        'rate_sum': float(df['rate'].sum()) if 'rate' in df else None,
        'checksum': hashlib.md5(open(rates_file, 'rb').read()).hexdigest()
    }


def summarize_clade_founder_nts(founder_nts_file):
    """Get summary of clade_founder_nts.csv."""
    if not os.path.exists(founder_nts_file):
        return None

    df = pd.read_csv(founder_nts_file)

    return {
        'num_rows': len(df),
        'columns': list(df.columns),
        'clades': sorted(df['clade'].unique().tolist()) if 'clade' in df else None,
        'num_sites': df['site'].nunique() if 'site' in df else None,
        'checksum': hashlib.md5(open(founder_nts_file, 'rb').read()).hexdigest()
    }


def summarize_spike_founders(spike_founders_dir):
    """Summarize gene-sliced (spike) founder FASTAs.

    Records per-file md5s so ``validate_test.py`` can compare against the
    canonical cross-version sentinel (see ``SPIKE_19A_SENTINEL_MD5`` in
    ``validate_test.py``).

    NOTE: This helper assumes the v2 single-gene layout where the gene-sliced
    founder directory is hardcoded to ``clade_founders/<SPIKE_FOUNDER_DIRNAME>/``.
    When v3 introduces multi-gene slicing, this should fan out over every
    configured ``gene_slice_outputs`` entry. See the comment block above
    ``GENE_FOUNDER_DIRNAMES`` in ``workflow/Snakefile`` for the wider context.
    """
    if not os.path.exists(spike_founders_dir):
        return None

    fasta_files = sorted(
        f for f in os.listdir(spike_founders_dir) if f.endswith(".fa")
    )
    summary = {
        "num_files": len(fasta_files),
        "files": fasta_files,
        "checksums": {},
    }
    for fname in fasta_files:
        clade = os.path.splitext(fname)[0]
        path = os.path.join(spike_founders_dir, fname)
        with open(path, "rb") as f:
            summary["checksums"][clade] = hashlib.md5(f.read()).hexdigest()
    return summary


def summarize_dasm_ready_pcps(pcps_dir):
    """Summarize the DASM-ready PCP outputs (#80).

    Records md5s for the aggregate and any per-clade dasm_ready files, plus
    a quick check on the ``in_train`` column. Asserts are deferred to
    ``validate_test.py``.
    """
    if not os.path.exists(pcps_dir):
        return None

    aggregate = os.path.join(pcps_dir, "all_pcps_dasm_ready.csv")
    if not os.path.exists(aggregate):
        return None

    agg_df = pd.read_csv(aggregate)
    per_clade_files = sorted(
        f for f in os.listdir(pcps_dir)
        if f.endswith("_pcps_dasm_ready.csv") and not f.startswith("all_")
    )

    # Compute n_train / n_val via bool-equality (not ~) so a hypothetical
    # NaN in ``in_train`` surfaces as a count mismatch rather than a
    # TypeError; the empty-input path in curate_pcps_for_dasm.py writes a
    # header-only file so this branch may see zero rows.
    has_in_train = "in_train" in agg_df.columns
    n_train = int((agg_df["in_train"] == True).sum()) if has_in_train else None
    n_val = int((agg_df["in_train"] == False).sum()) if has_in_train else None

    summary = {
        "aggregate": {
            "num_rows": len(agg_df),
            "columns": list(agg_df.columns),
            "has_in_train": has_in_train,
            "n_train": n_train,
            "n_val": n_val,
            "checksum": hashlib.md5(Path(aggregate).read_bytes()).hexdigest(),
        },
        "per_clade_files": per_clade_files,
        "per_clade_checksums": {
            os.path.splitext(f)[0]: hashlib.md5(
                Path(pcps_dir, f).read_bytes()
            ).hexdigest()
            for f in per_clade_files
        },
    }
    return summary


def summarize_aamut_fitness(aamut_fitness_dir):
    """Get summary of aamut_fitness outputs."""
    if not os.path.exists(aamut_fitness_dir):
        return None

    all_fitness_file = os.path.join(aamut_fitness_dir, 'all_aamut_fitness.csv')
    if not os.path.exists(all_fitness_file):
        return None

    df = pd.read_csv(all_fitness_file)

    # List per-clade files
    clade_files = [f for f in os.listdir(aamut_fitness_dir)
                   if f.endswith('_aamut_fitness.csv') and f != 'all_aamut_fitness.csv']

    summary = {
        'num_rows': len(df),
        'columns': list(df.columns),
        'clade_files': sorted(clade_files),
        'num_clade_files': len(clade_files),
    }

    # Add numerical summaries for key columns
    if 'expected_count' in df.columns:
        summary['expected_count_sum'] = float(df['expected_count'].sum())
    if 'actual_count' in df.columns:
        summary['actual_count_sum'] = float(df['actual_count'].sum())
    if 'predicted_count' in df.columns:
        summary['predicted_count_sum'] = float(df['predicted_count'].sum())
    if 'delta_fitness' in df.columns:
        summary['delta_fitness_sum'] = float(df['delta_fitness'].sum())

    # Add checksum for the aggregated file
    summary['checksum'] = hashlib.md5(open(all_fitness_file, 'rb').read()).hexdigest()

    return summary


def main(results_dir):
    """Generate test summary for the given results directory."""
    results_path = Path(results_dir)
    mat_path = results_path / "mat"
    pcps_path = results_path / "pcps"

    summary = {
        'results_directory': str(results_path),
        'files': {},
        'pcp_summary': None,
        'clade_counts_summary': None,
        'rates_table_summary': None,
        'clade_founder_nts_summary': None,
        'spike_founders_summary': None,
        'dasm_ready_pcps_summary': None,
        'aamut_fitness_summary': None
    }

    # Check key output files
    key_files = [
        pcps_path / "all_pcps.csv",
        mat_path / "clade_counts.csv",
        mat_path / "filtered_metadata.tsv.gz",
        mat_path / "filtered_tree_raw.pb"
    ]

    for filepath in key_files:
        summary['files'][str(filepath)] = get_file_stats(str(filepath))

    # Get PCP summary
    pcp_file = pcps_path / "all_pcps.csv"
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

    # Get rates table summary
    rates_file = results_path / "rates_table.csv"
    summary['rates_table_summary'] = summarize_rates_table(str(rates_file))

    # Get clade founder nucleotides summary
    founder_nts_file = results_path / "clade_founder_nts.csv"
    summary['clade_founder_nts_summary'] = summarize_clade_founder_nts(str(founder_nts_file))

    # Get gene-sliced (spike) founder summary
    spike_founders_dir = results_path / "clade_founders" / SPIKE_FOUNDER_DIRNAME
    summary['spike_founders_summary'] = summarize_spike_founders(str(spike_founders_dir))

    # Get DASM-ready PCP summary (#80)
    summary['dasm_ready_pcps_summary'] = summarize_dasm_ready_pcps(str(pcps_path))

    # Get amino acid fitness summary
    aamut_fitness_dir = results_path / "aamut_fitness"
    summary['aamut_fitness_summary'] = summarize_aamut_fitness(str(aamut_fitness_dir))

    # Print summary in readable format
    print(json.dumps(summary, indent=2))

    return summary


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python test_summary.py <results_directory>")
        sys.exit(1)
    
    main(sys.argv[1])