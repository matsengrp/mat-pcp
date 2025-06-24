#!/usr/bin/env python
"""Aggregate clade-specific PCP files into a single output file."""

import pandas as pd
import os

def main():
    """Aggregate all clade PCP files into a single file with clade annotations."""
    
    # Read all PCP files and concatenate them
    dfs = []
    for pcp_file in snakemake.input.pcps:
        df = pd.read_csv(pcp_file)
        # Add clade information
        clade = os.path.basename(pcp_file).replace('_pcps.csv', '')
        df['clade'] = clade
        dfs.append(df)
    
    # Concatenate all dataframes
    all_pcps_df = pd.concat(dfs, ignore_index=True)
    
    # Save the aggregated results
    all_pcps_df.to_csv(snakemake.output.all_pcps, index=False)
    
    print(f"Aggregated {len(dfs)} clade PCP files into {snakemake.output.all_pcps}")
    print(f"Total PCPs: {len(all_pcps_df)}")
    print(f"Clades processed: {sorted(all_pcps_df['clade'].unique())}")

if __name__ == "__main__":
    main()