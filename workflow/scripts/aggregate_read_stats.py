#!/usr/bin/env python3

import pandas as pd
import sys
from pathlib import Path


def main():
    input_files = snakemake.input
    output_file = snakemake.output[0]
    log_file = snakemake.log[0]
    
    log = open(log_file, 'w')
    sys.stderr = sys.stdout = log
    
    print(f"Aggregating {len(input_files)} samples")
    
    all_stats = []
    for file_path in input_files:
        if Path(file_path).exists():
            df = pd.read_csv(file_path, sep='\t')
            all_stats.append(df)
        else:
            print(f"WARNING: {file_path} not found")
    
    if not all_stats:
        print("ERROR: No valid stats files found")
        sys.exit(1)
    
    df_combined = pd.concat(all_stats, ignore_index=True)
    
    for col in ['total_reads', 'total_bases', 'mean_length', 'median_length', 
                'mean_quality', 'median_quality', 'n50']:
        if col in df_combined.columns:
            df_combined[col] = pd.to_numeric(df_combined[col], errors='coerce')
    
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    df_combined.to_csv(output_file, sep='\t', index=False)
    
    print(f"Aggregated {len(df_combined)} samples")
    if 'total_reads' in df_combined.columns:
        print(f"Total reads: {df_combined['total_reads'].sum():,.0f}")
    
    log.close()


if __name__ == "__main__":
    main()
