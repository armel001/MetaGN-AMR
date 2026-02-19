#!/usr/bin/env python3
"""
Generate taxonomy summary report for independent filtering approach
"""

import pandas as pd
from pathlib import Path
from datetime import datetime
import sys

# Snakemake inputs
kraken_reports = snakemake.input.kraken_reports
diversity_file = snakemake.input.diversity
matrix_files = snakemake.input.matrices
filtering_stats_files = snakemake.input.filtering_stats
output_file = snakemake.output.summary

def parse_kraken_report(report_path):
    """Parse Kraken2 report and extract classification statistics"""
    stats = {'total_reads': 0, 'classified': 0, 'unclassified': 0, 
             'classification_rate': 0, 'top_taxa': {}}
    
    try:
        with open(report_path, 'r') as f:
            first_line = True
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 6:
                    continue
                
                try:
                    num_clade = int(parts[1])
                    rank = parts[3].strip()
                    name = parts[5].strip()
                    
                    if first_line:
                        stats['unclassified'] = num_clade
                        first_line = False
                        continue
                    
                    if rank in ['R', 'R1', 'R2'] and 'root' in name.lower():
                        if stats['classified'] == 0:
                            stats['classified'] = num_clade
                            stats['total_reads'] = stats['classified'] + stats['unclassified']
                    
                    if rank in ['D', 'K', 'P'] and num_clade > 0:
                        if rank not in stats['top_taxa']:
                            stats['top_taxa'][rank] = {}
                        stats['top_taxa'][rank][name.strip()] = num_clade
                
                except (ValueError, IndexError):
                    continue
        
        if stats['total_reads'] > 0:
            stats['classification_rate'] = (stats['classified'] / stats['total_reads']) * 100
        elif stats['unclassified'] > 0:
            stats['total_reads'] = stats['unclassified']
    
    except Exception as e:
        print(f"Error parsing {report_path}: {e}", file=sys.stderr)
    
    return stats

def get_matrix_stats(matrix_path):
    """Extract statistics from abundance matrix"""
    try:
        df = pd.read_csv(matrix_path, sep='\t', index_col=0)
        abundance = df.drop('taxonomy_id', axis=1)
        
        return {
            'n_taxa': len(df),
            'n_samples': len(abundance.columns),
            'total_reads': int(abundance.sum().sum()),
            'mean_reads_per_sample': int(abundance.sum(axis=0).mean()),
            'mean_taxa_per_sample': int((abundance > 0).sum(axis=0).mean())
        }
    except Exception as e:
        print(f"Error reading {matrix_path}: {e}", file=sys.stderr)
        return None

# Collect Kraken2 statistics
kraken_stats = []
for report_path in kraken_reports:
    sample_name = Path(report_path).parent.parent.name
    stats = parse_kraken_report(report_path)
    stats['sample'] = sample_name
    kraken_stats.append(stats)

kraken_df = pd.DataFrame(kraken_stats)

# Load diversity
try:
    diversity_df = pd.read_csv(diversity_file, sep='\t')
except Exception as e:
    diversity_df = pd.DataFrame()

# Matrix statistics
matrix_stats = {}
level_names = {'S': 'Species', 'G': 'Genus', 'F': 'Family', 'P': 'Phylum'}

for matrix_path in matrix_files:
    filename = Path(matrix_path).stem
    parts = filename.split('_')
    level = parts[2] if len(parts) > 2 else 'unknown'
    stats = get_matrix_stats(matrix_path)
    if stats:
        matrix_stats[level] = stats

# Calculate totals
total_reads = kraken_df['total_reads'].sum()
total_classified = kraken_df['classified'].sum()
total_unclassified = kraken_df['unclassified'].sum()
avg_classification = (total_classified / total_reads * 100) if total_reads > 0 else 0

# Generate report
with open(output_file, 'w') as f:
    f.write("=" * 80 + "\n")
    f.write("TAXONOMIC ANALYSIS - INDEPENDENT SAMPLES APPROACH\n")
    f.write("=" * 80 + "\n")
    f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"Samples: {len(kraken_stats)}\n\n")
    
    f.write("FILTERING STRATEGY:\n")
    f.write("  Approach: Biologically heterogeneous/independent samples\n")
    f.write("  Principle: Preserve environment-specific biomarkers\n")
    f.write("  Criteria: Robustness in AT LEAST ONE sample\n\n")
    
    # Section 1: Classification
    f.write("=" * 80 + "\n")
    f.write("1. CLASSIFICATION STATISTICS (KRAKEN2)\n")
    f.write("=" * 80 + "\n\n")
    
    f.write("Overall:\n")
    f.write(f"  Total reads       : {total_reads:,}\n")
    f.write(f"  Classified        : {total_classified:,} ({avg_classification:.2f}%)\n")
    f.write(f"  Unclassified      : {total_unclassified:,}\n\n")
    
    f.write("Per sample:\n")
    f.write("-" * 80 + "\n")
    f.write(f"{'Sample':<25} {'Total reads':>15} {'Classified':>15} {'Rate (%)':>10}\n")
    f.write("-" * 80 + "\n")
    
    for _, row in kraken_df.iterrows():
        f.write(f"{row['sample']:<25} {row['total_reads']:>15,} "
                f"{row['classified']:>15,} {row['classification_rate']:>10.2f}\n")
    
    # Section 2: Diversity
    if not diversity_df.empty:
        f.write("\n\n" + "=" * 80 + "\n")
        f.write("2. ALPHA DIVERSITY (POST-FILTERING)\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("Summary statistics:\n")
        f.write("-" * 80 + "\n")
        div_stats = diversity_df[['shannon_diversity', 'simpson_diversity',
                                   'observed_richness', 'pielou_evenness']].describe()
        f.write(div_stats.to_string())
        f.write("\n\n")
        
        f.write("Per sample:\n")
        f.write("-" * 80 + "\n")
        f.write(f"{'Sample':<25} {'Shannon':>10} {'Simpson':>10} "
                f"{'Richness':>10} {'Evenness':>10}\n")
        f.write("-" * 80 + "\n")
        
        for _, row in diversity_df.iterrows():
            f.write(f"{row['sample']:<25} {row['shannon_diversity']:>10.3f} "
                    f"{row['simpson_diversity']:>10.3f} {row['observed_richness']:>10.0f} "
                    f"{row['pielou_evenness']:>10.3f}\n")
    
    # Section 3: Composition
    f.write("\n\n" + "=" * 80 + "\n")
    f.write("3. TAXONOMIC COMPOSITION (POST-FILTERING)\n")
    f.write("=" * 80 + "\n\n")
    
    if matrix_stats:
        f.write(f"{'Level':<20} {'Taxa':>10} {'Total reads':>15} {'Mean taxa/sample':>20}\n")
        f.write("-" * 80 + "\n")
        
        for level in ['P', 'F', 'G', 'S']:
            if level in matrix_stats:
                stats = matrix_stats[level]
                f.write(f"{level_names.get(level, level):<20} {stats['n_taxa']:>10,} "
                        f"{stats['total_reads']:>15,} {stats['mean_taxa_per_sample']:>20.1f}\n")
    
    # Section 4: Filtering details
    f.write("\n\n" + "=" * 80 + "\n")
    f.write("4. FILTERING STATISTICS BY TAXONOMIC LEVEL\n")
    f.write("=" * 80 + "\n\n")
    
    for stats_file in filtering_stats_files:
        level = Path(stats_file).stem.split('_')[2]
        f.write(f"\n--- {level_names.get(level, level)} level ---\n")
        try:
            with open(stats_file, 'r') as sf:
                f.write(sf.read())
        except Exception as e:
            f.write(f"Error: {e}\n")
    
    # Section 5: Output files
    f.write("\n\n" + "=" * 80 + "\n")
    f.write("5. OUTPUT FILES\n")
    f.write("=" * 80 + "\n\n")
    
    f.write("Filtered abundance matrices:\n")
    for level in ['P', 'F', 'G', 'S']:
        f.write(f"  results/taxonomy/abundance_matrix_{level}_filtered_independent.tsv\n")
    
    f.write("\nRelative abundance (%):\n")
    for level in ['P', 'F', 'G', 'S']:
        f.write(f"  results/taxonomy/abundance_matrix_{level}_relative_independent.tsv\n")
    
    f.write("\nCLR-transformed:\n")
    for level in ['P', 'F', 'G', 'S']:
        f.write(f"  results/taxonomy/abundance_matrix_{level}_clr_independent.tsv\n")
    
    f.write("\nDiversity:\n")
    f.write("  results/taxonomy/alpha_diversity_independent.tsv\n")
    
    f.write("\nFiltering statistics:\n")
    for level in ['P', 'F', 'G', 'S']:
        f.write(f"  results/taxonomy/filtering_stats_{level}_independent.txt\n")
    
    f.write("\n" + "=" * 80 + "\n")

print(f"Summary saved: {output_file}")
