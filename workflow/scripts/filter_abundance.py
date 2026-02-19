#!/usr/bin/env python3
"""
Filter low-abundance taxa based on multiple criteria (standard approach)
"""

import pandas as pd
import sys

# Snakemake parameters
input_matrix = snakemake.input.matrix
output_filtered = snakemake.output.filtered
output_removed = snakemake.output.removed
min_reads = snakemake.params.min_reads
min_samples = snakemake.params.min_samples
min_abundance = snakemake.params.min_abundance

print(f"\n{'='*80}")
print(f"STANDARD FILTERING (HOMOGENEOUS SAMPLES)")
print(f"{'='*80}\n")
print(f"Filtering parameters:")
print(f"  Min total reads          : {min_reads}")
print(f"  Min sample prevalence    : {min_samples}")
print(f"  Min relative abundance   : {min_abundance} ({min_abundance*100}%)")

# Load matrix
try:
    df = pd.read_csv(input_matrix, sep='\t', index_col=0)
except Exception as e:
    print(f"Error loading {input_matrix}: {e}", file=sys.stderr)
    sys.exit(1)

if 'taxonomy_id' not in df.columns:
    print("Error: 'taxonomy_id' column not found", file=sys.stderr)
    sys.exit(1)

taxonomy_col = df['taxonomy_id']
abundance_cols = df.drop('taxonomy_id', axis=1)

print(f"\nInput matrix: {len(df)} taxa × {len(abundance_cols.columns)} samples")

# Calculate sample totals for relative abundance
sample_totals = abundance_cols.sum(axis=0)

# Criterion 1: Total reads across all samples
total_reads_per_taxon = abundance_cols.sum(axis=1)
total_reads_filter = total_reads_per_taxon >= min_reads

print(f"\nCriterion 1 (total reads ≥ {min_reads}):")
print(f"  Taxa passing: {total_reads_filter.sum()}/{len(df)}")

# Criterion 2: Present in at least X samples
presence_per_taxon = (abundance_cols > 0).sum(axis=1)
presence_filter = presence_per_taxon >= min_samples

print(f"\nCriterion 2 (present in ≥ {min_samples} samples):")
print(f"  Taxa passing: {presence_filter.sum()}/{len(df)}")

# Criterion 3: Minimum relative abundance in at least one sample
relative_abundance = abundance_cols.div(sample_totals, axis=1)
max_rel_abundance_per_taxon = relative_abundance.max(axis=1)
abundance_filter = max_rel_abundance_per_taxon >= min_abundance

print(f"\nCriterion 3 (rel. abundance ≥ {min_abundance} in ≥1 sample):")
print(f"  Taxa passing: {abundance_filter.sum()}/{len(df)}")

# Combine filters (AND logic)
keep_filter = total_reads_filter & presence_filter & abundance_filter
df_filtered = df[keep_filter].copy()
df_removed = df[~keep_filter].copy()

# Save filtered matrix
df_filtered.to_csv(output_filtered, sep='\t')

reads_kept = df_filtered.drop('taxonomy_id', axis=1).sum().sum()
reads_removed = df_removed.drop('taxonomy_id', axis=1).sum().sum()

print(f"\n{'='*80}")
print(f"FILTERING RESULTS")
print(f"{'='*80}")
print(f"Taxa retained : {len(df_filtered)} ({len(df_filtered)/len(df)*100:.1f}%)")
print(f"Taxa removed  : {len(df_removed)} ({len(df_removed)/len(df)*100:.1f}%)")
print(f"Reads retained: {reads_kept:,}")
print(f"Reads removed : {reads_removed:,}")

# Save removed taxa with reasons
removed_info = []
for idx in df_removed.index:
    reasons = []
    
    if not total_reads_filter[idx]:
        reasons.append(f"total_reads={int(total_reads_per_taxon[idx])}<{min_reads}")
    if not presence_filter[idx]:
        reasons.append(f"present_in={int(presence_per_taxon[idx])}<{min_samples}_samples")
    if not abundance_filter[idx]:
        reasons.append(f"max_rel_abundance={max_rel_abundance_per_taxon[idx]:.4f}<{min_abundance}")
    
    removed_info.append({
        'taxon': idx,
        'taxonomy_id': int(df_removed.loc[idx, 'taxonomy_id']),
        'total_reads': int(total_reads_per_taxon[idx]),
        'n_samples_present': int(presence_per_taxon[idx]),
        'max_relative_abundance': f"{max_rel_abundance_per_taxon[idx]:.6f}",
        'filter_reasons': '; '.join(reasons)
    })

removed_df = pd.DataFrame(removed_info)
removed_df.to_csv(output_removed, sep='\t', index=False)

print(f"\nFiltered matrix: {output_filtered}")
print(f"Removed taxa   : {output_removed}")
print(f"{'='*80}\n")
