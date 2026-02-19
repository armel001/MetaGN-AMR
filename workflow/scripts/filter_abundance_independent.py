#!/usr/bin/env python3
"""
Taxonomic filtering for biologically heterogeneous/independent samples.
Preserves environment-specific biomarkers by retaining taxa robust in at least one sample.
"""

import pandas as pd
import sys

# Snakemake parameters
input_matrix = snakemake.input.matrix
output_filtered = snakemake.output.filtered
output_removed = snakemake.output.removed
output_stats = snakemake.output.stats
min_reads_per_sample = snakemake.params.min_reads_per_sample
min_rel_abundance_per_sample = snakemake.params.min_rel_abundance_per_sample
min_samples_present = snakemake.params.min_samples_present
absolute_min_reads = snakemake.params.absolute_min_reads

print(f"\n{'='*80}")
print(f"INDEPENDENT SAMPLES FILTERING")
print(f"{'='*80}\n")

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
n_samples = len(abundance_cols.columns)
n_taxa_initial = len(df)

print(f"Input matrix: {n_taxa_initial:,} taxa × {n_samples} samples")
print(f"\nFiltering parameters:")
print(f"  Min reads per sample     : {min_reads_per_sample}")
print(f"  Min relative abundance   : {min_rel_abundance_per_sample} ({min_rel_abundance_per_sample*100}%)")
print(f"  Min sample prevalence    : {min_samples_present}")
print(f"  Absolute min reads       : {absolute_min_reads}\n")

# Calculate sample totals
sample_totals = abundance_cols.sum(axis=0)
print(f"Sequencing depth per sample:")
for sample, total in sample_totals.items():
    print(f"  {sample:<20} : {total:>10,} reads")

# Criterion 1: Minimum reads in at least one sample
max_reads_per_taxon = abundance_cols.max(axis=1)
criterion1 = max_reads_per_taxon >= min_reads_per_sample

print(f"\n{'─'*80}")
print(f"CRITERION 1: ≥ {min_reads_per_sample} reads in ≥1 sample")
print(f"{'─'*80}")
print(f"Taxa passing: {criterion1.sum():,}/{n_taxa_initial:,} ({criterion1.sum()/n_taxa_initial*100:.1f}%)")

# Criterion 2: Minimum relative abundance in at least one sample
relative_abundance = abundance_cols.div(sample_totals, axis=1)
max_rel_abundance_per_taxon = relative_abundance.max(axis=1)
criterion2 = max_rel_abundance_per_taxon >= min_rel_abundance_per_sample

print(f"\n{'─'*80}")
print(f"CRITERION 2: ≥ {min_rel_abundance_per_sample*100}% relative abundance in ≥1 sample")
print(f"{'─'*80}")
print(f"Taxa passing: {criterion2.sum():,}/{n_taxa_initial:,} ({criterion2.sum()/n_taxa_initial*100:.1f}%)")

# Criterion 3: Minimum sample prevalence
presence_per_taxon = (abundance_cols > 0).sum(axis=1)
criterion3 = presence_per_taxon >= min_samples_present

print(f"\n{'─'*80}")
print(f"CRITERION 3: Present in ≥ {min_samples_present} sample(s)")
print(f"{'─'*80}")
print(f"Taxa passing: {criterion3.sum():,}/{n_taxa_initial:,} ({criterion3.sum()/n_taxa_initial*100:.1f}%)")

# Criterion 4: Absolute minimum reads across all samples
total_reads_per_taxon = abundance_cols.sum(axis=1)
criterion4 = total_reads_per_taxon >= absolute_min_reads

print(f"\n{'─'*80}")
print(f"CRITERION 4: ≥ {absolute_min_reads} reads total (all samples)")
print(f"{'─'*80}")
print(f"Taxa passing: {criterion4.sum():,}/{n_taxa_initial:,} ({criterion4.sum()/n_taxa_initial*100:.1f}%)")

# Combine criteria (AND logic)
keep_filter = criterion1 & criterion2 & criterion3 & criterion4
df_filtered = df[keep_filter].copy()
df_removed = df[~keep_filter].copy()

n_taxa_kept = len(df_filtered)
n_taxa_removed = len(df_removed)
reads_kept = df_filtered.drop('taxonomy_id', axis=1).sum().sum()
reads_removed = df_removed.drop('taxonomy_id', axis=1).sum().sum()
total_reads = reads_kept + reads_removed

print(f"\n{'='*80}")
print(f"FILTERING RESULTS")
print(f"{'='*80}")
print(f"Taxa retained    : {n_taxa_kept:>6,} ({n_taxa_kept/n_taxa_initial*100:>5.1f}%)")
print(f"Taxa removed     : {n_taxa_removed:>6,} ({n_taxa_removed/n_taxa_initial*100:>5.1f}%)")
print(f"Reads retained   : {reads_kept:>12,} ({reads_kept/total_reads*100:>5.1f}%)")
print(f"Reads removed    : {reads_removed:>12,} ({reads_removed/total_reads*100:>5.1f}%)")

# Analyze removed taxa
removed_info = []
for idx in df_removed.index:
    reasons = []
    
    if not criterion1[idx]:
        reasons.append(f"max_reads={int(max_reads_per_taxon[idx])}<{min_reads_per_sample}")
    if not criterion2[idx]:
        reasons.append(f"max_rel_abund={max_rel_abundance_per_taxon[idx]:.6f}<{min_rel_abundance_per_sample}")
    if not criterion3[idx]:
        reasons.append(f"prevalence={int(presence_per_taxon[idx])}<{min_samples_present}_samples")
    if not criterion4[idx]:
        reasons.append(f"total_reads={int(total_reads_per_taxon[idx])}<{absolute_min_reads}")
    
    best_sample = abundance_cols.loc[idx].idxmax()
    best_sample_reads = abundance_cols.loc[idx, best_sample]
    best_sample_pct = relative_abundance.loc[idx, best_sample]
    
    removed_info.append({
        'taxon': idx,
        'taxonomy_id': int(df_removed.loc[idx, 'taxonomy_id']),
        'total_reads_all_samples': int(total_reads_per_taxon[idx]),
        'n_samples_present': int(presence_per_taxon[idx]),
        'max_reads_any_sample': int(max_reads_per_taxon[idx]),
        'max_relative_abundance': f"{max_rel_abundance_per_taxon[idx]:.6f}",
        'best_sample': best_sample,
        'reads_in_best_sample': int(best_sample_reads),
        'abundance_in_best_sample': f"{best_sample_pct:.6f}",
        'filter_reasons': '; '.join(reasons)
    })

removed_df = pd.DataFrame(removed_info)
removed_df = removed_df.sort_values('max_reads_any_sample', ascending=False)
removed_df.to_csv(output_removed, sep='\t', index=False)

# Per-sample statistics
sample_stats = []
for sample in abundance_cols.columns:
    sample_data = abundance_cols[sample]
    sample_data_filtered = df_filtered.drop('taxonomy_id', axis=1)[sample]
    
    n_taxa_before = (sample_data > 0).sum()
    n_taxa_after = (sample_data_filtered > 0).sum()
    reads_before = sample_data.sum()
    reads_after = sample_data_filtered.sum()
    
    sample_stats.append({
        'sample': sample,
        'taxa_before': n_taxa_before,
        'taxa_after': n_taxa_after,
        'taxa_removed': n_taxa_before - n_taxa_after,
        'reads_before': int(reads_before),
        'reads_after': int(reads_after),
        'reads_removed': int(reads_before - reads_after),
        'pct_taxa_kept': f"{n_taxa_after/n_taxa_before*100:.1f}%" if n_taxa_before > 0 else "0.0%",
        'pct_reads_kept': f"{reads_after/reads_before*100:.1f}%" if reads_before > 0 else "0.0%"
    })

stats_df = pd.DataFrame(sample_stats)
print(f"\n{'='*80}")
print(f"PER-SAMPLE STATISTICS")
print(f"{'='*80}\n")
print(stats_df.to_string(index=False))

# Save outputs
df_filtered.to_csv(output_filtered, sep='\t')

with open(output_stats, 'w') as f:
    f.write("="*80 + "\n")
    f.write("FILTERING STATISTICS - INDEPENDENT SAMPLES\n")
    f.write("="*80 + "\n\n")
    
    f.write("PARAMETERS:\n")
    f.write(f"  min_reads_per_sample         : {min_reads_per_sample}\n")
    f.write(f"  min_rel_abundance_per_sample : {min_rel_abundance_per_sample} ({min_rel_abundance_per_sample*100}%)\n")
    f.write(f"  min_samples_present          : {min_samples_present}\n")
    f.write(f"  absolute_min_reads           : {absolute_min_reads}\n\n")
    
    f.write("OVERALL RESULTS:\n")
    f.write(f"  Taxa before filtering  : {n_taxa_initial:,}\n")
    f.write(f"  Taxa after filtering   : {n_taxa_kept:,} ({n_taxa_kept/n_taxa_initial*100:.1f}%)\n")
    f.write(f"  Taxa removed           : {n_taxa_removed:,} ({n_taxa_removed/n_taxa_initial*100:.1f}%)\n")
    f.write(f"  Reads retained         : {reads_kept:,} ({reads_kept/total_reads*100:.1f}%)\n")
    f.write(f"  Reads removed          : {reads_removed:,} ({reads_removed/total_reads*100:.1f}%)\n\n")
    
    f.write("PER-SAMPLE:\n")
    f.write(stats_df.to_string(index=False))
    f.write("\n\n")
    
    f.write("INDIVIDUAL CRITERIA:\n")
    f.write(f"  Criterion 1 (reads/sample)   : {criterion1.sum():,} taxa\n")
    f.write(f"  Criterion 2 (rel. abundance) : {criterion2.sum():,} taxa\n")
    f.write(f"  Criterion 3 (prevalence)     : {criterion3.sum():,} taxa\n")
    f.write(f"  Criterion 4 (absolute noise) : {criterion4.sum():,} taxa\n")

print(f"\n{'='*80}")
print(f"OUTPUT FILES")
print(f"{'='*80}")
print(f"Filtered matrix : {output_filtered}")
print(f"Removed taxa    : {output_removed}")
print(f"Statistics      : {output_stats}")
print(f"{'='*80}\n")
