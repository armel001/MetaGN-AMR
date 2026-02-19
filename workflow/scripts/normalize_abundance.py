#!/usr/bin/env python3
"""
Normalize abundance matrices with actual read length statistics
Normalization by gigabases (Gb) for Nanopore long-read data
"""

import pandas as pd
import numpy as np
import sys

# Snakemake inputs
input_matrix = snakemake.input.matrix
output_relative = snakemake.output.relative
output_clr = snakemake.output.clr

# Path to sequencing stats
STATS_FILE = "results/stats/sequencing_stats.tsv"

print(f"\n{'='*80}")
print(f"ABUNDANCE NORMALIZATION - NANOPORE LONG-READ")
print(f"{'='*80}\n")

# ============================================================================
# Load sequencing statistics
# ============================================================================

try:
    seq_stats = pd.read_csv(STATS_FILE, sep='\t', index_col='sample_id')
    print(f"Loaded sequencing statistics from: {STATS_FILE}")
    print(f"\nSequencing statistics:")
    print(seq_stats[['total_reads', 'total_bases', 'mean_length', 'n50']])
except Exception as e:
    print(f"Warning: Could not load {STATS_FILE}: {e}", file=sys.stderr)
    print("Falling back to fixed mean read length of 2500 bp")
    seq_stats = None

# ============================================================================
# Load abundance matrix
# ============================================================================

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

n_taxa = len(df)
n_samples = len(abundance_cols.columns)

print(f"\nInput matrix: {n_taxa} taxa × {n_samples} samples")

# ============================================================================
# Get mean read length per sample
# ============================================================================

sample_read_lengths = {}

for sample in abundance_cols.columns:
    if seq_stats is not None and sample in seq_stats.index:
        # Use actual mean length from stats
        sample_read_lengths[sample] = seq_stats.loc[sample, 'mean_length']
    else:
        # Fallback to default
        sample_read_lengths[sample] = 2500
        print(f"Warning: No stats for {sample}, using default 2500 bp")

print(f"\n{'─'*80}")
print(f"Read lengths per sample:")
print(f"{'─'*80}")
for sample, length in sample_read_lengths.items():
    print(f"  {sample:<20} : {length:>8.1f} bp")

# ============================================================================
# Normalize by gigabases
# ============================================================================

print(f"\n{'─'*80}")
print(f"NORMALIZATION BY GIGABASES")
print(f"{'─'*80}")

# Convert reads to bases for each sample
# (each sample can have different mean read length)
bases_matrix = pd.DataFrame(index=abundance_cols.index, 
                            columns=abundance_cols.columns)

for sample in abundance_cols.columns:
    read_length = sample_read_lengths[sample]
    bases_matrix[sample] = abundance_cols[sample] * read_length

# Total bases per sample (in bp)
bases_totals = bases_matrix.sum(axis=0).astype(float)

# Convert to Gb (1 Gb = 1e9 bp)
gb_totals = bases_totals / 1e9

print(f"\nSequencing output per sample:")
print(f"{'Sample':<20} {'Reads':>12} {'Total Gb':>12} {'Mean length':>12}")
print(f"{'─'*60}")

for sample in abundance_cols.columns:
    reads = abundance_cols[sample].sum()
    gb = gb_totals[sample]
    mean_len = sample_read_lengths[sample]
    print(f"{sample:<20} {reads:>12,} {gb:>12.4f} {mean_len:>12.1f} bp")

# Calculate relative abundance normalized by Gb
# Each taxon's proportion of total sequencing output
relative_abundance = bases_matrix.div(bases_totals, axis=1) * 100

# ============================================================================
# Verify relative abundance
# ============================================================================

sums = relative_abundance.sum(axis=0)
print(f"\n{'─'*80}")
print(f"VERIFICATION: Relative abundance sums (should be 100%)")
print(f"{'─'*80}")

all_correct = True
for sample, total in sums.items():
    status = "✓" if abs(total - 100.0) < 0.01 else "✗"
    if abs(total - 100.0) >= 0.01:
        all_correct = False
    print(f"  {sample:<20} : {total:>10.6f}% {status}")

if not all_correct:
    print("\n⚠️  Warning: Some sums not exactly 100% (check for rounding errors)")

# ============================================================================
# CLR transformation (Centered Log-Ratio)
# ============================================================================

print(f"\n{'─'*80}")
print(f"CLR TRANSFORMATION")
print(f"{'─'*80}")

# Use read counts for CLR (not Gb-normalized values)
# CLR is scale-independent, so read counts are appropriate
pseudocount = 1
count_pseudo = abundance_cols + pseudocount

print(f"Pseudocount added: {pseudocount}")
print(f"Transformation: CLR = log(count + {pseudocount}) - geometric_mean")

# Log-transformation
log_matrix = np.log(count_pseudo)

# Geometric mean per sample
geometric_means = log_matrix.mean(axis=0)

# CLR: subtract geometric mean
clr_matrix = log_matrix.sub(geometric_means, axis=1)

# Verification: CLR sums should be ≈ 0
clr_sums = clr_matrix.sum(axis=0)
print(f"\nCLR sums (should be ≈ 0):")

all_clr_correct = True
for sample, total in clr_sums.items():
    status = "✓" if abs(total) < 0.001 else "✗"
    if abs(total) >= 0.001:
        all_clr_correct = False
    print(f"  {sample:<20} : {total:>14.8f} {status}")

if not all_clr_correct:
    print("\n⚠️  Warning: CLR sums not close to 0 (numerical precision issue)")

# ============================================================================
# Save outputs
# ============================================================================

# Relative abundance (Gb-normalized, with taxonomy_id)
relative_output = pd.concat([taxonomy_col, relative_abundance], axis=1)
relative_output.to_csv(output_relative, sep='\t')

# CLR (with taxonomy_id)
clr_output = pd.concat([taxonomy_col, clr_matrix], axis=1)
clr_output.to_csv(output_clr, sep='\t')

print(f"\n{'='*80}")
print(f"OUTPUT FILES GENERATED")
print(f"{'='*80}")
print(f"✓ Relative abundance (Gb-normalized): {output_relative}")
print(f"✓ CLR-transformed                    : {output_clr}")
print(f"\nNormalization method: GIGABASE-BASED using actual read lengths")
print(f"{'='*80}\n")

# ============================================================================
# Summary statistics
# ============================================================================

print(f"SUMMARY:")
print(f"  Taxa           : {n_taxa}")
print(f"  Samples        : {n_samples}")
print(f"  Total Gb       : {gb_totals.sum():.4f} Gb")
print(f"  Mean Gb/sample : {gb_totals.mean():.4f} Gb")
print(f"  Read length    : {pd.Series(sample_read_lengths).mean():.1f} bp (mean across samples)")
print(f"\n{'='*80}\n")
