#!/usr/bin/env python3
"""
Filtre les taxons faiblement représentés selon plusieurs critères
"""

import pandas as pd
import numpy as np
import sys

# Paramètres Snakemake
input_matrix = snakemake.input.matrix
output_filtered = snakemake.output.filtered
output_removed = snakemake.output.removed
min_reads = snakemake.params.min_reads
min_samples = snakemake.params.min_samples
min_abundance = snakemake.params.min_abundance

print(f"\n=== Filtrage des taxons faiblement représentés ===")
print(f"Paramètres de filtrage:")
print(f"  - Minimum de reads total: {min_reads}")
print(f"  - Présence minimum dans N échantillons: {min_samples}")
print(f"  - Abondance relative minimum: {min_abundance} ({min_abundance*100}%)")

# Lire la matrice
try:
    df = pd.read_csv(input_matrix, sep='\t', index_col=0)
except Exception as e:
    print(f"ERREUR lors de la lecture de {input_matrix}: {e}", file=sys.stderr)
    sys.exit(1)

# Séparer taxonomy_id et les colonnes d'abondance
if 'taxonomy_id' not in df.columns:
    print("ERREUR: Colonne 'taxonomy_id' non trouvée", file=sys.stderr)
    sys.exit(1)

taxonomy_col = df['taxonomy_id']
abundance_cols = df.drop('taxonomy_id', axis=1)

print(f"\nMatrice d'entrée: {len(df)} taxons × {len(abundance_cols.columns)} échantillons")

# Calculer les totaux par échantillon pour les abondances relatives
sample_totals = abundance_cols.sum(axis=0)

# ============================================================================
# Critère 1: Total de reads sur tous les échantillons
# ============================================================================
total_reads_per_taxon = abundance_cols.sum(axis=1)
total_reads_filter = total_reads_per_taxon >= min_reads

print(f"\nCritère 1 (total reads ≥ {min_reads}):")
print(f"  Taxons conservés: {total_reads_filter.sum()}/{len(df)}")

# ============================================================================
# Critère 2: Présence dans au moins X échantillons
# ============================================================================
presence_per_taxon = (abundance_cols > 0).sum(axis=1)
presence_filter = presence_per_taxon >= min_samples

print(f"\nCritère 2 (présent dans ≥ {min_samples} échantillons):")
print(f"  Taxons conservés: {presence_filter.sum()}/{len(df)}")

# ============================================================================
# Critère 3: Abondance relative minimale dans au moins un échantillon
# ============================================================================
relative_abundance = abundance_cols.div(sample_totals, axis=1)
max_rel_abundance_per_taxon = relative_abundance.max(axis=1)
abundance_filter = max_rel_abundance_per_taxon >= min_abundance

print(f"\nCritère 3 (abondance relative ≥ {min_abundance} dans ≥1 échantillon):")
print(f"  Taxons conservés: {abundance_filter.sum()}/{len(df)}")

# ============================================================================
# Combiner les filtres (logique AND)
# ============================================================================
keep_filter = total_reads_filter & presence_filter & abundance_filter

df_filtered = df[keep_filter].copy()
df_removed = df[~keep_filter].copy()

# Sauvegarder la matrice filtrée
df_filtered.to_csv(output_filtered, sep='\t')

print(f"\n=== Résultats du filtrage ===")
print(f"Taxons conservés: {len(df_filtered)} ({len(df_filtered)/len(df)*100:.1f}%)")
print(f"Taxons retirés: {len(df_removed)} ({len(df_removed)/len(df)*100:.1f}%)")
print(f"Reads conservés: {df_filtered.drop('taxonomy_id', axis=1).sum().sum():,}")
print(f"Reads retirés: {df_removed.drop('taxonomy_id', axis=1).sum().sum():,}")

# ============================================================================
# Sauvegarder la liste des taxons retirés avec raisons
# ============================================================================
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
        'taxonomy_id': df_removed.loc[idx, 'taxonomy_id'],
        'total_reads': int(total_reads_per_taxon[idx]),
        'n_samples_present': int(presence_per_taxon[idx]),
        'max_relative_abundance': f"{max_rel_abundance_per_taxon[idx]:.6f}",
        'filter_reasons': '; '.join(reasons)
    })

removed_df = pd.DataFrame(removed_info)
removed_df.to_csv(output_removed, sep='\t', index=False)

print(f"\n✓ Matrice filtrée sauvegardée: {output_filtered}")
print(f"✓ Liste des taxons retirés: {output_removed}")
