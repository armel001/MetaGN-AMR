#!/usr/bin/env python3
"""
Normalise les abondances pour permettre les comparaisons entre échantillons
- Abondance relative (%)
- CLR (Centered Log-Ratio) pour analyses compositionnelles
"""

import pandas as pd
import numpy as np
import sys

# Paramètres Snakemake
input_matrix = snakemake.input.matrix
output_relative = snakemake.output.relative
output_clr = snakemake.output.clr

print(f"\n=== Normalisation des abondances ===")

# Lire la matrice
try:
    df = pd.read_csv(input_matrix, sep='\t', index_col=0)
except Exception as e:
    print(f"ERREUR lors de la lecture de {input_matrix}: {e}", file=sys.stderr)
    sys.exit(1)

# Séparer taxonomy_id et abondances
taxonomy_col = df['taxonomy_id']
abundance = df.drop('taxonomy_id', axis=1)

print(f"Matrice: {len(df)} taxons × {len(abundance.columns)} échantillons")

# ============================================================================
# 1. Abondance relative (pourcentage)
# ============================================================================
print("\n1. Calcul des abondances relatives...")

sample_totals = abundance.sum(axis=0)
relative_abundance = abundance.div(sample_totals, axis=1) * 100

# Vérification
for col in relative_abundance.columns:
    col_sum = relative_abundance[col].sum()
    if not np.isclose(col_sum, 100.0, atol=0.01):
        print(f"  ATTENTION: La somme pour {col} = {col_sum:.2f}% (devrait être 100%)", 
              file=sys.stderr)

relative_df = pd.concat([taxonomy_col, relative_abundance], axis=1)
relative_df.to_csv(output_relative, sep='\t')

print(f"  ✓ Abondances relatives calculées")
print(f"  ✓ Vérification: somme par échantillon ≈ 100%")

# ============================================================================
# 2. CLR transformation (Centered Log-Ratio)
# ============================================================================
print("\n2. Transformation CLR (Centered Log-Ratio)...")
print("   Note: CLR est recommandé pour les données compositionnelles")

# Ajouter une pseudocount pour éviter log(0)
# Utilisation d'une petite valeur (1 read)
pseudocount = 1
abundance_pseudo = abundance + pseudocount

# Calculer la moyenne géométrique par échantillon
# geom_mean = exp(mean(log(x)))
log_abundance = np.log(abundance_pseudo)
log_geom_means = log_abundance.mean(axis=0)
geom_means = np.exp(log_geom_means)

# CLR = log(x_i / geometric_mean)
clr_data = log_abundance.sub(log_geom_means, axis=1)

# Vérification: la somme des CLR par échantillon devrait être ~0
clr_sums = clr_data.sum(axis=0)
if not all(np.abs(clr_sums) < 1e-10):
    print(f"  ATTENTION: Somme CLR par échantillon non nulle (max: {clr_sums.abs().max():.2e})",
          file=sys.stderr)

clr_df = pd.concat([taxonomy_col, clr_data], axis=1)
clr_df.to_csv(output_clr, sep='\t')

print(f"  ✓ Transformation CLR appliquée")
print(f"  ✓ Vérification: somme CLR par échantillon ≈ 0")

# ============================================================================
# Résumé
# ============================================================================
print(f"\n=== Résumé ===")
print(f"✓ Abondances relatives: {output_relative}")
print(f"   - Range: [0, 100]%")
print(f"   - Usage: Visualisations, comparaisons simples")

print(f"\n✓ Transformation CLR: {output_clr}")
print(f"   - Range: [{clr_data.min().min():.2f}, {clr_data.max().max():.2f}]")
print(f"   - Usage: PCA, analyses multivariées, tests statistiques")
print(f"   - Pseudocount utilisé: {pseudocount}")

print("\nRecommandations:")
print("  - Pour barplots/heatmaps: utilisez les abondances relatives")
print("  - Pour PCoA/NMDS/PCA: utilisez les transformations CLR")
print("  - Pour analyses différentielles: utilisez les counts bruts")
