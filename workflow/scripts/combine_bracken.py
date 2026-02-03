#!/usr/bin/env python3
"""
Combine les outputs Bracken de tous les échantillons en une matrice d'abondance
"""

import pandas as pd
from pathlib import Path
import sys

# Paramètres Snakemake
input_files = snakemake.input
output_matrix = snakemake.output.abundance_matrix
output_metadata = snakemake.output.metadata
level = snakemake.wildcards.level

# Vérifier qu'on a des fichiers
if not input_files:
    print("ERREUR: Aucun fichier d'entrée fourni", file=sys.stderr)
    sys.exit(1)

# Liste pour stocker les données
all_data = []
sample_stats = []

print(f"\n=== Combinaison des échantillons au niveau {level} ===")
print(f"Nombre de fichiers à combiner: {len(input_files)}")

for filepath in input_files:
    # Extraire le nom de l'échantillon
    sample_name = Path(filepath).parent.parent.name
    
    try:
        # Lire le fichier Bracken
        df = pd.read_csv(filepath, sep='\t')
        
        # Vérifier que le fichier contient les colonnes attendues
        required_cols = ['name', 'taxonomy_id', 'new_est_reads', 'fraction_total_reads']
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            print(f"ATTENTION: Colonnes manquantes dans {filepath}: {missing_cols}", file=sys.stderr)
            continue
        
        # Garder les colonnes importantes
        df_subset = df[required_cols].copy()
        df_subset['sample'] = sample_name
        
        all_data.append(df_subset)
        
        # Statistiques par échantillon
        sample_stats.append({
            'sample': sample_name,
            'n_taxa': len(df),
            'total_reads': int(df['new_est_reads'].sum()),
            'taxonomic_level': level
        })
        
        print(f"  ✓ {sample_name}: {len(df)} taxons, {int(df['new_est_reads'].sum()):,} reads")
        
    except Exception as e:
        print(f"ERREUR lors de la lecture de {filepath}: {e}", file=sys.stderr)
        continue

# Vérifier qu'on a récupéré des données
if not all_data:
    print("ERREUR: Aucune donnée n'a pu être combinée", file=sys.stderr)
    sys.exit(1)

# Combiner tous les échantillons
combined = pd.concat(all_data, ignore_index=True)

print(f"\n=== Création de la matrice d'abondance ===")

# Créer la matrice d'abondance (taxons en lignes, échantillons en colonnes)
abundance_wide = combined.pivot_table(
    index='name',
    columns='sample',
    values='new_est_reads',
    fill_value=0
)

# Ajouter la colonne taxonomy_id
taxonomy_mapping = combined[['name', 'taxonomy_id']].drop_duplicates().set_index('name')
abundance_wide = abundance_wide.join(taxonomy_mapping)

# Réorganiser : taxonomy_id en première colonne
cols = ['taxonomy_id'] + [col for col in abundance_wide.columns if col != 'taxonomy_id']
abundance_wide = abundance_wide[cols]

# Trier par abondance totale (décroissant)
abundance_wide['_total'] = abundance_wide.drop('taxonomy_id', axis=1).sum(axis=1)
abundance_wide = abundance_wide.sort_values('_total', ascending=False)
abundance_wide = abundance_wide.drop('_total', axis=1)

# Sauvegarder la matrice
abundance_wide.to_csv(output_matrix, sep='\t')

print(f"  ✓ Matrice sauvegardée: {output_matrix}")
print(f"  ✓ Dimensions: {len(abundance_wide)} taxons × {len(abundance_wide.columns)-1} échantillons")

# Créer le fichier de métadonnées des échantillons
metadata_df = pd.DataFrame(sample_stats)
metadata_df = metadata_df.sort_values('sample')
metadata_df.to_csv(output_metadata, sep='\t', index=False)

print(f"  ✓ Métadonnées sauvegardées: {output_metadata}")

# Résumé final
print(f"\n=== Résumé ===")
print(f"Niveau taxonomique: {level}")
print(f"Nombre total de taxons uniques: {len(abundance_wide)}")
print(f"Nombre d'échantillons: {len(sample_stats)}")
print(f"Total de reads classifiés: {metadata_df['total_reads'].sum():,}")
print(f"Moyenne de reads par échantillon: {metadata_df['total_reads'].mean():.0f}")
print(f"Moyenne de taxons par échantillon: {metadata_df['n_taxa'].mean():.1f}")
