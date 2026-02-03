#!/usr/bin/env python3
"""
Génère un résumé statistique complet du projet de métagénomique
"""

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
import sys

# Paramètres Snakemake
kraken_reports = snakemake.input.kraken_reports
diversity_file = snakemake.input.diversity
matrix_files = snakemake.input.matrices
output_file = snakemake.output.summary

def parse_kraken_summary(report_path):
    """
    Parse un rapport Kraken2 pour extraire les statistiques
    
    Format Kraken2:
    %       num_clade   num_direct  rank    taxid   name
    93.73   1060163     1060163     U       0       unclassified
    6.27    70875       9959        R       1       root
    """
    stats = {
        'total_reads': 0,
        'classified': 0,
        'unclassified': 0,
        'classification_rate': 0,
        'top_taxa': {}
    }
    
    try:
        with open(report_path, 'r') as f:
            first_line = True
            
            for line in f:
                parts = line.strip().split('\t')
                
                if len(parts) < 6:
                    continue
                
                try:
                    pct = float(parts[0])
                    num_clade = int(parts[1])
                    num_direct = int(parts[2])
                    rank = parts[3].strip()
                    taxid = parts[4].strip()
                    name = parts[5].strip()
                    
                    # Première ligne = unclassified
                    if first_line:
                        stats['unclassified'] = num_clade
                        first_line = False
                        continue
                    
                    # Ligne "root" = total classifiés
                    if rank in ['R', 'R1', 'R2'] and 'root' in name.lower():
                        if stats['classified'] == 0:  # Prendre seulement la première ligne "root"
                            stats['classified'] = num_clade
                            stats['total_reads'] = stats['classified'] + stats['unclassified']
                    
                    # Capturer les rangs taxonomiques élevés (Domain, Kingdom, Phylum)
                    # D = Domain, K = Kingdom, P = Phylum
                    if rank in ['D', 'K', 'P'] and num_clade > 0:
                        # Nettoyer le nom (enlever l'indentation)
                        clean_name = name.strip()
                        if rank not in stats['top_taxa']:
                            stats['top_taxa'][rank] = {}
                        stats['top_taxa'][rank][clean_name] = num_clade
                
                except (ValueError, IndexError) as e:
                    continue
        
        # Calculer le taux de classification
        if stats['total_reads'] > 0:
            stats['classification_rate'] = (stats['classified'] / stats['total_reads']) * 100
        
        # Si total_reads est toujours 0, utiliser unclassified
        if stats['total_reads'] == 0 and stats['unclassified'] > 0:
            stats['total_reads'] = stats['unclassified']
    
    except Exception as e:
        print(f"ERREUR lors du parsing de {report_path}: {e}", file=sys.stderr)
    
    return stats

def get_matrix_stats(matrix_path):
    """Extrait les statistiques d'une matrice d'abondance"""
    try:
        df = pd.read_csv(matrix_path, sep='\t', index_col=0)
        abundance_data = df.drop('taxonomy_id', axis=1)
        
        return {
            'n_taxa': len(df),
            'n_samples': len(abundance_data.columns),
            'total_reads': int(abundance_data.sum().sum()),
            'mean_reads_per_sample': int(abundance_data.sum(axis=0).mean()),
            'mean_taxa_per_sample': int((abundance_data > 0).sum(axis=0).mean()),
            'min_abundance': int(abundance_data[abundance_data > 0].min().min()) if (abundance_data > 0).any().any() else 0,
            'max_abundance': int(abundance_data.max().max()) if abundance_data.max().max() > 0 else 0
        }
    except Exception as e:
        print(f"ERREUR lors de la lecture de {matrix_path}: {e}", file=sys.stderr)
        return None

print(f"\n=== Génération du résumé du projet ===")

# Collecter les statistiques Kraken2
kraken_stats = []
for report_path in kraken_reports:
    sample_name = Path(report_path).parent.parent.name
    stats = parse_kraken_summary(report_path)
    stats['sample'] = sample_name
    kraken_stats.append(stats)
    
    # Debug
    print(f"\n{sample_name}:")
    print(f"  Total reads: {stats['total_reads']:,}")
    print(f"  Classified: {stats['classified']:,} ({stats['classification_rate']:.2f}%)")
    print(f"  Unclassified: {stats['unclassified']:,}")

kraken_df = pd.DataFrame(kraken_stats)

# Charger les statistiques de diversité
try:
    diversity_df = pd.read_csv(diversity_file, sep='\t')
except Exception as e:
    print(f"ATTENTION: Impossible de lire {diversity_file}: {e}", file=sys.stderr)
    diversity_df = pd.DataFrame()

# Statistiques des matrices par niveau
matrix_stats = {}
level_names = {
    'S': 'Espèce (Species)',
    'G': 'Genre (Genus)',
    'F': 'Famille (Family)',
    'P': 'Phylum'
}

for matrix_path in matrix_files:
    filename = Path(matrix_path).stem
    parts = filename.split('_')
    level = parts[2] if len(parts) > 2 else 'unknown'
    
    stats = get_matrix_stats(matrix_path)
    if stats:
        matrix_stats[level] = stats

# ============================================================================
# Générer le rapport texte
# ============================================================================

with open(output_file, 'w') as f:
    f.write("=" * 80 + "\n")
    f.write("ANALYSE TAXONOMIQUE MÉTAGÉNOMIQUE - RÉSUMÉ COMPLET\n")
    f.write("=" * 80 + "\n")
    f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"Nombre d'échantillons analysés: {len(kraken_stats)}\n")
    f.write("\n")
    
    # === SECTION 1: CLASSIFICATION KRAKEN2 ===
    f.write("=" * 80 + "\n")
    f.write("1. STATISTIQUES DE CLASSIFICATION (KRAKEN2)\n")
    f.write("=" * 80 + "\n\n")
    
    total_reads = kraken_df['total_reads'].sum()
    total_classified = kraken_df['classified'].sum()
    total_unclassified = kraken_df['unclassified'].sum()
    avg_classification = (total_classified / total_reads * 100) if total_reads > 0 else 0
    
    f.write("Statistiques globales:\n")
    f.write(f"  Total de reads séquencés: {total_reads:,}\n")
    f.write(f"  Total de reads classifiés: {total_classified:,} ({avg_classification:.2f}%)\n")
    f.write(f"  Total de reads non classifiés: {total_unclassified:,} ({100-avg_classification:.2f}%)\n\n")
    
    f.write("Par échantillon:\n")
    f.write("-" * 80 + "\n")
    f.write(f"{'Échantillon':<25} {'Reads totaux':>15} {'Classifiés':>15} {'Non classifiés':>15} {'Taux (%)':>10}\n")
    f.write("-" * 80 + "\n")
    
    for _, row in kraken_df.iterrows():
        f.write(f"{row['sample']:<25} {row['total_reads']:>15,} "
                f"{row['classified']:>15,} {row['unclassified']:>15,} "
                f"{row['classification_rate']:>10.2f}\n")
    
    # Taxons de haut niveau détectés
    f.write("\n\nTaxons de haut niveau détectés (tous échantillons):\n")
    f.write("-" * 80 + "\n")
    
    # Agréger les taxons de tous les échantillons
    all_top_taxa = {}
    for stats in kraken_stats:
        for rank, taxa_dict in stats['top_taxa'].items():
            if rank not in all_top_taxa:
                all_top_taxa[rank] = {}
            for taxon, count in taxa_dict.items():
                if taxon not in all_top_taxa[rank]:
                    all_top_taxa[rank][taxon] = 0
                all_top_taxa[rank][taxon] += count
    
    rank_names = {'D': 'Domaines', 'K': 'Kingdoms', 'P': 'Phyla'}
    
    for rank in ['D', 'K', 'P']:
        if rank in all_top_taxa and all_top_taxa[rank]:
            f.write(f"\n{rank_names.get(rank, rank)}:\n")
            for taxon, count in sorted(all_top_taxa[rank].items(), key=lambda x: x[1], reverse=True):
                pct = (count / total_classified * 100) if total_classified > 0 else 0
                f.write(f"  {taxon:<50} {count:>12,} reads ({pct:>6.2f}%)\n")
    
    # === SECTION 2: DIVERSITÉ ALPHA ===
    if not diversity_df.empty:
        f.write("\n\n" + "=" * 80 + "\n")
        f.write("2. INDICES DE DIVERSITÉ ALPHA (NIVEAU ESPÈCE)\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("Statistiques descriptives:\n")
        f.write("-" * 80 + "\n")
        diversity_stats = diversity_df[['shannon_diversity', 'simpson_diversity', 
                                         'observed_richness', 'pielou_evenness']].describe()
        f.write(diversity_stats.to_string())
        f.write("\n\n")
        
        f.write("Par échantillon:\n")
        f.write("-" * 80 + "\n")
        f.write(f"{'Échantillon':<25} {'Shannon':>10} {'Simpson':>10} "
                f"{'Richesse':>10} {'Équitabilité':>13}\n")
        f.write("-" * 80 + "\n")
        
        for _, row in diversity_df.iterrows():
            f.write(f"{row['sample']:<25} {row['shannon_diversity']:>10.3f} "
                    f"{row['simpson_diversity']:>10.3f} {row['observed_richness']:>10.0f} "
                    f"{row['pielou_evenness']:>13.3f}\n")
    
    # === SECTION 3: COMPOSITION TAXONOMIQUE ===
    f.write("\n\n" + "=" * 80 + "\n")
    f.write("3. COMPOSITION TAXONOMIQUE PAR NIVEAU\n")
    f.write("=" * 80 + "\n\n")
    
    if matrix_stats:
        f.write(f"{'Niveau':<30} {'Taxons':>12} {'Reads totaux':>15} "
                f"{'Moy taxons/échantillon':>25}\n")
        f.write("-" * 80 + "\n")
        
        for level in ['P', 'F', 'G', 'S']:
            if level in matrix_stats:
                stats = matrix_stats[level]
                f.write(f"{level_names.get(level, level):<30} {stats['n_taxa']:>12,} "
                        f"{stats['total_reads']:>15,} "
                        f"{stats['mean_taxa_per_sample']:>25.1f}\n")
    
    # === SECTION 4: QUALITÉ DES DONNÉES ===
    f.write("\n\n" + "=" * 80 + "\n")
    f.write("4. ÉVALUATION DE LA QUALITÉ\n")
    f.write("=" * 80 + "\n\n")
    
    min_reads = kraken_df['total_reads'].min()
    max_reads = kraken_df['total_reads'].max()
    mean_reads = kraken_df['total_reads'].mean()
    read_ratio = max_reads / min_reads if min_reads > 0 else 0
    
    f.write("Profondeur de séquençage:\n")
    f.write(f"  Minimum: {min_reads:,} reads\n")
    f.write(f"  Maximum: {max_reads:,} reads\n")
    f.write(f"  Moyenne: {mean_reads:,.0f} reads\n")
    f.write(f"  Ratio max/min: {read_ratio:.2f}x\n")
    
    if read_ratio > 5:
        f.write(f"  ⚠ Profondeur très variable (>{5}x) → Normalisation recommandée\n")
    elif read_ratio > 2:
        f.write(f"  ℹ Profondeur modérément variable → Normalisation pour stats\n")
    else:
        f.write(f"  ✓ Profondeur homogène\n")
    
    f.write(f"\nTaux de classification:\n")
    f.write(f"  Moyen: {avg_classification:.2f}%\n")
    f.write(f"  Range: {kraken_df['classification_rate'].min():.2f}% - "
            f"{kraken_df['classification_rate'].max():.2f}%\n")
    
    if avg_classification < 10:
        f.write(f"  ⚠ Taux très faible pour environnements peu caractérisés\n")
    elif avg_classification < 30:
        f.write(f"  ℹ Taux normal pour échantillons environnementaux complexes\n")
    else:
        f.write(f"  ✓ Bon taux pour échantillons environnementaux\n")
    
    # === SECTION 5: FICHIERS ===
    f.write("\n\n" + "=" * 80 + "\n")
    f.write("5. FICHIERS PRÊTS POUR L'ANALYSE\n")
    f.write("=" * 80 + "\n\n")
    
    f.write("Matrices d'abondance filtrées:\n")
    for level in ['P', 'F', 'G', 'S']:
        f.write(f"  - results/taxonomy/abundance_matrix_{level}_filtered.tsv\n")
    
    f.write("\nMatrices normalisées (% relatif):\n")
    for level in ['P', 'F', 'G', 'S']:
        f.write(f"  - results/taxonomy/abundance_matrix_{level}_relative.tsv\n")
    
    f.write("\nMatrices CLR (pour PCA/PCoA):\n")
    for level in ['P', 'F', 'G', 'S']:
        f.write(f"  - results/taxonomy/abundance_matrix_{level}_clr.tsv\n")
    
    f.write("\nDiversité:\n")
    f.write("  - results/taxonomy/alpha_diversity.tsv\n")
    
    f.write("\n" + "=" * 80 + "\n")
    f.write("Fin du rapport\n")
    f.write("=" * 80 + "\n")

print(f"\n✓ Résumé sauvegardé: {output_file}")
print(f"\nRÉSUMÉ RAPIDE:")
print(f"  Échantillons: {len(kraken_stats)}")
print(f"  Total reads: {total_reads:,}")
print(f"  Taux classification: {avg_classification:.2f}%")
if not diversity_df.empty:
    print(f"  Shannon moyen: {diversity_df['shannon_diversity'].mean():.3f}")
