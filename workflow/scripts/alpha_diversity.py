#!/usr/bin/env python3
"""
Calcule les indices de diversité alpha pour chaque échantillon
- Shannon (H')
- Simpson (1-D)
- Richesse observée
- Équitabilité de Pielou (J')
- Estimation Chao1
"""

import pandas as pd
import numpy as np
import sys

def shannon_index(counts):
    """Indice de Shannon (H')"""
    counts = counts[counts > 0]
    if len(counts) == 0:
        return 0
    proportions = counts / counts.sum()
    return -np.sum(proportions * np.log(proportions))

def simpson_index(counts):
    """Indice de Simpson (1 - D)"""
    counts = counts[counts > 0]
    if len(counts) == 0:
        return 0
    proportions = counts / counts.sum()
    return 1 - np.sum(proportions ** 2)

def observed_richness(counts):
    """Nombre d'espèces observées (S)"""
    return np.sum(counts > 0)

def pielou_evenness(counts):
    """Équitabilité de Pielou (J' = H' / ln(S))"""
    richness = observed_richness(counts)
    if richness <= 1:
        return 0
    H = shannon_index(counts)
    return H / np.log(richness)

def chao1_richness(counts):
    """Estimation de richesse Chao1"""
    counts = counts[counts > 0]
    n = len(counts)  # Richesse observée
    
    f1 = np.sum(counts == 1)  # Nombre de singletons
    f2 = np.sum(counts == 2)  # Nombre de doubletons
    
    if f2 > 0:
        # Formule standard
        return n + (f1 * (f1 - 1)) / (2 * (f2 + 1))
    elif f1 > 0:
        # Cas spécial: pas de doubletons
        return n + (f1 * (f1 - 1)) / 2
    else:
        # Pas de singletons: la richesse est complètement échantillonnée
        return n

def inverse_simpson(counts):
    """Indice de Simpson inverse (1/D)"""
    counts = counts[counts > 0]
    if len(counts) == 0:
        return 0
    proportions = counts / counts.sum()
    D = np.sum(proportions ** 2)
    return 1 / D if D > 0 else 0

# Paramètres Snakemake
input_matrix = snakemake.input.matrix
output_file = snakemake.output.diversity

print(f"\n=== Calcul de la diversité alpha ===")

# Lire la matrice (niveau espèce)
try:
    df = pd.read_csv(input_matrix, sep='\t', index_col=0)
except Exception as e:
    print(f"ERREUR lors de la lecture de {input_matrix}: {e}", file=sys.stderr)
    sys.exit(1)

# Extraire les données d'abondance
abundance_data = df.drop('taxonomy_id', axis=1)

print(f"Matrice: {len(df)} espèces × {len(abundance_data.columns)} échantillons")
print(f"\nCalcul des indices de diversité pour chaque échantillon...")

# Calculer les indices pour chaque échantillon
results = []

for sample in abundance_data.columns:
    counts = abundance_data[sample].values
    
    # Calculer tous les indices
    obs_rich = observed_richness(counts)
    shannon = shannon_index(counts)
    simpson = simpson_index(counts)
    inv_simpson = inverse_simpson(counts)
    evenness = pielou_evenness(counts)
    chao1 = chao1_richness(counts)
    total_reads = int(counts.sum())
    
    results.append({
        'sample': sample,
        'observed_richness': int(obs_rich),
        'shannon_diversity': shannon,
        'simpson_diversity': simpson,
        'inverse_simpson': inv_simpson,
        'pielou_evenness': evenness,
        'chao1_richness': chao1,
        'total_reads': total_reads
    })
    
    print(f"  ✓ {sample}: S={obs_rich}, H'={shannon:.3f}, "
          f"Simpson={simpson:.3f}, J'={evenness:.3f}")

# Créer le DataFrame et sauvegarder
results_df = pd.DataFrame(results)
results_df = results_df.sort_values('sample')
results_df.to_csv(output_file, sep='\t', index=False, float_format='%.4f')

print(f"\n=== Résumé statistique ===")
print("\nStatistiques descriptives:")
print(results_df[['observed_richness', 'shannon_diversity', 'simpson_diversity', 
                   'pielou_evenness', 'chao1_richness']].describe().to_string())

print(f"\n✓ Indices de diversité sauvegardés: {output_file}")

# Interprétation rapide
print(f"\n=== Interprétation ===")
avg_shannon = results_df['shannon_diversity'].mean()
if avg_shannon < 2:
    print("Shannon moyen < 2: Diversité relativement faible")
elif avg_shannon < 3:
    print("Shannon moyen 2-3: Diversité modérée")
elif avg_shannon < 4:
    print("Shannon moyen 3-4: Diversité élevée")
else:
    print("Shannon moyen > 4: Très haute diversité")

avg_evenness = results_df['pielou_evenness'].mean()
if avg_evenness > 0.7:
    print(f"Équitabilité élevée ({avg_evenness:.2f}): distribution relativement uniforme")
elif avg_evenness > 0.5:
    print(f"Équitabilité modérée ({avg_evenness:.2f}): quelques espèces dominantes")
else:
    print(f"Équitabilité faible ({avg_evenness:.2f}): forte domination de quelques espèces")
