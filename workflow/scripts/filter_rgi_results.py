#!/usr/bin/env python3

import pandas as pd
import sys
from pathlib import Path


def main():
    input_file = snakemake.input.rgi
    output_file = snakemake.output.filtered
    min_identity = snakemake.params.min_identity
    min_coverage = snakemake.params.min_coverage
    model_types = snakemake.params.model_types
    log_file = snakemake.log[0]
    
    log = open(log_file, 'w')
    sys.stderr = sys.stdout = log
    
    print("=" * 70)
    print("RGI Results Quality Filtering")
    print("=" * 70)
    print(f"\nInput:  {input_file}")
    print(f"Output: {output_file}")
    print(f"\nParameters:")
    print(f"  - Min identity:  {min_identity}%")
    print(f"  - Min coverage:  {min_coverage}%")
    print(f"  - Model types:   {', '.join(model_types)}")
    print("\n" + "-" * 70)
    
    if not Path(input_file).exists():
        print(f"\nERROR: Input file not found: {input_file}")
        sys.exit(1)
    
    print(f"\nLoading RGI results...")
    df = pd.read_csv(input_file, sep='\t', low_memory=False)
    print(f"✓ Loaded {len(df):,} rows")
    print(f"✓ Found {len(df.columns)} columns")
    
    print(f"\nInitial statistics:")
    print(f"  - Total hits: {len(df):,}")
    
    if 'Best_Hit_ARO' in df.columns:
        print(f"  - Unique ARGs: {df['Best_Hit_ARO'].nunique():,}")
    
    if 'Best_Identities' in df.columns:
        print(f"  - Identity range: {df['Best_Identities'].min():.1f}% - {df['Best_Identities'].max():.1f}%")
        print(f"  - Identity mean: {df['Best_Identities'].mean():.1f}%")
    
    if 'Percentage Length of Reference Sequence' in df.columns:
        print(f"  - Coverage range: {df['Percentage Length of Reference Sequence'].min():.1f}% - {df['Percentage Length of Reference Sequence'].max():.1f}%")
        print(f"  - Coverage mean: {df['Percentage Length of Reference Sequence'].mean():.1f}%")
    
    if 'Model_type' in df.columns:
        print(f"\n  Model types found:")
        for model, count in df['Model_type'].value_counts().items():
            print(f"    • {model}: {count:,} ({100*count/len(df):.1f}%)")
    
    print(f"\n{'Applying filters:'}")
    df_filtered = df.copy()
    
    # Filtre 1: Identité
    if 'Best_Identities' in df.columns:
        before = len(df_filtered)
        df_filtered = df_filtered[df_filtered['Best_Identities'] >= min_identity]
        removed = before - len(df_filtered)
        print(f"  1. Identity ≥ {min_identity}%")
        print(f"     → Removed: {removed:,} hits ({100*removed/before:.1f}%)")
        print(f"     → Retained: {len(df_filtered):,} hits")
    
    # Filtre 2: Couverture
    if 'Percentage Length of Reference Sequence' in df.columns:
        before = len(df_filtered)
        if before > 0:  # ← FIX: Éviter division par zéro
            df_filtered = df_filtered[df_filtered['Percentage Length of Reference Sequence'] >= min_coverage]
            removed = before - len(df_filtered)
            print(f"  2. Coverage ≥ {min_coverage}%")
            print(f"     → Removed: {removed:,} hits ({100*removed/before:.1f}%)")
            print(f"     → Retained: {len(df_filtered):,} hits")
        else:
            print(f"  2. Coverage ≥ {min_coverage}%")
            print(f"     → Skipped (no hits remaining after previous filter)")
    
    # Filtre 3: Type de modèle
    if 'Model_type' in df.columns and model_types:
        before = len(df_filtered)
        if before > 0:  # ← FIX: Éviter division par zéro
            df_filtered = df_filtered[df_filtered['Model_type'].isin(model_types)]
            removed = before - len(df_filtered)
            print(f"  3. Model type filter")
            print(f"     → Removed: {removed:,} hits ({100*removed/before:.1f}%)")
            print(f"     → Retained: {len(df_filtered):,} hits")
        else:
            print(f"  3. Model type filter")
            print(f"     → Skipped (no hits remaining after previous filters)")
    
    print(f"\n{'Final statistics:'}")
    print(f"  - Total hits retained: {len(df_filtered):,}")
    
    if len(df) > 0:
        print(f"  - Percentage retained: {100*len(df_filtered)/len(df):.2f}%")
    
    if 'Best_Hit_ARO' in df_filtered.columns and len(df_filtered) > 0:
        print(f"  - Unique ARGs: {df_filtered['Best_Hit_ARO'].nunique():,}")
    
    # Nettoyer les noms de colonnes
    print(f"\nCleaning column names...")
    df_filtered.columns = (df_filtered.columns
                          .str.replace(' ', '_')
                          .str.replace('(', '')
                          .str.replace(')', '')
                          .str.replace('%', 'percent')
                          .str.lower())
    print(f"  ✓ Standardized {len(df_filtered.columns)} column names")
    
    # Créer le répertoire de sortie
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Sauvegarder (même si vide)
    print(f"\nSaving filtered results...")
    df_filtered.to_csv(output_file, sep='\t', index=False)
    file_size = output_path.stat().st_size / 1024
    print(f"  ✓ Saved to: {output_file}")
    print(f"  ✓ File size: {file_size:.2f} KB")
    
    # Warning si aucun résultat
    if len(df_filtered) == 0:
        print(f"\n{'!' * 70}")
        print(f"WARNING: No ARGs passed quality filters!")
        print(f"{'!' * 70}")
        print(f"\nConsider adjusting filter thresholds in config.yaml:")
        print(f"  - Current min_identity: {min_identity}%")
        print(f"  - Current min_coverage: {min_coverage}%")
        print(f"\nFor this sample, identities range from {df['Best_Identities'].min():.1f}% to {df['Best_Identities'].max():.1f}%")
        print(f"You might want to lower min_identity to ~40-50% for ONT data")
    
    print("\n" + "=" * 70)
    print("FILTERING COMPLETED")
    print("=" * 70)
    print(f"  Input:  {len(df):,} hits")
    print(f"  Output: {len(df_filtered):,} hits ({100*len(df_filtered)/len(df):.1f}% retained)" if len(df) > 0 else "  Output: 0 hits")
    print("=" * 70 + "\n")
    
    log.close()


if __name__ == "__main__":
    main()
