#!/usr/bin/env python3
"""
Agrege les resultats RGI de tous les echantillons en un seul fichier
"""

import pandas as pd
import sys
from pathlib import Path


def main():
    # Récuperer les parametres
    input_files = snakemake.input.selected
    output_file = snakemake.output.aggregated
    summary_file = snakemake.output.summary
    sample_names = snakemake.params.samples
    log_file = snakemake.log[0]
    
    # Configurer le logging
    log = open(log_file, 'w')
    sys.stderr = sys.stdout = log
    
    print("=" * 70)
    print("RGI Results Aggregation")
    print("=" * 70)
    print(f"\nNumber of samples to aggregate: {len(sample_names)}")
    print(f"Output file: {output_file}")
    print("\n" + "-" * 70)
    
    # Liste pour stocker les dataframes
    all_data = []
    sample_stats = []
    
    # Charger chaque fichier
    print(f"\nLoading individual sample files...")
    for i, (file_path, sample_name) in enumerate(zip(input_files, sample_names), 1):
        print(f"\n[{i}/{len(sample_names)}] Processing: {sample_name}")
        
        # Verifier que le fichier existe
        if not Path(file_path).exists():
            print(f"  ⚠ WARNING: File not found, skipping: {file_path}")
            continue
        
        try:
            # Charger le fichier
            df = pd.read_csv(file_path, sep='\t', low_memory=False)
            
            # Ajouter la colonne sample_id
            df['sample_id'] = sample_name
            
            # Statistiques du fichier
            n_rows = len(df)
            n_unique_args = df['best_hit_aro'].nunique() if 'best_hit_aro' in df.columns else 0
            
            print(f"  ✓ Loaded: {n_rows:,} rows, {n_unique_args} unique ARGs")
            
            # Stocker les stats
            sample_stats.append({
                'sample_id': sample_name,
                'n_args': n_rows,
                'n_unique_args': n_unique_args
            })
            
            # Ajouter au dataset global
            all_data.append(df)
            
        except Exception as e:
            print(f"  ✗ ERROR loading file: {str(e)}")
            continue
    
    # Verifier qu'on a au moins un échantillon
    if len(all_data) == 0:
        print(f"\nERROR: No valid sample files found!")
        sys.exit(1)
    
    print(f"\n" + "-" * 70)
    print(f"Successfully loaded {len(all_data)} sample(s)")
    
    # Concatener tous les dataframes
    print(f"\nCombining all samples...")
    df_combined = pd.concat(all_data, ignore_index=True)
    
    print(f"  ✓ Combined dataset created")
    print(f"  - Total rows: {len(df_combined):,}")
    print(f"  - Total columns: {len(df_combined.columns)}")
    
    # Reorganiser les colonnes pour mettre sample_id en premier
    cols = df_combined.columns.tolist()
    if 'sample_id' in cols:
        cols.remove('sample_id')
        cols = ['sample_id'] + cols
        df_combined = df_combined[cols]
    
    # Statistiques globales
    print(f"\n" + "-" * 70)
    print(f"Global Statistics:")
    print(f"  - Total ARG observations: {len(df_combined):,}")
    print(f"  - Unique samples: {df_combined['sample_id'].nunique()}")
    
    if 'best_hit_aro' in df_combined.columns:
        print(f"  - Unique ARGs detected: {df_combined['best_hit_aro'].nunique():,}")
        
        # Top 10 ARGs les plus fréquents
        print(f"\n  Top 10 most frequent ARGs:")
        top_args = df_combined['best_hit_aro'].value_counts().head(10)
        for j, (arg, count) in enumerate(top_args.items(), 1):
            percent = 100 * count / len(df_combined)
            print(f"    {j:2d}. {arg}: {count:,} ({percent:.1f}%)")
    
    if 'drug_class' in df_combined.columns:
        print(f"\n  Unique drug classes: {df_combined['drug_class'].nunique()}")
        
        # Top 10 classes d'antibiotiques
        print(f"\n  Top 10 drug classes:")
        # Gérer les classes multiples séparées par "; "
        drug_classes_expanded = df_combined['drug_class'].str.split('; ').explode()
        top_classes = drug_classes_expanded.value_counts().head(10)
        for j, (drug_class, count) in enumerate(top_classes.items(), 1):
            percent = 100 * count / len(drug_classes_expanded)
            print(f"    {j:2d}. {drug_class}: {count:,} ({percent:.1f}%)")
    
    if 'resistance_mechanism' in df_combined.columns:
        print(f"\n  Unique resistance mechanisms: {df_combined['resistance_mechanism'].nunique()}")
        
        print(f"\n  Resistance mechanisms distribution:")
        mech_counts = df_combined['resistance_mechanism'].value_counts().head(10)
        for mech, count in mech_counts.items():
            percent = 100 * count / len(df_combined)
            print(f"    • {mech}: {count:,} ({percent:.1f}%)")
    
    # Creer le repertoire de sortie
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Sauvegarder le fichier agrege
    print(f"\n" + "-" * 70)
    print(f"Saving aggregated data...")
    try:
        df_combined.to_csv(output_file, sep='\t', index=False)
        file_size = output_path.stat().st_size / 1024 / 1024  # MB
        print(f"  ✓ Saved to: {output_file}")
        print(f"  ✓ File size: {file_size:.2f} MB")
    except Exception as e:
        print(f"\nERROR: Failed to save aggregated file")
        print(f"  {str(e)}")
        sys.exit(1)
    
    # Creer le fichier de resume
    print(f"\nCreating summary file...")
    try:
        with open(summary_file, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("RGI Aggregation Summary\n")
            f.write("=" * 70 + "\n\n")
            
            f.write(f"Total samples processed: {len(all_data)}\n")
            f.write(f"Total ARG observations: {len(df_combined):,}\n")
            f.write(f"Unique ARGs detected: {df_combined['best_hit_aro'].nunique():,}\n\n")
            
            f.write("-" * 70 + "\n")
            f.write("Per-sample statistics:\n")
            f.write("-" * 70 + "\n")
            f.write(f"{'Sample':<20} {'N_ARGs':>10} {'Unique_ARGs':>15}\n")
            f.write("-" * 70 + "\n")
            
            for stat in sample_stats:
                f.write(f"{stat['sample_id']:<20} {stat['n_args']:>10,} {stat['n_unique_args']:>15}\n")
            
            f.write("\n" + "=" * 70 + "\n")
            
            # Statistiques par echantillon
            f.write("\nDetailed statistics by sample:\n")
            f.write("=" * 70 + "\n")
            for sample in sample_names:
                sample_data = df_combined[df_combined['sample_id'] == sample]
                f.write(f"\n{sample}:\n")
                f.write(f"  - Total ARGs: {len(sample_data):,}\n")
                
                if 'best_hit_aro' in sample_data.columns:
                    f.write(f"  - Unique ARGs: {sample_data['best_hit_aro'].nunique()}\n")
                
                if 'drug_class' in sample_data.columns:
                    drug_classes = sample_data['drug_class'].str.split('; ').explode()
                    top_3_classes = drug_classes.value_counts().head(3)
                    f.write(f"  - Top 3 drug classes:\n")
                    for drug_class, count in top_3_classes.items():
                        f.write(f"      • {drug_class}: {count}\n")
        
        print(f"  ✓ Summary saved to: {summary_file}")
        
    except Exception as e:
        print(f"\nWARNING: Failed to save summary file")
        print(f"  {str(e)}")
    
    # Resume final
    print("\n" + "=" * 70)
    print("AGGREGATION COMPLETED SUCCESSFULLY")
    print("=" * 70)
    print(f"  Samples:      {len(all_data)}")
    print(f"  Total ARGs:   {len(df_combined):,}")
    print(f"  Unique ARGs:  {df_combined['best_hit_aro'].nunique():,}")
    print(f"  Output size:  {file_size:.2f} MB")
    print("=" * 70 + "\n")
    
    log.close()


if __name__ == "__main__":
    main()
