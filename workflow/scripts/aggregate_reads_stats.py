#!/usr/bin/env python3
"""
Aggregate read statistics across all samples and steps
Computes data loss at each pipeline step
"""

import pandas as pd
import sys
from pathlib import Path


def main():
    # Snakemake parameters
    raw_files = snakemake.input.raw
    noh_files = snakemake.input.noh
    filtered_files = snakemake.input.filtered
    output_file = snakemake.output.aggregated
    summary_file = snakemake.output.summary
    sample_names = snakemake.params.samples
    log_file = snakemake.log[0]

    # Logger
    log = open(log_file, 'w')
    sys.stderr = sys.stdout = log

    print("=" * 70)
    print("Aggregate Read Statistics")
    print("=" * 70)

    # Load all stats files
    all_files = list(raw_files) + list(noh_files) + list(filtered_files)
    all_data = []

    for f in all_files:
        if Path(f).exists():
            df = pd.read_csv(f, sep='\t')
            all_data.append(df)
        else:
            print(f"  ⚠ WARNING: {f} not found")

    if not all_data:
        print("ERROR: No stats files found")
        sys.exit(1)

    # Combine all stats
    df_all = pd.concat(all_data, ignore_index=True)

    # Define step order
    step_order = ['raw', 'host_depleted', 'filtered']
    df_all['step'] = pd.Categorical(df_all['step'], categories=step_order, ordered=True)
    df_all = df_all.sort_values(['sample_id', 'step'])

    print(f"\n  Samples: {df_all['sample_id'].nunique()}")
    print(f"  Steps:   {df_all['step'].nunique()}")

    # Compute data loss between steps
    print(f"\nComputing data loss...")
    loss_rows = []

    for sample in sample_names:
        sample_data = df_all[df_all['sample_id'] == sample].set_index('step')

        for i, step in enumerate(step_order[1:], 1):
            prev_step = step_order[i - 1]

            if prev_step in sample_data.index and step in sample_data.index:
                reads_before = sample_data.loc[prev_step, 'n_reads']
                reads_after = sample_data.loc[step, 'n_reads']
                bases_before = sample_data.loc[prev_step, 'total_bases']
                bases_after = sample_data.loc[step, 'total_bases']

                reads_lost = reads_before - reads_after
                bases_lost = bases_before - bases_after
                pct_reads_retained = (reads_after / reads_before * 100) if reads_before > 0 else 0
                pct_bases_retained = (bases_after / bases_before * 100) if bases_before > 0 else 0

                loss_rows.append({
                    'sample_id': sample,
                    'step_from': prev_step,
                    'step_to': step,
                    'reads_before': reads_before,
                    'reads_after': reads_after,
                    'reads_lost': reads_lost,
                    'pct_reads_retained': round(pct_reads_retained, 2),
                    'pct_reads_lost': round(100 - pct_reads_retained, 2),
                    'bases_before': bases_before,
                    'bases_after': bases_after,
                    'bases_lost': bases_lost,
                    'pct_bases_retained': round(pct_bases_retained, 2),
                    'pct_bases_lost': round(100 - pct_bases_retained, 2)
                })

    df_loss = pd.DataFrame(loss_rows)

    # Save aggregated stats
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    df_all.to_csv(output_file, sep='\t', index=False)
    print(f"\n  ✓ Aggregated stats saved: {output_file}")

    # Save summary report
    with open(summary_file, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("READ STATISTICS SUMMARY\n")
        f.write("=" * 70 + "\n\n")

        f.write(f"{'Sample':<15} {'Step':<15} {'Reads':>12} {'Bases (Gb)':>12} {'Mean Len':>10} {'Quality':>10}\n")
        f.write("-" * 70 + "\n")

        for _, row in df_all.iterrows():
            f.write(f"{row['sample_id']:<15} {row['step']:<15} "
                   f"{row['n_reads']:>12,} "
                   f"{row['total_bases']/1e9:>12.3f} "
                   f"{row['mean_length']:>10.1f} "
                   f"{row['mean_quality']:>10.2f}\n")

        f.write("\n\n" + "=" * 70 + "\n")
        f.write("DATA LOSS SUMMARY\n")
        f.write("=" * 70 + "\n\n")

        f.write(f"{'Sample':<15} {'Transition':<30} {'Reads Lost':>12} {'% Retained':>12}\n")
        f.write("-" * 70 + "\n")

        for _, row in df_loss.iterrows():
            transition = f"{row['step_from']} → {row['step_to']}"
            f.write(f"{row['sample_id']:<15} {transition:<30} "
                   f"{row['reads_lost']:>12,} "
                   f"{row['pct_reads_retained']:>11.1f}%\n")

    print(f"  ✓ Summary saved: {summary_file}")

    # Print data loss to log
    print(f"\n{'=' * 70}")
    print(f"DATA LOSS REPORT")
    print(f"{'=' * 70}")

    for sample in sample_names:
        print(f"\n  {sample}:")
        sample_loss = df_loss[df_loss['sample_id'] == sample]

        for _, row in sample_loss.iterrows():
            print(f"    {row['step_from']} → {row['step_to']}:")
            print(f"      Reads:  {row['reads_before']:>10,} → {row['reads_after']:>10,} "
                  f"({row['pct_reads_retained']:.1f}% retained, "
                  f"{row['pct_reads_lost']:.1f}% lost)")
            print(f"      Bases:  {row['bases_before']/1e9:>10.3f} → "
                  f"{row['bases_after']/1e9:>10.3f} Gb "
                  f"({row['pct_bases_retained']:.1f}% retained)")

    print(f"\n{'=' * 70}")
    print("AGGREGATION COMPLETED")
    print(f"{'=' * 70}\n")

    log.close()


if __name__ == "__main__":
    main()
