#!/usr/bin/env python3
"""
Extract read statistics from a FASTQ file
Works for raw, host-depleted, and filtered reads
"""

import gzip
import sys
import re
from pathlib import Path


def get_step_name(file_path):
    """Infer pipeline step from file path"""
    if "mythesis-data" in file_path or "data" in file_path:
        return "raw"
    elif "noh" in file_path:
        return "host_depleted"
    elif "fp_trimmed" in file_path or "fastplong" in file_path:
        return "filtered"
    return "unknown"


def parse_fastq_stats(fastq_file):
    """
    Parse FASTQ file and compute statistics
    Returns dict with read stats
    """
    n_reads = 0
    total_bases = 0
    lengths = []
    qualities = []

    opener = gzip.open if str(fastq_file).endswith('.gz') else open
    mode = 'rt'

    with opener(fastq_file, mode) as f:
        while True:
            header = f.readline().strip()
            if not header:
                break

            seq    = f.readline().strip()
            plus   = f.readline().strip()
            qual   = f.readline().strip()

            if not header.startswith('@'):
                continue

            read_len = len(seq)
            n_reads += 1
            total_bases += read_len
            lengths.append(read_len)

            # Mean quality score (Phred)
            if qual:
                mean_qual = sum(ord(c) - 33 for c in qual) / len(qual)
                qualities.append(mean_qual)

    if n_reads == 0:
        return None

    lengths.sort()

    # N50
    n50 = 0
    cumsum = 0
    target = total_bases / 2
    for l in sorted(lengths, reverse=True):
        cumsum += l
        if cumsum >= target:
            n50 = l
            break

    return {
        'n_reads': n_reads,
        'total_bases': total_bases,
        'mean_length': round(total_bases / n_reads, 2),
        'median_length': lengths[len(lengths) // 2],
        'min_length': lengths[0],
        'max_length': lengths[-1],
        'n50': n50,
        'mean_quality': round(sum(qualities) / len(qualities), 2) if qualities else 0
    }


def main():
    # Snakemake parameters
    # Detect which input key is available (raw or fastq)
    if hasattr(snakemake.input, 'raw'):
        input_file = snakemake.input.raw
    else:
        input_file = snakemake.input.fastq

    output_file = snakemake.output.stats
    log_file = snakemake.log[0]

    # Extract sample name and step
    sample_id = Path(output_file).parts[1]
    step = get_step_name(input_file)

    # Logger
    log = open(log_file, 'w')
    sys.stderr = sys.stdout = log

    print("=" * 70)
    print(f"Read Statistics Extraction")
    print("=" * 70)
    print(f"\n  Sample:  {sample_id}")
    print(f"  Step:    {step}")
    print(f"  Input:   {input_file}")
    print(f"  Output:  {output_file}")
    print("\n" + "-" * 70)

    # Check input file
    if not Path(input_file).exists():
        print(f"\nERROR: File not found: {input_file}")
        sys.exit(1)

    file_size = Path(input_file).stat().st_size / 1024 / 1024
    print(f"\n  File size: {file_size:.2f} MB")

    # Parse FASTQ
    print(f"\n  Parsing FASTQ...")
    stats = parse_fastq_stats(input_file)

    if stats is None:
        print(f"\n  WARNING: No reads found in {input_file}")
        stats = {k: 0 for k in ['n_reads', 'total_bases', 'mean_length', 
                                  'median_length', 'min_length', 'max_length', 
                                  'n50', 'mean_quality']}

    # Print stats
    print(f"\n  Results:")
    print(f"    • Reads:          {stats['n_reads']:,}")
    print(f"    • Total bases:    {stats['total_bases']:,} ({stats['total_bases']/1e9:.3f} Gb)")
    print(f"    • Mean length:    {stats['mean_length']:,.1f} bp")
    print(f"    • Median length:  {stats['median_length']:,} bp")
    print(f"    • Min length:     {stats['min_length']:,} bp")
    print(f"    • Max length:     {stats['max_length']:,} bp")
    print(f"    • N50:            {stats['n50']:,} bp")
    print(f"    • Mean quality:   {stats['mean_quality']:.2f}")

    # Save output
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)

    headers = ['sample_id', 'step', 'n_reads', 'total_bases', 'mean_length',
               'median_length', 'min_length', 'max_length', 'n50', 'mean_quality']

    values = [sample_id, step,
              stats['n_reads'], stats['total_bases'],
              stats['mean_length'], stats['median_length'],
              stats['min_length'], stats['max_length'],
              stats['n50'], stats['mean_quality']]

    with open(output_file, 'w') as f:
        f.write('\t'.join(headers) + '\n')
        f.write('\t'.join(str(v) for v in values) + '\n')

    print(f"\n  ✓ Saved: {output_file}")
    print("\n" + "=" * 70)
    print("COMPLETED")
    print("=" * 70 + "\n")

    log.close()


if __name__ == "__main__":
    main()
