#!/usr/bin/env python3

import re
import sys
from pathlib import Path


def parse_nanostats(file_path):
    stats = {}
    with open(file_path, 'r') as f:
        for line in f:
            if ':' not in line:
                continue
            key, value = line.split(':', 1)
            key = key.strip()
            value = value.strip().replace(',', '')
            
            match = re.search(r'([\d\.]+)', value)
            if match:
                num = match.group(1)
                stats[key] = float(num) if '.' in num else int(num)
    return stats


def extract_key_stats(stats, sample_id):
    mappings = {
        'total_reads': ['Number of reads', 'Active channels'],
        'total_bases': ['Total bases', 'Number of bases'],
        'mean_length': ['Mean read length'],
        'median_length': ['Median read length'],
        'mean_quality': ['Mean read quality'],
        'median_quality': ['Median read quality'],
        'n50': ['Read length N50', 'N50']
    }
    
    result = {'sample_id': sample_id}
    for key, possible_names in mappings.items():
        result[key] = next((stats[name] for name in possible_names if name in stats), 'NA')
    
    return result


def main():
    input_file = snakemake.input.stats
    output_file = snakemake.output[0]
    log_file = snakemake.log[0]
    
    sample_id = Path(input_file).parts[1]
    
    log = open(log_file, 'w')
    sys.stderr = sys.stdout = log
    
    print(f"Extracting stats for {sample_id}")
    
    if not Path(input_file).exists():
        print(f"ERROR: {input_file} not found")
        sys.exit(1)
    
    all_stats = parse_nanostats(input_file)
    key_stats = extract_key_stats(all_stats, sample_id)
    
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w') as f:
        headers = ['sample_id', 'total_reads', 'total_bases', 'mean_length', 
                  'median_length', 'mean_quality', 'median_quality', 'n50']
        f.write('\t'.join(headers) + '\n')
        f.write('\t'.join(str(key_stats.get(h, 'NA')) for h in headers) + '\n')
    
    print(f"Stats extracted: {key_stats['total_reads']} reads")
    log.close()


if __name__ == "__main__":
    main()
