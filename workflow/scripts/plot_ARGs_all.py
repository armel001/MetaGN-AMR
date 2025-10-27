# scripts/plot_ARGs_all.py

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# sys.argv[0] = script name
# sys.argv[1] = output plot file
# sys.argv[2] = output TSV file
# sys.argv[3:] = input files

output_file = sys.argv[1]
tsv_file = sys.argv[2]
input_files = sys.argv[3:]

print("Generating stacked barplot...")
print("Output PNG:", output_file)
print("Output TSV:", tsv_file)
print("Input files:", input_files)

dfs = []

for summary_file in input_files:
    sample_name = Path(summary_file).stem.split("_")[0]
    try:
        df = pd.read_csv(summary_file, sep="\t")
        df.columns = df.columns.str.strip().str.upper()
        if "GENE" not in df.columns:
            print(f"Skipping {summary_file}: no GENE column")
            continue
        df["Sample"] = sample_name
        dfs.append(df[["GENE", "Sample"]])
    except Exception as e:
        print(f"Error reading {summary_file}: {e}")

if not dfs:
    print("No valid input found.")
    sys.exit(1)

all_data = pd.concat(dfs)

# Count gene occurrences per sample
gene_counts = all_data.groupby(["Sample", "GENE"]).size().reset_index(name="Count")

# Create pivot table
pivot_df = gene_counts.pivot(index="Sample", columns="GENE", values="Count").fillna(0)

# Plot
pivot_df.plot(kind="bar", stacked=True, figsize=(20, 16), colormap="tab20")
plt.title("Distribution of resistance genes across samples")
plt.ylabel("Gene Count")
plt.xlabel("Sample")
plt.xticks(rotation=45, ha='right')
plt.legend(title="Resistance Genes", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.savefig(output_file, bbox_inches="tight")

# Save table
print("Saving TSV to:", tsv_file)
pivot_df.to_csv(tsv_file, sep="\t")
