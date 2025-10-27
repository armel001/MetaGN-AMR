# scripts/plot_ARGs_by_gene.py

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

abricate_summary_file = sys.argv[1]
output_plot_file = sys.argv[2]
sample = sys.argv[3]

# Charger les données
df = pd.read_csv(abricate_summary_file, sep="\t")

# Filtrer et compter les gènes
gene_counts = df['GENE'].value_counts().reset_index()
gene_counts.columns = ['Gene', 'Count']

# Plot
plt.figure(figsize=(14, 10))
sns.barplot(data=gene_counts, x='Count', y='Gene', palette="viridis")
plt.title(f"Distribution of resistance genes in {sample}")
plt.xlabel("Count")
plt.ylabel("Gene")
plt.tight_layout()

# Sauvegarde
plt.savefig(output_plot_file)
