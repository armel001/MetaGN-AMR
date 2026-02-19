#!/usr/bin/env python3
"""
Generate ARG visualizations with professional journal-ready formatting
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from scipy.cluster.hierarchy import dendrogram, linkage
import json
import sys
from pathlib import Path

# Logger
log_file = snakemake.log[0]
log = open(log_file, 'w')
sys.stderr = sys.stdout = log

print("=" * 70)
print("ARG Visualization Generation - Journal Ready")
print("=" * 70)

# Create output directory
Path("results/figures").mkdir(parents=True, exist_ok=True)

# ============================================================================
# PROFESSIONAL STYLE CONFIGURATION
# ============================================================================

# Set style: 'ticks' is cleaner than 'whitegrid' for journals
sns.set_style("ticks")

# Define MASTER PALETTE (tab20 has 20 distinct colors)
master_palette = sns.color_palette("tab20", 25)

# High-resolution output settings
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 12
plt.rcParams['font.weight'] = 'normal'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'

print("\n✓ Professional style configured")
print("  - Style: ticks")
print("  - Palette: tab20 (25 colors)")
print("  - DPI: 300")

# ============================================================================
# LOAD DATA
# ============================================================================

print("\n[1/8] Loading data...")

arg_matrix_norm = pd.read_csv(snakemake.input.arg_matrix_normalized, sep='\t', index_col=0)
arg_matrix_counts = pd.read_csv(snakemake.input.arg_matrix_counts, sep='\t', index_col=0)
arg_matrix_presence = pd.read_csv(snakemake.input.arg_matrix_presence, sep='\t', index_col=0)
drug_class = pd.read_csv(snakemake.input.drug_class, sep='\t')
mechanism = pd.read_csv(snakemake.input.mechanism, sep='\t')
family = pd.read_csv(snakemake.input.family, sep='\t')
arg_counts = pd.read_csv(snakemake.input.arg_counts, sep='\t')

print(f"  ✓ Normalized matrix: {arg_matrix_norm.shape}")
print(f"  ✓ Counts matrix: {arg_matrix_counts.shape}")
print(f"  ✓ Drug classes: {len(drug_class)} rows")

# ============================================================================
# FIGURE 1: DRUG CLASSES DISTRIBUTION
# ============================================================================

print("\n[2/8] Generating Figure 1: Drug Classes Distribution...")

fig, ax = plt.subplots(figsize=(12, 6))

# Pivot for stacked bar
drug_pivot = drug_class.pivot(index='sample_id', columns='drug_class', values='relative_abundance').fillna(0)

# Use master palette
drug_pivot.plot(kind='bar', stacked=True, ax=ax, color=master_palette[:len(drug_pivot.columns)])

# Professional formatting
ax.set_title("Antibiotic Resistance Classes Distribution", fontsize=16, pad=20, weight='bold')
ax.set_ylabel("Relative Abundance (%)", fontsize=14, weight='bold')
ax.set_xlabel("")  # Remove x-label
ax.tick_params(axis='x', rotation=0, labelsize=12)  # Horizontal labels
ax.tick_params(axis='y', labelsize=12)
ax.margins(x=0.01)  # Remove white space

# Legend outside the plot
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left',
         title='Drug Class', frameon=False, fontsize=10)

# Remove top and right spines
sns.despine()

plt.tight_layout()
plt.savefig(snakemake.output.fig1, dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved: 1_drug_classes_distribution.png")

# ============================================================================
# FIGURE 2: ALPHA DIVERSITY BY SAMPLE (BOXPLOTS)
# ============================================================================

print("\n[3/8] Generating Figure 2: Alpha Diversity by Sample (Boxplots)...")

from scipy.stats import entropy

# Calculate Shannon diversity
def shannon_diversity(row):
    row = row[row > 0]
    if len(row) == 0:
        return 0
    proportions = row / row.sum()
    return entropy(proportions, base=np.e)

diversity_shannon = arg_matrix_norm.apply(shannon_diversity, axis=1)
richness = (arg_matrix_norm > 0).sum(axis=1)

# Create DataFrame
diversity_df = pd.DataFrame({
    'sample_id': diversity_shannon.index,
    'shannon': diversity_shannon.values,
    'richness': richness.values
})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# ============================================================================
# SHANNON INDEX - Boxplot per sample
# ============================================================================

# Prepare data for boxplot (each sample = one boxplot)
shannon_data = [[val] for val in diversity_df['shannon'].values]  # Each sample as a list

# Create boxplots
positions = np.arange(1, len(diversity_df) + 1)
bp1 = ax1.boxplot(shannon_data,
                  positions=positions,
                  widths=0.6,
                  patch_artist=True,
                  boxprops=dict(linewidth=1.5),
                  medianprops=dict(color='darkred', linewidth=2.5),
                  whiskerprops=dict(color='black', linewidth=1.5),
                  capprops=dict(color='black', linewidth=1.5),
                  showfliers=False)  # No outliers since we have single values

# Color each boxplot with master palette
for i, (patch, sample) in enumerate(zip(bp1['boxes'], diversity_df['sample_id'])):
    color = master_palette[i % len(master_palette)]
    patch.set_facecolor(color)
    patch.set_edgecolor('black')
    patch.set_alpha(0.7)

# Add individual points on top
ax1.scatter(positions, diversity_df['shannon'], 
           s=100, c=[master_palette[i % len(master_palette)] for i in range(len(diversity_df))],
           edgecolors='black', linewidth=1.5, zorder=3, alpha=0.9)

# Professional formatting
ax1.set_title("Alpha Diversity: Shannon Index", fontsize=16, pad=20, weight='bold')
ax1.set_ylabel("Shannon Index", fontsize=14, weight='bold')
ax1.set_xlabel("")
ax1.set_xticks(positions)
ax1.set_xticklabels(diversity_df['sample_id'], rotation=45, ha='right', fontsize=11)
ax1.tick_params(axis='y', labelsize=12)
ax1.margins(x=0.02, y=0.1)

# Add mean line
mean_shannon = diversity_df['shannon'].mean()
ax1.axhline(y=mean_shannon, color='gray', linestyle='--', linewidth=2, 
           alpha=0.5, label=f'Mean = {mean_shannon:.2f}')
ax1.legend(frameon=False, fontsize=10, loc='best')

# ============================================================================
# RICHNESS - Boxplot per sample
# ============================================================================

# Prepare data for boxplot
richness_data = [[val] for val in diversity_df['richness'].values]

# Create boxplots
bp2 = ax2.boxplot(richness_data,
                  positions=positions,
                  widths=0.6,
                  patch_artist=True,
                  boxprops=dict(linewidth=1.5),
                  medianprops=dict(color='darkred', linewidth=2.5),
                  whiskerprops=dict(color='black', linewidth=1.5),
                  capprops=dict(color='black', linewidth=1.5),
                  showfliers=False)

# Color each boxplot
for i, patch in enumerate(bp2['boxes']):
    color = master_palette[i % len(master_palette)]
    patch.set_facecolor(color)
    patch.set_edgecolor('black')
    patch.set_alpha(0.7)

# Add individual points
ax2.scatter(positions, diversity_df['richness'], 
           s=100, c=[master_palette[i % len(master_palette)] for i in range(len(diversity_df))],
           edgecolors='black', linewidth=1.5, zorder=3, alpha=0.9)

# Professional formatting
ax2.set_title("ARG Richness", fontsize=16, pad=20, weight='bold')
ax2.set_ylabel("Number of Unique ARGs", fontsize=14, weight='bold')
ax2.set_xlabel("")
ax2.set_xticks(positions)
ax2.set_xticklabels(diversity_df['sample_id'], rotation=45, ha='right', fontsize=11)
ax2.tick_params(axis='y', labelsize=12)
ax2.margins(x=0.02, y=0.1)

# Add mean line
mean_richness = diversity_df['richness'].mean()
ax2.axhline(y=mean_richness, color='gray', linestyle='--', linewidth=2, 
           alpha=0.5, label=f'Mean = {mean_richness:.0f}')
ax2.legend(frameon=False, fontsize=10, loc='best')

# Remove top and right spines
sns.despine(fig=fig)

plt.tight_layout()
plt.savefig(snakemake.output.fig2, dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved: 2_alpha_diversity.png")


# ============================================================================
# FIGURE 3: HEATMAP TOP 30 ARGs
# ============================================================================

print("\n[4/8] Generating Figure 3: Heatmap Top 30 ARGs...")

# Select top 30 most abundant ARGs
top_args = arg_matrix_norm.sum(axis=0).nlargest(30).index
arg_matrix_top = arg_matrix_norm[top_args]

# Log transform
arg_matrix_log = np.log10(arg_matrix_top + 1)

# Create clustermap
g = sns.clustermap(
    arg_matrix_log,
    cmap='YlOrRd',
    standard_scale=1,
    figsize=(14, 10),
    cbar_kws={'label': 'log₁₀(Normalized Copy Number + 1)'},
    xticklabels=True,
    yticklabels=True,
    linewidths=0.3,
    cbar_pos=(0.02, 0.83, 0.03, 0.15),
    dendrogram_ratio=(0.1, 0.1)
)

# Title
g.fig.suptitle('Top 30 ARGs Heatmap (Hierarchical Clustering)', 
               fontsize=16, weight='bold', x=0.5, y=0.98)

# Adjust labels
g.ax_heatmap.tick_params(labelsize=9)

plt.savefig(snakemake.output.fig3, dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved: 3_heatmap_top30.png")

# ============================================================================
# FIGURE 4: RESISTANCE MECHANISMS
# ============================================================================

print("\n[5/8] Generating Figure 4: Resistance Mechanisms...")

mechanism_total = mechanism.groupby('resistance_mechanism')['count'].sum().sort_values(ascending=True)

fig, ax = plt.subplots(figsize=(10, 7))

# Use master palette
colors_mech = [master_palette[i % len(master_palette)] for i in range(len(mechanism_total))]
bars = ax.barh(mechanism_total.index, mechanism_total.values, color=colors_mech, edgecolor='black', linewidth=0.7)

# Professional formatting
ax.set_title("Distribution of Antibiotic Resistance Mechanisms", fontsize=16, pad=20, weight='bold')
ax.set_xlabel("Number of ARG Observations", fontsize=14, weight='bold')
ax.set_ylabel("")  # Remove y-label for cleaner look
ax.tick_params(axis='both', labelsize=12)
ax.margins(y=0.01)

# Remove top and right spines
sns.despine()

plt.tight_layout()
plt.savefig(snakemake.output.fig4, dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved: 4_resistance_mechanisms.png")

# ============================================================================
# FIGURE 5: RAREFACTION CURVES
# ============================================================================

print("\n[6/8] Generating Figure 5: Rarefaction Curves...")

def rarefaction_curve(row, step=10):
    """Calculate rarefaction curve for a sample"""
    total = int(row.sum())
    
    # Edge case : pas assez de données
    if total == 0:
        return [0], [0]
    
    if total < step:
        return [total], [int((row > 0).sum())]
    
    counts = row[row > 0].values
    
    # Créer array d'observations
    observations = []
    for i, count in enumerate(counts):
        observations.extend([i] * int(count))
    
    observations = np.array(observations)
    
    # Générer les steps
    steps = list(range(step, total, step))
    
    # Toujours inclure le total
    if not steps or steps[-1] != total:
        steps.append(total)
    
    # Calculer la richesse pour chaque step
    richness_values = []
    for n in steps:
        np.random.seed(42)
        subsample = np.random.choice(observations, size=n, replace=False)
        richness_values.append(len(np.unique(subsample)))
    
    return steps, richness_values

fig, ax = plt.subplots(figsize=(10, 7))

for idx, (sample, row) in enumerate(arg_matrix_counts.iterrows()):
    try:
        steps, richness = rarefaction_curve(row, step=10)
        
        if len(steps) > 1:
            ax.plot(steps, richness, label=sample, 
                   color=master_palette[idx % len(master_palette)], 
                   linewidth=2.5, alpha=0.8)
        else:
            # Afficher comme point si un seul step
            ax.scatter(steps, richness, label=sample, 
                      color=master_palette[idx % len(master_palette)], 
                      s=100, zorder=3)
    
    except Exception as e:
        print(f"  ⚠ Warning: Could not compute rarefaction for {sample}: {e}")
        continue

# Professional formatting
ax.set_title("ARG Rarefaction Curves", fontsize=16, pad=20, weight='bold')
ax.set_xlabel("Number of ARG Observations", fontsize=14, weight='bold')
ax.set_ylabel("Number of Unique ARGs", fontsize=14, weight='bold')
ax.tick_params(axis='both', labelsize=12)
ax.margins(x=0.01)
ax.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)

# Legend outside
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left',
         frameon=False, fontsize=10)

sns.despine()
plt.tight_layout()
plt.savefig(snakemake.output.fig5, dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved: 5_rarefaction_curves.png")


# ============================================================================
# FIGURE 6: PCoA BETA DIVERSITY
# ============================================================================

print("\n[7/8] Generating Figure 6: PCoA Beta Diversity...")

from sklearn.metrics import pairwise_distances

# Calculate Bray-Curtis distances
def bray_curtis_distance(u, v):
    return np.sum(np.abs(u - v)) / np.sum(u + v)

distances = pairwise_distances(arg_matrix_norm.values, metric=bray_curtis_distance)

# MDS (equivalent to PCoA)
mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
coords = mds.fit_transform(distances)

fig, ax = plt.subplots(figsize=(9, 7))

# Use master palette
for idx, sample in enumerate(arg_matrix_norm.index):
    ax.scatter(coords[idx, 0], coords[idx, 1], s=150, 
              c=[master_palette[idx % len(master_palette)]], 
              edgecolors='black', linewidth=1.5, alpha=0.8, label=sample)

# Add sample labels
for i, sample in enumerate(arg_matrix_norm.index):
    ax.annotate(sample, (coords[i, 0], coords[i, 1]), 
                xytext=(8, 8), textcoords='offset points', fontsize=10, weight='bold')

# Professional formatting
ax.set_title("Principal Coordinates Analysis (PCoA)\nBray-Curtis Distance", 
            fontsize=16, pad=20, weight='bold')
ax.set_xlabel("PCo1", fontsize=14, weight='bold')
ax.set_ylabel("PCo2", fontsize=14, weight='bold')
ax.tick_params(axis='both', labelsize=12)
ax.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)

# Remove top and right spines
sns.despine()

plt.tight_layout()
plt.savefig(snakemake.output.fig6, dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved: 6_pcoa_beta_diversity.png")

# ============================================================================
# FIGURE 7: TOP 20 ARGs
# ============================================================================

print("\n[8/8] Generating Figure 7: Top 20 ARGs...")

top20_args = arg_counts.groupby('best_hit_aro')['count'].sum().nlargest(20).sort_values(ascending=True)

fig, ax = plt.subplots(figsize=(10, 8))

# Use master palette
colors_args = [master_palette[i % len(master_palette)] for i in range(len(top20_args))]
bars = ax.barh(top20_args.index, top20_args.values, color=colors_args, 
               edgecolor='black', linewidth=0.7)

# Professional formatting
ax.set_title("Top 20 Most Abundant ARGs", fontsize=16, pad=20, weight='bold')
ax.set_xlabel("Number of Observations", fontsize=14, weight='bold')
ax.set_ylabel("")  # Remove y-label
ax.tick_params(axis='both', labelsize=11)
ax.margins(y=0.01)

# Remove top and right spines
sns.despine()

plt.tight_layout()
plt.savefig(snakemake.output.fig7, dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved: 7_top20_args.png")

# ============================================================================
# GENERATE GOOGLE COLAB NOTEBOOK
# ============================================================================

print("\n[9/9] Generating Google Colab notebook...")

notebook_content = {
    "nbformat": 4,
    "nbformat_minor": 0,
    "metadata": {
        "colab": {
            "provenance": [],
            "toc_visible": True
        },
        "kernelspec": {
            "name": "python3",
            "display_name": "Python 3"
        },
        "language_info": {
            "name": "python"
        }
    },
    "cells": [
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "# ARG Visualization - Journal-Ready Figures\n",
                "\n",
                "**Pipeline**: MetagenAMR  \n",
                "**Analysis**: Antibiotic Resistance Genes  \n",
                "**Normalization**: Normalized copy number per Gb  \n",
                "**Style**: Professional journal-ready formatting\n",
                "\n",
                "---"
            ]
        },
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "## 1. Setup - Professional Style Configuration"
            ]
        },
        {
            "cell_type": "code",
            "metadata": {},
            "source": [
                "# Install required packages\n",
                "!pip install pandas numpy matplotlib seaborn scikit-learn scipy -q\n",
                "\n",
                "import pandas as pd\n",
                "import numpy as np\n",
                "import matplotlib.pyplot as plt\n",
                "import seaborn as sns\n",
                "from sklearn.manifold import MDS\n",
                "from sklearn.metrics import pairwise_distances\n",
                "from scipy.stats import entropy\n",
                "import warnings\n",
                "warnings.filterwarnings('ignore')\n",
                "\n",
                "# PROFESSIONAL STYLE CONFIGURATION\n",
                "# 'ticks' is cleaner than 'whitegrid' for journals\n",
                "sns.set_style('ticks')\n",
                "\n",
                "# Define MASTER PALETTE (tab20 has 20 distinct colors)\n",
                "master_palette = sns.color_palette('tab20', 25)\n",
                "\n",
                "# High-resolution output\n",
                "%config InlineBackend.figure_format = 'retina'\n",
                "plt.rcParams['figure.dpi'] = 150\n",
                "plt.rcParams['savefig.dpi'] = 300\n",
                "plt.rcParams['font.size'] = 12\n",
                "plt.rcParams['axes.labelweight'] = 'bold'\n",
                "plt.rcParams['axes.titleweight'] = 'bold'\n",
                "\n",
                "print('✓ Professional style configured!')\n",
                "print('  - Style: ticks')\n",
                "print('  - Palette: tab20 (25 colors)')\n",
                "print('  - DPI: 300 for saving')"
            ],
            "execution_count": None,
            "outputs": []
        },
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "## 2. Upload Data Files\n",
                "\n",
                "Upload these TSV files from your pipeline:\n",
                "- `arg_matrix_normalized.tsv`\n",
                "- `arg_matrix_counts.tsv`\n",
                "- `drug_class_abundance.tsv`\n",
                "- `mechanism_abundance.tsv`\n",
                "- `arg_counts.tsv`"
            ]
        },
        {
            "cell_type": "code",
            "metadata": {},
            "source": [
                "from google.colab import files\n",
                "uploaded = files.upload()\n",
                "print('\\n✓ Files uploaded successfully!')"
            ],
            "execution_count": None,
            "outputs": []
        },
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "## 3. Load Data"
            ]
        },
        {
            "cell_type": "code",
            "metadata": {},
            "source": [
                "arg_matrix_norm = pd.read_csv('arg_matrix_normalized.tsv', sep='\\t', index_col=0)\n",
                "arg_matrix_counts = pd.read_csv('arg_matrix_counts.tsv', sep='\\t', index_col=0)\n",
                "drug_class = pd.read_csv('drug_class_abundance.tsv', sep='\\t')\n",
                "mechanism = pd.read_csv('mechanism_abundance.tsv', sep='\\t')\n",
                "arg_counts = pd.read_csv('arg_counts.tsv', sep='\\t')\n",
                "\n",
                "print('✓ Data loaded!')\n",
                "print(f'Samples: {len(arg_matrix_norm)}')\n",
                "print(f'ARGs: {len(arg_matrix_norm.columns)}')"
            ],
            "execution_count": None,
            "outputs": []
        },
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "## 4. Figure 1: Drug Classes Distribution"
            ]
        },
        {
            "cell_type": "code",
            "metadata": {},
            "source": [
                "fig, ax = plt.subplots(figsize=(12, 6))\n",
                "\n",
                "drug_pivot = drug_class.pivot(index='sample_id', columns='drug_class', \n",
                "                               values='relative_abundance').fillna(0)\n",
                "drug_pivot.plot(kind='bar', stacked=True, ax=ax, \n",
                "                color=master_palette[:len(drug_pivot.columns)])\n",
                "\n",
                "# Professional formatting\n",
                "ax.set_title('Antibiotic Resistance Classes Distribution', \n",
                "            fontsize=16, pad=20, weight='bold')\n",
                "ax.set_ylabel('Relative Abundance (%)', fontsize=14, weight='bold')\n",
                "ax.set_xlabel('')\n",
                "ax.tick_params(axis='x', rotation=0, labelsize=12)\n",
                "ax.tick_params(axis='y', labelsize=12)\n",
                "ax.margins(x=0.01)\n",
                "\n",
                "# Legend outside\n",
                "ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left',\n",
                "         title='Drug Class', frameon=False, fontsize=10)\n",
                "\n",
                "sns.despine()\n",
                "plt.tight_layout()\n",
                "plt.show()"
            ],
            "execution_count": None,
            "outputs": []
        },
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "## 5. Figure 2: Alpha Diversity"
            ]
        },
        {
            "cell_type": "code",
            "metadata": {},
            "source": [
                "def shannon_diversity(row):\n",
                "    row = row[row > 0]\n",
                "    if len(row) == 0:\n",
                "        return 0\n",
                "    proportions = row / row.sum()\n",
                "    return entropy(proportions, base=np.e)\n",
                "\n",
                "diversity_shannon = arg_matrix_norm.apply(shannon_diversity, axis=1)\n",
                "richness = (arg_matrix_norm > 0).sum(axis=1)\n",
                "\n",
                "fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))\n",
                "\n",
                "# Shannon\n",
                "ax1.bar(diversity_shannon.index, diversity_shannon.values, \n",
                "       color=master_palette[0], edgecolor='black', linewidth=0.7)\n",
                "ax1.set_title('Alpha Diversity: Shannon Index', fontsize=16, pad=20, weight='bold')\n",
                "ax1.set_ylabel('Shannon Index', fontsize=14, weight='bold')\n",
                "ax1.set_xlabel('')\n",
                "ax1.tick_params(axis='x', rotation=45, labelsize=12)\n",
                "ax1.margins(x=0.01)\n",
                "\n",
                "# Richness\n",
                "ax2.bar(richness.index, richness.values, \n",
                "       color=master_palette[1], edgecolor='black', linewidth=0.7)\n",
                "ax2.set_title('ARG Richness', fontsize=16, pad=20, weight='bold')\n",
                "ax2.set_ylabel('Number of Unique ARGs', fontsize=14, weight='bold')\n",
                "ax2.set_xlabel('')\n",
                "ax2.tick_params(axis='x', rotation=45, labelsize=12)\n",
                "ax2.margins(x=0.01)\n",
                "\n",
                "sns.despine(fig=fig)\n",
                "plt.tight_layout()\n",
                "plt.show()"
            ],
            "execution_count": None,
            "outputs": []
        },
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "## 6. Continue with remaining figures...\n",
                "\n",
                "*(Add cells for Figures 3-7 following the same professional formatting pattern)*"
            ]
        }
    ]
}

with open(snakemake.output.colab_notebook, 'w') as f:
    json.dump(notebook_content, f, indent=2)

print(f"  ✓ Saved: ARG_Visualization_Colab.ipynb")

# ============================================================================
# SUMMARY
# ============================================================================

print("\n" + "=" * 70)
print("VISUALIZATION GENERATION COMPLETED")
print("=" * 70)
print("\nGenerated files:")
print("  ✓ 7 high-resolution figures (300 dpi)")
print("  ✓ 1 interactive Colab notebook")
print("\nStyle configuration:")
print("  - Theme: seaborn 'ticks' (journal-ready)")
print("  - Palette: tab20 (25 distinct colors)")
print("  - Legends: positioned outside plots")
print("  - Spines: top and right removed")
print("  - Labels: bold, professional sizing")
print("=" * 70 + "\n")

log.close()
