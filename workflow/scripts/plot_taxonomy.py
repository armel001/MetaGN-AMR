#!/usr/bin/env python3
"""
Generate all taxonomic figures in a single script:
  1. Stacked barplot - Genus level composition
  2. Stacked barplot - Species level composition
  3. Multi-panel barplot - Alpha diversity metrics
  4. Clustered heatmap - Top 30 species
  5. Grouped barplot - Priority pathogens
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend (for cluster/server)
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
import os
import sys

# ── Snakemake params ──────────────────────────────────────────────────────────
genus_matrix_path   = snakemake.input.genus_matrix
species_matrix_path = snakemake.input.species_matrix
alpha_div_path      = snakemake.input.alpha_div

out_composition_genus   = snakemake.output.composition_genus
out_composition_species = snakemake.output.composition_species
out_alpha_diversity     = snakemake.output.alpha_diversity
out_heatmap             = snakemake.output.heatmap_species
out_pathogens           = snakemake.output.pathogens

top_genera  = snakemake.params.get("top_genera",  15)
top_species = snakemake.params.get("top_species", 20)
top_heatmap = snakemake.params.get("top_heatmap", 30)
dpi         = snakemake.params.get("dpi",         300)
pathogens   = snakemake.params.get("pathogens",   [])

# ── Ensure output directory exists ───────────────────────────────────────────
os.makedirs(os.path.dirname(out_composition_genus), exist_ok=True)

# ── Color palettes ────────────────────────────────────────────────────────────
# Colorblind-friendly - categorical
PALETTE = [
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
    "#8CD17D", "#B6992D", "#499894", "#86BCB6", "#D37295",
    "#FABFD2", "#A0CBE8", "#FFBE7D", "#F1CE63", "#D4E6F1",
    "#C7C7C7"  # Others (always grey)
]

# Per-sample colors (consistent across all figures)
SAMPLE_COLORS = ["#4E79A7", "#F28E2B", "#59A14F", "#E15759"]

# Pathogen categories
PRIORITY_PATHOGENS = [
    "Escherichia coli", "Klebsiella pneumoniae",
    "Pseudomonas aeruginosa", "Clostridioides difficile",
    "Clostridium perfringens", "Acinetobacter baumannii",
    "Enterococcus faecium", "Bacteroides fragilis"
]
BENEFICIAL_SPECIES = [
    "Faecalibacterium prausnitzii",
    "Agathobacter rectalis",
    "Roseburia intestinalis"
]

# ── Helper functions ──────────────────────────────────────────────────────────

def load_matrix(path):
    """Load abundance matrix, remove taxonomy_id column"""
    df = pd.read_csv(path, sep='\t', index_col=0)
    if 'taxonomy_id' in df.columns:
        df = df.drop(columns=['taxonomy_id'])
    return df


def build_top_n_df(df, top_n, others_label_pct=None):
    """
    Keep top N taxa by mean abundance.
    Aggregate remaining taxa as 'Others'.
    Returns transposed DataFrame (samples as rows).
    """
    mean_abund   = df.mean(axis=1)
    top_taxa     = mean_abund.nlargest(top_n).index.tolist()
    other_taxa   = [t for t in df.index if t not in top_taxa]

    plot_df = df.loc[top_taxa].copy()

    if other_taxa:
        label = f"Others (<{100 / top_n:.0f}%)" if others_label_pct is None else others_label_pct
        others_row      = df.loc[other_taxa].sum(axis=0)
        others_row.name = label
        plot_df = pd.concat([plot_df, others_row.to_frame().T])

    return plot_df.T   # samples × taxa


def style_ax(ax):
    """Apply common axis styling"""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.grid(True, linestyle='--', alpha=0.35, color='grey')
    ax.set_axisbelow(True)


def stacked_barplot(df_matrix, top_n, title, out_path, dpi):
    """
    Stacked barplot of relative abundances.

    Parameters
    ----------
    df_matrix : pd.DataFrame
        Taxon × Sample, values in %
    top_n     : int
        Number of top taxa to show individually
    title     : str
    out_path  : str
    dpi       : int
    """
    print(f"\n  → Generating: {out_path}")

    plot_df = build_top_n_df(df_matrix, top_n)  # samples × taxa
    taxa    = plot_df.columns.tolist()
    colors  = PALETTE[:len(taxa) - 1] + ["#C7C7C7"]  # Last = Others (grey)

    fig, ax = plt.subplots(figsize=(10, 7))
    bottom  = np.zeros(len(plot_df))

    for i, taxon in enumerate(taxa):
        values = plot_df[taxon].values
        ax.bar(
            plot_df.index,
            values,
            bottom      = bottom,
            color       = colors[i],
            label       = taxon,
            width       = 0.60,
            edgecolor   = 'white',
            linewidth   = 0.5
        )
        bottom += values

    # Formatting
    ax.set_xlabel("Sample",                  fontsize=12, labelpad=10)
    ax.set_ylabel("Relative Abundance (%)",  fontsize=12, labelpad=10)
    ax.set_title(title,                      fontsize=14, fontweight='bold', pad=15)
    ax.set_ylim(0, 108)
    ax.set_xticks(range(len(plot_df)))
    ax.set_xticklabels(plot_df.index,        fontsize=11)
    ax.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, _: f"{x:.0f}%")
    )
    style_ax(ax)

    # Legend
    handles = [
        mpatches.Patch(color=colors[i], label=t)
        for i, t in enumerate(taxa)
    ]
    legend = ax.legend(
        handles        = handles,
        title          = "Taxa",
        title_fontsize = 10,
        fontsize       = 8,
        loc            = 'upper left',
        bbox_to_anchor = (1.01, 1),
        borderaxespad  = 0,
        frameon        = True,
        framealpha     = 0.9,
        edgecolor      = 'grey'
    )
    # Italicize species names
    for text in legend.get_texts():
        if not text.get_text().startswith("Others"):
            text.set_style('italic')

    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"     ✓ Saved ({dpi} dpi)")


# ── FIGURE 1 & 2 : Stacked barplots ─────────────────────────────────────────

print("\n" + "="*70)
print("FIGURE 1 - Composition: Genus Level")
print("="*70)

genus_df = load_matrix(genus_matrix_path)
print(f"  Matrix: {genus_df.shape[0]} genera × {genus_df.shape[1]} samples")

stacked_barplot(
    df_matrix = genus_df,
    top_n     = top_genera,
    title     = "Bacterial Composition: Genus Level",
    out_path  = out_composition_genus,
    dpi       = dpi
)

print("\n" + "="*70)
print("FIGURE 2 - Composition: Species Level")
print("="*70)

species_df = load_matrix(species_matrix_path)
print(f"  Matrix: {species_df.shape[0]} species × {species_df.shape[1]} samples")

stacked_barplot(
    df_matrix = species_df,
    top_n     = top_species,
    title     = "Bacterial Composition: Species Level",
    out_path  = out_composition_species,
    dpi       = dpi
)

# ── FIGURE 3 : Alpha diversity ───────────────────────────────────────────────

print("\n" + "="*70)
print("FIGURE 3 - Alpha Diversity")
print("="*70)
print(f"  → Generating: {out_alpha_diversity}")

alpha_df = pd.read_csv(alpha_div_path, sep='\t')

required = ['sample', 'shannon_diversity', 'simpson_diversity',
            'observed_richness', 'pielou_evenness']
missing_cols = [c for c in required if c not in alpha_df.columns]
if missing_cols:
    print(f"Error: Missing columns {missing_cols}", file=sys.stderr)
    sys.exit(1)

print(f"  Loaded: {len(alpha_df)} samples")

metrics = [
    {
        "col"   : "shannon_diversity",
        "title" : "Shannon Diversity (H')",
        "ylabel": "H'",
        "color" : "#4E79A7",
        "fmt"   : ".3f"
    },
    {
        "col"   : "simpson_diversity",
        "title" : "Simpson Diversity (1-D)",
        "ylabel": "1-D",
        "color" : "#F28E2B",
        "fmt"   : ".3f"
    },
    {
        "col"   : "observed_richness",
        "title" : "Observed Richness (S)",
        "ylabel": "Number of species",
        "color" : "#59A14F",
        "fmt"   : ".0f"
    },
    {
        "col"   : "pielou_evenness",
        "title" : "Pielou's Evenness (J')",
        "ylabel": "J'",
        "color" : "#E15759",
        "fmt"   : ".3f"
    }
]

fig = plt.figure(figsize=(13, 9))
fig.suptitle("Alpha Diversity Metrics - Hospital Wastewater",
             fontsize=15, fontweight='bold', y=1.01)
gs = gridspec.GridSpec(2, 2, hspace=0.50, wspace=0.38)

for idx, m in enumerate(metrics):
    ax      = fig.add_subplot(gs[idx // 2, idx % 2])
    values  = alpha_df[m["col"]].values
    samples = alpha_df['sample'].values
    ylim_max = max(values) * 1.25

    bars = ax.bar(
        samples,
        values,
        color     = m["color"],
        alpha     = 0.85,
        width     = 0.55,
        edgecolor = 'white',
        linewidth = 0.8
    )

    # Value labels
    for bar, val in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + ylim_max * 0.02,
            format(val, m["fmt"]),
            ha='center', va='bottom',
            fontsize=9, fontweight='bold'
        )

    ax.set_title(m["title"],  fontsize=12, fontweight='bold', pad=10)
    ax.set_ylabel(m["ylabel"], fontsize=10)
    ax.set_xlabel("Sample",   fontsize=10)
    ax.set_ylim(0, ylim_max)
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, fontsize=9, rotation=15, ha='right')
    style_ax(ax)

plt.tight_layout()
plt.savefig(out_alpha_diversity, dpi=dpi, bbox_inches='tight', facecolor='white')
plt.close()
print(f"     ✓ Saved ({dpi} dpi)")

# ── FIGURE 4 : Heatmap ───────────────────────────────────────────────────────

print("\n" + "="*70)
print("FIGURE 4 - Heatmap Top Species")
print("="*70)
print(f"  → Generating: {out_heatmap}")

# Select top N
mean_abund    = species_df.mean(axis=1)
top_taxa      = mean_abund.nlargest(top_heatmap).index.tolist()
heatmap_df    = species_df.loc[top_taxa]

print(f"  Top {top_heatmap} species selected from {species_df.shape[0]} total")

# Log10 transformation
log_df = np.log10(heatmap_df + 0.01)

# Hierarchical clustering on rows
row_link  = linkage(pdist(log_df.values, metric='euclidean'), method='ward')
row_order = leaves_list(row_link)
log_df_ordered    = log_df.iloc[row_order]
top_taxa_ordered  = [top_taxa[i] for i in row_order]

fig_h = max(10, top_heatmap * 0.38)
fig, ax = plt.subplots(figsize=(8, fig_h))

# Heatmap image
img = ax.imshow(
    log_df_ordered.values,
    aspect      = 'auto',
    cmap        = 'YlOrRd',
    interpolation = 'nearest'
)

# Colorbar
cbar = fig.colorbar(img, ax=ax, shrink=0.35, pad=0.02)
cbar.set_label('log₁₀(Relative abundance % + 0.01)', fontsize=9)
cbar.ax.tick_params(labelsize=8)

# X-axis (top)
ax.set_xticks(range(len(log_df_ordered.columns)))
ax.set_xticklabels(log_df_ordered.columns, fontsize=11)
ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_position('top')

# Y-axis with color annotation
ax.set_yticks(range(len(top_taxa_ordered)))
ytick_labels = ax.set_yticklabels(top_taxa_ordered, fontsize=8, style='italic')

for label, species in zip(ytick_labels, top_taxa_ordered):
    if any(p.lower() in species.lower() for p in PRIORITY_PATHOGENS):
        label.set_color('#CC0000')
        label.set_fontweight('bold')
    elif any(b.lower() in species.lower() for b in BENEFICIAL_SPECIES):
        label.set_color('#1a7a1a')
        label.set_fontweight('bold')

# Cell grid
ax.set_xticks(np.arange(-0.5, len(log_df_ordered.columns), 1), minor=True)
ax.set_yticks(np.arange(-0.5, len(top_taxa_ordered), 1),       minor=True)
ax.grid(which='minor', color='white', linewidth=0.6)
ax.tick_params(which='minor', bottom=False, left=False)

ax.set_title("Top Species Abundance Heatmap",
             fontsize=13, fontweight='bold', pad=40)

# Category legend
legend_elements = [
    Line2D([0], [0], color='#CC0000', lw=2.5, label='Priority pathogen'),
    Line2D([0], [0], color='#1a7a1a', lw=2.5, label='Beneficial species'),
    Line2D([0], [0], color='black',   lw=2.5, label='Commensal')
]
ax.legend(
    handles        = legend_elements,
    loc            = 'lower right',
    bbox_to_anchor = (1.0, -0.06),
    fontsize       = 8,
    frameon        = True,
    framealpha     = 0.9
)

plt.tight_layout()
plt.savefig(out_heatmap, dpi=dpi, bbox_inches='tight', facecolor='white')
plt.close()
print(f"     ✓ Saved ({dpi} dpi)")

# ── FIGURE 5 : Priority pathogens ────────────────────────────────────────────

print("\n" + "="*70)
print("FIGURE 5 - Priority Pathogens")
print("="*70)
print(f"  → Generating: {out_pathogens}")

found   = [p for p in pathogens if p in species_df.index]
missing = [p for p in pathogens if p not in species_df.index]

if missing:
    print(f"  ⚠️  Below filter threshold (not plotted):")
    for p in missing:
        print(f"       - {p}")

if not found:
    print("  Error: No priority pathogens found in matrix", file=sys.stderr)
    sys.exit(1)

print(f"  ✓ Plotting {len(found)} pathogens")

path_df   = species_df.loc[found]
samples   = species_df.columns.tolist()
n_samples = len(samples)
n_path    = len(found)

bar_width = 0.75 / n_samples
offsets   = np.linspace(
    -(n_samples - 1) / 2,
     (n_samples - 1) / 2,
    n_samples
) * bar_width

x = np.arange(n_path)

fig, ax = plt.subplots(figsize=(max(10, n_path * 1.6), 7))

for i, (sample, color) in enumerate(zip(samples, SAMPLE_COLORS[:n_samples])):
    values = path_df[sample].values
    bars   = ax.bar(
        x + offsets[i],
        values,
        width     = bar_width,
        color     = color,
        alpha     = 0.85,
        label     = sample,
        edgecolor = 'white',
        linewidth = 0.5
    )

    # Labels only if value > 0.1%
    for bar, val in zip(bars, values):
        if val > 0.1:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.08,
                f"{val:.1f}%",
                ha       = 'center',
                va       = 'bottom',
                fontsize = 7,
                rotation = 90
            )

# Alert threshold
max_val = path_df.values.max()
ax.axhline(y=1.0, color='red', linestyle='--', alpha=0.5, linewidth=1.2)
ax.text(n_path - 0.5, 1.05, 'Alert threshold (1%)',
        color='red', fontsize=8, alpha=0.7, ha='right')

# Formatting
ax.set_xlabel("Priority Pathogens",               fontsize=12, labelpad=10)
ax.set_ylabel("Relative Abundance (%)\n(Gb-normalized)", fontsize=12, labelpad=10)
ax.set_title("Priority Pathogen Abundance",        fontsize=14, fontweight='bold', pad=15)
ax.set_xticks(x)
ax.set_xticklabels(found, fontsize=9, style='italic', rotation=30, ha='right')
ax.set_ylim(0, max(max_val * 1.3, 2.0))
ax.legend(title="Hospital", fontsize=10, title_fontsize=10,
          loc='upper right', framealpha=0.9)
style_ax(ax)

plt.tight_layout()
plt.savefig(out_pathogens, dpi=dpi, bbox_inches='tight', facecolor='white')
plt.close()
print(f"     ✓ Saved ({dpi} dpi)")

# ── Final summary ─────────────────────────────────────────────────────────────

print("\n" + "="*70)
print("ALL FIGURES GENERATED SUCCESSFULLY")
print("="*70)
print(f"  1. {out_composition_genus}")
print(f"  2. {out_composition_species}")
print(f"  3. {out_alpha_diversity}")
print(f"  4. {out_heatmap}")
print(f"  5. {out_pathogens}")
print(f"\n  DPI    : {dpi}")
print(f"  Format : PNG")
print("="*70 + "\n")
