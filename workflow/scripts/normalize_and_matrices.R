#!/usr/bin/env Rscript

library(tidyverse)

# Paramètres Snakemake
input_rgi <- snakemake@input$rgi_aggregated
input_stats <- snakemake@input$sequencing_stats
log_file <- snakemake@log[[1]]

# Outputs
output_counts <- snakemake@output$arg_counts
output_relative <- snakemake@output$arg_relative
output_normalized <- snakemake@output$arg_normalized  # ← Changé de rpm à normalized
output_presence <- snakemake@output$arg_presence
output_matrix_counts <- snakemake@output$arg_matrix_counts
output_matrix_relative <- snakemake@output$arg_matrix_relative
output_matrix_normalized <- snakemake@output$arg_matrix_normalized  # ← Changé
output_matrix_presence <- snakemake@output$arg_matrix_presence
output_drug_class <- snakemake@output$drug_class_abundance
output_mechanism <- snakemake@output$mechanism_abundance
output_family <- snakemake@output$family_abundance

# Logger
log <- file(log_file, open = "wt")
sink(log)
sink(log, type = "message")

cat(rep("=", 70), "\n", sep = "")
cat("R Analysis: Normalization and Matrices\n")
cat(rep("=", 70), "\n", sep = "")

# ----------------------------------------------------------------------------
# 1. CHARGEMENT DES DONNÉES
# ----------------------------------------------------------------------------

cat("\n[1/7] Loading data...\n")

rgi_all <- read_tsv(input_rgi, col_types = cols(), show_col_types = FALSE)
sequencing_stats <- read_tsv(input_stats, col_types = cols(), show_col_types = FALSE)

cat("  ✓ RGI: ", nrow(rgi_all), " observations\n", sep = "")
cat("  ✓ Stats: ", nrow(sequencing_stats), " samples\n", sep = "")

# ----------------------------------------------------------------------------
# 2. COMPTAGES BRUTS
# ----------------------------------------------------------------------------

cat("\n[2/7] Computing raw counts...\n")

arg_counts <- rgi_all %>%
  count(sample_id, best_hit_aro, name = "count")

write_tsv(arg_counts, output_counts)
cat("  ✓ Saved: ", basename(output_counts), " (", nrow(arg_counts), " rows)\n", sep = "")

# ----------------------------------------------------------------------------
# 3. ABONDANCES RELATIVES
# ----------------------------------------------------------------------------

cat("\n[3/7] Computing relative abundances...\n")

arg_relative <- arg_counts %>%
  group_by(sample_id) %>%
  mutate(
    total_args = sum(count),
    relative_abundance = (count / total_args) * 100
  ) %>%
  ungroup()

write_tsv(arg_relative, output_relative)
cat("  ✓ Saved: ", basename(output_relative), "\n", sep = "")

# ----------------------------------------------------------------------------
# 4. NORMALIZED COPY NUMBER (par gigabase)
# ----------------------------------------------------------------------------

cat("\n[4/7] Computing normalized copy number (per Gb)...\n")

# Calculer la taille du dataset en Gb
sequencing_stats_gb <- sequencing_stats %>%
  mutate(total_bases_gb = total_bases / 1e9)

cat("  Dataset sizes (Gb):\n")
for (i in 1:nrow(sequencing_stats_gb)) {
  cat(sprintf("    %s: %.3f Gb\n", 
              sequencing_stats_gb$sample_id[i], 
              sequencing_stats_gb$total_bases_gb[i]))
}

# Normalisation : ARG count / Total bases (Gb)
arg_normalized <- arg_counts %>%
  left_join(sequencing_stats_gb %>% select(sample_id, total_bases_gb), 
            by = "sample_id") %>%
  mutate(
    normalized_copy_number = count / total_bases_gb
  )

write_tsv(arg_normalized, output_normalized)
cat("  ✓ Saved: ", basename(output_normalized), "\n", sep = "")
cat("  ✓ Normalized copy number = ARG count / Dataset size (Gb)\n")

# ----------------------------------------------------------------------------
# 5. PRÉSENCE/ABSENCE
# ----------------------------------------------------------------------------

cat("\n[5/7] Creating presence/absence table...\n")

arg_presence <- arg_counts %>%
  mutate(presence = 1) %>%
  select(sample_id, best_hit_aro, presence)

write_tsv(arg_presence, output_presence)
cat("  ✓ Saved: ", basename(output_presence), "\n", sep = "")

# ----------------------------------------------------------------------------
# 6. MATRICES
# ----------------------------------------------------------------------------

cat("\n[6/7] Creating matrices...\n")

# Matrice comptages
arg_matrix_counts <- arg_counts %>%
  pivot_wider(names_from = best_hit_aro, values_from = count, values_fill = 0)

write_tsv(arg_matrix_counts, output_matrix_counts)
cat("  ✓ arg_matrix_counts: ", nrow(arg_matrix_counts), " × ", 
    ncol(arg_matrix_counts) - 1, "\n", sep = "")

# Matrice relative
arg_matrix_relative <- arg_relative %>%
  select(sample_id, best_hit_aro, relative_abundance) %>%
  pivot_wider(names_from = best_hit_aro, values_from = relative_abundance, values_fill = 0)

write_tsv(arg_matrix_relative, output_matrix_relative)
cat("  ✓ arg_matrix_relative\n")

# Matrice normalized copy number
arg_matrix_normalized <- arg_normalized %>%
  select(sample_id, best_hit_aro, normalized_copy_number) %>%
  pivot_wider(names_from = best_hit_aro, values_from = normalized_copy_number, values_fill = 0)

write_tsv(arg_matrix_normalized, output_matrix_normalized)
cat("  ✓ arg_matrix_normalized (copy number per Gb)\n")

# Matrice présence
arg_matrix_presence <- arg_presence %>%
  pivot_wider(names_from = best_hit_aro, values_from = presence, values_fill = 0)

write_tsv(arg_matrix_presence, output_matrix_presence)
cat("  ✓ arg_matrix_presence\n")

# ----------------------------------------------------------------------------
# 7. AGRÉGATIONS FONCTIONNELLES
# ----------------------------------------------------------------------------

cat("\n[7/7] Creating functional aggregations...\n")

# Par classe d'antibiotique
drug_class_abundance <- rgi_all %>%
  separate_rows(drug_class, sep = "; ") %>%
  mutate(drug_class = str_trim(drug_class)) %>%
  count(sample_id, drug_class, name = "count") %>%
  group_by(sample_id) %>%
  mutate(
    total = sum(count),
    relative_abundance = (count / total) * 100
  ) %>%
  ungroup()

write_tsv(drug_class_abundance, output_drug_class)
cat("  ✓ drug_class_abundance (", n_distinct(drug_class_abundance$drug_class), 
    " classes)\n", sep = "")

# Par mécanisme
mechanism_abundance <- rgi_all %>%
  count(sample_id, resistance_mechanism, name = "count") %>%
  group_by(sample_id) %>%
  mutate(
    total = sum(count),
    relative_abundance = (count / total) * 100
  ) %>%
  ungroup()

write_tsv(mechanism_abundance, output_mechanism)
cat("  ✓ mechanism_abundance (", n_distinct(mechanism_abundance$resistance_mechanism), 
    " mechanisms)\n", sep = "")

# Par famille
family_abundance <- rgi_all %>%
  count(sample_id, amr_gene_family, name = "count") %>%
  group_by(sample_id) %>%
  mutate(
    total = sum(count),
    relative_abundance = (count / total) * 100
  ) %>%
  ungroup()

write_tsv(family_abundance, output_family)
cat("  ✓ family_abundance (", n_distinct(family_abundance$amr_gene_family), 
    " families)\n", sep = "")

# ----------------------------------------------------------------------------
# RÉSUMÉ
# ----------------------------------------------------------------------------

cat("\n", rep("=", 70), "\n", sep = "")
cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
cat(rep("=", 70), "\n", sep = "")
cat("Normalization method: Normalized copy number per Gb\n")
cat("  Formula: ARG count / (Total bases / 10⁹)\n")
cat("\nFiles created in results/r_analysis/:\n")
cat("  • Counts and normalization: 4 files\n")
cat("  • Matrices: 4 files\n")
cat("  • Functional aggregations: 3 files\n")
cat("  • Total: 11 files\n")
cat(rep("=", 70), "\n\n", sep = "")

sink(type = "message")
sink()
close(log)
