#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  out_plot <- if (length(args) >= 1) args[1] else "results/summary/r_all_ARGs.png"
  out_tsv  <- if (length(args) >= 2) args[2] else "results/summary/r_all_ARGs_counts.tsv"
  dir.create(dirname(out_plot), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(out_tsv),  recursive = TRUE, showWarnings = FALSE)
  png(out_plot, width = 1200, height = 800, res = 150)
  plot.new(); text(0.5, 0.5, "No data: insufficient arguments", cex = 1.2)
  dev.off()
  write_tsv(tibble(Sample=character(), GENE=character()), out_tsv)
  quit(save="no", status=0)
}

out_plot   <- args[1]
out_tsv    <- args[2]
input_files <- args[-c(1,2)]

dir.create(dirname(out_plot), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_tsv),  recursive = TRUE, showWarnings = FALSE)

dfs <- list()

safe_size <- function(fp) {
  if (!file.exists(fp)) return(0L)
  sz <- tryCatch(file.info(fp)$size, error = function(e) 0L)
  ifelse(is.na(sz), 0L, as.integer(sz))
}

for (file in input_files) {
  if (safe_size(file) <= 0L) {
    cat("[r_plot_ARGs_all] Skip (missing/empty):", file, "\n")
    next
  }
  sample <- strsplit(basename(file), "_")[[1]][1]

  df <- tryCatch(
    read_tsv(file, show_col_types = FALSE),
    error = function(e) {
      cat("[r_plot_ARGs_all] Read error for", file, ":", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(df) || nrow(df) == 0) {
    cat("[r_plot_ARGs_all] Skip (no rows):", file, "\n")
    next
  }

  # normalise les noms (utile si variantes)
  names(df) <- toupper(trimws(names(df)))
  if (!"GENE" %in% names(df)) {
    cat("[r_plot_ARGs_all] Skip (no 'GENE' col):", file, "Found:", paste(names(df), collapse=", "), "\n")
    next
  }

  df <- df %>% transmute(Sample = sample, GENE = as.character(GENE))
  dfs[[length(dfs) + 1]] <- df
}

if (length(dfs) == 0) {
  # Aucun fichier exploitable → placeholders
  png(out_plot, width = 1200, height = 800, res = 150)
  plot.new(); text(0.5, 0.5, "No data: no valid ARG summaries", cex = 1.2)
  dev.off()
  write_tsv(tibble(Sample=character(), GENE=character()), out_tsv)
  quit(save="no", status=0)
}

all_data <- bind_rows(dfs) %>% filter(!is.na(GENE) & nzchar(GENE))

if (nrow(all_data) == 0) {
  png(out_plot, width = 1200, height = 800, res = 150)
  plot.new(); text(0.5, 0.5, "No data: empty after filtering", cex = 1.2)
  dev.off()
  write_tsv(tibble(Sample=character(), GENE=character()), out_tsv)
  quit(save="no", status=0)
}

counts <- all_data %>%
  group_by(Sample, GENE) %>%
  summarise(Count = n(), .groups = "drop")

# Tableau pivot large (comme ton script)
pivot_table <- counts %>%
  mutate(Sample = as.character(Sample), GENE = as.character(GENE)) %>%
  pivot_wider(names_from = GENE, values_from = Count, values_fill = 0)

write_tsv(pivot_table, out_tsv)

# Pour éviter un plot illisible si trop de gènes, tu peux limiter aux top N (ex: 30)
# Commenter/décommenter selon besoin :
# topN <- 30
# top_genes <- counts %>% group_by(GENE) %>% summarise(Total = sum(Count), .groups="drop") %>%
#   arrange(desc(Total)) %>% slice_head(n = topN) %>% pull(GENE)
# counts <- counts %>% filter(GENE %in% top_genes)

p <- ggplot(counts, aes(x = Sample, y = Count, fill = GENE)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "Distribution of resistance genes across samples",
    x = "Sample", y = "Gene Count", fill = "Gene"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width  = unit(0.8, "lines"))

ggsave(out_plot, plot = p, width = 12, height = 8, dpi = 150, device = "png")
quit(save="no", status=0)
