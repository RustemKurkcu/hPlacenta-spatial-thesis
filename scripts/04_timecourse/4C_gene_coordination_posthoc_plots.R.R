# ======================================================================
# scripts/04_timecourse/04C_gene_coordination_posthoc_plots.R
# Lightweight post-hoc plotting for 04C using existing CSV outputs.
#
# Purpose:
# - Avoid rerunning heavy COM/null computations in 04C.
# - Create publication-friendly summary plots from saved results.
#
# Inputs:
# - output/tables/gene_coordination_scores.csv
#
# Outputs:
# - output/tables/gene_coordination_scores_with_stats.csv
# - output/tables/gene_coordination_significant_only.csv
# - output/figures/gene_coordination_dotplot_all.png
# - output/figures/gene_coordination_dotplot_significant_only.png
# - output/figures/gene_coordination_top_hits_barplot.png
# ======================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

source("config/config.R")
source("scripts/R/utils.R")

ensure_dir(DIR_TABLES); ensure_dir(DIR_FIGURES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "04C_gene_coordination_posthoc_plots.log")

infile <- file.path(DIR_TABLES, "gene_coordination_scores.csv")
if (!file.exists(infile)) stop("Missing input: ", infile, ". Run 04C first.")

res <- read.csv(infile, stringsAsFactors = FALSE)
if (nrow(res) == 0) stop("Input table is empty: ", infile)

if (!("coordination_z" %in% colnames(res))) {
  stop("coordination_z column not found in ", infile)
}

if (!("p_value" %in% colnames(res))) {
  res$p_value <- 2 * pnorm(-abs(res$coordination_z))
}
res$p_adj <- p.adjust(res$p_value, method = "BH")
res$significant <- !is.na(res$p_adj) & res$p_adj < 0.05

out_all <- file.path(DIR_TABLES, "gene_coordination_scores_with_stats.csv")
out_sig <- file.path(DIR_TABLES, "gene_coordination_significant_only.csv")
write.csv(res, out_all, row.names = FALSE)
write.csv(res %>% filter(significant), out_sig, row.names = FALSE)

# Biological ordering helper
order_celltypes_bio <- function(ct) {
  ct <- unique(as.character(ct))
  grp <- dplyr::case_when(
    grepl("EVT|CTB|STB|Troph", ct, ignore.case = TRUE) ~ "01_trophoblast",
    grepl("NK|Tcell|T cell|Bcell|B cell|mac|mono|myelo|immune", ct, ignore.case = TRUE) ~ "02_immune",
    grepl("Endo|vascular|vessel", ct, ignore.case = TRUE) ~ "03_endothelial",
    grepl("Fib|strom|mesench", ct, ignore.case = TRUE) ~ "04_stromal",
    grepl("gland|epith", ct, ignore.case = TRUE) ~ "05_epithelial",
    TRUE ~ "99_other"
  )
  tibble::tibble(celltype = ct, grp = grp) %>%
    arrange(grp, celltype) %>%
    pull(celltype)
}

if ("geneset" %in% colnames(res) && "celltype" %in% colnames(res)) {
  geneset_order <- unique(res$geneset)
  celltype_order <- order_celltypes_bio(res$celltype)
  
  # 1) Dot plot (all)
  p_all <- res %>%
    mutate(
      geneset = factor(geneset, levels = geneset_order),
      celltype = factor(celltype, levels = rev(celltype_order)),
      sig_label = ifelse(significant, "FDR<0.05", "NS")
    ) %>%
    ggplot(aes(x = geneset, y = celltype)) +
    geom_point(aes(size = pmax(1, n_genes), fill = coordination_z, shape = sig_label), color = "black", alpha = 0.9) +
    scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick4", midpoint = 0) +
    scale_shape_manual(values = c("FDR<0.05" = 21, "NS" = 24)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Gene coordination summary (all results)",
      subtitle = "Dot size = n_genes; fill = coordination z; shape = significance",
      x = "Gene set", y = "Cell type", fill = "coordination z", size = "n_genes", shape = "Trend"
    )
  save_plot(p_all, file.path(DIR_FIGURES, "gene_coordination_dotplot_all.png"), w = 11, h = 9)
  
  # 2) Significant-only dot plot
  res_sig <- res %>% filter(significant)
  if (nrow(res_sig) > 0) {
    p_sig <- res_sig %>%
      mutate(
        geneset = factor(geneset, levels = geneset_order),
        celltype = factor(celltype, levels = rev(celltype_order))
      ) %>%
      ggplot(aes(x = geneset, y = celltype)) +
      geom_point(aes(size = pmax(1, n_genes), fill = coordination_z), shape = 21, color = "black", alpha = 0.95) +
      scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick4", midpoint = 0) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = "Gene coordination summary (significant only)",
        subtitle = "Filtered to BH FDR<0.05",
        x = "Gene set", y = "Cell type", fill = "coordination z", size = "n_genes"
      )
    save_plot(p_sig, file.path(DIR_FIGURES, "gene_coordination_dotplot_significant_only.png"), w = 11, h = 9)
  }
  
  # 3) Top hits bar plot
  n_top <- min(20L, nrow(res))
  top_hits <- res %>%
    arrange(p_adj, desc(abs(coordination_z))) %>%
    slice_head(n = n_top) %>%
    mutate(label = paste(celltype, geneset, sep = " | "))
  
  p_top <- ggplot(top_hits, aes(x = reorder(label, coordination_z), y = coordination_z, fill = coordination_z)) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick4", midpoint = 0) +
    theme_classic() +
    labs(title = "Top coordination hits", subtitle = "Ranked by FDR then |z|", x = "Cell type | Gene set", y = "coordination z")
  save_plot(p_top, file.path(DIR_FIGURES, "gene_coordination_top_hits_barplot.png"), w = 11, h = 7)
}

log_msg("Done post-hoc coordination plotting.", logfile)