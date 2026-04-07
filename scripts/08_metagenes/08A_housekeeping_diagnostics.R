# ======================================================================
# scripts/08_metagenes/08A_housekeeping_diagnostics.R
# Housekeeping diagnostics: detection-frequency vs mean expression
# 
# PURPOSE:
#   Identify stable genes and sample shifts; QC not normalization replacement.
#   Based on Bo/Nelson Chapter 3 logic.
#
# METHOD:
#   - Calculate detection frequency (% cells expressing each gene)
#   - Calculate mean expression per gene
#   - Identify housekeeping candidates (high detection + moderate expression)
#   - Generate QC plots for technical stability assessment
#
# OUTPUT:
#   - Per-object detection vs mean expression CSVs
#   - Top 200 housekeeping candidate genes
#   - QC diagnostic plots
#
# USAGE:
#   source("scripts/08_metagenes/08A_housekeeping_diagnostics.R")
#
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
})

source("config/config.R")
source("scripts/R/utils.R")

ensure_dir(DIR_TABLES); ensure_dir(DIR_FIGURES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "08A_housekeeping_diagnostics.log")

# Create metagenes subdirectories
dir.create(file.path(DIR_TABLES, "metagenes"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(DIR_FIGURES, "metagenes"), showWarnings = FALSE, recursive = TRUE)

log_msg("Starting housekeeping diagnostics...", logfile)

# ---- Choose object(s) to QC ----
paths <- c(
  file.path(DIR_OBJS, "multiome_reference_processed.rds"),
  file.path(DIR_OBJS, "slidetags_mapped_to_multiome.rds"),
  file.path(DIR_OBJS, "starmap_mapped_to_multiome.rds")
)

paths <- paths[file.exists(paths)]
if (length(paths) == 0) stop("No processed objects found in ", DIR_OBJS)

log_msg(paste0("Running housekeeping diagnostics on: ", 
               paste(basename(paths), collapse=", ")), logfile)

# ---- Helper functions ----

#' Calculate housekeeping statistics
#' 
#' @param obj Seurat object
#' @param assay Assay to use
#' @param slot Slot to use (counts or data)
#' 
#' @return Tibble with gene, detect_frac, mean_expr
#' 
hk_stats <- function(obj, assay="RNA", slot="counts") {
  DefaultAssay(obj) <- assay
  m <- GetAssayData(obj, slot = slot)
  
  # Convert to sparse matrix if needed
  if (!inherits(m, "dgCMatrix")) m <- as(m, "dgCMatrix")
  
  # Calculate statistics
  detect_frac <- Matrix::rowSums(m > 0) / ncol(m)
  mean_expr   <- Matrix::rowMeans(m)
  
  tibble(
    gene = rownames(m),
    detect_frac = as.numeric(detect_frac),
    mean_expr = as.numeric(mean_expr)
  ) %>%
    mutate(mean_expr_log1p = log1p(mean_expr))
}

#' Pick housekeeping candidates
#' 
#' @param df Data frame from hk_stats
#' @param min_detect Minimum detection frequency (default: 0.90)
#' @param expr_quantile Expression quantile threshold (default: 0.50)
#' 
#' @return Filtered tibble with housekeeping candidates
#' 
pick_housekeeping <- function(df, min_detect=0.90, expr_quantile=0.50) {
  # Stable: high detect + not just extreme high-expression ribosomal
  expr_cut <- quantile(df$mean_expr_log1p, probs = expr_quantile, na.rm = TRUE)
  
  df %>%
    filter(detect_frac >= min_detect, mean_expr_log1p >= expr_cut) %>%
    arrange(desc(detect_frac), desc(mean_expr_log1p))
}

# ---- Process each object ----
for (p in paths) {
  nm <- tools::file_path_sans_ext(basename(p))
  log_msg(paste0("Processing: ", nm), logfile)
  
  obj <- readRDS(p)
  stopifnot(inherits(obj, "Seurat"))
  
  # Select assay (prefer RNA for raw counts)
  assay_use <- if ("RNA" %in% names(obj@assays)) "RNA" else DefaultAssay(obj)
  slot_use  <- "counts"
  
  log_msg(paste0("  Using assay: ", assay_use, ", slot: ", slot_use), logfile)
  
  # Calculate statistics
  df <- hk_stats(obj, assay=assay_use, slot=slot_use)
  
  # Pick housekeeping candidates
  hk <- pick_housekeeping(df, min_detect=0.90, expr_quantile=0.50) %>% 
    head(200)
  
  log_msg(sprintf("  Found %d housekeeping candidates (top 200)", nrow(hk)), logfile)
  
  # Save results
  write.csv(df, 
            file.path(DIR_TABLES, "metagenes", paste0(nm, "_detect_vs_meanexpr.csv")), 
            row.names = FALSE)
  write.csv(hk, 
            file.path(DIR_TABLES, "metagenes", paste0(nm, "_housekeeping_candidates_top200.csv")), 
            row.names = FALSE)
  
  # Plot: detection vs mean expression
  p1 <- ggplot(df, aes(x = mean_expr_log1p, y = detect_frac)) +
    geom_point(size=0.4, alpha=0.4, color="gray60") +
    geom_point(data=hk, aes(x=mean_expr_log1p, y=detect_frac), 
               size=0.8, color="red", alpha=0.6) +
    geom_hline(yintercept=0.90, linetype="dashed", color="blue", alpha=0.5) +
    labs(
      title=paste0("Housekeeping Diagnostic: ", nm),
      subtitle=sprintf("Red points: Top 200 housekeeping candidates (detect >= 90%%, n=%d)", nrow(hk)),
      x="log1p(mean expression)",
      y="Detection fraction (% cells expressing)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  save_plot(p1, 
            file.path(DIR_FIGURES, "metagenes", paste0(nm, "_housekeeping_detect_vs_mean.png")), 
            w=8, h=6)
  
  log_msg(paste0("  Saved figures and tables for ", nm), logfile)
}

log_msg("08A complete.", logfile)

cat("\n")
cat(strrep("=", 70), "\n")
cat("HOUSEKEEPING DIAGNOSTICS COMPLETE\n")
cat(strrep("=", 70), "\n")
cat(sprintf("Tables: %s/metagenes/*_detect_vs_meanexpr.csv\n", DIR_TABLES))
cat(sprintf("        %s/metagenes/*_housekeeping_candidates_top200.csv\n", DIR_TABLES))
cat(sprintf("Figures: %s/metagenes/*_housekeeping_detect_vs_mean.png\n", DIR_FIGURES))
cat(strrep("=", 70), "\n\n")

cat("INTERPRETATION:\n")
cat("- Red points: Stable housekeeping candidates (high detection, moderate expression)\n")
cat("- Use for: Technical QC, sample quality assessment, normalization validation\n")
cat("- NOT for: Replacing standard normalization methods\n\n")

# ======================================================================
# END OF SCRIPT
# ======================================================================