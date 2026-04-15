#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(jsonlite)
})

source("R/spatial_color_themes.R")
source("R/celltype_dictionary.R")

PIPELINE_NAME <- "02c_plot_spatial_misi"
PIPELINE_VERSION <- "1.0.0"
RUN_TIMESTAMP <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
RUN_SEED <- 42L
set.seed(RUN_SEED)

OUT_ROOT <- "output"
DIR_LOGS <- file.path(OUT_ROOT, "logs")
DIR_OBJECTS <- file.path(OUT_ROOT, "objects")
DIR_FIGURES <- file.path(OUT_ROOT, "figures")
DIR_REPORTS <- file.path(OUT_ROOT, "reports")
for (d in c(DIR_LOGS, DIR_OBJECTS, DIR_FIGURES, DIR_REPORTS)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
LOG_FILE <- file.path(DIR_LOGS, "02c_plot_spatial_misi.log")

log_msg <- function(..., .level = "INFO") {
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  txt <- paste0(..., collapse = "")
  line <- paste0(stamp, " [", .level, "] ", txt)
  cat(line, "\n")
  cat(line, "\n", file = LOG_FILE, append = TRUE)
}

resolve_object_path <- function(path_rds) {
  path_qs <- sub("\\.rds$", ".qs", path_rds)
  if (file.exists(path_qs)) return(path_qs)
  if (file.exists(path_rds)) return(path_rds)
  NA_character_
}

read_object <- function(path) {
  if (grepl("\\.qs$", path)) {
    if (!requireNamespace("qs", quietly = TRUE)) stop("Need package 'qs' to read: ", path)
    return(qs::qread(path))
  }
  readRDS(path)
}

record_artifact_manifest <- function(
  manifest_path,
  pipeline,
  version,
  run_timestamp,
  seed,
  source_data,
  compute_script,
  plotting_script,
  figures_output
) {
  manifest <- list(
    pipeline = pipeline,
    version = version,
    run_timestamp = run_timestamp,
    seed = seed,
    source_data = source_data,
    compute_script = compute_script,
    plotting_script = plotting_script,
    figures_output = figures_output
  )
  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
  invisible(manifest)
}

input_obj <- file.path(DIR_OBJECTS, "02_scored_misi_ido1.rds")
input_obj_found <- resolve_object_path(input_obj)
if (is.na(input_obj_found)) {
  stop("Missing scored object: ", input_obj,
       "\nRun scripts/01_active_pipeline/02_spatial_misi_ido1_scoring.R first.")
}

log_msg("Loading scored object: ", input_obj_found)
seu <- read_object(input_obj_found)

features_to_plot <- c(
  "MMP_ECM_Remodeling",
  "Immune_Exhaustion_Tolerogenic",
  "MegL_Sulfide_Vulnerability",
  "IDO1_Tolerogenic_Shield",
  "MISI_Vulnerability"
)

missing_features <- setdiff(features_to_plot, colnames(seu@meta.data))
if (length(missing_features) > 0) {
  stop("Missing required score columns: ", paste(missing_features, collapse = ", "))
}

reductions <- c("umap_harmony", "tsne_harmony")
missing_reductions <- setdiff(reductions, names(seu@reductions))
if (length(missing_reductions) > 0) {
  stop("Missing required reductions: ", paste(missing_reductions, collapse = ", "))
}

fig_paths <- character(0)
celltype_col <- pick_celltype_source_column(seu@meta.data)

for (red in reductions) {
  for (feat in features_to_plot) {
    p <- FeaturePlot(
      object = seu,
      features = feat,
      reduction = red,
      cols = c("grey92", "#CB181D"),
      pt.size = 0.25,
      order = TRUE
    ) +
      ggtitle(paste0(feat, " on ", red)) +
      theme_thesis_spatial()

    base <- file.path(DIR_FIGURES, paste0("02c_", feat, "_", red))
    pdf_path <- paste0(base, ".pdf")
    png_path <- paste0(base, ".png")

    ggsave(pdf_path, p, width = 8, height = 6)
    ggsave(png_path, p, width = 8, height = 6, dpi = 300)

    fig_paths <- c(fig_paths, pdf_path, png_path)
  }

  if (!is.na(celltype_col)) {
    cell_levels <- sort(unique(as.character(seu@meta.data[[celltype_col]])))
    cell_cols <- get_universal_colors(cell_levels)
    p_ct <- DimPlot(
      object = seu,
      reduction = red,
      group.by = celltype_col,
      label = TRUE,
      pt.size = 0.2
    ) +
      scale_color_manual(values = cell_cols, drop = FALSE) +
      ggtitle(paste0("Cell type overlay on ", red, " (", celltype_col, ")")) +
      theme_thesis_spatial()

    base_ct <- file.path(DIR_FIGURES, paste0("02c_celltype_overlay_", red))
    pdf_ct <- paste0(base_ct, ".pdf")
    png_ct <- paste0(base_ct, ".png")
    ggsave(pdf_ct, p_ct, width = 10, height = 8)
    ggsave(png_ct, p_ct, width = 10, height = 8, dpi = 300)
    fig_paths <- c(fig_paths, pdf_ct, png_ct)
  }
}

# Spatial heat maps in physical coordinates for each score
if (all(c("x_um", "y_um", "week") %in% colnames(seu@meta.data))) {
  spatial_df <- seu@meta.data
  spatial_df$week <- as.character(spatial_df$week)
  for (feat in features_to_plot) {
    d <- spatial_df[, c("x_um", "y_um", "week", feat), drop = FALSE]
    names(d)[4] <- "score"
    d <- d[is.finite(d$x_um) & is.finite(d$y_um) & is.finite(d$score), , drop = FALSE]
    if (nrow(d) == 0) next

    p_red <- ggplot(d, aes(x = x_um, y = y_um, color = score)) +
      geom_point(size = 0.15, alpha = 0.8) +
      facet_wrap(~week, scales = "free") +
      scale_color_gradient(low = "white", high = "#8B0000") +
      coord_equal() +
      ggtitle(paste0(feat, " spatial heatmap (white → dark red)")) +
      theme_thesis_spatial()

    p_div <- ggplot(d, aes(x = x_um, y = y_um, color = score)) +
      geom_point(size = 0.15, alpha = 0.8) +
      facet_wrap(~week, scales = "free") +
      scale_color_gradient2(low = "#0B1F5E", mid = "white", high = "#CB181D", midpoint = median(d$score, na.rm = TRUE)) +
      coord_equal() +
      ggtitle(paste0(feat, " spatial heatmap (navy → white → red)")) +
      theme_thesis_spatial()

    base_red <- file.path(DIR_FIGURES, paste0("02c_spatial_", feat, "_white_to_red"))
    base_div <- file.path(DIR_FIGURES, paste0("02c_spatial_", feat, "_navy_white_red"))
    ggsave(paste0(base_red, ".pdf"), p_red, width = 12, height = 8)
    ggsave(paste0(base_red, ".png"), p_red, width = 12, height = 8, dpi = 320)
    ggsave(paste0(base_div, ".pdf"), p_div, width = 12, height = 8)
    ggsave(paste0(base_div, ".png"), p_div, width = 12, height = 8, dpi = 320)
    fig_paths <- c(
      fig_paths,
      paste0(base_red, ".pdf"), paste0(base_red, ".png"),
      paste0(base_div, ".pdf"), paste0(base_div, ".png")
    )
  }
} else {
  log_msg("Skipping spatial heat maps: need x_um, y_um, and week columns in metadata.", .level = "WARN")
}

manifest_path <- file.path(DIR_REPORTS, "02c_plot_spatial_misi_manifest.json")
record_artifact_manifest(
  manifest_path = manifest_path,
  pipeline = PIPELINE_NAME,
  version = PIPELINE_VERSION,
  run_timestamp = RUN_TIMESTAMP,
  seed = RUN_SEED,
  source_data = input_obj_found,
  compute_script = "scripts/01_active_pipeline/02_spatial_misi_ido1_scoring.R",
  plotting_script = "scripts/01_active_pipeline/02c_plot_spatial_misi.R",
  figures_output = fig_paths
)

log_msg("Saved MISI plots and manifest: ", manifest_path)
