#!/usr/bin/env Rscript

# =============================================================================
# Script: 01c_plot_spatial_embeddings.R
# Purpose: Decoupled visualization-only script for integrated STARmap object.
#          Loads precomputed object and generates embedding + tissue maps.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})

source("R/spatial_color_themes.R")
source("R/celltype_dictionary.R")

OUT_ROOT <- "output"
DIR_OBJECTS <- file.path(OUT_ROOT, "objects")
DIR_FIGURES <- file.path(OUT_ROOT, "figures")
DIR_REPORTS <- file.path(OUT_ROOT, "reports")
for (d in c(DIR_FIGURES, DIR_REPORTS)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "[INFO]", paste0(..., collapse = ""), "\n")
}

resolve_object_path <- function(path_rds) {
  path_qs <- sub("\\.rds$", ".qs", path_rds)
  if (file.exists(path_qs)) return(path_qs)
  if (file.exists(path_rds)) return(path_rds)
  stop("Missing object: ", path_rds)
}

read_object <- function(path) {
  if (grepl("\\.qs$", path)) {
    if (!requireNamespace("qs", quietly = TRUE)) stop("Need package 'qs' to read: ", path)
    return(qs::qread(path))
  }
  readRDS(path)
}

record_plot_manifest_01c <- function(manifest_path, source_data, plotting_script, figures_output) {
  manifest <- list(
    source_data = source_data,
    plotting_script = plotting_script,
    figures_output = figures_output,
    run_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
}

obj_path <- resolve_object_path(file.path(DIR_OBJECTS, "01_integrated_harmony_sct_4weeks.rds"))
log_msg("Loading integrated object: ", obj_path)
seu <- read_object(obj_path)

# Directive 1: force sequential developmental ordering in all plots.
seu$week <- factor(as.character(seu$week), levels = c("W7", "W8-2", "W9", "W11"))

# Prefer predicted.celltype if present for tissue map coloring.
celltype_col <- if ("predicted.celltype" %in% colnames(seu@meta.data)) "predicted.celltype" else pick_celltype_source_column(seu@meta.data)
if (is.na(celltype_col)) celltype_col <- "seurat_clusters"
seu$celltype_plot <- rename_with_true_names(as.character(seu@meta.data[[celltype_col]]))
cell_levels <- sort(unique(as.character(seu$celltype_plot)))
cell_cols <- get_universal_colors(cell_levels)

# -----------------------------------------------------------------------------
# A) Standard embedding plots by week and cluster
# -----------------------------------------------------------------------------
p_umap_week <- DimPlot(seu, reduction = "umap_harmony", group.by = "week", pt.size = 0.2) +
  scale_color_manual(values = misi_week_colors, drop = FALSE) +
  ggtitle("UMAP (Harmony) by Week") +
  theme_thesis_spatial()

p_tsne_week <- DimPlot(seu, reduction = "tsne_harmony", group.by = "week", pt.size = 0.2) +
  scale_color_manual(values = misi_week_colors, drop = FALSE) +
  ggtitle("t-SNE (Harmony) by Week") +
  theme_thesis_spatial()

p_umap_cluster <- DimPlot(seu, reduction = "umap_harmony", group.by = "seurat_clusters", label = TRUE, pt.size = 0.2) +
  ggtitle("UMAP (Harmony) by Cluster") +
  theme_thesis_spatial()

p_tsne_cluster <- DimPlot(seu, reduction = "tsne_harmony", group.by = "seurat_clusters", label = TRUE, pt.size = 0.2) +
  ggtitle("t-SNE (Harmony) by Cluster") +
  theme_thesis_spatial()

# -----------------------------------------------------------------------------
# B) Tissue map (physical reconstruction): x_um vs y_um, colored by celltype
# -----------------------------------------------------------------------------
spatial_df <- seu@meta.data
p_tissue_celltype <- ggplot(spatial_df, aes(x = x_um, y = y_um, color = celltype_plot)) +
  geom_point(size = 0.25, alpha = 0.85) +
  facet_wrap(~week, ncol = 4) +
  scale_color_manual(values = cell_cols, drop = FALSE) +
  coord_fixed() +
  labs(
    title = paste0("Physical Tissue Reconstruction by Cell Type (", celltype_col, " -> canonical)"),
    x = "x_um", y = "y_um", color = "Cell Type"
  ) +
  theme_thesis_spatial()

out_files <- c(
  file.path(DIR_FIGURES, "01c_embeddings_week_umap_tsne.pdf"),
  file.path(DIR_FIGURES, "01c_embeddings_week_umap_tsne.png"),
  file.path(DIR_FIGURES, "01c_embeddings_clusters_umap_tsne.pdf"),
  file.path(DIR_FIGURES, "01c_embeddings_clusters_umap_tsne.png"),
  file.path(DIR_FIGURES, "01c_tissue_map_celltype_by_week.pdf"),
  file.path(DIR_FIGURES, "01c_tissue_map_celltype_by_week.png")
)

ggsave(out_files[1], p_umap_week + p_tsne_week, width = 14, height = 6)
ggsave(out_files[2], p_umap_week + p_tsne_week, width = 14, height = 6, dpi = 320)
ggsave(out_files[3], p_umap_cluster + p_tsne_cluster, width = 14, height = 6)
ggsave(out_files[4], p_umap_cluster + p_tsne_cluster, width = 14, height = 6, dpi = 320)
ggsave(out_files[5], p_tissue_celltype, width = 16, height = 7)
ggsave(out_files[6], p_tissue_celltype, width = 16, height = 7, dpi = 320)

plot_manifest_path <- file.path(DIR_REPORTS, "01c_plot_spatial_embeddings_manifest.json")
record_plot_manifest_01c(
  manifest_path = plot_manifest_path,
  source_data = obj_path,
  plotting_script = "scripts/01_active_pipeline/01c_plot_spatial_embeddings.R",
  figures_output = out_files
)

log_msg("Saved plotting artifacts and manifest: ", plot_manifest_path)
