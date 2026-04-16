#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

source("R/spatial_color_themes.R")
source("R/celltype_dictionary.R")

OUT_ROOT <- "output"
DIR_OBJECTS <- file.path(OUT_ROOT, "objects")
DIR_FIGURES <- file.path(OUT_ROOT, "figures")
DIR_REPORTS <- file.path(OUT_ROOT, "reports")
if (!dir.exists(DIR_FIGURES)) dir.create(DIR_FIGURES, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "[INFO]", paste0(..., collapse = ""), "\n")
}

resolve_object_path <- function(path_rds) {
  path_with_ct_qs <- sub("\\.rds$", "_with_celltypes.qs", path_rds)
  path_with_ct_rds <- sub("\\.rds$", "_with_celltypes.rds", path_rds)
  path_qs <- sub("\\.rds$", ".qs", path_rds)
  if (file.exists(path_with_ct_qs)) return(path_with_ct_qs)
  if (file.exists(path_with_ct_rds)) return(path_with_ct_rds)
  if (file.exists(path_qs)) return(path_qs)
  if (file.exists(path_rds)) return(path_rds)
  NA_character_
}

read_object <- function(path) {
  if (grepl("\\.qs$", path) || (!grepl("\\.rds$", path) && requireNamespace("qs", quietly = TRUE))) {
    if (!requireNamespace("qs", quietly = TRUE)) stop("Need package 'qs' to read: ", path)
    return(qs::qread(path))
  }
  readRDS(path)
}

record_plot_manifest_01c <- function(
  manifest_path,
  source_data,
  compute_script,
  plotting_script,
  figures_output
) {
  manifest <- list(
    source_data = source_data,
    compute_script = compute_script,
    plotting_script = plotting_script,
    figures_output = figures_output,
    run_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
}

choose_celltype_col <- function(md) {
  for (cand in c("predicted.celltype", "cell_label_display", "celltype_corrected", "celltype_true_name", "celltype_original", "seurat_clusters")) {
    if (cand %in% colnames(md)) return(cand)
  }
  NA_character_
}

obj_path <- file.path(DIR_OBJECTS, "01_integrated_harmony_sct_4weeks.rds")
obj_path_found <- resolve_object_path(obj_path)
if (is.na(obj_path_found)) {
  stop("Missing expected compute artifact: ", obj_path,
       "\nRun scripts/01_active_pipeline/01_preprocess_harmony_embeddings.R first.")
}

seu <- read_object(obj_path_found)

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

out_files <- c(
  file.path(DIR_FIGURES, "01c_embeddings_week_umap_tsne.pdf"),
  file.path(DIR_FIGURES, "01c_embeddings_week_umap_tsne.png"),
  file.path(DIR_FIGURES, "01c_embeddings_clusters_umap_tsne.pdf"),
  file.path(DIR_FIGURES, "01c_embeddings_clusters_umap_tsne.png")
)

week_patch <- p_umap_week + p_tsne_week
cluster_patch <- p_umap_cluster + p_tsne_cluster

ggsave(out_files[1], week_patch, width = 14, height = 6)
ggsave(out_files[2], week_patch, width = 14, height = 6, dpi = 320)
ggsave(out_files[3], cluster_patch, width = 14, height = 6)
ggsave(out_files[4], cluster_patch, width = 14, height = 6, dpi = 320)

celltype_col <- choose_celltype_col(seu@meta.data)
if (!is.na(celltype_col)) {
  log_msg("Using celltype column for overlays: ", celltype_col)
  seu$celltype_plot <- rename_with_true_names(as.character(seu@meta.data[[celltype_col]]))
  unk_frac <- mean(seu$celltype_plot %in% c("Unknown Cell Type", "unknown1", "unknown2", "unknown3"), na.rm = TRUE)
  log_msg("Canonical celltype coverage (non-unknown): ", round(100 * (1 - unk_frac), 2), "%")
  cell_levels <- sort(unique(as.character(seu$celltype_plot)))
  cell_cols <- get_universal_colors(cell_levels)
  p_umap_celltype <- DimPlot(seu, reduction = "umap_harmony", group.by = "celltype_plot", label = TRUE, pt.size = 0.2) +
    scale_color_manual(values = cell_cols, drop = FALSE) +
    ggtitle(paste0("UMAP (Harmony) by Cell Type (", celltype_col, " → canonical)")) +
    theme_thesis_spatial()
  p_tsne_celltype <- DimPlot(seu, reduction = "tsne_harmony", group.by = "celltype_plot", label = TRUE, pt.size = 0.2) +
    scale_color_manual(values = cell_cols, drop = FALSE) +
    ggtitle(paste0("t-SNE (Harmony) by Cell Type (", celltype_col, " → canonical)")) +
    theme_thesis_spatial()
  p_ct <- p_umap_celltype + p_tsne_celltype
  ct_pdf <- file.path(DIR_FIGURES, "01c_embeddings_celltype_umap_tsne.pdf")
  ct_png <- file.path(DIR_FIGURES, "01c_embeddings_celltype_umap_tsne.png")
  ggsave(ct_pdf, p_ct, width = 14, height = 6)
  ggsave(ct_png, p_ct, width = 14, height = 6, dpi = 320)
  out_files <- c(out_files, ct_pdf, ct_png)
}

plot_manifest_path <- file.path(DIR_REPORTS, "01c_plot_spatial_embeddings_manifest.json")
record_plot_manifest_01c(
  manifest_path = plot_manifest_path,
  source_data = obj_path_found,
  compute_script = "scripts/01_active_pipeline/01_preprocess_harmony_embeddings.R",
  plotting_script = "scripts/01_active_pipeline/01c_plot_spatial_embeddings.R",
  figures_output = out_files
)

log_msg("Saved plotting artifacts and manifest: ", plot_manifest_path)
