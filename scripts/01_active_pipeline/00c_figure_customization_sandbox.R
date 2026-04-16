#!/usr/bin/env Rscript

# =============================================================================
# Script: 00c_figure_customization_sandbox.R
# Purpose: Master visualization sandbox / tutorial for thesis figure customization.
#
# DESIGN INTENT (READ THIS FIRST)
# - This script is intentionally comment-heavy.
# - It demonstrates "world-class" thesis-oriented figures from a pre-scored object.
# - It does NOT recompute SCTransform/Harmony/UMAP/scores. It only loads and plots.
# - Every figure is tracked in a manifest to guarantee reproducibility.
#
# HOW TO USE THIS SANDBOX
# 1) Run as-is to generate baseline figures.
# 2) Change one aesthetic at a time (color hex, point size, alpha, facet layout).
# 3) Re-run and compare outputs in output/figures.
# 4) Keep winning variants and note rationale in thesis methods notes.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(jsonlite)
  library(forcats)
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

record_artifact_manifest <- function(
  manifest_path,
  source_data,
  plotting_script,
  figures_output,
  notes,
  hypothesis,
  methods_blurb,
  thesis_aim
) {
  manifest <- list(
    source_data = source_data,
    plotting_script = plotting_script,
    figures_output = figures_output,
    notes = notes,
    hypothesis = hypothesis,
    methods_blurb = methods_blurb,
    thesis_aim = thesis_aim,
    run_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
}

# -----------------------------------------------------------------------------
# Load scored object and enforce temporal ordering
# -----------------------------------------------------------------------------
obj_path <- resolve_object_path(file.path(DIR_OBJECTS, "02_scored_misi_ido1.rds"))
log_msg("Loading scored object: ", obj_path)
seu <- read_object(obj_path)

seu$week <- factor(as.character(seu$week), levels = c("W7", "W8-2", "W9", "W11"))
seu$celltype_plot <- rename_with_true_names(as.character(seu$predicted.celltype))

# -----------------------------------------------------------------------------
# Plot 1: World-class Figure Idea A
# "Temporal burden map": which cell types have highest/lowest MISI by week?
# -----------------------------------------------------------------------------
# WHY INFORMATIVE:
# - Directly answers "which biological compartments are most vulnerable over time?"
# - Useful as a central thesis panel or summary figure.
#
# TWEAK GUIDE:
# - Change fill gradient with scale_fill_gradient2(low=...,mid=...,high=...)
# - Change label text size in geom_text(size=...)
score_tbl <- seu@meta.data %>%
  dplyr::group_by(week, celltype_plot) %>%
  dplyr::summarise(
    mean_misi = mean(MISI_Vulnerability, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(n >= 30)

p_heat_celltype_week <- ggplot(score_tbl, aes(x = week, y = fct_reorder(celltype_plot, mean_misi, .fun = max), fill = mean_misi)) +
  geom_tile() +
  scale_fill_gradient2(low = "#0B1F5E", mid = "white", high = "#CB181D", midpoint = median(score_tbl$mean_misi, na.rm = TRUE)) +
  labs(title = "Mean MISI Vulnerability by Cell Type and Week", x = "Week", y = "Cell Type", fill = "Mean MISI") +
  theme_thesis_spatial()

# -----------------------------------------------------------------------------
# Plot 2: World-class Figure Idea B
# "Top-vs-bottom celltype ranking": strongest contrasts for narrative framing
# -----------------------------------------------------------------------------
# WHY INFORMATIVE:
# - Gives a clean "leaderboard" of vulnerable vs protected niches.
# - Easy to interpret in slides and defense Q&A.
#
# TWEAK GUIDE:
# - Change bar color with fill = "#..."
# - Change number of shown groups by head(..., n = 10)
rank_tbl <- seu@meta.data %>%
  dplyr::group_by(celltype_plot) %>%
  dplyr::summarise(mean_misi = mean(MISI_Vulnerability, na.rm = TRUE), n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n >= 30) %>%
  dplyr::arrange(dplyr::desc(mean_misi))

p_top <- ggplot(head(rank_tbl, 10), aes(x = fct_reorder(celltype_plot, mean_misi), y = mean_misi)) +
  geom_col(fill = "#CB181D") +
  coord_flip() +
  labs(title = "Top 10 Cell Types by Mean MISI Vulnerability", x = "Cell Type", y = "Mean MISI") +
  theme_thesis_spatial()

p_bottom <- ggplot(tail(rank_tbl, 10), aes(x = fct_reorder(celltype_plot, mean_misi), y = mean_misi)) +
  geom_col(fill = "#0B1F5E") +
  coord_flip() +
  labs(title = "Bottom 10 Cell Types by Mean MISI Vulnerability", x = "Cell Type", y = "Mean MISI") +
  theme_thesis_spatial()

# -----------------------------------------------------------------------------
# Plot 3: World-class Figure Idea C
# Physical tissue map: top quartile vulnerability overlaid on full tissue
# -----------------------------------------------------------------------------
# WHY INFORMATIVE:
# - Grounds abstract score into real anatomical context.
# - Highlights putative "toxic-switch" microdomains for biological interpretation.
#
# TWEAK GUIDE:
# - Change highlight color by editing "#FF0000" below.
# - Increase highlighted dot size by changing size=0.35.
q75 <- stats::quantile(seu$MISI_Vulnerability, probs = 0.75, na.rm = TRUE)
spatial_df <- seu@meta.data %>%
  dplyr::mutate(is_top = MISI_Vulnerability >= q75)

p_spatial_niche <- ggplot(spatial_df, aes(x = x_um, y = y_um)) +
  geom_point(color = "grey85", size = 0.12, alpha = 0.6) +
  geom_point(data = subset(spatial_df, is_top), color = "#FF0000", size = 0.35, alpha = 0.95) +
  facet_wrap(~week, ncol = 4) +
  coord_fixed() +
  labs(
    title = "Top Quartile MISI Vulnerability in Physical Tissue Space",
    subtitle = "Gray = all cells, Red = top 25% vulnerable cells",
    x = "x_um", y = "y_um"
  ) +
  theme_thesis_spatial()

# -----------------------------------------------------------------------------
# Plot 4: Cluster vulnerability profiling (violin by Seurat cluster)
# -----------------------------------------------------------------------------
cluster_tbl <- seu@meta.data %>%
  dplyr::mutate(seurat_clusters = as.factor(seurat_clusters))

p_vln_cluster_misi <- ggplot(cluster_tbl, aes(x = seurat_clusters, y = MISI_Vulnerability, fill = seurat_clusters)) +
  geom_violin(scale = "width", color = "grey20", alpha = 0.9, linewidth = 0.15) +
  geom_boxplot(width = 0.10, outlier.size = 0.2, alpha = 0.55, color = "grey15") +
  labs(
    title = "Cluster Vulnerability Profiling",
    subtitle = "MISI_Vulnerability distribution by Seurat cluster",
    x = "Seurat Cluster", y = "MISI_Vulnerability"
  ) +
  theme_thesis_spatial() +
  theme(legend.position = "none")

# -----------------------------------------------------------------------------
# Plot 5: Spatial niche mapping for top-2 vulnerable Seurat clusters
# -----------------------------------------------------------------------------
top2_clusters <- cluster_tbl %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::summarise(mean_misi = mean(MISI_Vulnerability, na.rm = TRUE), n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n >= 20) %>%
  dplyr::arrange(dplyr::desc(mean_misi), dplyr::desc(n)) %>%
  dplyr::slice_head(n = 2) %>%
  dplyr::pull(seurat_clusters)

spatial_cluster_df <- cluster_tbl %>%
  dplyr::mutate(highlight_cluster = ifelse(seurat_clusters %in% top2_clusters, as.character(seurat_clusters), "Other")) %>%
  dplyr::mutate(highlight_cluster = factor(highlight_cluster, levels = c(as.character(top2_clusters), "Other")))

p_spatial_top2_clusters <- ggplot(spatial_cluster_df, aes(x = x_um, y = y_um)) +
  geom_point(data = subset(spatial_cluster_df, highlight_cluster == "Other"), color = "grey87", size = 0.12, alpha = 0.55) +
  geom_point(data = subset(spatial_cluster_df, highlight_cluster != "Other"), aes(color = highlight_cluster), size = 0.28, alpha = 0.95) +
  facet_wrap(~week, ncol = 4) +
  scale_color_manual(values = c("#CB181D", "#6A51A3")) +
  coord_fixed() +
  labs(
    title = "Spatial Niche Mapping of Top-2 Vulnerable Seurat Clusters",
    subtitle = "Background = all other clusters; color = highest-mean-MISI clusters",
    x = "x_um", y = "y_um", color = "Top Vulnerable\nClusters"
  ) +
  theme_thesis_spatial()

# -----------------------------------------------------------------------------
# Save all sandbox figures
# -----------------------------------------------------------------------------
out_files <- c(
  file.path(DIR_FIGURES, "00c_heatmap_celltype_week_misi.pdf"),
  file.path(DIR_FIGURES, "00c_heatmap_celltype_week_misi.png"),
  file.path(DIR_FIGURES, "00c_bar_top10_misi_celltypes.pdf"),
  file.path(DIR_FIGURES, "00c_bar_top10_misi_celltypes.png"),
  file.path(DIR_FIGURES, "00c_bar_bottom10_misi_celltypes.pdf"),
  file.path(DIR_FIGURES, "00c_bar_bottom10_misi_celltypes.png"),
  file.path(DIR_FIGURES, "00c_spatial_top_quartile_misi.pdf"),
  file.path(DIR_FIGURES, "00c_spatial_top_quartile_misi.png"),
  file.path(DIR_FIGURES, "00c_violin_misi_by_seurat_cluster.pdf"),
  file.path(DIR_FIGURES, "00c_violin_misi_by_seurat_cluster.png"),
  file.path(DIR_FIGURES, "00c_spatial_top2_vulnerable_clusters.pdf"),
  file.path(DIR_FIGURES, "00c_spatial_top2_vulnerable_clusters.png")
)

ggsave(out_files[1], p_heat_celltype_week, width = 10, height = 8)
ggsave(out_files[2], p_heat_celltype_week, width = 10, height = 8, dpi = 320)
ggsave(out_files[3], p_top, width = 8, height = 6)
ggsave(out_files[4], p_top, width = 8, height = 6, dpi = 320)
ggsave(out_files[5], p_bottom, width = 8, height = 6)
ggsave(out_files[6], p_bottom, width = 8, height = 6, dpi = 320)
ggsave(out_files[7], p_spatial_niche, width = 16, height = 6)
ggsave(out_files[8], p_spatial_niche, width = 16, height = 6, dpi = 320)
ggsave(out_files[9], p_vln_cluster_misi, width = 11, height = 6)
ggsave(out_files[10], p_vln_cluster_misi, width = 11, height = 6, dpi = 320)
ggsave(out_files[11], p_spatial_top2_clusters, width = 16, height = 6)
ggsave(out_files[12], p_spatial_top2_clusters, width = 16, height = 6, dpi = 320)

manifest_path <- file.path(DIR_REPORTS, "00c_figure_customization_sandbox_manifest.json")
record_artifact_manifest(
  manifest_path = manifest_path,
  source_data = obj_path,
  plotting_script = "scripts/01_active_pipeline/00c_figure_customization_sandbox.R",
  figures_output = out_files,
  notes = c(
    "Decoupled visualization-only sandbox",
    "World-class thesis plots: temporal burden heatmap, top/bottom rankings, vulnerability niche overlay",
    "Cluster vulnerability profiling and top-2 vulnerable cluster spatial mapping",
    "Tune colors/sizes directly in this script for final figure polishing"
  ),
  hypothesis = "MISI vulnerability is not random but concentrated in specific cell-state clusters.",
  methods_blurb = "Extract top quantiles of MISI scores mapped to physical tissue coordinates.",
  thesis_aim = "Identify and visualize spatially localized, high-risk trophoblast-immune niches that support the toxic-switch hypothesis."
)

log_msg("Sandbox figures + manifest saved: ", manifest_path)
