#!/usr/bin/env Rscript

# =============================================================================
# Script: 02c_plot_spatial_misi.R
# Purpose: Decoupled visualization-only script for precomputed MISI scores.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(jsonlite)
  library(dplyr)
  library(patchwork)
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
for (d in c(DIR_LOGS, DIR_OBJECTS, DIR_FIGURES, DIR_REPORTS)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
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
  stop("Missing object: ", path_rds)
}

read_object <- function(path) {
  if (grepl("\\.qs$", path)) {
    if (!requireNamespace("qs", quietly = TRUE)) stop("Need package 'qs' to read: ", path)
    return(qs::qread(path))
  }
  readRDS(path)
}

record_plot_manifest_02c <- function(
  manifest_path,
  pipeline,
  version,
  run_timestamp,
  seed,
  source_data,
  plotting_script,
  figures_output
) {
  manifest <- list(
    pipeline = pipeline,
    version = version,
    run_timestamp = run_timestamp,
    seed = seed,
    source_data = source_data,
    plotting_script = plotting_script,
    figures_output = figures_output
  )
  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
}

input_obj <- resolve_object_path(file.path(DIR_OBJECTS, "02_scored_misi_ido1.rds"))
log_msg("Loading scored object: ", input_obj)
seu <- read_object(input_obj)

# Directive 1: force developmental order.
seu$week <- factor(as.character(seu$week), levels = c("W7", "W8-2", "W9", "W11"))

features_to_plot <- c(
  "MMP_ECM_Remodeling",
  "Immune_Exhaustion_Tolerogenic",
  "MegL_Sulfide_Vulnerability",
  "IDO1_Tolerogenic_Shield",
  "MISI_Vulnerability"
)

# -----------------------------------------------------------------------------
# A) Existing embedding feature plots (for consistency)
# -----------------------------------------------------------------------------
fig_paths <- character(0)
celltype_col <- if ("predicted.celltype" %in% colnames(seu@meta.data)) "predicted.celltype" else pick_celltype_source_column(seu@meta.data)
if (is.na(celltype_col)) celltype_col <- "seurat_clusters"
seu$celltype_plot <- rename_with_true_names(as.character(seu@meta.data[[celltype_col]]))

for (red in c("umap_harmony", "tsne_harmony")) {
  p_celltype <- DimPlot(
    object = seu,
    reduction = red,
    group.by = "celltype_plot",
    label = TRUE,
    repel = TRUE,
    pt.size = 0.2
  ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(
      title = paste0(toupper(sub("_harmony", "", red)), " Cell Type Map"),
      subtitle = "Reference map for interpreting signature hotspots",
      color = "Cell Type"
    ) +
    theme_thesis_spatial() +
    theme(legend.position = "bottom", legend.key.height = grid::unit(0.35, "cm"))

  for (feat in features_to_plot) {
    p_feat <- FeaturePlot(
      object = seu,
      features = feat,
      reduction = red,
      cols = c("grey92", "#CB181D"),
      pt.size = 0.25,
      order = TRUE
    ) +
      ggtitle(paste0(feat, " on ", red, " (score overlay)")) +
      theme_thesis_spatial()

    p_side_by_side <- p_celltype + p_feat + patchwork::plot_layout(widths = c(1, 1))

    base_single <- file.path(DIR_FIGURES, paste0("02c_", feat, "_", red))
    base_pair <- file.path(DIR_FIGURES, paste0("02c_", feat, "_", red, "_with_celltype"))

    ggsave(paste0(base_single, ".pdf"), p_feat, width = 8, height = 6)
    ggsave(paste0(base_single, ".png"), p_feat, width = 8, height = 6, dpi = 320)
    ggsave(paste0(base_pair, ".pdf"), p_side_by_side, width = 16, height = 6.5)
    ggsave(paste0(base_pair, ".png"), p_side_by_side, width = 16, height = 6.5, dpi = 320)

    fig_paths <- c(
      fig_paths,
      paste0(base_single, ".pdf"), paste0(base_single, ".png"),
      paste0(base_pair, ".pdf"), paste0(base_pair, ".png")
    )
  }
}

# -----------------------------------------------------------------------------
# B) Directive 3: Vulnerability niche combo plot in physical space
# -----------------------------------------------------------------------------
spatial_df <- seu@meta.data %>%
  dplyr::mutate(
    week = factor(as.character(week), levels = c("W7", "W8-2", "W9", "W11")),
    celltype_plot = rename_with_true_names(as.character(.data[[celltype_col]]))
  )

q75 <- stats::quantile(spatial_df$MISI_Vulnerability, probs = 0.75, na.rm = TRUE)
spatial_df$is_vulnerable_top25 <- spatial_df$MISI_Vulnerability >= q75

# Overlay plot: all cells in light gray; top-25% vulnerability in bright red.
p_vuln_overlay <- ggplot(spatial_df, aes(x = x_um, y = y_um)) +
  geom_point(color = "grey85", size = 0.15, alpha = 0.6) +
  geom_point(
    data = subset(spatial_df, is_vulnerable_top25),
    color = "#CB181D", size = 0.30, alpha = 0.95
  ) +
  facet_wrap(~week, ncol = 4) +
  coord_fixed() +
  labs(
    title = "Vulnerability Niche Overlay (Top 25% MISI_Vulnerability)",
    subtitle = "Background: all cells (gray); Highlight: top quartile vulnerability (red)",
    x = "x_um", y = "y_um"
  ) +
  theme_thesis_spatial()

# Secondary layer: reduce clutter by restricting to top 3 vulnerable cell types.
top3_celltypes <- spatial_df %>%
  dplyr::filter(is_vulnerable_top25) %>%
  dplyr::group_by(celltype_plot) %>%
  dplyr::summarise(mean_vuln = mean(MISI_Vulnerability, na.rm = TRUE), n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n >= 20) %>%
  dplyr::arrange(dplyr::desc(mean_vuln), dplyr::desc(n)) %>%
  dplyr::slice_head(n = 3) %>%
  dplyr::pull(celltype_plot)

vuln_top3_df <- spatial_df %>%
  dplyr::filter(is_vulnerable_top25, celltype_plot %in% top3_celltypes)

p_vuln_celltypes <- ggplot(vuln_top3_df, aes(x = x_um, y = y_um, color = celltype_plot)) +
  geom_point(size = 0.22, alpha = 0.9) +
  facet_grid(celltype_plot ~ week, scales = "fixed") +
  scale_color_manual(values = get_universal_colors(sort(unique(as.character(vuln_top3_df$celltype_plot))))) +
  coord_fixed() +
  labs(
    title = "Top-Quartile Vulnerable Cells: Top 3 Cell Types by Week",
    subtitle = "Decluttered niche view focused on the strongest vulnerable cell-type compartments",
    x = "x_um", y = "y_um", color = "Cell Type"
  ) +
  theme_thesis_spatial() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.text.y = element_text(size = 8),
    strip.text.x = element_text(size = 9)
  )

base_overlay <- file.path(DIR_FIGURES, "02c_vulnerability_niche_overlay")
base_combo <- file.path(DIR_FIGURES, "02c_vulnerability_niche_celltype_combo")
ggsave(paste0(base_overlay, ".pdf"), p_vuln_overlay, width = 16, height = 6)
ggsave(paste0(base_overlay, ".png"), p_vuln_overlay, width = 16, height = 6, dpi = 320)
ggsave(paste0(base_combo, ".pdf"), p_vuln_celltypes, width = 18, height = 14)
ggsave(paste0(base_combo, ".png"), p_vuln_celltypes, width = 18, height = 14, dpi = 320)
fig_paths <- c(fig_paths,
               paste0(base_overlay, ".pdf"), paste0(base_overlay, ".png"),
               paste0(base_combo, ".pdf"), paste0(base_combo, ".png"))

manifest_path <- file.path(DIR_REPORTS, "02c_plot_spatial_misi_manifest.json")
record_plot_manifest_02c(
  manifest_path = manifest_path,
  pipeline = PIPELINE_NAME,
  version = PIPELINE_VERSION,
  run_timestamp = RUN_TIMESTAMP,
  seed = RUN_SEED,
  source_data = input_obj,
  plotting_script = "scripts/01_active_pipeline/02c_plot_spatial_misi.R",
  figures_output = fig_paths
)

log_msg("Saved MISI plots and manifest: ", manifest_path)
