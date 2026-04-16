#!/usr/bin/env Rscript

# =============================================================================
# Script: 03c_plot_spatial_cellchat.R
# Purpose: Visualization-only script for spatial CellChat outputs.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(jsonlite)
})

source("R/spatial_color_themes.R")
source("R/celltype_dictionary.R")

PIPELINE_NAME <- "03c_plot_spatial_cellchat"
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
LOG_FILE <- file.path(DIR_LOGS, "03c_plot_spatial_cellchat.log")

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

record_artifact_manifest <- function(
  manifest_path,
  pipeline,
  version,
  run_timestamp,
  seed,
  source_data,
  plotting_script,
  figures_output,
  notes,
  hypothesis,
  methods_blurb,
  thesis_aim
) {
  manifest <- list(
    pipeline = pipeline,
    version = version,
    run_timestamp = run_timestamp,
    seed = seed,
    source_data = source_data,
    plotting_script = plotting_script,
    figures_output = figures_output,
    notes = notes,
    hypothesis = hypothesis,
    methods_blurb = methods_blurb,
    thesis_aim = thesis_aim
  )
  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
}

pathway_activity_from_cellchat <- function(cellchat, pathways) {
  arr <- cellchat@netP$prob
  if (is.null(arr) || length(dim(arr)) != 3) return(list())

  pathway_names <- dimnames(arr)[[3]]
  group_names <- dimnames(arr)[[1]]
  out <- vector("list", length(pathways))
  names(out) <- pathways

  for (p in pathways) {
    idx <- which(tolower(pathway_names) == tolower(p))
    if (length(idx) == 0) {
      out[[p]] <- NULL
      next
    }
    mat <- arr[, , idx[1], drop = FALSE][, , 1]
    sender <- rowSums(mat, na.rm = TRUE)
    receiver <- colSums(mat, na.rm = TRUE)
    total <- sender + receiver
    names(total) <- group_names
    out[[p]] <- total
  }
  out
}

input_obj <- resolve_object_path(file.path(DIR_OBJECTS, "03_spatial_cellchat.rds"))
manifest_path <- file.path(DIR_REPORTS, "03c_plot_spatial_cellchat_manifest.json")

log_msg("Loading spatial CellChat bundle: ", input_obj)
bundle <- readRDS(input_obj)
cellchat <- bundle$cellchat
spatial_df <- bundle$spatial_meta
spatial_df$week <- factor(as.character(spatial_df$week), levels = c("W7", "W8-2", "W9", "W11"))

if (requireNamespace("SpatialCellChat", quietly = TRUE)) {
  net_visual_circle_fn <- SpatialCellChat::netVisual_circle
  log_msg("Using SpatialCellChat visualization backend.")
} else if (requireNamespace("CellChat", quietly = TRUE)) {
  net_visual_circle_fn <- CellChat::netVisual_circle
  log_msg("Using CellChat visualization backend.", .level = "WARN")
} else {
  stop("Neither 'SpatialCellChat' nor 'CellChat' is installed for visualization.")
}

# -----------------------------------------------------------------------------
# Cell-type reference spatial map
# -----------------------------------------------------------------------------
cell_levels <- sort(unique(as.character(spatial_df$celltype_plot)))
cell_cols <- get_universal_colors(cell_levels)

p_ref <- ggplot(spatial_df, aes(x = x_um, y = y_um, color = celltype_plot)) +
  geom_point(size = 0.18, alpha = 0.85) +
  facet_wrap(~week, ncol = 4) +
  scale_color_manual(values = cell_cols, drop = FALSE) +
  coord_fixed() +
  labs(title = "Cell Type Reference (Physical Tissue Space)", x = "x_um", y = "y_um", color = "Cell Type") +
  theme_thesis_spatial() +
  theme(legend.position = "bottom")

# -----------------------------------------------------------------------------
# Global communication circle plot (export only; base graphics)
# -----------------------------------------------------------------------------
fig_paths <- character(0)
group_size <- as.numeric(table(cellchat@idents))
names(group_size) <- names(table(cellchat@idents))

circle_pdf <- file.path(DIR_FIGURES, "03c_global_communication_circle.pdf")
circle_png <- file.path(DIR_FIGURES, "03c_global_communication_circle.png")

pdf(circle_pdf, width = 10, height = 10)
net_visual_circle_fn(
  cellchat@net$weight,
  vertex.weight = group_size,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Global Spatial Communication Network"
)
dev.off()

png(circle_png, width = 3000, height = 3000, res = 300)
net_visual_circle_fn(
  cellchat@net$weight,
  vertex.weight = group_size,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Global Spatial Communication Network"
)
dev.off()
fig_paths <- c(fig_paths, circle_pdf, circle_png)

# -----------------------------------------------------------------------------
# Side-by-side pathway overlays (reference + pathway activity proxy)
# -----------------------------------------------------------------------------
target_pathways <- c("HLA", "SPP1", "VEGF", "TGFb")
pathway_scores_by_group <- pathway_activity_from_cellchat(
  cellchat = cellchat,
  pathways = target_pathways
)

for (pw in target_pathways) {
  group_score_map <- pathway_scores_by_group[[pw]]
  if (is.null(group_score_map)) {
    log_msg("Pathway not found in CellChat netP: ", pw, .level = "WARN")
    next
  }

  overlay_df <- spatial_df %>%
    dplyr::mutate(pathway_score = unname(group_score_map[as.character(celltype_plot)])) %>%
    dplyr::mutate(pathway_score = ifelse(is.na(pathway_score), 0, pathway_score))

  p_overlay <- ggplot(overlay_df, aes(x = x_um, y = y_um, color = pathway_score)) +
    geom_point(size = 0.20, alpha = 0.9) +
    facet_wrap(~week, ncol = 4) +
    coord_fixed() +
    scale_color_gradient(low = "grey93", high = "#CB181D") +
    labs(
      title = paste0("Spatial Pathway Overlay: ", pw),
      subtitle = "Higher color intensity indicates stronger pathway communication context",
      x = "x_um", y = "y_um", color = "Pathway\nScore"
    ) +
    theme_thesis_spatial()

  p_pair <- p_ref + p_overlay + patchwork::plot_layout(widths = c(1.15, 1))

  out_pdf <- file.path(DIR_FIGURES, paste0("03c_", pw, "_reference_plus_overlay.pdf"))
  out_png <- file.path(DIR_FIGURES, paste0("03c_", pw, "_reference_plus_overlay.png"))
  ggsave(out_pdf, p_pair, width = 20, height = 7)
  ggsave(out_png, p_pair, width = 20, height = 7, dpi = 300)
  fig_paths <- c(fig_paths, out_pdf, out_png)
}

record_artifact_manifest(
  manifest_path = manifest_path,
  pipeline = PIPELINE_NAME,
  version = PIPELINE_VERSION,
  run_timestamp = RUN_TIMESTAMP,
  seed = RUN_SEED,
  source_data = input_obj,
  plotting_script = "scripts/01_active_pipeline/03c_plot_spatial_cellchat.R",
  figures_output = fig_paths,
  notes = c(
    "Decoupled visualization-only script for precomputed spatial CellChat object",
    "Global communication circle plot exported as PDF + PNG",
    "Patchwork side-by-side panels: cell-type reference + pathway spatial overlays"
  ),
  hypothesis = "Vulnerable spatial niches exhibit enhanced tolerogenic/exhaustion signaling to break down the IDO1 shield.",
  methods_blurb = "Spatial CellChat utilizing physical distance matrices to penalize long-range interactions.",
  thesis_aim = "Visualize niche-aware communication programs and localize pathogenic pathway activity in high-vulnerability tissue microdomains."
)

log_msg("Saved CellChat figures and manifest: ", manifest_path)
