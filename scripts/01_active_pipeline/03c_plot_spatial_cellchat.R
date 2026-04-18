#!/usr/bin/env Rscript

# =============================================================================
# Script: 03c_plot_spatial_cellchat.R
# Purpose: Plot week-wise + merged Spatial CellChat outputs from Script 03.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(jsonlite)
})

source("R/celltype_dictionary.R")

PIPELINE_NAME <- "03c_plot_spatial_cellchat"
PIPELINE_VERSION <- "3.0.0"
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
  NA_character_
}

read_object <- function(path) {
  if (is.na(path)) return(NULL)
  if (grepl("\\.qs$", path)) {
    if (!requireNamespace("qs", quietly = TRUE)) {
      log_msg("qs file found but package 'qs' not installed: ", path, .level = "WARN")
      return(NULL)
    }
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
  plotting_script,
  figures_output,
  notes
) {
  manifest <- list(
    pipeline = pipeline,
    version = version,
    run_timestamp = run_timestamp,
    seed = seed,
    source_data = source_data,
    plotting_script = plotting_script,
    figures_output = figures_output,
    notes = notes
  )
  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
}

if (requireNamespace("SpatialCellChat", quietly = TRUE)) {
  net_circle_fn <- SpatialCellChat::netVisual_circle
  net_bubble_fn <- SpatialCellChat::netVisual_bubble
  log_msg("Using SpatialCellChat plotting backend.")
} else if (requireNamespace("CellChat", quietly = TRUE)) {
  net_circle_fn <- CellChat::netVisual_circle
  net_bubble_fn <- CellChat::netVisual_bubble
  log_msg("Using CellChat plotting backend.", .level = "WARN")
} else {
  stop("Neither SpatialCellChat nor CellChat is installed.")
}

base_candidates <- c(
  file.path(DIR_OBJECTS, "03_spatial_cellchat_W7.rds"),
  file.path(DIR_OBJECTS, "03_spatial_cellchat_W8-2.rds"),
  file.path(DIR_OBJECTS, "03_spatial_cellchat_W9.rds"),
  file.path(DIR_OBJECTS, "03_spatial_cellchat_W11.rds"),
  file.path(DIR_OBJECTS, "03_spatial_cellchat_merged.rds"),
  file.path(DIR_OBJECTS, "03b_spatial_cellchat_full_W7.rds"),
  file.path(DIR_OBJECTS, "03b_spatial_cellchat_full_W8-2.rds"),
  file.path(DIR_OBJECTS, "03b_spatial_cellchat_full_W9.rds"),
  file.path(DIR_OBJECTS, "03b_spatial_cellchat_full_W11.rds"),
  file.path(DIR_OBJECTS, "03b_spatial_cellchat_full_merged.rds")
)

resolved <- vapply(base_candidates, resolve_object_path, character(1))
resolved <- resolved[!is.na(resolved)]
if (length(resolved) == 0) stop("No week/merged CellChat objects found in output/objects.")

fig_paths <- character(0)

for (obj_path in resolved) {
  run_tag <- sub("\\.rds$|\\.qs$", "", basename(obj_path))
  log_msg("Plotting object: ", run_tag)

  cellchat <- read_object(obj_path)
  if (is.null(cellchat)) next

  idents <- levels(cellchat@idents)
  idents_lower <- tolower(idents)

  source_idx <- which(grepl("evt|extravillous", idents_lower))
  if (length(source_idx) == 0) source_idx <- 1

  target_idx <- which(grepl("macroph|immune|lymph|t cell|b cell|nk|monocyte|maternal", idents_lower))
  if (length(target_idx) == 0) {
    target_idx <- setdiff(seq_along(idents), source_idx)
    target_idx <- head(target_idx, 3)
  }

  net_weight <- cellchat@net$weight
  if (is.null(net_weight) || nrow(net_weight) == 0 || ncol(net_weight) == 0) {
    log_msg("Skipping ", run_tag, ": empty net$weight matrix.", .level = "WARN")
    next
  }

  group_size_full <- table(cellchat@idents)
  active_nodes <- rownames(net_weight)
  if (is.null(active_nodes) || length(active_nodes) == 0) {
    log_msg("Skipping ", run_tag, ": net$weight has no rownames.", .level = "WARN")
    next
  }
  # Align vertex weights exactly to active graph nodes to avoid igraph length mismatch.
  group_size <- as.numeric(group_size_full[active_nodes])
  group_size[is.na(group_size)] <- 1

  active_colors <- get_universal_colors(active_nodes)

  circle_pdf <- file.path(DIR_FIGURES, paste0(run_tag, "_circle.pdf"))
  circle_png <- file.path(DIR_FIGURES, paste0(run_tag, "_circle.png"))

  pdf(circle_pdf, width = 10, height = 10)
  net_circle_fn(
    net_weight,
    vertex.weight = group_size,
    color.use = active_colors,
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = paste0("Global Network: ", run_tag)
  )
  dev.off()

  png(circle_png, width = 3000, height = 3000, res = 300)
  net_circle_fn(
    net_weight,
    vertex.weight = group_size,
    color.use = active_colors,
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = paste0("Global Network: ", run_tag)
  )
  dev.off()

  bubble_pdf <- file.path(DIR_FIGURES, paste0(run_tag, "_evt_to_immune_bubble.pdf"))
  bubble_png <- file.path(DIR_FIGURES, paste0(run_tag, "_evt_to_immune_bubble.png"))

  pdf(bubble_pdf, width = 14, height = 8)
  net_bubble_fn(
    object = cellchat,
    sources.use = source_idx,
    targets.use = target_idx,
    remove.isolate = FALSE,
    title.name = paste0("EVT to Immune Signaling: ", run_tag)
  )
  dev.off()

  png(bubble_png, width = 4200, height = 2400, res = 300)
  net_bubble_fn(
    object = cellchat,
    sources.use = source_idx,
    targets.use = target_idx,
    remove.isolate = FALSE,
    title.name = paste0("EVT to Immune Signaling: ", run_tag)
  )
  dev.off()

  fig_paths <- c(fig_paths, circle_pdf, circle_png, bubble_pdf, bubble_png)
}

manifest_path <- file.path(DIR_REPORTS, "03c_plot_spatial_cellchat_manifest.json")
record_artifact_manifest(
  manifest_path = manifest_path,
  pipeline = PIPELINE_NAME,
  version = PIPELINE_VERSION,
  run_timestamp = RUN_TIMESTAMP,
  seed = RUN_SEED,
  source_data = resolved,
  plotting_script = "scripts/01_active_pipeline/03c_plot_spatial_cellchat.R",
  figures_output = fig_paths,
  notes = c(
    "Plots week-wise and merged CellChat objects produced by Script 03",
    "Supports both .rds and .qs object inputs",
    "Exports circle and EVT->immune bubble plots as PDF+PNG"
  )
)

log_msg("Saved CellChat figures + manifest: ", manifest_path)
