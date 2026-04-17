#!/usr/bin/env Rscript

# =============================================================================
# Script: 03c_plot_spatial_cellchat.R
# Purpose: Decoupled plotting for dual-architecture Spatial CellChat outputs.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(jsonlite)
})

source("R/spatial_color_themes.R")

PIPELINE_NAME <- "03c_plot_spatial_cellchat"
PIPELINE_VERSION <- "2.0.0"
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

bundle_paths <- c(
  file.path(DIR_OBJECTS, "03_spatial_juxtacrine_fast.rds"),
  file.path(DIR_OBJECTS, "03_spatial_paracrine_fast.rds"),
  file.path(DIR_OBJECTS, "03_spatial_juxtacrine_full.rds"),
  file.path(DIR_OBJECTS, "03_spatial_paracrine_full.rds")
)
resolved_paths <- vapply(bundle_paths, resolve_object_path, character(1))

fig_paths <- character(0)

for (obj_path in resolved_paths) {
  run_tag <- sub("\\.rds$", "", basename(obj_path))
  log_msg("Plotting bundle: ", run_tag)

  bundle <- readRDS(obj_path)
  cellchat <- bundle$cellchat

  idents <- levels(cellchat@idents)
  idents_lower <- tolower(idents)

  source_idx <- which(grepl("evt|extravillous", idents_lower))
  if (length(source_idx) == 0) source_idx <- 1

  target_idx <- which(grepl("macroph|immune|lymph|t cell|b cell|nk|monocyte|maternal", idents_lower))
  if (length(target_idx) == 0) {
    target_idx <- setdiff(seq_along(idents), source_idx)
    target_idx <- head(target_idx, 3)
  }

  group_size <- as.numeric(table(cellchat@idents))
  names(group_size) <- names(table(cellchat@idents))

  # Circle plot
  circle_pdf <- file.path(DIR_FIGURES, paste0(run_tag, "_circle.pdf"))
  circle_png <- file.path(DIR_FIGURES, paste0(run_tag, "_circle.png"))

  pdf(circle_pdf, width = 10, height = 10)
  net_circle_fn(
    cellchat@net$weight,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = paste0("Global Network: ", run_tag)
  )
  dev.off()

  png(circle_png, width = 3000, height = 3000, res = 300)
  net_circle_fn(
    cellchat@net$weight,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = paste0("Global Network: ", run_tag)
  )
  dev.off()

  # Bubble plot (EVT -> immune targets)
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
  source_data = resolved_paths,
  plotting_script = "scripts/01_active_pipeline/03c_plot_spatial_cellchat.R",
  figures_output = fig_paths,
  notes = c(
    "Decoupled visualization for 4 dual-architecture CellChat bundles",
    "Each bundle exports global network circle plot + EVT-to-immune bubble plot",
    "All figures saved as both PDF and PNG"
  ),
  hypothesis = "Vulnerable niches show architecture-specific communication programs with distinct contact and diffusion signatures.",
  methods_blurb = "Sequential plotting of juxtacrine/paracrine fast/full SpatialCellChat bundles with focused EVT-to-immune visualization.",
  thesis_aim = "Compare synaptic-local vs paracrine-long-range signaling around vulnerable maternal-fetal interface niches."
)

log_msg("Saved dual-architecture CellChat figures + manifest: ", manifest_path)
