#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(jsonlite)
})

source("R/spatial_color_themes.R")

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
if (!file.exists(input_obj)) {
  stop("Missing scored object: ", input_obj,
       "\nRun scripts/01_active_pipeline/02_spatial_misi_ido1_scoring.R first.")
}

log_msg("Loading scored object: ", input_obj)
seu <- readRDS(input_obj)

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
}

manifest_path <- file.path(DIR_REPORTS, "02c_plot_spatial_misi_manifest.json")
record_artifact_manifest(
  manifest_path = manifest_path,
  pipeline = PIPELINE_NAME,
  version = PIPELINE_VERSION,
  run_timestamp = RUN_TIMESTAMP,
  seed = RUN_SEED,
  source_data = input_obj,
  compute_script = "scripts/01_active_pipeline/02_spatial_misi_ido1_scoring.R",
  plotting_script = "scripts/01_active_pipeline/02c_plot_spatial_misi.R",
  figures_output = fig_paths
)

log_msg("Saved MISI plots and manifest: ", manifest_path)
