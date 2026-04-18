#!/usr/bin/env Rscript

# =============================================================================
# Script: 03d_convert_rds_to_qs.R
# Purpose: Convert one or more .rds objects to .qs for faster I/O and smaller size.
# Usage:   Rscript scripts/01_active_pipeline/03d_convert_rds_to_qs.R [optional paths...]
#          If no paths are provided, converts common Script03 outputs in output/objects.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (!requireNamespace("qs", quietly = TRUE)) stop("Package 'qs' is required for conversion.")

default_targets <- c(
  "output/objects/03_spatial_cellchat_W7.rds",
  "output/objects/03_spatial_cellchat_W8-2.rds",
  "output/objects/03_spatial_cellchat_W9.rds",
  "output/objects/03_spatial_cellchat_W11.rds",
  "output/objects/03_spatial_cellchat_merged.rds",
  "output/objects/03b_spatial_cellchat_full_W7.rds",
  "output/objects/03b_spatial_cellchat_full_W8-2.rds",
  "output/objects/03b_spatial_cellchat_full_W9.rds",
  "output/objects/03b_spatial_cellchat_full_W11.rds",
  "output/objects/03b_spatial_cellchat_full_merged.rds"
)

targets <- if (length(args) > 0) args else default_targets
targets <- targets[file.exists(targets)]
if (length(targets) == 0) stop("No existing .rds files found to convert.")

for (p in targets) {
  obj <- readRDS(p)
  out_qs <- sub("\\.rds$", ".qs", p)
  qs::qsave(obj, out_qs, preset = "high")
  cat("Converted:", p, "->", out_qs, "\n")
}

cat("Done. Converted", length(targets), "file(s).\n")
