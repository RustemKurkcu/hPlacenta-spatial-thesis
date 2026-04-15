#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

source("R/celltype_dictionary.R")

OUT_ROOT <- "output"
DIR_OBJECTS <- file.path(OUT_ROOT, "objects")
if (!dir.exists(DIR_OBJECTS)) dir.create(DIR_OBJECTS, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "[INFO]", paste0(..., collapse = ""), "\n")
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

write_object <- function(object, path_rds) {
  if (requireNamespace("qs", quietly = TRUE)) {
    path_qs <- sub("\\.rds$", ".qs", path_rds)
    qs::qsave(object, path_qs, preset = "high")
    return(path_qs)
  }
  saveRDS(object, path_rds)
  path_rds
}

# Fast post-hoc metadata augmentation (no SCTransform/PCA/Harmony recomputation)
input_candidates <- c(
  file.path(DIR_OBJECTS, "02_scored_misi_ido1.rds"),
  file.path(DIR_OBJECTS, "01_integrated_harmony_sct_4weeks.rds")
)
input_obj <- NA_character_
for (cand in input_candidates) {
  found <- resolve_object_path(cand)
  if (!is.na(found)) {
    input_obj <- found
    break
  }
}
if (is.na(input_obj)) {
  stop("No input object found. Expected one of: ", paste(input_candidates, collapse = ", "))
}

log_msg("Loading object: ", input_obj)
seu <- read_object(input_obj)
n_before <- ncol(seu@meta.data)

before_pred <- if ("predicted.celltype" %in% colnames(seu@meta.data)) sum(!is.na(seu@meta.data$predicted.celltype)) else 0
seu <- augment_seurat_with_predicted_celltypes(seu)
after_pred <- if ("predicted.celltype" %in% colnames(seu@meta.data)) sum(!is.na(seu@meta.data$predicted.celltype)) else 0
log_msg("Predicted celltype coverage: ", before_pred, " -> ", after_pred)

seu <- add_true_celltype_metadata(seu)
n_after <- ncol(seu@meta.data)
log_msg("Added/updated celltype metadata columns. meta.data cols: ", n_before, " -> ", n_after)

output_obj <- sub("\\.(rds|qs)$", "_with_celltypes.rds", input_obj)
saved <- write_object(seu, output_obj)
log_msg("Saved object with celltype metadata: ", saved)
log_msg("Done. This step only updates metadata labels and does not recompute embeddings/scores.")
