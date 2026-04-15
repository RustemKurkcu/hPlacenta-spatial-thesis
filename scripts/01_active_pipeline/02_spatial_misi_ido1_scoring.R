#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(jsonlite)
})
source("R/celltype_dictionary.R")

if (!requireNamespace("UCell", quietly = TRUE)) {
  stop("Package 'UCell' is required. Install with: remotes::install_github('carmonalab/UCell')")
}

PIPELINE_NAME <- "02_spatial_misi_ido1_scoring"
PIPELINE_VERSION <- "1.0.0"
RUN_TIMESTAMP <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
RUN_SEED <- 42L
set.seed(RUN_SEED)

OUT_ROOT <- "output"
DIR_LOGS <- file.path(OUT_ROOT, "logs")
DIR_OBJECTS <- file.path(OUT_ROOT, "objects")
DIR_REPORTS <- file.path(OUT_ROOT, "reports")
for (d in c(DIR_LOGS, DIR_OBJECTS, DIR_REPORTS)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
LOG_FILE <- file.path(DIR_LOGS, "02_spatial_misi_ido1_scoring.log")

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

record_artifact_manifest <- function(
  manifest_path,
  pipeline,
  version,
  run_timestamp,
  seed,
  source_data,
  compute_script,
  plotting_script,
  object_output,
  reports_output,
  module_definitions,
  assay_used
) {
  manifest <- list(
    pipeline = pipeline,
    version = version,
    run_timestamp = run_timestamp,
    seed = seed,
    source_data = source_data,
    compute_script = compute_script,
    plotting_script = plotting_script,
    object_output = object_output,
    reports_output = reports_output,
    module_definitions = module_definitions,
    assay_used = assay_used
  )
  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
  invisible(manifest)
}

input_obj <- file.path(DIR_OBJECTS, "01_integrated_harmony_sct_4weeks.rds")
input_obj_found <- resolve_object_path(input_obj)
if (is.na(input_obj_found)) {
  stop("Missing input object: ", input_obj,
       "\nRun scripts/01_active_pipeline/01_preprocess_harmony_embeddings.R first.")
}

log_msg("Loading input object: ", input_obj_found)
seu <- read_object(input_obj_found)
seu <- add_true_celltype_metadata(seu)
if ("celltype_true_name" %in% colnames(seu@meta.data)) {
  seu$celltype_corrected <- seu$celltype_true_name
}
log_msg("Added celltype columns: celltype_original, celltype_true_name, cell_label_display, celltype_corrected")

if (!"RNA" %in% names(seu@assays)) {
  stop("RNA assay not found. Available assays: ", paste(names(seu@assays), collapse = ", "))
}
DefaultAssay(seu) <- "RNA"
log_msg("Default assay set to RNA for interpretable biological scoring.")

signatures <- list(
  MMP_ECM_Remodeling = c("MMP2", "MMP7", "MMP14", "CTSB", "TWIST1", "TWIST2", "VIM", "FN1"),
  Immune_Exhaustion_Tolerogenic = c("TIGIT", "PDCD1", "LAG3", "HAVCR2", "IDO1", "S100A9", "S100A8", "IL10", "ARG1"),
  MegL_Sulfide_Vulnerability = c("PTGS2", "WNT2", "HMOX1", "NFE2L2", "CBS", "CDKN1A", "GPX2"),
  IDO1_Tolerogenic_Shield = c("IDO1", "TGFB1", "IL10", "CD274", "PDCD1LG2", "HAVCR2", "LGALS9", "CD80", "CD86", "ENTPD1", "NT5E", "FOXP3")
)

log_msg("Scoring UCell signatures...")
if (exists("AddModuleScore_UCell", where = asNamespace("UCell"), inherits = FALSE)) {
  seu <- UCell::AddModuleScore_UCell(seu, features = signatures)
} else {
  seu <- UCell::ScoreSignatures_UCell(seu, features = signatures)
}

for (nm in names(signatures)) {
  uc_col <- paste0(nm, "_UCell")
  if (!uc_col %in% colnames(seu@meta.data)) {
    stop("Expected UCell output column missing: ", uc_col)
  }
  seu[[nm]] <- seu@meta.data[[uc_col]]
}

seu$MISI_Vulnerability <- with(
  seu@meta.data,
  MMP_ECM_Remodeling + Immune_Exhaustion_Tolerogenic + MegL_Sulfide_Vulnerability - IDO1_Tolerogenic_Shield
)

output_obj <- file.path(DIR_OBJECTS, "02_scored_misi_ido1.rds")
output_obj_saved <- write_object(seu, output_obj)
log_msg("Saved scored object: ", output_obj_saved)

manifest_path <- file.path(DIR_REPORTS, "02_spatial_misi_ido1_manifest.json")
record_artifact_manifest(
  manifest_path = manifest_path,
  pipeline = PIPELINE_NAME,
  version = PIPELINE_VERSION,
  run_timestamp = RUN_TIMESTAMP,
  seed = RUN_SEED,
  source_data = input_obj_found,
  compute_script = "scripts/01_active_pipeline/02_spatial_misi_ido1_scoring.R",
  plotting_script = "scripts/01_active_pipeline/02c_plot_spatial_misi.R",
  object_output = output_obj_saved,
  reports_output = manifest_path,
  module_definitions = signatures,
  assay_used = "RNA"
)
log_msg("Saved manifest: ", manifest_path)
