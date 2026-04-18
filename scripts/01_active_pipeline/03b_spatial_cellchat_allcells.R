#!/usr/bin/env Rscript

# =============================================================================
# Script: 03b_spatial_cellchat_allcells.R
# Purpose: Run Spatial CellChat across all cells with conservative memory settings.
# Notes:   One-worker profile for high-memory hosts and stable long runs.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(jsonlite)
})

source("R/celltype_dictionary.R")

PIPELINE_NAME <- "03b_spatial_cellchat_allcells"
PIPELINE_VERSION <- "1.0.0"
RUN_TIMESTAMP <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
RUN_SEED <- 42L
set.seed(RUN_SEED)

OUT_ROOT <- "output"
DIR_LOGS <- file.path(OUT_ROOT, "logs")
DIR_OBJECTS <- file.path(OUT_ROOT, "objects")
DIR_REPORTS <- file.path(OUT_ROOT, "reports")
for (d in c(DIR_LOGS, DIR_OBJECTS, DIR_REPORTS)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
LOG_FILE <- file.path(DIR_LOGS, "03b_spatial_cellchat_allcells.log")

log_msg <- function(..., .level = "INFO") {
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  txt <- paste0(..., collapse = "")
  line <- paste0(stamp, " [", .level, "] ", txt)
  cat(line, "\n")
  cat(line, "\n", file = LOG_FILE, append = TRUE)
}

configure_future_runtime <- function(max_size_gb = 80, workers = 1) {
  old_plan <- NULL
  old_max <- getOption("future.globals.maxSize")
  if (requireNamespace("future", quietly = TRUE)) {
    old_plan <- future::plan()
    options(future.globals.maxSize = as.numeric(max_size_gb) * 1024^3)
    future::plan(future::sequential)
    log_msg("Configured future: plan=sequential, requested_workers=", workers, ", maxSizeGiB=", max_size_gb)
  }
  list(old_plan = old_plan, old_max = old_max)
}

restore_future_runtime <- function(state) {
  if (!requireNamespace("future", quietly = TRUE)) return(invisible(NULL))
  if (!is.null(state$old_max)) options(future.globals.maxSize = state$old_max)
  if (!is.null(state$old_plan)) future::plan(state$old_plan)
  invisible(NULL)
}

compute_min_nonzero_distance <- function(coords_mat) {
  xy <- as.matrix(coords_mat[, c("x_cent", "y_cent"), drop = FALSE])
  nn <- RANN::nn2(data = xy, query = xy, k = 2)
  d <- nn$nn.dists[, 2]
  d <- d[is.finite(d) & d > 0]
  if (length(d) == 0) stop("Could not compute min non-zero distance.")
  min(d)
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

get_normalized_assay_data <- function(seu, assay = "RNA", preferred_layers = c("data", "lognorm", "counts")) {
  is_valid_matrix <- function(x) !is.null(x) && (is.matrix(x) || inherits(x, "Matrix"))
  has_shape <- function(x) is_valid_matrix(x) && nrow(x) > 0 && ncol(x) > 0

  for (ly in preferred_layers) {
    mat <- tryCatch(SeuratObject::GetAssayData(seu, assay = assay, layer = ly), error = function(e) NULL)
    if (has_shape(mat)) return(list(mat = mat, source = paste0("layer:", ly), assay = assay))
  }

  for (sl in unique(c(preferred_layers, "data", "counts"))) {
    mat <- tryCatch(SeuratObject::GetAssayData(seu, assay = assay, slot = sl), error = function(e) NULL)
    if (has_shape(mat)) return(list(mat = mat, source = paste0("slot:", sl), assay = assay))
  }

  stop("Unable to extract assay data from assay '", assay, "'.")
}

extract_cellchat_input_matrix <- function(seu) {
  if ("SCT" %in% Assays(seu)) {
    log_msg("Detected SCT assay. Using SCT normalized data.")
    DefaultAssay(seu) <- "SCT"
    payload <- get_normalized_assay_data(seu, assay = "SCT", preferred_layers = c("data"))
    data_input <- payload$mat
  } else {
    DefaultAssay(seu) <- "RNA"
    payload <- get_normalized_assay_data(seu, assay = "RNA")
    data_input <- payload$mat
    if (!grepl("data|lognorm", payload$source)) {
      if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package required.")
      lib_size <- Matrix::colSums(data_input)
      lib_size[lib_size == 0] <- 1
      data_input <- Matrix::t(Matrix::t(data_input) / lib_size) * 1e4
      data_input <- log1p(data_input)
    }
  }
  list(seu = seu, data_input = data_input)
}

record_artifact_manifest <- function(manifest_path, source_data, output_object, notes) {
  manifest <- list(
    pipeline = PIPELINE_NAME,
    version = PIPELINE_VERSION,
    run_timestamp = RUN_TIMESTAMP,
    seed = RUN_SEED,
    source_data = source_data,
    output_object = output_object,
    script_path = "scripts/01_active_pipeline/03b_spatial_cellchat_allcells.R",
    notes = notes
  )
  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
}

input_obj <- resolve_object_path(file.path(DIR_OBJECTS, "02_scored_misi_ido1.rds"))
seu <- read_object(input_obj)

celltype_col <- if ("predicted.celltype" %in% colnames(seu@meta.data)) "predicted.celltype" else pick_celltype_source_column(seu@meta.data)
if (is.na(celltype_col)) celltype_col <- "seurat_clusters"
seu$celltype_plot <- rename_with_true_names(as.character(seu@meta.data[[celltype_col]]))
Idents(seu) <- "celltype_plot"

ext <- extract_cellchat_input_matrix(seu)
seu <- ext$seu
data_input <- ext$data_input

coords <- as.matrix(seu@meta.data[, c("x_um", "y_um"), drop = FALSE])
colnames(coords) <- c("x_cent", "y_cent")
rownames(coords) <- colnames(seu)
meta <- seu@meta.data %>% dplyr::mutate(group = as.character(celltype_plot))
rownames(meta) <- colnames(seu)

using_spatial_cellchat <- requireNamespace("SpatialCellChat", quietly = TRUE)
db_human <- tryCatch(SpatialCellChat::CellChatDB.human, error = function(e) CellChat::CellChatDB.human)

if (using_spatial_cellchat) {
  create_fn <- SpatialCellChat::createSpatialCellChat
  subset_fn <- SpatialCellChat::subsetData
  over_gene_fn <- SpatialCellChat::identifyOverExpressedGenes
  over_inter_fn <- SpatialCellChat::identifyOverExpressedInteractions
  commprob_fn <- SpatialCellChat::computeCommunProb
  filter_fn <- SpatialCellChat::filterCommunication
  pathway_fn <- SpatialCellChat::computeCommunProbPathway
  aggregate_fn <- SpatialCellChat::aggregateNet
  subset_db_fn <- SpatialCellChat::subsetDB
} else {
  create_fn <- CellChat::createCellChat
  subset_fn <- CellChat::subsetData
  over_gene_fn <- CellChat::identifyOverExpressedGenes
  over_inter_fn <- CellChat::identifyOverExpressedInteractions
  commprob_fn <- CellChat::computeCommunProb
  filter_fn <- CellChat::filterCommunication
  pathway_fn <- CellChat::computeCommunProbPathway
  aggregate_fn <- CellChat::aggregateNet
  subset_db_fn <- CellChat::subsetDB
}

create_args <- list(object = data_input, meta = meta, group.by = "group", datatype = "spatial", coordinates = coords)
create_formals <- names(formals(create_fn))
if ("spatial.factors" %in% create_formals) create_args$spatial.factors <- list(ratio = 1, tol = 0)
if ("scale.factors" %in% create_formals) create_args$scale.factors <- list(spot = 1, spot.diameter = 1)

cellchat_base <- do.call(create_fn, create_args)
cellchat_base@DB <- db_human
comm_formals <- names(formals(commprob_fn))
dynamic_scale_distance <- 1 / compute_min_nonzero_distance(coords)

future_state <- configure_future_runtime(max_size_gb = 80, workers = 1)
on.exit(restore_future_runtime(future_state), add = TRUE)

run_all_cells <- function(architecture_name, db_search, interaction_range, output_object, manifest_out) {
  cellchat <- cellchat_base
  cellchat@DB <- subset_db_fn(db_human, search = db_search)
  if (nrow(cellchat@DB$interaction) == 0) stop("Filtered DB is empty.")
  cellchat@DB$interaction$annotation <- "Secreted Signaling"

  cellchat <- subset_fn(cellchat)
  cellchat <- over_gene_fn(cellchat)
  cellchat <- over_inter_fn(cellchat)

  comm_args <- list(
    object = cellchat,
    distance.use = TRUE,
    interaction.range = interaction_range,
    contact.dependent = FALSE
  )
  if ("scale.distance" %in% comm_formals) comm_args$scale.distance <- dynamic_scale_distance
  if ("type" %in% comm_formals) comm_args$type <- "truncatedMean"
  if ("trim" %in% comm_formals) comm_args$trim <- 0

  cellchat <- do.call(commprob_fn, comm_args)
  gc(verbose = FALSE)

  tryCatch({
    cellchat <- filter_fn(cellchat, min.cells = 8)
    cellchat <- pathway_fn(cellchat)
    cellchat <- aggregate_fn(cellchat)
  }, error = function(e) log_msg("WARN during aggregation: ", conditionMessage(e)))
  gc(verbose = FALSE)

  saveRDS(list(architecture = architecture_name, cellchat = cellchat), output_object)
  record_artifact_manifest(manifest_out, input_obj, output_object,
      c(paste0("architecture=", architecture_name), paste0("workers=1"), paste0("interaction.range=", interaction_range), "pathways=ALL"))
}

run_all_cells("all_cells_juxtacrine_full", c("Cell-Cell Contact", "ECM-Receptor"), 10, file.path(DIR_OBJECTS, "03b_spatial_allcells_juxtacrine_full.rds"), file.path(DIR_REPORTS, "03b_manifest_allcells_juxt_full.json"))
run_all_cells("all_cells_paracrine_full", "Secreted Signaling", 100, file.path(DIR_OBJECTS, "03b_spatial_allcells_paracrine_full.rds"), file.path(DIR_REPORTS, "03b_manifest_allcells_para_full.json"))

log_msg("All-cells run complete.")
