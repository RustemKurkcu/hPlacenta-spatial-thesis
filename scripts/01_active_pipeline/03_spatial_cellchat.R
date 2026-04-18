#!/usr/bin/env Rscript

# =============================================================================
# Script: 03_spatial_cellchat.R
# Purpose: Spatial CellChat list-and-merge architecture by week (W7/W8-2/W9/W11)
# Notes:   Computes each slide/week independently, then merges CellChat objects.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(jsonlite)
})

source("R/celltype_dictionary.R")

PIPELINE_NAME <- "03_spatial_cellchat"
PIPELINE_VERSION <- "3.0.0"
RUN_TIMESTAMP <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
RUN_SEED <- 42L
set.seed(RUN_SEED)

OUT_ROOT <- "output"
DIR_LOGS <- file.path(OUT_ROOT, "logs")
DIR_OBJECTS <- file.path(OUT_ROOT, "objects")
DIR_REPORTS <- file.path(OUT_ROOT, "reports")
for (d in c(DIR_LOGS, DIR_OBJECTS, DIR_REPORTS)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
LOG_FILE <- file.path(DIR_LOGS, "03_spatial_cellchat.log")

log_msg <- function(..., .level = "INFO") {
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  txt <- paste0(..., collapse = "")
  line <- paste0(stamp, " [", .level, "] ", txt)
  cat(line, "\n")
  cat(line, "\n", file = LOG_FILE, append = TRUE)
}

configure_future_runtime <- function(max_size_gb = 80, workers = 2) {
  old_plan <- NULL
  old_max <- getOption("future.globals.maxSize")
  if (requireNamespace("future", quietly = TRUE)) {
    old_plan <- future::plan()
    options(future.globals.maxSize = as.numeric(max_size_gb) * 1024^3)
    future::plan(future::multisession, workers = as.integer(workers))
    log_msg(
      "Configured future runtime: plan=multisession (workers=", workers,
      "), future.globals.maxSize=", round(getOption("future.globals.maxSize") / 1024^3, 2), " GiB."
    )
  } else {
    log_msg("Package 'future' not installed; skipping future runtime configuration.", .level = "WARN")
  }
  list(old_plan = old_plan, old_max = old_max)
}

restore_future_runtime <- function(state) {
  if (!requireNamespace("future", quietly = TRUE)) return(invisible(NULL))
  if (!is.null(state$old_max)) options(future.globals.maxSize = state$old_max)
  if (!is.null(state$old_plan)) future::plan(state$old_plan)
  invisible(NULL)
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
    log_msg("Detected SCT assay. Using SCT normalized data for CellChat.")
    DefaultAssay(seu) <- "SCT"
    payload <- get_normalized_assay_data(seu, assay = "SCT", preferred_layers = c("data"))
    data_input <- payload$mat
  } else {
    DefaultAssay(seu) <- "RNA"
    payload <- get_normalized_assay_data(seu, assay = "RNA")
    data_input <- payload$mat
    if (!grepl("data|lognorm", payload$source)) {
      log_msg("RNA layer:data is empty. Applying manual CPM + log1p normalization directly.", .level = "WARN")
      if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package required.")
      lib_size <- Matrix::colSums(data_input)
      lib_size[lib_size == 0] <- 1
      data_input <- Matrix::t(Matrix::t(data_input) / lib_size) * 1e4
      data_input <- log1p(data_input)
    }
  }
  list(seu = seu, data_input = data_input, payload = payload)
}

compute_min_nonzero_distance <- function(coords_mat) {
  xy <- as.matrix(coords_mat[, c("x_cent", "y_cent"), drop = FALSE])
  if (nrow(xy) < 2) stop("Need at least two cells to compute min distance.")

  if (requireNamespace("RANN", quietly = TRUE)) {
    nn <- RANN::nn2(data = xy, query = xy, k = 2)
    d <- nn$nn.dists[, 2]
    d <- d[is.finite(d) & d > 0]
    if (length(d) > 0) return(min(d))
  }

  dmat <- as.matrix(stats::dist(xy))
  dvec <- dmat[upper.tri(dmat)]
  dvec <- dvec[is.finite(dvec) & dvec > 0]
  if (length(dvec) == 0) stop("Could not compute min non-zero distance.")
  min(dvec)
}

record_artifact_manifest <- function(manifest_path, source_data, output_objects, notes = NULL) {
  manifest <- list(
    pipeline = PIPELINE_NAME,
    version = PIPELINE_VERSION,
    run_timestamp = RUN_TIMESTAMP,
    seed = RUN_SEED,
    source_data = source_data,
    output_objects = output_objects,
    script_path = "scripts/01_active_pipeline/03_spatial_cellchat.R",
    notes = notes
  )
  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
}

interaction_range_um <- 100

input_obj <- resolve_object_path(file.path(DIR_OBJECTS, "02_scored_misi_ido1.rds"))
log_msg("Loading scored object: ", input_obj)
seu <- read_object(input_obj)

if (!all(c("week", "x_um", "y_um") %in% colnames(seu@meta.data))) {
  stop("Seurat object must include week, x_um, y_um metadata.")
}

seu$week <- factor(as.character(seu$week), levels = c("W7", "W8-2", "W9", "W11"))
celltype_col <- if ("predicted.celltype" %in% colnames(seu@meta.data)) "predicted.celltype" else pick_celltype_source_column(seu@meta.data)
if (is.na(celltype_col)) celltype_col <- "seurat_clusters"
seu$celltype_plot <- rename_with_true_names(as.character(seu@meta.data[[celltype_col]]))
Idents(seu) <- "celltype_plot"

week_list <- SplitObject(seu, split.by = "week")
week_list <- week_list[names(week_list) %in% c("W7", "W8-2", "W9", "W11")]
if (length(week_list) == 0) stop("No week-specific subsets were created.")

using_spatial_cellchat <- requireNamespace("SpatialCellChat", quietly = TRUE)
db_human <- tryCatch(SpatialCellChat::CellChatDB.human, error = function(e) CellChat::CellChatDB.human)

if (using_spatial_cellchat) {
  log_msg("Detected SpatialCellChat package; using SpatialCellChat API (CellChat v3).")
  create_fn <- SpatialCellChat::createSpatialCellChat
  subset_fn <- SpatialCellChat::subsetData
  over_gene_fn <- SpatialCellChat::identifyOverExpressedGenes
  over_inter_fn <- SpatialCellChat::identifyOverExpressedInteractions
  commprob_fn <- SpatialCellChat::computeCommunProb
  filter_fn <- SpatialCellChat::filterCommunication
  pathway_fn <- SpatialCellChat::computeCommunProbPathway
  aggregate_fn <- SpatialCellChat::aggregateNet
  subset_db_fn <- SpatialCellChat::subsetDB
  # SpatialCellChat does not re-export mergeCellChat; use CellChat namespace.
  merge_fn <- CellChat::mergeCellChat
} else {
  if (!requireNamespace("CellChat", quietly = TRUE)) stop("Neither SpatialCellChat nor CellChat is installed.")
  log_msg("Using CellChat fallback API.", .level = "WARN")
  create_fn <- CellChat::createCellChat
  subset_fn <- CellChat::subsetData
  over_gene_fn <- CellChat::identifyOverExpressedGenes
  over_inter_fn <- CellChat::identifyOverExpressedInteractions
  commprob_fn <- CellChat::computeCommunProb
  filter_fn <- CellChat::filterCommunication
  pathway_fn <- CellChat::computeCommunProbPathway
  aggregate_fn <- CellChat::aggregateNet
  subset_db_fn <- CellChat::subsetDB
  merge_fn <- CellChat::mergeCellChat
}

create_formals <- names(formals(create_fn))
comm_formals <- names(formals(commprob_fn))
thesis_pathways <- c("MMP", "WNT", "TGFb", "CXCL", "CCL", "SPP1", "NOTCH", "FN1", "MHC-I", "MHC-II", "CD45", "TIGIT", "PD-L1")

future_state <- configure_future_runtime(max_size_gb = 80, workers = 2)
on.exit(restore_future_runtime(future_state), add = TRUE)

process_one_week <- function(wk, seu_week) {
  log_msg("--- Processing week: ", wk, " (cells=", ncol(seu_week), ") ---")
  if (ncol(seu_week) < 50) {
    stop("Too few cells in week ", wk, " for stable CellChat run.")
  }

  ext <- extract_cellchat_input_matrix(seu_week)
  seu_week <- ext$seu
  data_input <- ext$data_input

  coords <- as.matrix(seu_week@meta.data[, c("x_um", "y_um"), drop = FALSE])
  colnames(coords) <- c("x_cent", "y_cent")
  rownames(coords) <- colnames(seu_week)

  meta <- seu_week@meta.data %>% dplyr::mutate(group = as.character(celltype_plot))
  rownames(meta) <- colnames(seu_week)

  create_args <- list(
    object = data_input,
    meta = meta,
    group.by = "group",
    datatype = "spatial",
    coordinates = coords
  )
  if ("spatial.factors" %in% create_formals) create_args$spatial.factors <- list(ratio = 1, tol = 0)
  if ("scale.factors" %in% create_formals) create_args$scale.factors <- list(spot = 1, spot.diameter = 1)

  cellchat <- do.call(create_fn, create_args)
  cellchat@DB <- subset_db_fn(db_human, search = "Secreted Signaling")
  cellchat@DB$interaction <- cellchat@DB$interaction[cellchat@DB$interaction$pathway_name %in% thesis_pathways, ]
  if (nrow(cellchat@DB$interaction) == 0) stop("Filtered DB is empty for week ", wk)

  # Force diffusion kernel path for continuous coordinates
  cellchat@DB$interaction$annotation <- "Secreted Signaling"

  cellchat <- subset_fn(cellchat)
  cellchat <- over_gene_fn(cellchat)
  cellchat <- over_inter_fn(cellchat)

  min_dist <- compute_min_nonzero_distance(coords)
  dynamic_scale_distance <- min(0.99, 1 / min_dist)
  log_msg("Week ", wk, " dynamic scale.distance (capped) = ", signif(dynamic_scale_distance, 6))

  comm_args <- list(
    object = cellchat,
    distance.use = TRUE,
    interaction.range = interaction_range_um,
    contact.dependent = FALSE
  )
  if ("scale.distance" %in% comm_formals) comm_args$scale.distance <- dynamic_scale_distance
  if ("type" %in% comm_formals) comm_args$type <- "truncatedMean"
  if ("trim" %in% comm_formals) comm_args$trim <- 0

  cellchat <- do.call(commprob_fn, comm_args)
  cellchat <- filter_fn(cellchat, min.cells = 8)
  cellchat <- pathway_fn(cellchat)

  # Keep aggregate optional to reduce crash risk in constrained environments
  cellchat <- tryCatch(
    aggregate_fn(cellchat),
    error = function(e) {
      log_msg("Week ", wk, " WARN during aggregateNet: ", conditionMessage(e), .level = "WARN")
      cellchat
    }
  )

  out_week <- file.path(DIR_OBJECTS, paste0("03_spatial_cellchat_", wk, ".rds"))
  saveRDS(cellchat, out_week)
  log_msg("Saved weekly CellChat object: ", out_week)
  gc(verbose = FALSE)

  list(cellchat = cellchat, output = out_week, min_dist = min_dist)
}

week_names <- names(week_list)
week_results <- setNames(vector("list", length(week_names)), week_names)
cellchat_list <- setNames(vector("list", length(week_names)), week_names)
for (wk in week_names) {
  res <- process_one_week(wk, week_list[[wk]])
  week_results[[wk]] <- list(output = res$output, min_dist = res$min_dist)
  cellchat_list[[wk]] <- res$cellchat
  rm(res)
  gc(verbose = FALSE)
}

cellchat_merged <- merge_fn(cellchat_list, add.names = names(cellchat_list))

merged_out <- file.path(DIR_OBJECTS, "03_spatial_cellchat_merged.rds")
saveRDS(cellchat_merged, merged_out)
log_msg("Saved merged CellChat object: ", merged_out)

weekly_outputs <- unname(vapply(week_results, `[[`, character(1), "output"))
manifest_path <- file.path(DIR_REPORTS, "03_spatial_cellchat_manifest.json")
record_artifact_manifest(
  manifest_path = manifest_path,
  source_data = input_obj,
  output_objects = c(weekly_outputs, merged_out),
  notes = c(
    "Architecture: week-wise split then mergeCellChat (thesis pathways only)",
    paste0("weeks_processed=", paste(names(week_results), collapse = ",")),
    paste0("pathways=", paste(thesis_pathways, collapse = ",")),
    paste0("interaction.range_um=", interaction_range_um),
    "distance.use=TRUE, contact.dependent=FALSE, min.cells=8",
    "workers=2, future.globals.maxSize=80 GiB"
  )
)

log_msg("Spatial CellChat list-and-merge run complete.")
