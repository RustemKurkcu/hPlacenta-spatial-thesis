#!/usr/bin/env Rscript

# =============================================================================
# Script: 03b_spatial_cellchat_allcells.R
# Purpose: Full-pathway Spatial CellChat list-and-merge architecture by week.
# Notes:   Runs ALL pathways (no thesis pathway filter), saves week + merged outputs
#          as both .rds and .qs, and uses sequential future plan for lower RAM.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(jsonlite)
})

source("R/celltype_dictionary.R")

PIPELINE_NAME <- "03b_spatial_cellchat_allcells"
PIPELINE_VERSION <- "2.0.0"
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

as_gib <- function(bytes) {
  as.numeric(bytes) / 1024^3
}

estimate_prob_cell_tensor_gib <- function(n_cells, n_lr, bytes_per_value = 8) {
  as_gib(as.numeric(n_cells) * as.numeric(n_cells) * as.numeric(n_lr) * as.numeric(bytes_per_value))
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

save_dual_object <- function(object, path_rds) {
  saveRDS(object, path_rds)
  out <- c(path_rds)
  if (requireNamespace("qs", quietly = TRUE)) {
    path_qs <- sub("\\.rds$", ".qs", path_rds)
    qs::qsave(object, path_qs, preset = "high")
    out <- c(out, path_qs)
  } else {
    log_msg("Package 'qs' not installed; skipping .qs export for ", basename(path_rds), .level = "WARN")
  }
  out
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
  list(seu = seu, data_input = data_input)
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
    script_path = "scripts/01_active_pipeline/03b_spatial_cellchat_allcells.R",
    notes = notes
  )
  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
}

interaction_range_um <- 100
pathway_tensor_soft_limit_gib <- as.numeric(Sys.getenv("CELLCHAT_PATHWAY_TENSOR_SOFT_LIMIT_GIB", "90"))

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

future_state <- configure_future_runtime(max_size_gb = 80, workers = 1)
on.exit(restore_future_runtime(future_state), add = TRUE)

process_one_week <- function(wk, seu_week) {
  log_msg("--- Processing week (FULL pathways): ", wk, " (cells=", ncol(seu_week), ") ---")
  if (ncol(seu_week) < 50) stop("Too few cells in week ", wk)

  ext <- extract_cellchat_input_matrix(seu_week)
  seu_week <- ext$seu
  data_input <- ext$data_input

  coords <- as.matrix(seu_week@meta.data[, c("x_um", "y_um"), drop = FALSE])
  colnames(coords) <- c("x_cent", "y_cent")
  rownames(coords) <- colnames(seu_week)

  meta <- seu_week@meta.data %>% dplyr::mutate(group = as.character(celltype_plot))
  rownames(meta) <- colnames(seu_week)

  create_args <- list(object = data_input, meta = meta, group.by = "group", datatype = "spatial", coordinates = coords)
  if ("spatial.factors" %in% create_formals) create_args$spatial.factors <- list(ratio = 1, tol = 0)
  if ("scale.factors" %in% create_formals) create_args$scale.factors <- list(spot = 1, spot.diameter = 1)

  cellchat <- do.call(create_fn, create_args)
  cellchat@DB <- subset_db_fn(db_human, search = "Secreted Signaling")
  if (nrow(cellchat@DB$interaction) == 0) stop("Filtered DB is empty for week ", wk)
  cellchat@DB$interaction$annotation <- "Secreted Signaling"

  cellchat <- subset_fn(cellchat)
  cellchat <- over_gene_fn(cellchat)
  cellchat <- over_inter_fn(cellchat)

  min_dist <- compute_min_nonzero_distance(coords)
  dynamic_scale_distance <- min(0.99, 1 / min_dist)

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
  gc(verbose = FALSE)
  cellchat <- filter_fn(cellchat, min.cells = 8)
  n_lr <- tryCatch(nrow(cellchat@LR$LRsig), error = function(e) NA_integer_)
  est_tensor_gib <- if (!is.na(n_lr)) estimate_prob_cell_tensor_gib(ncol(seu_week), n_lr) else NA_real_
  if (!is.na(est_tensor_gib) && est_tensor_gib > pathway_tensor_soft_limit_gib) {
    log_msg(
      "Week ", wk, " WARN estimated Prob.cell tensor size is ~",
      format(round(est_tensor_gib, 1), nsmall = 1), " GiB (cells=", ncol(seu_week),
      ", LR=", n_lr, "). This is above soft limit ",
      pathway_tensor_soft_limit_gib, " GiB and may OOM during computeCommunProbPathway.",
      .level = "WARN"
    )
  }

  pathway_ok <- TRUE
  pathway_error <- NULL
  cellchat <- tryCatch(
    pathway_fn(cellchat),
    error = function(e) {
      pathway_ok <<- FALSE
      pathway_error <<- conditionMessage(e)
      log_msg("Week ", wk, " WARN during computeCommunProbPathway: ", pathway_error, .level = "WARN")
      cellchat
    }
  )

  if (pathway_ok) {
    cellchat <- tryCatch(
      aggregate_fn(cellchat),
      error = function(e) {
        log_msg("Week ", wk, " WARN during aggregateNet: ", conditionMessage(e), .level = "WARN")
        cellchat
      }
    )
  } else {
    log_msg(
      "Week ", wk, " pathway-level outputs were skipped due to error; saving object with inferred LR-level communication intact.",
      .level = "WARN"
    )
  }

  out_week <- file.path(DIR_OBJECTS, paste0("03b_spatial_cellchat_full_", wk, ".rds"))
  out_week_paths <- save_dual_object(cellchat, out_week)
  log_msg("Saved full-pathway weekly object(s): ", paste(out_week_paths, collapse = ", "))
  gc(verbose = FALSE)

  list(
    cellchat = cellchat,
    output = out_week_paths,
    min_dist = min_dist,
    pathway_ok = pathway_ok,
    pathway_error = pathway_error,
    est_tensor_gib = est_tensor_gib
  )
}

week_names <- names(week_list)
week_results <- setNames(vector("list", length(week_names)), week_names)
cellchat_list <- setNames(vector("list", length(week_names)), week_names)
merge_candidates <- character(0)
for (wk in week_names) {
  res <- process_one_week(wk, week_list[[wk]])
  week_results[[wk]] <- list(
    output = res$output,
    min_dist = res$min_dist,
    pathway_ok = res$pathway_ok,
    pathway_error = res$pathway_error,
    est_tensor_gib = res$est_tensor_gib
  )
  cellchat_list[[wk]] <- res$cellchat
  if (isTRUE(res$pathway_ok)) merge_candidates <- c(merge_candidates, wk)
  rm(res)
  gc(verbose = FALSE)
}

merged_out_paths <- character(0)
if (length(merge_candidates) == 0) {
  log_msg("No week completed computeCommunProbPathway successfully; skipping merge.", .level = "WARN")
} else {
  merge_list <- cellchat_list[merge_candidates]
  cellchat_merged <- merge_fn(merge_list, add.names = names(merge_list))
  merged_out <- file.path(DIR_OBJECTS, "03b_spatial_cellchat_full_merged.rds")
  merged_out_paths <- save_dual_object(cellchat_merged, merged_out)
  log_msg("Saved full-pathway merged object(s): ", paste(merged_out_paths, collapse = ", "))
}

weekly_outputs <- unlist(lapply(week_results, `[[`, "output"), use.names = FALSE)
manifest_path <- file.path(DIR_REPORTS, "03b_spatial_cellchat_allcells_manifest.json")
record_artifact_manifest(
  manifest_path = manifest_path,
  source_data = input_obj,
  output_objects = c(weekly_outputs, merged_out_paths),
  notes = c(
    "Architecture: week-wise split then mergeCellChat (FULL pathways)",
    paste0("weeks_processed=", paste(names(week_results), collapse = ",")),
    paste0("weeks_pathway_ok=", paste(merge_candidates, collapse = ",")),
    paste0("interaction.range_um=", interaction_range_um),
    "distance.use=TRUE, contact.dependent=FALSE, min.cells=8",
    paste0("pathway_tensor_soft_limit_gib=", pathway_tensor_soft_limit_gib),
    "workers=1 (sequential), future.globals.maxSize=80 GiB"
  )
)

log_msg("03b full-pathway week-wise + merged run complete.")
