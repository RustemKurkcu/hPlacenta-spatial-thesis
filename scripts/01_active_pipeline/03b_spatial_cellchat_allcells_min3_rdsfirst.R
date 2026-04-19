#!/usr/bin/env Rscript

# =============================================================================
# Script: 03b_spatial_cellchat_allcells_min3_rdsfirst.R
# Purpose: Thesis-safe Spatial CellChat / CellChat run across weeks with
#          RDS-first input resolution, min.cells=3, and organized outputs.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(jsonlite)
})

source("R/celltype_dictionary.R")

PIPELINE_NAME <- "03b_spatial_cellchat_allcells_min3_rdsfirst"
PIPELINE_VERSION <- "2.1.0"
RUN_TIMESTAMP <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
RUN_SEED <- 42L
set.seed(RUN_SEED)

OUT_ROOT <- Sys.getenv("CELLCHAT_OUT_ROOT", file.path("output", "cellchat", "allcells_min3"))
INPUT_OBJECT_DIR <- Sys.getenv("CELLCHAT_INPUT_DIR", file.path("output", "objects"))
INPUT_BASENAME <- Sys.getenv("CELLCHAT_INPUT_BASENAME", "02_scored_misi_ido1")
MIN_CELLS_FILTER <- as.integer(Sys.getenv("CELLCHAT_MIN_CELLS", "3"))
INTERACTION_RANGE_UM <- as.numeric(Sys.getenv("CELLCHAT_INTERACTION_RANGE_UM", "100"))
PATHWAY_TENSOR_SOFT_LIMIT_GIB <- as.numeric(Sys.getenv("CELLCHAT_PATHWAY_TENSOR_SOFT_LIMIT_GIB", "90"))
DB_SEARCH <- Sys.getenv("CELLCHAT_DB_SEARCH", "Secreted Signaling")

DIR_LOGS <- file.path(OUT_ROOT, "logs")
DIR_OBJECTS <- file.path(OUT_ROOT, "objects")
DIR_REPORTS <- file.path(OUT_ROOT, "reports")
DIR_TABLES <- file.path(OUT_ROOT, "tables")
DIR_TABLES_WEEKLY <- file.path(DIR_TABLES, "weekly")
DIR_TABLES_QC <- file.path(DIR_TABLES, "qc")

for (d in c(DIR_LOGS, DIR_OBJECTS, DIR_REPORTS, DIR_TABLES, DIR_TABLES_WEEKLY, DIR_TABLES_QC)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
LOG_FILE <- file.path(DIR_LOGS, paste0(PIPELINE_NAME, ".log"))

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

as_gib <- function(bytes) as.numeric(bytes) / 1024^3
estimate_prob_cell_tensor_gib <- function(n_cells, n_lr, bytes_per_value = 8) {
  as_gib(as.numeric(n_cells) * as.numeric(n_cells) * as.numeric(n_lr) * as.numeric(bytes_per_value))
}

resolve_object_path <- function(path_rds) {
  path_qs <- sub("\\.rds$", ".qs", path_rds)
  if (file.exists(path_rds)) return(path_rds)
  if (file.exists(path_qs)) return(path_qs)
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
    ok <- tryCatch({ qs::qsave(object, path_qs, preset = "high"); TRUE }, error = function(e) {
      log_msg("Package 'qs' available but qsave failed for ", basename(path_rds), ": ", conditionMessage(e), .level = "WARN")
      FALSE
    })
    if (isTRUE(ok)) out <- c(out, path_qs)
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
  min(as.matrix(stats::dist(xy))[upper.tri(as.matrix(stats::dist(xy)))], na.rm = TRUE)
}

record_artifact_manifest <- function(manifest_path, source_data, output_objects, notes = NULL) {
  manifest <- list(
    pipeline = PIPELINE_NAME,
    version = PIPELINE_VERSION,
    run_timestamp = RUN_TIMESTAMP,
    seed = RUN_SEED,
    source_data = source_data,
    output_objects = output_objects,
    script_path = "scripts/01_active_pipeline/03b_spatial_cellchat_allcells_min3_rdsfirst.R",
    notes = notes
  )
  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
}

save_celltype_counts <- function(seu_week, wk) {
  counts <- seu_week@meta.data |>
    dplyr::count(celltype_plot, name = "n_cells", sort = TRUE) |>
    dplyr::mutate(week = wk, .before = 1)
  out <- file.path(DIR_TABLES_WEEKLY, paste0("celltype_counts_", wk, ".tsv"))
  write.table(counts, out, sep = "\t", row.names = FALSE, quote = FALSE)
  out
}

input_obj <- resolve_object_path(file.path(INPUT_OBJECT_DIR, paste0(INPUT_BASENAME, ".rds")))
log_msg("Loading scored object: ", input_obj)
seu <- read_object(input_obj)

seu$week <- factor(as.character(seu$week), levels = c("W7", "W8-2", "W9", "W11"))
celltype_col <- if ("predicted.celltype" %in% colnames(seu@meta.data)) "predicted.celltype" else pick_celltype_source_column(seu@meta.data)
if (is.na(celltype_col)) celltype_col <- "seurat_clusters"
seu$celltype_plot <- rename_with_true_names(as.character(seu@meta.data[[celltype_col]]))
Idents(seu) <- "celltype_plot"

week_list <- SplitObject(seu, split.by = "week")
week_list <- week_list[names(week_list) %in% c("W7", "W8-2", "W9", "W11")]

using_spatial_cellchat <- requireNamespace("SpatialCellChat", quietly = TRUE)
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
  merge_fn <- CellChat::mergeCellChat
  db_human <- SpatialCellChat::CellChatDB.human
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
  merge_fn <- CellChat::mergeCellChat
  db_human <- CellChat::CellChatDB.human
}

create_formals <- names(formals(create_fn))
comm_formals <- names(formals(commprob_fn))
future_state <- configure_future_runtime(max_size_gb = 80, workers = 1)
on.exit(restore_future_runtime(future_state), add = TRUE)

process_one_week <- function(wk, seu_week) {
  celltype_count_file <- save_celltype_counts(seu_week, wk)
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
  cellchat@DB <- subset_db_fn(db_human, search = DB_SEARCH)
  cellchat <- subset_fn(cellchat)
  cellchat <- over_gene_fn(cellchat)
  cellchat <- over_inter_fn(cellchat)

  min_dist <- compute_min_nonzero_distance(coords)
  comm_args <- list(object = cellchat, distance.use = TRUE, interaction.range = INTERACTION_RANGE_UM, contact.dependent = FALSE)
  if ("scale.distance" %in% comm_formals) comm_args$scale.distance <- min(0.99, 1 / min_dist)
  if ("type" %in% comm_formals) comm_args$type <- "truncatedMean"
  if ("trim" %in% comm_formals) comm_args$trim <- 0

  cellchat <- do.call(commprob_fn, comm_args)
  cellchat <- filter_fn(cellchat, min.cells = MIN_CELLS_FILTER)

  pathway_ok <- TRUE
  pathway_error <- NULL
  cellchat <- tryCatch(pathway_fn(cellchat), error = function(e) { pathway_ok <<- FALSE; pathway_error <<- conditionMessage(e); cellchat })
  if (pathway_ok) {
    cellchat <- tryCatch(aggregate_fn(cellchat), error = function(e) cellchat)
  }

  out_week <- file.path(DIR_OBJECTS, paste0("spatial_cellchat_full_", wk, ".rds"))
  out_week_paths <- save_dual_object(cellchat, out_week)
  list(cellchat = cellchat, output = out_week_paths, min_dist = min_dist, pathway_ok = pathway_ok, pathway_error = pathway_error, celltype_count_file = celltype_count_file)
}

week_results <- list()
cellchat_list <- list()
merge_candidates <- character(0)
for (wk in names(week_list)) {
  res <- process_one_week(wk, week_list[[wk]])
  week_results[[wk]] <- res
  cellchat_list[[wk]] <- res$cellchat
  if (isTRUE(res$pathway_ok)) merge_candidates <- c(merge_candidates, wk)
}

week_summary_df <- do.call(rbind, lapply(names(week_results), function(wk) {
  x <- week_results[[wk]]
  data.frame(week = wk, min_dist_um = x$min_dist, pathway_ok = x$pathway_ok, pathway_error = ifelse(is.null(x$pathway_error), "", x$pathway_error), output_rds = paste(x$output, collapse = ";"), celltype_count_file = x$celltype_count_file, stringsAsFactors = FALSE)
}))
week_summary_file <- file.path(DIR_TABLES_QC, "week_summary.tsv")
write.table(week_summary_df, week_summary_file, sep = "\t", quote = FALSE, row.names = FALSE)

merged_out_paths <- character(0)
if (length(merge_candidates) >= 2) {
  cellchat_merged <- merge_fn(cellchat_list[merge_candidates], add.names = merge_candidates)
  merged_out_paths <- save_dual_object(cellchat_merged, file.path(DIR_OBJECTS, "spatial_cellchat_full_merged.rds"))
}

weekly_outputs <- unlist(lapply(week_results, function(x) x$output), use.names = FALSE)
record_artifact_manifest(file.path(DIR_REPORTS, "artifact_manifest.json"), input_obj, c(weekly_outputs, merged_out_paths, week_summary_file),
  notes = c(paste0("min.cells=", MIN_CELLS_FILTER), paste0("db_search=", DB_SEARCH), paste0("pathway_tensor_soft_limit_gib=", PATHWAY_TENSOR_SOFT_LIMIT_GIB)))

log_msg("Run complete.")
