#!/usr/bin/env Rscript

# =============================================================================
# Script: 03f_force_pathway_aggregation_and_merge.R
# Purpose: RAM-safe, cell-level pathway aggregation for full-pathway weekly
#          SpatialCellChat objects, then 4-week merge.
# Notes:   Preserves n_cells x n_cells spatial geometry (does NOT collapse to
#          group-level arrays).
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(jsonlite)
})

if (requireNamespace("SpatialCellChat", quietly = TRUE)) {
  merge_fn <- CellChat::mergeCellChat
} else if (requireNamespace("CellChat", quietly = TRUE)) {
  merge_fn <- CellChat::mergeCellChat
} else {
  stop("Neither SpatialCellChat nor CellChat is installed.")
}

PIPELINE_NAME <- "03f_force_pathway_aggregation_and_merge"
PIPELINE_VERSION <- "2.0.1"
RUN_TIMESTAMP <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

OUT_ROOT <- "output"
DIR_LOGS <- file.path(OUT_ROOT, "logs")
DIR_OBJECTS <- file.path(OUT_ROOT, "objects")
DIR_REPORTS <- file.path(OUT_ROOT, "reports")
for (d in c(DIR_LOGS, DIR_OBJECTS, DIR_REPORTS)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
LOG_FILE <- file.path(DIR_LOGS, "03f_force_pathway_aggregation_and_merge.log")

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
  }
  out
}

record_manifest <- function(manifest_path, inputs, outputs, notes = NULL) {
  manifest <- list(
    pipeline = PIPELINE_NAME,
    version = PIPELINE_VERSION,
    run_timestamp = RUN_TIMESTAMP,
    input_objects = inputs,
    output_objects = outputs,
    script_path = "scripts/01_active_pipeline/03f_force_pathway_aggregation_and_merge.R",
    notes = notes
  )
  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
}

get_sparse3d_values <- function(s3d) {
  for (nm in c("x", "v", "value", "values")) {
    if (!is.null(s3d[[nm]])) return(list(name = nm, values = s3d[[nm]]))
  }
  stop("Could not find sparse3D values field in sparse3Darray.")
}

get_sparse3d_dim <- function(s3d) {
  if (!is.null(s3d$dim)) return(as.integer(s3d$dim))
  d <- attr(s3d, "dim")
  if (!is.null(d)) return(as.integer(d))
  stop("Could not determine sparse3Darray dimensions.")
}

get_sparse3d_dimnames <- function(s3d) {
  if (!is.null(s3d$dimnames)) return(s3d$dimnames)
  dn <- attr(s3d, "dimnames")
  if (!is.null(dn)) return(dn)
  NULL
}

build_pathway_tensor_cell_level <- function(cellchat, week_label = "unknown") {
  if (is.null(cellchat@net$prob.cell)) {
    stop("Week ", week_label, ": missing net$prob.cell.")
  }
  if (is.null(cellchat@LR$LRsig) || nrow(cellchat@LR$LRsig) == 0) {
    stop("Week ", week_label, ": missing LRsig.")
  }

  prob_cell <- cellchat@net$prob.cell
  dim_pc <- get_sparse3d_dim(prob_cell)
  dn_pc <- get_sparse3d_dimnames(prob_cell)
  val_info <- get_sparse3d_values(prob_cell)

  edge_i <- as.integer(prob_cell$i)
  edge_j <- as.integer(prob_cell$j)
  edge_k <- as.integer(prob_cell$k)
  edge_v <- as.numeric(val_info$values)

  if (length(edge_i) == 0) stop("Week ", week_label, ": prob.cell has zero edges.")

  lr_df <- cellchat@LR$LRsig
  pathways <- unique(as.character(lr_df$pathway_name))
  pathways <- pathways[!is.na(pathways) & nzchar(pathways)]
  if (length(pathways) == 0) stop("Week ", week_label, ": no valid pathway_name values in LRsig.")

  lr_path <- as.character(lr_df$pathway_name)
  pathway_chunks <- vector("list", length(pathways))

  log_msg("Week ", week_label, ": chunking ", length(pathways), " pathways across ", length(edge_v), " sparse edges.")

  for (p_idx in seq_along(pathways)) {
    p <- pathways[p_idx]
    lr_ids <- which(lr_path == p)
    sel <- which(edge_k %in% lr_ids)

    if (length(sel) == 0) {
      pathway_chunks[[p_idx]] <- NULL
      if (p_idx %% 10 == 0) gc(verbose = FALSE)
      next
    }

    mat <- Matrix::sparseMatrix(
      i = edge_i[sel],
      j = edge_j[sel],
      x = edge_v[sel],
      dims = c(dim_pc[1], dim_pc[2])
    )
    # Matrix::summary() can fail on some local Matrix builds due to class
    # registration issues (e.g., gTMatrix not defined). Prefer direct slots.
    if (inherits(mat, "dgCMatrix")) {
      dp <- diff(mat@p)
      j_out <- rep.int(seq_along(dp), dp)
      i_out <- mat@i + 1L  # convert 0-based CSC row index to 1-based
      v_out <- mat@x
      pathway_chunks[[p_idx]] <- list(
        i = as.integer(i_out),
        j = as.integer(j_out),
        k = rep.int(as.integer(p_idx), length(v_out)),
        v = as.numeric(v_out)
      )
      rm(dp, j_out, i_out, v_out)
    } else {
      sm <- Matrix::summary(mat)
      pathway_chunks[[p_idx]] <- list(
        i = as.integer(sm$i),
        j = as.integer(sm$j),
        k = rep.int(as.integer(p_idx), nrow(sm)),
        v = as.numeric(sm$x)
      )
      rm(sm)
    }

    rm(mat)
    if (p_idx %% 10 == 0) gc(verbose = FALSE)
  }

  keep <- vapply(pathway_chunks, Negate(is.null), logical(1))
  pathway_chunks <- pathway_chunks[keep]

  if (length(pathway_chunks) == 0) {
    stop("Week ", week_label, ": no non-zero pathway chunks could be constructed.")
  }

  i_all <- unlist(lapply(pathway_chunks, `[[`, "i"), use.names = FALSE)
  j_all <- unlist(lapply(pathway_chunks, `[[`, "j"), use.names = FALSE)
  k_all <- unlist(lapply(pathway_chunks, `[[`, "k"), use.names = FALSE)
  v_all <- unlist(lapply(pathway_chunks, `[[`, "v"), use.names = FALSE)

  dimnames_path <- NULL
  if (!is.null(dn_pc) && length(dn_pc) >= 2) {
    dimnames_path <- list(dn_pc[[1]], dn_pc[[2]], pathways)
  } else {
    dimnames_path <- list(NULL, NULL, pathways)
  }

  pathway_tensor <- list(
    i = as.integer(i_all),
    j = as.integer(j_all),
    k = as.integer(k_all),
    dim = c(as.integer(dim_pc[1]), as.integer(dim_pc[2]), as.integer(length(pathways))),
    dimnames = dimnames_path
  )
  pathway_tensor[[val_info$name]] <- as.numeric(v_all)
  attr(pathway_tensor, "class") <- "sparse3Darray"

  list(pathway_tensor = pathway_tensor, pathways = pathways)
}

has_complete_pathway_tensor <- function(cellchat) {
  if (is.null(cellchat@net$prob)) return(FALSE)
  d <- tryCatch(get_sparse3d_dim(cellchat@net$prob), error = function(e) NULL)
  if (is.null(d) || length(d) < 3) return(FALSE)
  d[3] > 0
}

weeks <- c("W7", "W8-2", "W9", "W11")
in_week_paths <- file.path(DIR_OBJECTS, paste0("03b_spatial_cellchat_full_", weeks, ".rds"))
in_week_paths <- vapply(in_week_paths, resolve_object_path, FUN.VALUE = character(1))

log_msg("Starting cell-level pathway repair for: ", paste(basename(in_week_paths), collapse = ", "))

cellchat_list <- setNames(vector("list", length(weeks)), weeks)
out_paths <- character(0)

for (wk in weeks) {
  in_path <- resolve_object_path(file.path(DIR_OBJECTS, paste0("03b_spatial_cellchat_full_", wk, ".rds")))
  log_msg("Loading ", in_path)
  obj <- read_object(in_path)

  if (has_complete_pathway_tensor(obj)) {
    log_msg("Week ", wk, ": net$prob already present; preserving existing tensor.")
  } else {
    log_msg("Week ", wk, ": net$prob missing; constructing cell-level pathway tensor from net$prob.cell.", .level = "WARN")
    rebuilt <- build_pathway_tensor_cell_level(obj, week_label = wk)
    obj@net$prob <- rebuilt$pathway_tensor
    if (is.null(obj@netP)) obj@netP <- list()
    obj@netP$pathways <- rebuilt$pathways
    obj@netP$prob <- rebuilt$pathway_tensor
    if (is.null(obj@netP$pval)) {
      obj@netP$pval <- rebuilt$pathway_tensor
      vinfo <- get_sparse3d_values(obj@netP$pval)
      obj@netP$pval[[vinfo$name]] <- rep(0, length(vinfo$values))
    }
  }

  out_week <- file.path(DIR_OBJECTS, paste0("03b_spatial_cellchat_full_", wk, "_completed.rds"))
  out_week_paths <- save_dual_object(obj, out_week)
  out_paths <- c(out_paths, out_week_paths)
  log_msg("Saved completed weekly object(s): ", paste(out_week_paths, collapse = ", "))

  cellchat_list[[wk]] <- obj
  rm(obj)
  gc(verbose = FALSE)
}

log_msg("Merging 4 completed weekly objects...")
cellchat_merged <- merge_fn(cellchat_list, add.names = names(cellchat_list))
if (!is.null(cellchat_merged@meta$datasets)) {
  cellchat_merged@meta$datasets <- factor(cellchat_merged@meta$datasets, levels = weeks)
}

merged_out <- file.path(DIR_OBJECTS, "03b_spatial_cellchat_full_merged_complete.rds")
merged_out_paths <- save_dual_object(cellchat_merged, merged_out)
out_paths <- c(out_paths, merged_out_paths)
log_msg("Saved merged completed object(s): ", paste(merged_out_paths, collapse = ", "))

manifest_path <- file.path(DIR_REPORTS, "03f_force_pathway_aggregation_and_merge_manifest.json")
record_manifest(
  manifest_path,
  inputs = unname(in_week_paths),
  outputs = out_paths,
  notes = c(
    "Cell-level (n_cells x n_cells x n_pathways) pathway tensor reconstruction from net$prob.cell",
    "Preserves spatial coordinate resolution for spatial pathway maps",
    "Uses dgCMatrix slot extraction fallback to bypass Matrix::summary class issues",
    "Writes *_completed weekly objects and merged_complete object"
  )
)

log_msg("03f force-pathway aggregation + merge completed.")
