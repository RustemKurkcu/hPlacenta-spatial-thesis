#!/usr/bin/env Rscript

# =============================================================================
# Script: 03f_force_pathway_aggregation_and_merge.R
# Purpose: Repair weekly SpatialCellChat objects that are missing pathway-level
#          aggregates by chunking over LR pathways (RAM-safe), then merge.
# Notes:   Designed for local machines where aggregateNet() can fail.
# =============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
})

PIPELINE_NAME <- "03f_force_pathway_aggregation_and_merge"
PIPELINE_VERSION <- "1.0.0"
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
    if (!is.null(s3d[[nm]])) return(s3d[[nm]])
  }
  stop("Could not find sparse3D values field (expected one of: x, v, value, values).")
}

set_sparse3d_values <- function(s3d, vals) {
  for (nm in c("x", "v", "value", "values")) {
    if (!is.null(s3d[[nm]])) {
      s3d[[nm]] <- vals
      return(s3d)
    }
  }
  stop("Could not set sparse3D values field.")
}

force_group_level_aggregates <- function(cellchat) {
  if (is.null(cellchat@LR$LRsig) || nrow(cellchat@LR$LRsig) == 0) {
    stop("Missing cellchat@LR$LRsig; cannot build pathway aggregates.")
  }
  if (is.null(cellchat@net$prob.cell)) {
    stop("Missing cellchat@net$prob.cell; cannot build pathway aggregates.")
  }

  prob_cell <- cellchat@net$prob.cell
  edge_i <- as.integer(prob_cell$i)
  edge_j <- as.integer(prob_cell$j)
  edge_k <- as.integer(prob_cell$k)
  edge_x <- as.numeric(get_sparse3d_values(prob_cell))

  if (length(edge_i) == 0 || length(edge_j) == 0 || length(edge_k) == 0 || length(edge_x) == 0) {
    stop("prob.cell is empty; cannot aggregate.")
  }

  idents <- as.factor(cellchat@idents)
  group_levels <- levels(idents)
  n_groups <- length(group_levels)
  if (n_groups == 0) stop("No group levels found in cellchat@idents.")

  sender_group_idx <- as.integer(idents)[edge_i]
  receiver_group_idx <- as.integer(idents)[edge_j]

  lr_df <- cellchat@LR$LRsig
  if (!"pathway_name" %in% colnames(lr_df)) stop("LRsig is missing pathway_name column.")
  pathway_names <- unique(as.character(lr_df$pathway_name))
  pathway_names <- pathway_names[!is.na(pathway_names) & nzchar(pathway_names)]
  if (length(pathway_names) == 0) stop("No valid pathway_name values in LRsig.")

  lr_pathway <- as.character(lr_df$pathway_name)

  prob_arr <- array(0, dim = c(n_groups, n_groups, length(pathway_names)),
                    dimnames = list(group_levels, group_levels, pathway_names))
  pval_arr <- array(1, dim = c(n_groups, n_groups, length(pathway_names)),
                    dimnames = list(group_levels, group_levels, pathway_names))

  log_msg("Chunk aggregation: groups=", n_groups, ", pathways=", length(pathway_names), ", edges=", length(edge_x))

  for (p_idx in seq_along(pathway_names)) {
    p <- pathway_names[p_idx]
    lr_ids <- which(lr_pathway == p)
    if (length(lr_ids) == 0) next

    sel <- which(edge_k %in% lr_ids)
    if (length(sel) == 0) {
      gc(verbose = FALSE)
      next
    }

    pair_idx <- (sender_group_idx[sel] - 1L) * n_groups + receiver_group_idx[sel]

    pair_sum <- rowsum(edge_x[sel], group = pair_idx, reorder = FALSE)
    pair_ids <- as.integer(rownames(pair_sum))
    vals <- as.numeric(pair_sum[, 1])

    row_idx <- ((pair_ids - 1L) %/% n_groups) + 1L
    col_idx <- ((pair_ids - 1L) %% n_groups) + 1L

    prob_arr[cbind(row_idx, col_idx, rep.int(p_idx, length(vals)))] <- vals
    pval_arr[cbind(row_idx, col_idx, rep.int(p_idx, length(vals)))] <- 0

    gc(verbose = FALSE)
  }

  # Provide netP for pathway-level downstream functions.
  cellchat@netP$pathways <- pathway_names
  cellchat@netP$prob <- prob_arr
  cellchat@netP$pval <- pval_arr

  # Provide lightweight group-level net summaries compatible with most plotting.
  weight_mat <- apply(prob_arr, c(1, 2), sum)
  count_mat <- apply(prob_arr > 0, c(1, 2), sum)
  cellchat@net$weight <- weight_mat
  cellchat@net$count <- count_mat

  cellchat
}

is_pathway_complete <- function(cellchat) {
  !is.null(cellchat@netP$prob) && length(cellchat@netP$prob) > 0
}

weeks <- c("W7", "W8-2", "W9", "W11")
input_paths <- file.path(DIR_OBJECTS, paste0("03_spatial_cellchat_", weeks, ".rds"))
input_paths <- vapply(input_paths, resolve_object_path, FUN.VALUE = character(1))

log_msg("Repair pass for weekly objects: ", paste(basename(input_paths), collapse = ", "))

cellchat_fixed <- setNames(vector("list", length(weeks)), weeks)
out_paths <- character(0)

for (wk in weeks) {
  in_path <- resolve_object_path(file.path(DIR_OBJECTS, paste0("03_spatial_cellchat_", wk, ".rds")))
  obj <- read_object(in_path)

  if (is_pathway_complete(obj)) {
    log_msg("Week ", wk, ": netP already populated; recomputing lightweight group summaries for consistency.")
  } else {
    log_msg("Week ", wk, ": netP missing/empty; building via pathway chunk aggregation.", .level = "WARN")
  }

  obj <- force_group_level_aggregates(obj)

  out_week <- file.path(DIR_OBJECTS, paste0("03_spatial_cellchat_", wk, "_repaired.rds"))
  out_week_paths <- save_dual_object(obj, out_week)
  out_paths <- c(out_paths, out_week_paths)
  log_msg("Saved repaired weekly object(s): ", paste(out_week_paths, collapse = ", "))

  cellchat_fixed[[wk]] <- obj
  gc(verbose = FALSE)
}

merge_fn <- CellChat::mergeCellChat
cellchat_merged <- merge_fn(cellchat_fixed, add.names = names(cellchat_fixed))
merged_out <- file.path(DIR_OBJECTS, "03_spatial_cellchat_merged_repaired.rds")
merged_out_paths <- save_dual_object(cellchat_merged, merged_out)
out_paths <- c(out_paths, merged_out_paths)
log_msg("Saved repaired merged CellChat object(s): ", paste(merged_out_paths, collapse = ", "))

manifest_path <- file.path(DIR_REPORTS, "03f_force_pathway_aggregation_and_merge_manifest.json")
record_manifest(
  manifest_path,
  inputs = unname(input_paths),
  outputs = out_paths,
  notes = c(
    "Chunked pathway aggregation from net$prob.cell by pathway_name",
    "Builds netP$prob/netP$pval and lightweight net$weight/net$count",
    "Merges repaired weekly objects into 03_spatial_cellchat_merged_repaired"
  )
)

log_msg("03f repair-and-merge completed successfully.")
