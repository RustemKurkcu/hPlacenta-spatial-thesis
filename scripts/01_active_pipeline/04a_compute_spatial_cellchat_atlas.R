#!/usr/bin/env Rscript

# =============================================================================
# Script: 04a_compute_spatial_cellchat_atlas.R
# Purpose: Heavy-compute stage for the 5-Phase Spatial CellChat Atlas.
#          - Build master LR and thesis pathway tables (Phase 5)
#          - Compute MISI vulnerability split and re-aggregated communication
#            summaries for target micro-niches (Phase 4 compute)
# Outputs: output/tables/*.csv, output/objects/04a_atlas_compute_cache.rds,
#          output/reports/04_methods_and_provenance.md
# =============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
})

OUT_ROOT <- "output"
DIR_OBJECTS <- file.path(OUT_ROOT, "objects")
DIR_TABLES <- file.path(OUT_ROOT, "tables")
DIR_REPORTS <- file.path(OUT_ROOT, "reports")
for (d in c(DIR_OBJECTS, DIR_TABLES, DIR_REPORTS)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

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

append_methods_log <- function(header, lines) {
  methods_path <- file.path(DIR_REPORTS, "04_methods_and_provenance.md")
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  cat(
    paste0("\n## ", stamp, " — ", header, "\n\n"),
    file = methods_path,
    append = TRUE
  )
  for (ln in lines) cat(paste0("- ", ln, "\n"), file = methods_path, append = TRUE)
}

safe_subset_communication <- function(obj) {
  if (requireNamespace("SpatialCellChat", quietly = TRUE)) {
    fn <- SpatialCellChat::subsetCommunication
  } else if (requireNamespace("CellChat", quietly = TRUE)) {
    fn <- CellChat::subsetCommunication
  } else {
    stop("Neither SpatialCellChat nor CellChat is installed.")
  }
  out <- tryCatch(fn(obj), error = function(e) NULL)
  if (is.null(out) || nrow(out) == 0) {
    out <- tryCatch(fn(obj, slot.name = "netP"), error = function(e) NULL)
  }
  out
}

get_sparse3d_values <- function(s3d) {
  for (nm in c("x", "v", "value", "values")) {
    if (!is.null(s3d[[nm]])) return(s3d[[nm]])
  }
  stop("Could not find sparse3D values field in sparse3Darray.")
}

aggregate_probcell_to_groups <- function(obj, group_labels) {
  prob_cell <- obj@net$prob.cell
  lr_df <- obj@LR$LRsig

  edge_i <- as.integer(prob_cell$i)
  edge_j <- as.integer(prob_cell$j)
  edge_k <- as.integer(prob_cell$k)
  edge_v <- as.numeric(get_sparse3d_values(prob_cell))

  if (length(edge_v) == 0) return(data.frame())

  gl <- as.character(group_labels)
  src <- gl[edge_i]
  tgt <- gl[edge_j]
  lr_name <- if ("interaction_name" %in% colnames(lr_df)) as.character(lr_df$interaction_name)[edge_k] else paste0("LR_", edge_k)
  pathway <- if ("pathway_name" %in% colnames(lr_df)) as.character(lr_df$pathway_name)[edge_k] else NA_character_

  df <- data.frame(
    source = src,
    target = tgt,
    interaction_name = lr_name,
    pathway_name = pathway,
    prob_value = edge_v,
    stringsAsFactors = FALSE
  )

  df |>
    dplyr::group_by(source, target, interaction_name, pathway_name) |>
    dplyr::summarise(prob_value = sum(prob_value, na.rm = TRUE), .groups = "drop")
}

weeks <- c("W7", "W8-2", "W9", "W11")
week_paths <- setNames(
  vapply(file.path(DIR_OBJECTS, paste0("03b_spatial_cellchat_full_", weeks, "_completed.rds")), resolve_object_path, FUN.VALUE = character(1)),
  weeks
)

target_celltype <- Sys.getenv("ATLAS_TARGET_CELLTYPE", unset = "EVT")
low_q <- as.numeric(Sys.getenv("ATLAS_MISI_LOW_Q", unset = "0.33"))
high_q <- as.numeric(Sys.getenv("ATLAS_MISI_HIGH_Q", unset = "0.67"))

master_lr <- list()
thesis_lr <- list()
vulnerability_tables <- list()
vulnerability_meta <- list()

for (wk in weeks) {
  obj <- read_object(week_paths[[wk]])

  comm_df <- safe_subset_communication(obj)
  if (!is.null(comm_df) && nrow(comm_df) > 0) {
    comm_df$week <- wk
    master_lr[[wk]] <- comm_df

    pathway_col <- c("pathway_name", "pathway", "signaling")[c("pathway_name", "pathway", "signaling") %in% colnames(comm_df)][1]
    if (!is.na(pathway_col)) {
      thesis_lr[[wk]] <- comm_df[as.character(comm_df[[pathway_col]]) %in% c("MMP", "IDO1", "TGFb"), , drop = FALSE]
    }
  }

  md <- obj@meta
  misi_col <- c("MISI_Vulnerability", "MISI", "MISI_score")[c("MISI_Vulnerability", "MISI", "MISI_score") %in% colnames(md)][1]
  if (is.na(misi_col)) {
    warning("Week ", wk, ": no MISI column found; skipping vulnerability split.")
    next
  }

  id_vec <- as.character(obj@idents)
  target_idx <- which(id_vec == target_celltype)
  if (length(target_idx) < 20) {
    warning("Week ", wk, ": too few target cells for split (", target_celltype, ").")
    next
  }

  v <- as.numeric(md[[misi_col]])
  q_low <- stats::quantile(v[target_idx], probs = low_q, na.rm = TRUE)
  q_high <- stats::quantile(v[target_idx], probs = high_q, na.rm = TRUE)

  split_group <- rep("Other", length(v))
  split_group[target_idx[v[target_idx] <= q_low]] <- paste0(target_celltype, "_Exposed")
  split_group[target_idx[v[target_idx] >= q_high]] <- paste0(target_celltype, "_Shielded")
  split_group <- factor(split_group, levels = c(paste0(target_celltype, "_Exposed"), paste0(target_celltype, "_Shielded"), "Other"))

  md$vulnerability_group <- split_group
  vulnerability_meta[[wk]] <- md[, c("vulnerability_group", misi_col, "x_cent", "y_cent", "x_um", "y_um")[c("vulnerability_group", misi_col, "x_cent", "y_cent", "x_um", "y_um") %in% colnames(md)], drop = FALSE]

  agg_df <- aggregate_probcell_to_groups(obj, split_group)
  agg_df$week <- wk
  vulnerability_tables[[wk]] <- agg_df

  gc(verbose = FALSE)
}

master_lr_df <- dplyr::bind_rows(master_lr)
thesis_lr_df <- dplyr::bind_rows(thesis_lr)
vulnerability_df <- dplyr::bind_rows(vulnerability_tables)

master_path <- file.path(DIR_TABLES, "04_master_ligand_receptor_table.csv")
thesis_path <- file.path(DIR_TABLES, "04_thesis_pathway_isolation_table.csv")
vuln_path <- file.path(DIR_TABLES, "04_vulnerability_split_lr_table.csv")
utils::write.csv(master_lr_df, master_path, row.names = FALSE)
utils::write.csv(thesis_lr_df, thesis_path, row.names = FALSE)
utils::write.csv(vulnerability_df, vuln_path, row.names = FALSE)

cache <- list(
  week_paths = week_paths,
  merged_path = resolve_object_path(file.path(DIR_OBJECTS, "03b_spatial_cellchat_full_merged_complete.rds")),
  target_celltype = target_celltype,
  low_q = low_q,
  high_q = high_q,
  vulnerability_meta = vulnerability_meta,
  tables = list(
    master = master_path,
    thesis = thesis_path,
    vulnerability = vuln_path
  )
)

cache_path <- file.path(DIR_OBJECTS, "04a_atlas_compute_cache.rds")
saveRDS(cache, cache_path)

append_methods_log(
  "Phase 4-5 compute",
  c(
    "Loaded weekly completed full-pathway CellChat objects and extracted communication tables with subsetCommunication().",
    "Generated master LR table and thesis-pathway isolation table (MMP, IDO1, TGFb) as CSV exports.",
    paste0("Split target cell type '", target_celltype, "' into Exposed/Shielded using MISI quantiles: low<=", low_q, ", high>=", high_q, "."),
    "Re-aggregated LR communication from net$prob.cell by vulnerability groups (source/target) and exported as a dedicated table.",
    "Saved compute cache for plotting stage at output/objects/04a_atlas_compute_cache.rds."
  )
)

cat("Compute stage complete.\n")
cat("- ", master_path, "\n", sep = "")
cat("- ", thesis_path, "\n", sep = "")
cat("- ", vuln_path, "\n", sep = "")
cat("- ", cache_path, "\n", sep = "")
