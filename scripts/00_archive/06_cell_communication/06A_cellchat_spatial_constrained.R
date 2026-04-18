<<<<<<< HEAD
# ======================================================================
# scripts/06_cell_communication/06A_cellchat_spatial_constrained.R
# OPTIONAL: CellChat analysis (ligand-receptor inference).
#
# This step is intentionally optional because CellChat can be heavy and
# requires additional dependencies.
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

source("config/config.R")
source("scripts/R/utils.R")

# --- Pre-run Checks ---

ensure_dir(DIR_OBJS); ensure_dir(DIR_TABLES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "06A_cellchat_spatial_constrained.log")
log_msg("Starting CellChat module 06A.", logfile)

if (!isTRUE(RUN_OPTIONAL_HEAVY)) {
  log_msg("RUN_OPTIONAL_HEAVY=FALSE. Skipping CellChat.", logfile)
  if (interactive()) {
    stop("RUN_OPTIONAL_HEAVY=FALSE so 06A is intentionally skipped. Set RUN_OPTIONAL_HEAVY <- TRUE to run CellChat.", call. = FALSE)
  }
  quit(save = "no")
}

if (!requireNamespace("CellChat", quietly = TRUE)) {
  stop("CellChat is not installed. Install it or set RUN_OPTIONAL_HEAVY=FALSE in config.")
}

suppressPackageStartupMessages({
  library(CellChat)
})

# --- Helper Functions ---

safe_read_rds_with_fallback <- function(candidates, label, logfile = NULL) {
  for (pth in candidates) {
    if (is.null(pth) || is.na(pth) || !nzchar(pth) || !file.exists(pth)) next
    obj <- tryCatch(readRDS(pth), error = function(e) e)
    if (!inherits(obj, "error")) {
      log_msg(paste0("[06A] Loaded ", label, " object: ", pth), logfile)
      return(obj)
    }
  }
  stop("[06A] Unable to load ", label, " object from candidate paths.")
}

pick_cellchat_group_col <- function(md) {
  cand <- unique(c("celltype_final_refined", "celltype_final_conservative", COL_PRED_CELLTYPE,
                   "celltype_author", "celltype", "cluster"))
  cand <- cand[!is.na(cand) & nzchar(cand)]
  hit <- cand[cand %in% colnames(md)]
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

resolve_cellchat_human_db <- function(logfile = NULL) {
  db <- tryCatch(get("CellChatDB.human", envir = asNamespace("CellChat")), error = function(e) NULL)
  if (is.null(db) && requireNamespace("CellChatDB", quietly = TRUE)) {
    db <- tryCatch(get("CellChatDB.human", envir = asNamespace("CellChatDB")), error = function(e) NULL)
  }
  if (is.null(db)) {
    env_db <- new.env(parent = emptyenv())
    try(utils::data("CellChatDB.human", package = "CellChat", envir = env_db), silent = TRUE)
    if (exists("CellChatDB.human", envir = env_db, inherits = FALSE)) {
      db <- get("CellChatDB.human", envir = env_db, inherits = FALSE)
    }
  }
  if (is.null(db)) stop("Could not resolve CellChatDB.human for 06A.")
  db
}

# --- Main Logic ---

# Load Slide-tags (transcriptome-wide data is best for CellChat)
obj <- safe_read_rds_with_fallback(
  candidates = unique(c(
    file.path(DIR_OBJS, "slidetags_harmonized.rds"),
    file.path(DIR_OBJS, "slidetags_mapped_to_multiome.rds"),
    file.path(DIR_OBJS, "slidetags_mapped.rds")
  )),
  label = "SlideTags", logfile = logfile
)

group_col <- pick_cellchat_group_col(obj@meta.data)
if (is.null(group_col)) stop("No usable cell-group column found for CellChat.")
obj$celltype_use <- as.character(obj@meta.data[[group_col]])

log_msg(paste0("[06A] Using group column: ", group_col), logfile)
log_msg(paste0("[06A] Unique groups: ", length(unique(obj$celltype_use))), logfile)

# Choose assay & ensure normalized data layer exists
assay_use <- if ("SCT" %in% names(obj@assays)) "SCT" else "RNA"
DefaultAssay(obj) <- assay_use
obj <- safe_join_layers(obj, assay = assay_use)
if (!has_data_layer(obj, assay_use)) obj <- NormalizeData(obj, verbose = FALSE)

# Prepare CellChat object
data_input <- get_assay_matrix(obj, assay = assay_use, layer = "data")
meta <- data.frame(labels = obj$celltype_use, row.names = colnames(obj))

cellchat <- createCellChat(object = data_input, meta = meta, group.by = "labels")

# Database: human (robust across versions)
cellchat@DB <- resolve_cellchat_human_db(logfile = logfile)

# Run standard workflow
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Save result
out_rds <- file.path(DIR_OBJS, "cellchat_slidetags.rds")
saveRDS(cellchat, out_rds)
log_msg(paste0("Saved CellChat object: ", out_rds), logfile)

# Export results table
df_net <- subsetCommunication(cellchat)
write.csv(df_net, file.path(DIR_TABLES, "cellchat_slidetags_interactions.csv"), row.names = FALSE)

# Save summary
summary_06a <- data.frame(
  module = "06A",
  dataset = "SlideTags",
  group_col = group_col,
  assay = assay_use,
  n_cells = ncol(data_input),
  n_genes = nrow(data_input),
  n_groups = length(unique(obj$celltype_use)),
  n_edges = nrow(df_net),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  stringsAsFactors = FALSE
)
write.csv(summary_06a, file.path(DIR_TABLES, "06A_cellchat_run_summary.csv"), row.names = FALSE)

=======
# ======================================================================
# scripts/06_cell_communication/06A_cellchat_spatial_constrained.R
# OPTIONAL: CellChat analysis (ligand-receptor inference).
#
# This step is intentionally optional because CellChat can be heavy and
# requires additional dependencies.
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

source("config/config.R")
source("scripts/R/utils.R")

# --- Pre-run Checks ---

ensure_dir(DIR_OBJS); ensure_dir(DIR_TABLES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "06A_cellchat_spatial_constrained.log")
log_msg("Starting CellChat module 06A.", logfile)

if (!isTRUE(RUN_OPTIONAL_HEAVY)) {
  log_msg("RUN_OPTIONAL_HEAVY=FALSE. Skipping CellChat.", logfile)
  if (interactive()) {
    stop("RUN_OPTIONAL_HEAVY=FALSE so 06A is intentionally skipped. Set RUN_OPTIONAL_HEAVY <- TRUE to run CellChat.", call. = FALSE)
  }
  quit(save = "no")
}

if (!requireNamespace("CellChat", quietly = TRUE)) {
  stop("CellChat is not installed. Install it or set RUN_OPTIONAL_HEAVY=FALSE in config.")
}

suppressPackageStartupMessages({
  library(CellChat)
})

# --- Helper Functions ---

safe_read_rds_with_fallback <- function(candidates, label, logfile = NULL) {
  for (pth in candidates) {
    if (is.null(pth) || is.na(pth) || !nzchar(pth) || !file.exists(pth)) next
    obj <- tryCatch(readRDS(pth), error = function(e) e)
    if (!inherits(obj, "error")) {
      log_msg(paste0("[06A] Loaded ", label, " object: ", pth), logfile)
      return(obj)
    }
  }
  stop("[06A] Unable to load ", label, " object from candidate paths.")
}

pick_cellchat_group_col <- function(md) {
  cand <- unique(c("celltype_final_refined", "celltype_final_conservative", COL_PRED_CELLTYPE,
                   "celltype_author", "celltype", "cluster"))
  cand <- cand[!is.na(cand) & nzchar(cand)]
  hit <- cand[cand %in% colnames(md)]
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

resolve_cellchat_human_db <- function(logfile = NULL) {
  db <- tryCatch(get("CellChatDB.human", envir = asNamespace("CellChat")), error = function(e) NULL)
  if (is.null(db) && requireNamespace("CellChatDB", quietly = TRUE)) {
    db <- tryCatch(get("CellChatDB.human", envir = asNamespace("CellChatDB")), error = function(e) NULL)
  }
  if (is.null(db)) {
    env_db <- new.env(parent = emptyenv())
    try(utils::data("CellChatDB.human", package = "CellChat", envir = env_db), silent = TRUE)
    if (exists("CellChatDB.human", envir = env_db, inherits = FALSE)) {
      db <- get("CellChatDB.human", envir = env_db, inherits = FALSE)
    }
  }
  if (is.null(db)) stop("Could not resolve CellChatDB.human for 06A.")
  db
}

# --- Main Logic ---

# Load Slide-tags (transcriptome-wide data is best for CellChat)
obj <- safe_read_rds_with_fallback(
  candidates = unique(c(
    file.path(DIR_OBJS, "slidetags_harmonized.rds"),
    file.path(DIR_OBJS, "slidetags_mapped_to_multiome.rds"),
    file.path(DIR_OBJS, "slidetags_mapped.rds")
  )),
  label = "SlideTags", logfile = logfile
)

group_col <- pick_cellchat_group_col(obj@meta.data)
if (is.null(group_col)) stop("No usable cell-group column found for CellChat.")
obj$celltype_use <- as.character(obj@meta.data[[group_col]])

log_msg(paste0("[06A] Using group column: ", group_col), logfile)
log_msg(paste0("[06A] Unique groups: ", length(unique(obj$celltype_use))), logfile)

# Choose assay & ensure normalized data layer exists
assay_use <- if ("SCT" %in% names(obj@assays)) "SCT" else "RNA"
DefaultAssay(obj) <- assay_use
obj <- safe_join_layers(obj, assay = assay_use)
if (!has_data_layer(obj, assay_use)) obj <- NormalizeData(obj, verbose = FALSE)

# Prepare CellChat object
data_input <- get_assay_matrix(obj, assay = assay_use, layer = "data")
meta <- data.frame(labels = obj$celltype_use, row.names = colnames(obj))

cellchat <- createCellChat(object = data_input, meta = meta, group.by = "labels")

# Database: human (robust across versions)
cellchat@DB <- resolve_cellchat_human_db(logfile = logfile)

# Run standard workflow
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Save result
out_rds <- file.path(DIR_OBJS, "cellchat_slidetags.rds")
saveRDS(cellchat, out_rds)
log_msg(paste0("Saved CellChat object: ", out_rds), logfile)

# Export results table
df_net <- subsetCommunication(cellchat)
write.csv(df_net, file.path(DIR_TABLES, "cellchat_slidetags_interactions.csv"), row.names = FALSE)

# Save summary
summary_06a <- data.frame(
  module = "06A",
  dataset = "SlideTags",
  group_col = group_col,
  assay = assay_use,
  n_cells = ncol(data_input),
  n_genes = nrow(data_input),
  n_groups = length(unique(obj$celltype_use)),
  n_edges = nrow(df_net),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  stringsAsFactors = FALSE
)
write.csv(summary_06a, file.path(DIR_TABLES, "06A_cellchat_run_summary.csv"), row.names = FALSE)

>>>>>>> 54c1a5c675ca223d183141f8a75de090237ce902
log_msg("06A CellChat module complete.", logfile)