# ======================================================================
# scripts/06_cell_communication/06C_cellchat_optional.R
#
# ADVANCED MODULE (optional): CellChat analysis.
#
# This script is fully optional.
# It will *skip* gracefully unless CellChat is installed.
#
# References:
#   * Jin et al. CellChat (Nat Commun 2021) and the Nat Protocols guide.
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

source("config/config.R")
source("scripts/R/utils.R")

# --- Pre-run Checks ---

if (!isTRUE(RUN_OPTIONAL_HEAVY)) {
  if (interactive()) {
    stop("RUN_OPTIONAL_HEAVY=FALSE so 06C is intentionally skipped. Set RUN_OPTIONAL_HEAVY <- TRUE in config to run.", call. = FALSE)
  }
  message("RUN_OPTIONAL_HEAVY=FALSE; skipping 06C.")
  quit(save = "no")
}

if (!requireNamespace("CellChat", quietly = TRUE)) {
  stop("CellChat is not installed. Install it or disable 06C by setting RUN_OPTIONAL_HEAVY=FALSE in config.")
}

suppressPackageStartupMessages({
  library(CellChat)
})

ensure_dir(DIR_OBJS); ensure_dir(DIR_FIGURES); ensure_dir(DIR_TABLES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "06C_cellchat_optional.log")
log_msg("Starting CellChat optional module.", logfile)

# --- Helper Functions ---

resolve_cellchat_human_db <- function() {
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
  if (is.null(db)) stop("Could not resolve CellChatDB.human for 06C.")
  db
}

safe_read_rds_with_fallback <- function(candidates, label) {
  for (pth in candidates) {
    if (is.null(pth) || is.na(pth) || !nzchar(pth) || !file.exists(pth)) next
    obj <- tryCatch(readRDS(pth), error = function(e) e)
    if (!inherits(obj, "error")) {
      log_msg(paste0("[06C] Loaded ", label, " object: ", pth), logfile)
      return(obj)
    }
  }
  stop("[06C] Could not load ", label, " object from candidate paths.")
}

pick_group_col <- function(md) {
  cand <- unique(c("celltype_final_refined", "celltype_final_conservative", COL_PRED_CELLTYPE,
                   "celltype_author", "celltype", "cluster"))
  cand <- cand[!is.na(cand) & nzchar(cand)]
  hit <- cand[cand %in% colnames(md)]
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

# --- Main Logic ---

obj_paths <- list(
  SlideTags = unique(c(file.path(DIR_OBJS, "slidetags_harmonized.rds"), 
                       file.path(DIR_OBJS, "slidetags_mapped_to_multiome.rds"))),
  STARmap   = unique(c(file.path(DIR_OBJS, "starmap_harmonized.rds"), 
                       file.path(DIR_OBJS, "starmap_mapped_to_multiome.rds")))
)

run_one <- function(path_candidates, name) {
  obj <- safe_read_rds_with_fallback(path_candidates, name)
  
  group_col <- pick_group_col(obj@meta.data)
  if (is.null(group_col)) {
    log_msg(paste0(name, ": no usable group column found; skipping."), logfile)
    return(NULL)
  }
  log_msg(paste0(name, ": using group column ", group_col), logfile)
  
  # Prepare Assay Data
  assay_use <- if ("SCT" %in% names(obj@assays)) "SCT" else "RNA"
  DefaultAssay(obj) <- assay_use
  obj <- safe_join_layers(obj, assay = assay_use)
  if (!has_data_layer(obj, assay_use)) obj <- NormalizeData(obj, verbose = FALSE)
  
  meta <- obj@meta.data
  meta$cellgroup <- as.character(meta[[group_col]])
  # Ensure week column exists
  obj <- ensure_week_column(obj, COL_WEEK_CANDIDATES)
  meta <- obj@meta.data
  meta$week <- as.character(meta$week %||% "all")
  
  out_list <- list()
  week_rows <- list()
  
  # Resolve Database Once
  CellChatDB <- resolve_cellchat_human_db()
  
  for (w in sort(unique(meta$week))) {
    cells_use <- rownames(meta)[meta$week == w]
    if (length(cells_use) < 200) {
      log_msg(paste0("[06C] ", name, " week ", w, ": insufficient cells (", length(cells_use), "); skipping."), logfile)
      next
    }
    
    x <- subset(obj, cells = cells_use)
    data.input <- get_assay_matrix(x, assay = assay_use, layer = "data")
    if (is.null(data.input)) next
    
    meta.input <- x@meta.data
    
    # Run CellChat Pipeline
    cellchat <- createCellChat(object = data.input, meta = meta.input, group.by = "cellgroup")
    cellchat@DB <- CellChatDB
    
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    
    out_list[[w]] <- cellchat
    
    # Export edges
    df <- subsetCommunication(cellchat)
    write.csv(df, file.path(DIR_TABLES, paste0(name, "_CellChat_week_", w, "_edges.csv")), row.names = FALSE)
    
    week_rows[[length(week_rows) + 1]] <- data.frame(
      dataset = name,
      week = as.character(w),
      n_cells = ncol(data.input),
      n_edges = nrow(df),
      stringsAsFactors = FALSE
    )
    log_msg(paste0("[06C] ", name, " week ", w, ": processed, edges=", nrow(df)), logfile)
  }
  
  if (length(out_list) > 0) {
    saveRDS(out_list, file.path(DIR_OBJS, paste0(name, "_cellchat_by_week.rds")))
    summary_df <- bind_rows(week_rows) %>% mutate(timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    write.csv(summary_df, file.path(DIR_TABLES, paste0(name, "_06C_cellchat_run_summary.csv")), row.names = FALSE)
  }
  
  invisible(out_list)
}

# --- Execution ---

for (nm in names(obj_paths)) {
  log_msg(paste0("Starting analysis for: ", nm), logfile)
  tryCatch({
    run_one(obj_paths[[nm]], nm)
  }, error = function(e) {
    log_msg(paste0("[06C] ", nm, " failed: ", e$message), logfile)
  })
}

log_msg("06C CellChat module complete.", logfile)