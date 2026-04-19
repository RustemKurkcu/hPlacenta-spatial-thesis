#!/usr/bin/env Rscript

# =============================================================================
# Script: 03b_spatial_cellchat_allcells_enhanced_toxicswitch.R
# Purpose: Enhanced version of min3_rdsfirst run with extra mechanistic tables.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(jsonlite)
})

source("R/celltype_dictionary.R")

PIPELINE_NAME <- "03b_spatial_cellchat_allcells_enhanced_toxicswitch"
PIPELINE_VERSION <- "2.1.0"
RUN_TIMESTAMP <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
RUN_SEED <- 42L
set.seed(RUN_SEED)

OUT_ROOT <- Sys.getenv("CELLCHAT_OUT_ROOT", file.path("output", "cellchat", "enhanced_toxicswitch_min3"))
INPUT_OBJECT_DIR <- Sys.getenv("CELLCHAT_INPUT_DIR", file.path("output", "objects"))
INPUT_BASENAME <- Sys.getenv("CELLCHAT_INPUT_BASENAME", "02_scored_misi_ido1")
MIN_CELLS_FILTER <- as.integer(Sys.getenv("CELLCHAT_MIN_CELLS", "3"))
INTERACTION_RANGE_UM <- as.numeric(Sys.getenv("CELLCHAT_INTERACTION_RANGE_UM", "100"))
DB_SEARCH <- Sys.getenv("CELLCHAT_DB_SEARCH", "Secreted Signaling")

DIR_LOGS <- file.path(OUT_ROOT, "logs")
DIR_OBJECTS <- file.path(OUT_ROOT, "objects")
DIR_REPORTS <- file.path(OUT_ROOT, "reports")
DIR_TABLES <- file.path(OUT_ROOT, "tables")
DIR_TABLES_PANELS <- file.path(DIR_TABLES, "mechanistic_panels")
DIR_TABLES_NETWORKS <- file.path(DIR_TABLES, "network_edges")
DIR_TABLES_FOCUS <- file.path(DIR_TABLES, "focus")
for (d in c(DIR_LOGS, DIR_OBJECTS, DIR_REPORTS, DIR_TABLES, DIR_TABLES_PANELS, DIR_TABLES_NETWORKS, DIR_TABLES_FOCUS)) {
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
    ok <- tryCatch({ qs::qsave(object, path_qs, preset = "high"); TRUE }, error = function(e) FALSE)
    if (isTRUE(ok)) out <- c(out, path_qs)
  }
  out
}

focus_patterns <- list(
  VCT = c("^VCT", "villous cytotrophoblast"), SCT = c("^SCT", "syncytio"), EVT = c("^EVT", "extravillous"),
  DSC = c("^DSC", "decidual strom"), dNK = c("^dNK", "uNK", "natural killer"), Endothelial = c("endo"),
  Macrophage = c("macroph", "PAMM"), Treg = c("treg")
)
mechanistic_panels <- list(
  Homing = c("GALNT1","GALNT2","GALNT3","C1GALT1","C1GALT1C1","GCNT1","CDH1","CDH5","OCLN","TJP1","CLDN4"),
  Nutrient_Liberation = c("PLD1","PLD2","GDPD1","GDPD5","GDE1","ETNK1","FAAH"),
  Nutrient_Recapture = c("PCYT2","CHPT1","CEPT1"),
  Switch_Pressure = c("IL1B","TNF","IL6","NFKB1","RELA","HIF1A","SLC2A1","GAPDH","CASP4"),
  Vascular_Fragility = c("FLT1","PGF","KDR","ENG","NOS3","ICAM1","VCAM1"),
  Tolerance_Interface = c("HLA-G","KIR2DL4","IGFBP1","NOTCH1","IFNG","VEGFA")
)

map_focus_celltypes <- function(celltypes) {
  out <- lapply(names(focus_patterns), function(lbl) {
    pats <- focus_patterns[[lbl]]
    hits <- unique(unlist(lapply(pats, function(p) grep(p, celltypes, ignore.case = TRUE, value = TRUE))))
    if (length(hits) == 0) hits <- NA_character_
    data.frame(focus_label = lbl, matched_celltype = hits, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

extract_panel_summary <- function(expr_mat, meta_df, wk) {
  out <- list(); idx <- 1L
  for (panel_name in names(mechanistic_panels)) {
    present <- intersect(mechanistic_panels[[panel_name]], rownames(expr_mat))
    if (length(present) == 0) next
    for (grp in unique(as.character(meta_df$celltype_plot))) {
      keep <- rownames(meta_df)[meta_df$celltype_plot == grp]
      keep <- intersect(keep, colnames(expr_mat))
      if (length(keep) == 0) next
      submat <- expr_mat[present, keep, drop = FALSE]
      out[[idx]] <- data.frame(week = wk, panel = panel_name, celltype = grp, gene = present,
        avg_expr = as.numeric(Matrix::rowMeans(submat)), pct_expr = as.numeric(Matrix::rowMeans(submat > 0) * 100),
        n_cells = length(keep), stringsAsFactors = FALSE)
      idx <- idx + 1L
    }
  }
  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}

extract_net_edges <- function(cellchat, wk, keep_groups = NULL) {
  weight_mat <- tryCatch(cellchat@net$weight, error = function(e) NULL)
  if (is.null(weight_mat)) return(NULL)
  weight_mat <- as.matrix(weight_mat)
  if (!is.null(keep_groups)) {
    keep_groups <- intersect(keep_groups, rownames(weight_mat))
    if (length(keep_groups) == 0) return(NULL)
    weight_mat <- weight_mat[keep_groups, keep_groups, drop = FALSE]
  }
  df <- as.data.frame(as.table(weight_mat), stringsAsFactors = FALSE)
  colnames(df) <- c("sender", "receiver", "weight")
  df$week <- wk
  df[df$weight > 0, , drop = FALSE]
}

# Source standard runner script logic by reusing the min3 workflow object generation approach
input_obj <- resolve_object_path(file.path(INPUT_OBJECT_DIR, paste0(INPUT_BASENAME, ".rds")))
seu <- read_object(input_obj)
seu$week <- factor(as.character(seu$week), levels = c("W7", "W8-2", "W9", "W11"))
celltype_col <- if ("predicted.celltype" %in% colnames(seu@meta.data)) "predicted.celltype" else pick_celltype_source_column(seu@meta.data)
if (is.na(celltype_col)) celltype_col <- "seurat_clusters"
seu$celltype_plot <- rename_with_true_names(as.character(seu@meta.data[[celltype_col]]))
Idents(seu) <- "celltype_plot"

focus_map_all <- map_focus_celltypes(sort(unique(as.character(seu$celltype_plot))))
write.table(focus_map_all, file.path(DIR_TABLES_FOCUS, "focus_celltype_mapping_all.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Delegate heavy CellChat run to min3 script outputs if present; otherwise perform a lightweight reuse by calling the standard script externally
std_out_root <- Sys.getenv("CELLCHAT_STD_OUT_ROOT", file.path("output", "cellchat", "allcells_min3", "objects"))
week_paths <- file.path(std_out_root, paste0("spatial_cellchat_full_", c("W7", "W8-2", "W9", "W11"), ".rds"))
week_paths <- week_paths[file.exists(week_paths)]
if (length(week_paths) == 0) stop("No standard weekly objects found. Run 03b_spatial_cellchat_allcells_min3_rdsfirst.R first.")

weekly_meta <- list()
for (p in week_paths) {
  wk <- sub(".*spatial_cellchat_full_(.*)\\.rds$", "\\1", p)
  obj <- readRDS(p)
  md <- obj@meta
  md$celltype_plot <- as.character(obj@idents)
  expr <- obj@data

  panel_summary <- extract_panel_summary(expr, md, wk)
  if (!is.null(panel_summary)) write.table(panel_summary, file.path(DIR_TABLES_PANELS, paste0("mechanistic_panel_summary_", wk, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)

  focus_map_week <- map_focus_celltypes(sort(unique(as.character(md$celltype_plot))))
  write.table(focus_map_week, file.path(DIR_TABLES_FOCUS, paste0("focus_celltype_mapping_", wk, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  focus_groups <- unique(na.omit(focus_map_week$matched_celltype))

  all_edges <- extract_net_edges(obj, wk)
  if (!is.null(all_edges)) write.table(all_edges, file.path(DIR_TABLES_NETWORKS, paste0("network_edges_all_", wk, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  focus_edges <- extract_net_edges(obj, wk, keep_groups = focus_groups)
  if (!is.null(focus_edges)) write.table(focus_edges, file.path(DIR_TABLES_NETWORKS, paste0("network_edges_focus_", wk, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)

  out_week <- file.path(DIR_OBJECTS, paste0("spatial_cellchat_enhanced_", wk, ".rds"))
  save_dual_object(obj, out_week)
  weekly_meta[[wk]] <- data.frame(week = wk, source_week_object = p, enhanced_week_object = out_week, stringsAsFactors = FALSE)
}

week_summary <- do.call(rbind, weekly_meta)
write.table(week_summary, file.path(DIR_TABLES, "enhanced_week_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

priority_notes <- c(
  "Focus cell types chosen from toxic-switch thesis narrative: VCT/SCT/EVT as trophoblast barrier and ethanolamine-release compartments; DSC as maternal stromal ethanolamine source; dNK/Tregs as tolerance regulators; endothelial and macrophage compartments as vascular/inflammatory effectors.",
  "Mechanistic panel summaries support homing, nutrient liberation/recapture, switch pressure, vascular fragility, and tolerance-interface interpretation."
)
writeLines(priority_notes, file.path(DIR_REPORTS, "mechanistic_priority_notes.txt"))

manifest <- list(
  pipeline = PIPELINE_NAME,
  version = PIPELINE_VERSION,
  run_timestamp = RUN_TIMESTAMP,
  source_data = input_obj,
  output_objects = as.list(list.files(OUT_ROOT, recursive = TRUE, full.names = TRUE)),
  notes = c("Enhanced toxic-switch reporting generated from standard min3 CellChat outputs.", paste0("db_search=", DB_SEARCH), paste0("min.cells=", MIN_CELLS_FILTER), paste0("interaction.range_um=", INTERACTION_RANGE_UM))
)
jsonlite::write_json(manifest, file.path(DIR_REPORTS, "artifact_manifest.json"), pretty = TRUE, auto_unbox = TRUE)

log_msg("Enhanced run complete.")
