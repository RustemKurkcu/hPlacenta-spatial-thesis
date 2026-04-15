source("R/config.R")
source("R/helpers_io.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(qs)
  library(Matrix)
})

log_file <- file.path(CFG$dirs$logs, "pipeline.log")
log_msg("18_cellphonedb_export starting.", log_file = log_file)

ensure_dir(CFG$dirs$tables)
ensure_dir(CFG$dirs$objects)

safe_write_table <- function(df_or_vec, path, ...) {
  ensure_dir(dirname(path))
  ok <- FALSE
  last_err <- NULL
  for (i in 1:3) {
    tryCatch({
      write.table(df_or_vec, path, ...)
      ok <- TRUE
    }, error = function(e) {
      last_err <<- e
      Sys.sleep(0.5 * i)
    })
    if (ok) break
  }
  if (!ok) {
    stop(
      "Failed to write file after retries: ", path,
      ". Possible causes: disk full/quota, file lock, or permission issue. Last error: ",
      conditionMessage(last_err)
    )
  }
}

get_expr_data <- function(seu_obj, assay = "RNA", layer_name = "data") {
  m <- tryCatch(SeuratObject::LayerData(seu_obj, assay = assay, layer = layer_name), error = function(e) NULL)
  if (!is.null(m) && ncol(m) > 0) return(m)
  m <- tryCatch(SeuratObject::GetAssayData(seu_obj, assay = assay, layer = layer_name), error = function(e) NULL)
  if (!is.null(m) && ncol(m) > 0) return(m)
  m <- tryCatch(SeuratObject::GetAssayData(seu_obj, assay = assay, slot = layer_name), error = function(e) NULL)
  if (is.null(m) || ncol(m) == 0) stop("Could not retrieve expression layer '", layer_name, "' from assay '", assay, "'.")
  m
}

obj_candidates <- c(
  file.path(CFG$dirs$objects, "seu_with_misi.qs"),
  file.path(CFG$dirs$objects, "seu_with_architecture_transfer_lite.qs"),
  file.path(CFG$dirs$objects, "seu_with_architecture_transfer.qs"),
  file.path(CFG$dirs$objects, "seu_with_scores.qs"),
  file.path(CFG$dirs$objects, "seu_clean.qs")
)
obj_candidates <- obj_candidates[file.exists(obj_candidates)]
if (length(obj_candidates) == 0) stop("No Seurat object found for script 18_cellphonedb_export")

seu <- NULL
obj_path <- NULL
for (cand in obj_candidates) {
  seu_try <- tryCatch(qs::qread(cand, nthreads = 1), error = function(e) NULL)
  if (!is.null(seu_try) && inherits(seu_try, "Seurat")) {
    seu <- seu_try
    obj_path <- cand
    break
  }
}
if (is.null(seu)) stop("Could not load any candidate Seurat object")
log_msg("Using object: ", obj_path, log_file = log_file)

if (!"condition" %in% colnames(seu@meta.data)) {
  if (all(c(CFG$cols$infection, CFG$cols$hpi) %in% colnames(seu@meta.data))) {
    seu$condition <- paste0(seu[[CFG$cols$infection]][, 1], "_", seu[[CFG$cols$hpi]][, 1])
  } else {
    stop("condition (or infection+hpi) metadata is required")
  }
}

cell_type_col <- NULL
if (!is.null(CFG$cols) && !is.null(CFG$cols$cell_type) && length(CFG$cols$cell_type) == 1 && nzchar(CFG$cols$cell_type)) {
  cell_type_col <- CFG$cols$cell_type
}
cell_type_candidates <- unique(c(cell_type_col, "celltype_refined", "celltype_author", "cell_type", "celltype"))
cell_type_candidates <- cell_type_candidates[!is.na(cell_type_candidates) & nzchar(cell_type_candidates)]
cell_type_col <- cell_type_candidates[cell_type_candidates %in% colnames(seu@meta.data)][1]
if (is.na(cell_type_col) || is.null(cell_type_col)) {
  stop(
    "No cell-type column available for cpdb metadata. Tried: ",
    paste(cell_type_candidates, collapse = ", ")
  )
}

all_conditions <- unique(as.character(seu$condition))
cond_ui <- all_conditions[grepl("^UI", all_conditions, ignore.case = TRUE) & grepl("24", all_conditions)]

cond_pathogen_24 <- all_conditions[
  grepl("24", all_conditions) &
    grepl("lm|l\\.monocytogenes|listeria|pf|plasmodium|tg|toxoplasma", all_conditions, ignore.case = TRUE)
]

if (length(cond_ui) == 0 || length(cond_pathogen_24) == 0) {
  log_msg("Could not auto-detect UI_24h + pathogen 24h conditions; exporting all conditions instead.", log_file = log_file)
  cond_keep <- all_conditions
} else {
  cond_keep <- unique(c(cond_ui[1], cond_pathogen_24))
}

cells_keep <- rownames(seu@meta.data)[as.character(seu$condition) %in% cond_keep]
if (length(cells_keep) < 100) stop("Too few cells after condition selection for CellPhoneDB export")
seu_sub <- subset(seu, cells = cells_keep)

export_prefix <- paste0("cpdb_", gsub("[^A-Za-z0-9_.-]", "_", paste(cond_keep, collapse = "_vs_")))

# metadata
meta <- seu_sub@meta.data
cpdb_meta <- data.frame(
  Cell = rownames(meta),
  cell_type = as.character(meta[[cell_type_col]]),
  stringsAsFactors = FALSE
)
safe_write_table(cpdb_meta, file.path(CFG$dirs$tables, paste0(export_prefix, "_meta.txt")), sep = "\t", quote = FALSE, row.names = FALSE)

# counts (sparse mtx)
expr <- get_expr_data(seu_sub, assay = "RNA", layer_name = "data")
if (!inherits(expr, "dgCMatrix")) expr <- as(expr, "dgCMatrix")

Matrix::writeMM(expr, file.path(CFG$dirs$tables, paste0(export_prefix, "_counts.mtx")))
safe_write_table(rownames(expr), file.path(CFG$dirs$tables, paste0(export_prefix, "_genes.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
safe_write_table(colnames(expr), file.path(CFG$dirs$tables, paste0(export_prefix, "_barcodes.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# DEG export loop: cell-type-specific DEGs for each test condition vs UI baseline.
seu_sub$joint_id <- paste0(as.character(seu_sub$condition), "_", as.character(seu_sub[[cell_type_col]]))
Idents(seu_sub) <- "joint_id"
ui_ref <- cond_keep[grepl("^UI", cond_keep, ignore.case = TRUE)][1]
if (is.na(ui_ref) || is.null(ui_ref)) {
  ui_ref <- cond_keep[1]
}
test_conditions <- cond_keep[cond_keep != ui_ref]
cell_types <- unique(as.character(seu_sub[[cell_type_col]]))

deg_status <- lapply(test_conditions, function(test_cond) {
  cond_deg_list <- list()

  for (ct in cell_types) {
    ident_test <- paste0(test_cond, "_", ct)
    ident_ref <- paste0(ui_ref, "_", ct)

    n_test <- sum(as.character(seu_sub$joint_id) == ident_test, na.rm = TRUE)
    n_ref <- sum(as.character(seu_sub$joint_id) == ident_ref, na.rm = TRUE)
    if (n_test < 10 || n_ref < 10) next

    deg <- tryCatch(
      FindMarkers(seu_sub, ident.1 = ident_test, ident.2 = ident_ref, assay = "RNA", logfc.threshold = 0.1, min.pct = 0.05),
      error = function(e) data.frame()
    )
    if (nrow(deg) == 0) next

    deg$gene <- rownames(deg)
    pcol <- if ("p_val_adj" %in% colnames(deg)) "p_val_adj" else if ("p_val" %in% colnames(deg)) "p_val" else NULL
    if (is.null(pcol)) next

    logcol <- if ("avg_log2FC" %in% colnames(deg)) "avg_log2FC" else if ("avg_logFC" %in% colnames(deg)) "avg_logFC" else NULL
    ct_deg <- deg %>%
      filter(.data[[pcol]] < 0.05) %>%
      transmute(
        cluster = ct,
        gene = gene,
        pvalue = .data[[pcol]],
        log2fc = if (!is.null(logcol)) .data[[logcol]] else NA_real_
      )
    if (nrow(ct_deg) > 0) cond_deg_list[[ct]] <- ct_deg
  }

  out_name <- paste0(
    export_prefix, "_DEGs_",
    gsub("[^A-Za-z0-9_.-]", "_", test_cond),
    "_vs_",
    gsub("[^A-Za-z0-9_.-]", "_", ui_ref),
    ".txt"
  )
  out_file <- file.path(CFG$dirs$tables, out_name)

  if (length(cond_deg_list) == 0) {
    return(data.frame(reference = ui_ref, test = test_cond, n_deg = 0, outfile = out_name, stringsAsFactors = FALSE))
  }

  final_cond_degs <- bind_rows(cond_deg_list)
  safe_write_table(final_cond_degs, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  data.frame(reference = ui_ref, test = test_cond, n_deg = nrow(final_cond_degs), outfile = out_name, stringsAsFactors = FALSE)
})
deg_status <- bind_rows(deg_status)
if (nrow(deg_status) > 0) {
  write.csv(deg_status, file.path(CFG$dirs$tables, paste0(export_prefix, "_deg_status.csv")), row.names = FALSE)
}

status <- data.frame(
  object_used = obj_path,
  export_prefix = export_prefix,
  condition_reference = ui_ref,
  condition_test = if (length(test_conditions) > 0) paste(test_conditions, collapse = ";") else NA_character_,
  n_cells_exported = length(cells_keep),
  stringsAsFactors = FALSE
)
write.csv(status, file.path(CFG$dirs$tables, paste0(export_prefix, "_status.csv")), row.names = FALSE)

log_msg("18_cellphonedb_export done.", log_file = log_file)
