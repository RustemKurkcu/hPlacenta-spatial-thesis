# ======================================================================
# scripts/03_mapping/03B_map_starmap_to_multiome.R
# Map STARmap-ISS spatial data to the Multiome reference taxonomy.
#
# KEY FIXES:
#   1. JoinLayers for split layers (counts.W8-2, counts.W9, etc.)
#   2. Smart assay selection for mapping (gene overlap with reference)
#   3. Robust UMAP assay fallback when imputed assay is not finite/usable
#
# Output:
#   output/objects/starmap_mapped.rds
#   output/figures/starmap_*.png
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

source("config/config.R")
source("scripts/R/utils.R")
if (file.exists("Slide-tags/seurat_compat_utils.R")) source("Slide-tags/seurat_compat_utils.R")

ensure_dir(DIR_OBJS); ensure_dir(DIR_FIGURES); ensure_dir(DIR_TABLES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "03B_map_starmap_to_multiome.log")

# ---- Load reference ----
ref_path <- file.path(DIR_OBJS, "multiome_reference_processed.rds")
if (!file.exists(ref_path)) stop("Missing processed reference: ", ref_path, " (run 02A first)")
ref <- readRDS(ref_path)

norm_method <- ref@misc$norm_method_use %||%
  (if ("SCT" %in% names(ref@assays)) "SCT" else "LogNormalize")
log_msg(paste0("Reference normalization = ", norm_method), logfile)

# ---- Load STARmap ----
if (!file.exists(PATH_STARMAP_RDS)) stop("Missing STARmap RDS: ", PATH_STARMAP_RDS)
log_msg(paste0("Loading STARmap: ", PATH_STARMAP_RDS), logfile)
star <- readRDS(PATH_STARMAP_RDS)

print_seurat_diagnostic(star, "STARmap (as loaded)")

star <- fix_metadata_reduction_clash(star)
star <- ensure_week_column(star, COL_WEEK_CANDIDATES)
star <- ensure_spatial_coords(
  star,
  c("x_um", "spatial_x", "spatial_x_um", "X", "x", "center_x"),
  c("y_um", "spatial_y", "spatial_y_um", "Y", "y", "center_y")
)

# Author-provided labels
author_col <- pick_first_present(star@meta.data, c("predicted.celltype", COL_AUTHOR_CELLTYPE))
if (!is.null(author_col)) {
  star$celltype_author <- as.character(star@meta.data[[author_col]])
} else {
  star$celltype_author <- NA_character_
}

# ---- Determine best assay for mapping ----
assays_available <- names(star@assays)
log_msg(paste0("Available assays: ", paste(assays_available, collapse = ", ")), logfile)

best_assay <- NULL
best_overlap <- 0L
ref_genes <- rownames(ref[["RNA"]])

for (a in assays_available) {
  star <- safe_join_layers(star, a)
  a_genes <- rownames(star[[a]])
  overlap <- length(intersect(a_genes, ref_genes))
  log_msg(sprintf("  Assay '%s': %d genes, %d overlap with reference", a, length(a_genes), overlap), logfile)
  if (overlap > best_overlap) {
    best_overlap <- overlap
    best_assay <- a
  }
}

if (is.null(best_assay)) stop("No valid assay found in STARmap object")
log_msg(sprintf("Best assay for mapping: '%s' (%d gene overlap)", best_assay, best_overlap), logfile)

# ---- Check if already mapped ----
has_pred <- "predicted.celltype" %in% colnames(star@meta.data)
has_pred_id <- COL_PRED_CELLTYPE %in% colnames(star@meta.data)

if (has_pred && !has_pred_id) {
  star$predicted.id <- star$predicted.celltype
  log_msg("Copied predicted.celltype -> predicted.id", logfile)
  has_pred_id <- TRUE
}

if (has_pred_id) {
  log_msg("STARmap already has predicted labels. Skipping transfer.", logfile)
} else if (best_overlap < 100) {
  log_msg("WARNING: Gene overlap too low for reliable transfer. Skipping.", logfile)
  star$predicted.id <- NA_character_
} else {
  log_msg("Running label transfer...", logfile)
  DefaultAssay(star) <- best_assay
  
  if (norm_method == "LogNormalize") {
    star <- NormalizeData(star, verbose = FALSE)
    star <- FindVariableFeatures(star, verbose = FALSE)
    star <- ScaleData(star, features = VariableFeatures(star), verbose = FALSE)
    star <- RunPCA(star, verbose = FALSE)
    
    anchors <- FindTransferAnchors(
      reference = ref, query = star,
      normalization.method = "LogNormalize", dims = DEFAULT_DIMS
    )
  } else {
    star <- SCTransform(star, verbose = FALSE)
    anchors <- FindTransferAnchors(
      reference = ref, query = star,
      normalization.method = "SCT", dims = DEFAULT_DIMS
    )
  }
  
  pred <- TransferData(anchorset = anchors, refdata = ref$celltype_ref, dims = DEFAULT_DIMS)
  star <- AddMetaData(star, pred)
}

# prediction.score.max
if (!("predicted.celltype.score" %in% colnames(star@meta.data)) &&
    !(COL_PRED_SCORE_MAX %in% colnames(star@meta.data))) {
  score_cols <- grep("^prediction\\.score\\.", colnames(star@meta.data), value = TRUE)
  if (length(score_cols) > 0) {
    star@meta.data[[COL_PRED_SCORE_MAX]] <- apply(
      star@meta.data[, score_cols, drop = FALSE], 1, max, na.rm = TRUE
    )
  }
}

# ---- Compute UMAP ----
log_msg("Computing UMAP for STARmap...", logfile)

# Prefer imputed when available, but guard against non-finite/problematic assay values.
umap_candidates <- unique(c("imputed", best_assay, "RNA_raw", "RNA", names(star@assays)))
umap_candidates <- umap_candidates[umap_candidates %in% names(star@assays)]

assay_is_usable <- function(obj, assay_name) {
  if (!(assay_name %in% names(obj@assays))) return(FALSE)
  if (nrow(obj[[assay_name]]) < 50) return(FALSE)
  
  if (exists("assay_all_finite", mode = "function")) {
    finite_ok <- assay_all_finite(obj, assay = assay_name, layer = "data")
    if (!finite_ok) finite_ok <- assay_all_finite(obj, assay = assay_name, layer = "counts")
    return(isTRUE(finite_ok))
  }
  
  TRUE
}

vis_assay <- NULL
for (cand in umap_candidates) {
  if (assay_is_usable(star, cand)) {
    vis_assay <- cand
    break
  }
}
if (is.null(vis_assay)) vis_assay <- best_assay

log_msg(paste0("  Visualization assay selected: ", vis_assay), logfile)

vis_dims <- 1:min(30, ncol(star) - 1)
if (length(vis_dims) < 2) stop("Not enough cells to compute PCA/UMAP")

umap_ok <- FALSE
for (cand in unique(c(vis_assay, umap_candidates))) {
  if (!assay_is_usable(star, cand)) next
  log_msg(paste0("  Trying UMAP with assay: ", cand), logfile)
  star_try <- tryCatch(
    ensure_pca_umap(star, assay = cand, dims = vis_dims,
                    umap_return_model = FALSE, force = TRUE),
    error = function(e) {
      log_msg(paste0("  UMAP failed for assay '", cand, "': ", conditionMessage(e)), logfile)
      NULL
    }
  )
  
  if (!is.null(star_try) && "umap" %in% names(star_try@reductions)) {
    star <- star_try
    vis_assay <- cand
    umap_ok <- TRUE
    break
  }
}

if (!umap_ok) stop("Failed to compute UMAP for STARmap with all assay candidates")

# Optional tSNE
if (isTRUE(DO_TSNE) && !("tsne" %in% names(star@reductions))) {
  log_msg("Computing tSNE for STARmap...", logfile)
  perplex <- min(30, max(5, floor((ncol(star) - 1) / 3)))
  star <- RunTSNE(star, dims = vis_dims, reduction = "pca",
                  perplexity = perplex, verbose = FALSE)
}

# ---- Save ----
out_rds <- file.path(DIR_OBJS, "starmap_mapped.rds")
saveRDS(star, out_rds)
log_msg(paste0("Saved: ", out_rds), logfile)

# ---- Plots ----
log_msg("Generating STARmap QC plots...", logfile)

celltype_col <- if ("predicted.id" %in% colnames(star@meta.data)) "predicted.id" else "celltype_author"

p1 <- DimPlot(star, group.by = celltype_col, label = TRUE, repel = TRUE) +
  ggtitle(paste0("STARmap UMAP (", celltype_col, ")"))
save_plot(p1, file.path(DIR_FIGURES, "starmap_umap_celltype.png"), w = 10, h = 8)

if (!all(is.na(star$spatial_x_use))) {
  p2 <- ggplot(star@meta.data, aes(x = spatial_x_use, y = spatial_y_use,
                                   color = .data[[celltype_col]])) +
    geom_point(size = 0.3) + coord_fixed() + theme_minimal() +
    ggtitle("STARmap: Spatial Cell Types")
  save_plot(p2, file.path(DIR_FIGURES, "starmap_spatial_celltype.png"), w = 12, h = 8)
}

if ("week" %in% colnames(star@meta.data) && !all(is.na(star$week))) {
  p3 <- DimPlot(star, group.by = "week", label = TRUE) +
    ggtitle("STARmap UMAP by Gestational Week")
  save_plot(p3, file.path(DIR_FIGURES, "starmap_umap_week.png"), w = 10, h = 8)
}

log_msg("03B complete.", logfile)
print_seurat_diagnostic(star, "STARmap (mapped)")

cat("\n", strrep("=", 70), "\n")
cat(sprintf("STARmap cells: %d\n", ncol(star)))
cat(sprintf("Cell types: %d\n", length(unique(star@meta.data[[celltype_col]]))))
cat(sprintf("Output: %s\n", out_rds))
cat(strrep("=", 70), "\n\n")
