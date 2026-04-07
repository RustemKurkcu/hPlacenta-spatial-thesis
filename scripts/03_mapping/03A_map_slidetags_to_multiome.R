# ======================================================================
# scripts/03_mapping/03A_map_slidetags_to_multiome.R
# Map Slide-tags spatial RNA to the Multiome reference taxonomy.
#
# KEY FIXES:
#   1. Renames umap_1/umap_2 metadata columns BEFORE computing UMAP
#      (these clash with the reduction and cause FeaturePlot to crash)
#   2. JoinLayers before normalization (Seurat v5 split layers)
#   3. Detects active assay correctly (SCT vs RNA)
#
# Output:
#   output/objects/slidetags_mapped.rds
#   output/figures/slidetags_*.png
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

source("config/config.R")
source("scripts/R/utils.R")

ensure_dir(DIR_OBJS); ensure_dir(DIR_FIGURES); ensure_dir(DIR_TABLES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "03A_map_slidetags_to_multiome.log")

# ---- Load reference ----
ref_path <- file.path(DIR_OBJS, "multiome_reference_processed.rds")
if (!file.exists(ref_path)) stop("Missing processed reference: ", ref_path, " (run 02A first)")
ref <- readRDS(ref_path)

norm_method <- ref@misc$norm_method_use %||%
  (if ("SCT" %in% names(ref@assays)) "SCT" else "LogNormalize")
log_msg(paste0("Using normalization method for mapping = ", norm_method), logfile)

# ---- Load query ----
qry_path <- if (file.exists(PATH_SLIDETAGS_RDS)) PATH_SLIDETAGS_RDS else PATH_SLIDETAGS_RAW
log_msg(paste0("Loading Slide-tags query: ", qry_path), logfile)
qry <- readRDS(qry_path)

print_seurat_diagnostic(qry, "Slide-tags (as loaded)")

# ---- CRITICAL: Fix metadata columns that clash with reductions ----
# The slidetags object has umap_1, umap_2 as metadata columns from a
# previous analysis. These MUST be renamed before computing a new UMAP
# or Seurat's FeaturePlot/DimPlot will crash with:
#   "non-numeric argument to mathematical function"
qry <- fix_metadata_reduction_clash(qry)

# ---- Standard metadata enrichment ----
qry <- ensure_week_column(qry, COL_WEEK_CANDIDATES)
qry <- ensure_spatial_coords(qry, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)

# Author-provided labels (keep; do not overwrite)
author_col <- pick_first_present(qry@meta.data, COL_AUTHOR_CELLTYPE)
if (!is.null(author_col)) {
  qry$celltype_author <- as.character(qry@meta.data[[author_col]])
} else {
  qry$celltype_author <- NA_character_
}

# ---- Label transfer (skip if already done) ----
already_mapped <- COL_PRED_CELLTYPE %in% colnames(qry@meta.data)
log_msg(paste0("Slide-tags already has predicted labels? ", already_mapped), logfile)

if (!already_mapped) {
  # Switch to RNA assay and join layers
  DefaultAssay(qry) <- "RNA"
  qry <- safe_join_layers(qry, "RNA")

  if (norm_method == "SCT") {
    log_msg("Running SCTransform on query.", logfile)
    qry <- SCTransform(qry, verbose = FALSE)
    log_msg("Finding transfer anchors (SCT).", logfile)
    anchors <- FindTransferAnchors(
      reference = ref, query = qry,
      normalization.method = "SCT", dims = DEFAULT_DIMS)
  } else {
    log_msg("Running LogNormalize on query.", logfile)
    qry <- NormalizeData(qry, verbose = FALSE)
    qry <- FindVariableFeatures(qry, verbose = FALSE)
    qry <- ScaleData(qry, features = VariableFeatures(qry), verbose = FALSE)
    qry <- RunPCA(qry, verbose = FALSE)
    log_msg("Finding transfer anchors (LogNormalize).", logfile)
    anchors <- FindTransferAnchors(
      reference = ref, query = qry,
      normalization.method = "LogNormalize", dims = DEFAULT_DIMS)
  }

  log_msg("Transferring cell types.", logfile)
  pred <- TransferData(anchorset = anchors, refdata = ref$celltype_ref, dims = DEFAULT_DIMS)
  qry <- AddMetaData(qry, pred)
} else {
  log_msg("Skipping transfer (predicted.id already present).", logfile)
}

qry <- set_scalar_meta(qry, "mapping_norm_method", norm_method)

# Compute prediction.score.max if missing
if (!(COL_PRED_SCORE_MAX %in% colnames(qry@meta.data))) {
  score_cols <- grep("^prediction\\.score\\.", colnames(qry@meta.data), value = TRUE)
  if (length(score_cols) > 0) {
    qry@meta.data[[COL_PRED_SCORE_MAX]] <- apply(
      qry@meta.data[, score_cols, drop = FALSE], 1, max, na.rm = TRUE)
  }
}

# ---- Compute UMAP for visualization ----
# We always recompute because the original object had no reduction slots.
# The active assay may be SCT (from the original build). We need to pick
# the right assay and ensure its layers are joined.
if (!("umap" %in% names(qry@reductions)) || TRUE) {
  log_msg("Computing UMAP for visualization...", logfile)

  # Determine best assay to use
  if ("SCT" %in% names(qry@assays)) {
    assay_use <- "SCT"
  } else {
    assay_use <- "RNA"
  }
  log_msg(paste0("  Using assay: ", assay_use), logfile)

  # ensure_pca_umap handles JoinLayers, metadata clash, normalization
  qry <- ensure_pca_umap(qry, assay = assay_use, dims = DEFAULT_DIMS,
                         umap_return_model = FALSE, force = TRUE)
}

# Optional tSNE
if (isTRUE(DO_TSNE) && !("tsne" %in% names(qry@reductions))) {
  log_msg("Computing tSNE for query...", logfile)
  perplex <- min(30, max(5, floor((ncol(qry) - 1) / 3)))
  qry <- RunTSNE(qry, dims = DEFAULT_DIMS, reduction = "pca",
                 perplexity = perplex, verbose = FALSE)
}

# ---- Save ----
out_rds <- file.path(DIR_OBJS, "slidetags_mapped.rds")
saveRDS(qry, out_rds)
log_msg(paste0("Saved mapped Slide-tags: ", out_rds), logfile)

# ---- Generate plots ----
log_msg("Generating mapping QC plots...", logfile)

p1 <- DimPlot(qry, group.by = "predicted.id", label = TRUE, repel = TRUE) +
  ggtitle("Slide-tags: Predicted Cell Types")
save_plot(p1, file.path(DIR_FIGURES, "slidetags_umap_predicted.png"), w = 10, h = 8)

if (!all(is.na(qry$celltype_author))) {
  p2 <- DimPlot(qry, group.by = "celltype_author", label = TRUE, repel = TRUE) +
    ggtitle("Slide-tags: Author Cell Types")
  save_plot(p2, file.path(DIR_FIGURES, "slidetags_umap_author.png"), w = 10, h = 8)
}

if (COL_PRED_SCORE_MAX %in% colnames(qry@meta.data)) {
  # Explicitly specify reduction to avoid Seurat v5 confusion when metadata
  # contains legacy columns that resemble embedding names.
  p3 <- FeaturePlot(qry, features = COL_PRED_SCORE_MAX, reduction = "umap") +
    ggtitle("Slide-tags: Prediction Confidence")
  save_plot(p3, file.path(DIR_FIGURES, "slidetags_umap_confidence.png"), w = 10, h = 8)
}

# Spatial scatter if coordinates available
if (!all(is.na(qry$spatial_x_use))) {
  p4 <- ggplot(qry@meta.data, aes(x = spatial_x_use, y = spatial_y_use,
                                    color = predicted.id)) +
    geom_point(size = 0.5) + coord_fixed() + theme_minimal() +
    ggtitle("Slide-tags: Spatial Cell Types")
  save_plot(p4, file.path(DIR_FIGURES, "slidetags_spatial_celltype.png"), w = 10, h = 8)
}

log_msg("03A complete.", logfile)
print_seurat_diagnostic(qry, "Slide-tags (mapped)")

cat("\n", strrep("=", 70), "\n")
cat(sprintf("Mapped cells: %d\n", ncol(qry)))
cat(sprintf("Cell types: %d\n", length(unique(qry$predicted.id))))
cat(sprintf("Output: %s\n", out_rds))
cat(strrep("=", 70), "\n\n")