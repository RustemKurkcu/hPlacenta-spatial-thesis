<<<<<<< HEAD
# ======================================================================
# scripts/02_preprocess/02A_preprocess_multiome_reference.R
# Preprocess the 10x multiome RNA reference for label transfer.
#
# Output:
#   output/objects/multiome_reference_processed.rds
#   output/figures/multiome_*.png
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

source("config/config.R")
source("scripts/R/utils.R")

ensure_dir(DIR_OBJS); ensure_dir(DIR_FIGURES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "02A_preprocess_multiome_reference.log")

log_msg("Starting multiome reference preprocessing...", logfile)

if (!file.exists(PATH_MULTIOME_RDS)) {
  stop("Missing multiome RDS: ", PATH_MULTIOME_RDS,
       "\nPlease build the object first or update PATH_MULTIOME_RDS in config/config.R")
}

log_msg(paste0("Loading multiome reference: ", PATH_MULTIOME_RDS), logfile)
ref <- readRDS(PATH_MULTIOME_RDS)
print_seurat_diagnostic(ref, "Multiome (raw)")

# Ensure week column exists
ref <- ensure_week_column(ref, COL_WEEK_CANDIDATES)

# Join layers if split (Seurat v5)
ref <- safe_join_layers(ref, "RNA")

# Heuristic: are the values in the counts layer integer-ish?
ctype <- infer_counts_type(ref, assay = "RNA", layer = "counts")
log_msg(paste0("Reference count-type heuristic = ", ctype), logfile)

if (ctype == "UMI_like") {
  log_msg("Running SCTransform on reference (UMI-like counts detected).", logfile)
  ref <- SCTransform(ref, verbose = FALSE)
  ref@misc$norm_method_use <- "SCT"

  log_msg("Running PCA/UMAP + clustering on reference (SCT).", logfile)
  ref <- RunPCA(ref, verbose = FALSE)
  ref <- RunUMAP(ref, dims = DEFAULT_DIMS, return.model = TRUE, verbose = FALSE)
  ref <- FindNeighbors(ref, dims = DEFAULT_DIMS, verbose = FALSE)
  ref <- FindClusters(ref, resolution = 0.5, verbose = FALSE)

} else {
  log_msg("Counts do NOT look UMI-like. Using LogNormalize pipeline.", logfile)
  ref <- NormalizeData(ref, verbose = FALSE)
  ref <- FindVariableFeatures(ref, verbose = FALSE)
  ref <- ScaleData(ref, features = VariableFeatures(ref), verbose = FALSE)

  ref <- RunPCA(ref, verbose = FALSE)
  ref <- RunUMAP(ref, dims = DEFAULT_DIMS, return.model = TRUE, verbose = FALSE)
  ref <- FindNeighbors(ref, dims = DEFAULT_DIMS, verbose = FALSE)
  ref <- FindClusters(ref, resolution = 0.5, verbose = FALSE)

  ref@misc$norm_method_use <- "LogNormalize"
}

# Optional tSNE
if (isTRUE(DO_TSNE) && !("tsne" %in% names(ref@reductions))) {
  log_msg("Computing tSNE embedding for reference...", logfile)
  perplex <- min(30, max(5, floor((ncol(ref) - 1) / 3)))
  ref <- RunTSNE(ref, dims = DEFAULT_DIMS, reduction = "pca",
                 perplexity = perplex, verbose = FALSE)
  log_msg("  tSNE complete.", logfile)
}

# Preserve the provided annotation column as a stable reference taxonomy
ref$celltype_ref <- ref@meta.data[[COL_REF_CELLTYPE]]
ref@misc$ref_celltype_col <- COL_REF_CELLTYPE

log_msg(sprintf("Reference has %d cells across %d cell types",
                ncol(ref), length(unique(ref$celltype_ref))), logfile)

# Quick QC plots
log_msg("Generating QC plots...", logfile)

p1 <- VlnPlot(ref, features = c("nCount_RNA", "nFeature_RNA"), pt.size = 0.05) +
  ggtitle("Multiome QC")
save_plot(p1, file.path(DIR_FIGURES, "multiome_QC_violin.png"), w = 10, h = 4)

p2 <- DimPlot(ref, group.by = "celltype_ref", label = TRUE) +
  ggtitle("Multiome reference UMAP (celltype_ref)")
save_plot(p2, file.path(DIR_FIGURES, "multiome_reference_umap_celltype.png"), w = 10, h = 8)

if ("tsne" %in% names(ref@reductions)) {
  log_msg("Generating tSNE plot...", logfile)
  p3 <- DimPlot(ref, reduction = "tsne", group.by = "celltype_ref", label = TRUE) +
    ggtitle("Multiome reference tSNE (celltype_ref)")
  save_plot(p3, file.path(DIR_FIGURES, "multiome_reference_tsne_celltype.png"), w = 10, h = 8)
}

# Save processed reference
out_rds <- file.path(DIR_OBJS, "multiome_reference_processed.rds")
saveRDS(ref, out_rds)
log_msg(paste0("Saved processed reference: ", out_rds), logfile)
log_msg("02A complete.", logfile)

print_seurat_diagnostic(ref, "Multiome (processed)")

cat("\n", strrep("=", 70), "\n")
cat(sprintf("Normalization method: %s\n", ref@misc$norm_method_use))
cat(sprintf("Number of cells: %d\n", ncol(ref)))
cat(sprintf("Number of cell types: %d\n", length(unique(ref$celltype_ref))))
cat(sprintf("Reductions: %s\n", paste(names(ref@reductions), collapse = ", ")))
cat(sprintf("Output: %s\n", out_rds))
=======
# ======================================================================
# scripts/02_preprocess/02A_preprocess_multiome_reference.R
# Preprocess the 10x multiome RNA reference for label transfer.
#
# Output:
#   output/objects/multiome_reference_processed.rds
#   output/figures/multiome_*.png
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

source("config/config.R")
source("scripts/R/utils.R")

ensure_dir(DIR_OBJS); ensure_dir(DIR_FIGURES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "02A_preprocess_multiome_reference.log")

log_msg("Starting multiome reference preprocessing...", logfile)

if (!file.exists(PATH_MULTIOME_RDS)) {
  stop("Missing multiome RDS: ", PATH_MULTIOME_RDS,
       "\nPlease build the object first or update PATH_MULTIOME_RDS in config/config.R")
}

log_msg(paste0("Loading multiome reference: ", PATH_MULTIOME_RDS), logfile)
ref <- readRDS(PATH_MULTIOME_RDS)
print_seurat_diagnostic(ref, "Multiome (raw)")

# Ensure week column exists
ref <- ensure_week_column(ref, COL_WEEK_CANDIDATES)

# Join layers if split (Seurat v5)
ref <- safe_join_layers(ref, "RNA")

# Heuristic: are the values in the counts layer integer-ish?
ctype <- infer_counts_type(ref, assay = "RNA", layer = "counts")
log_msg(paste0("Reference count-type heuristic = ", ctype), logfile)

if (ctype == "UMI_like") {
  log_msg("Running SCTransform on reference (UMI-like counts detected).", logfile)
  ref <- SCTransform(ref, verbose = FALSE)
  ref@misc$norm_method_use <- "SCT"

  log_msg("Running PCA/UMAP + clustering on reference (SCT).", logfile)
  ref <- RunPCA(ref, verbose = FALSE)
  ref <- RunUMAP(ref, dims = DEFAULT_DIMS, return.model = TRUE, verbose = FALSE)
  ref <- FindNeighbors(ref, dims = DEFAULT_DIMS, verbose = FALSE)
  ref <- FindClusters(ref, resolution = 0.5, verbose = FALSE)

} else {
  log_msg("Counts do NOT look UMI-like. Using LogNormalize pipeline.", logfile)
  ref <- NormalizeData(ref, verbose = FALSE)
  ref <- FindVariableFeatures(ref, verbose = FALSE)
  ref <- ScaleData(ref, features = VariableFeatures(ref), verbose = FALSE)

  ref <- RunPCA(ref, verbose = FALSE)
  ref <- RunUMAP(ref, dims = DEFAULT_DIMS, return.model = TRUE, verbose = FALSE)
  ref <- FindNeighbors(ref, dims = DEFAULT_DIMS, verbose = FALSE)
  ref <- FindClusters(ref, resolution = 0.5, verbose = FALSE)

  ref@misc$norm_method_use <- "LogNormalize"
}

# Optional tSNE
if (isTRUE(DO_TSNE) && !("tsne" %in% names(ref@reductions))) {
  log_msg("Computing tSNE embedding for reference...", logfile)
  perplex <- min(30, max(5, floor((ncol(ref) - 1) / 3)))
  ref <- RunTSNE(ref, dims = DEFAULT_DIMS, reduction = "pca",
                 perplexity = perplex, verbose = FALSE)
  log_msg("  tSNE complete.", logfile)
}

# Preserve the provided annotation column as a stable reference taxonomy
ref$celltype_ref <- ref@meta.data[[COL_REF_CELLTYPE]]
ref@misc$ref_celltype_col <- COL_REF_CELLTYPE

log_msg(sprintf("Reference has %d cells across %d cell types",
                ncol(ref), length(unique(ref$celltype_ref))), logfile)

# Quick QC plots
log_msg("Generating QC plots...", logfile)

p1 <- VlnPlot(ref, features = c("nCount_RNA", "nFeature_RNA"), pt.size = 0.05) +
  ggtitle("Multiome QC")
save_plot(p1, file.path(DIR_FIGURES, "multiome_QC_violin.png"), w = 10, h = 4)

p2 <- DimPlot(ref, group.by = "celltype_ref", label = TRUE) +
  ggtitle("Multiome reference UMAP (celltype_ref)")
save_plot(p2, file.path(DIR_FIGURES, "multiome_reference_umap_celltype.png"), w = 10, h = 8)

if ("tsne" %in% names(ref@reductions)) {
  log_msg("Generating tSNE plot...", logfile)
  p3 <- DimPlot(ref, reduction = "tsne", group.by = "celltype_ref", label = TRUE) +
    ggtitle("Multiome reference tSNE (celltype_ref)")
  save_plot(p3, file.path(DIR_FIGURES, "multiome_reference_tsne_celltype.png"), w = 10, h = 8)
}

# Save processed reference
out_rds <- file.path(DIR_OBJS, "multiome_reference_processed.rds")
saveRDS(ref, out_rds)
log_msg(paste0("Saved processed reference: ", out_rds), logfile)
log_msg("02A complete.", logfile)

print_seurat_diagnostic(ref, "Multiome (processed)")

cat("\n", strrep("=", 70), "\n")
cat(sprintf("Normalization method: %s\n", ref@misc$norm_method_use))
cat(sprintf("Number of cells: %d\n", ncol(ref)))
cat(sprintf("Number of cell types: %d\n", length(unique(ref$celltype_ref))))
cat(sprintf("Reductions: %s\n", paste(names(ref@reductions), collapse = ", ")))
cat(sprintf("Output: %s\n", out_rds))
>>>>>>> 54c1a5c675ca223d183141f8a75de090237ce902
cat(strrep("=", 70), "\n\n")