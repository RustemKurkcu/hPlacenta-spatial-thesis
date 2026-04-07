# ======================================================================
# scripts/08_metagene/08B_metagene_nmf.R
#
# ADVANCED MODULE (optional): metagene discovery via NMF.
#
# This is an “opt-in” program discovery layer:
#   * If you want fast, stable results: use curated gene signatures
#     (config/config_01_gene_signatures.R).
#   * If you want discovery: use NMF or WGCNA to learn modules directly
#     from your data.
#
# The script will run only if the `NMF` package is installed.
# Output:
#   * output/tables/nmf_metagenes_{dataset}.csv
#   * output/figures/nmf_metagene_scores_{dataset}.png
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

source("config/config.R")
source("scripts/R/utils.R")

ensure_dir(DIR_TABLES); ensure_dir(DIR_FIGURES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "08B_metagene_nmf.log")

if (!requireNamespace("NMF", quietly = TRUE)) {
  log_msg("NMF package not installed; skipping 08B. Install with install.packages('NMF').", logfile)
  quit(save = "no")
}

log_msg("Starting metagene discovery (NMF)...", logfile)

paths <- list(
  Multiome = file.path(DIR_OBJS, "multiome_reference_processed.rds"),
  SlideTags = file.path(DIR_OBJS, "slidetags_mapped_to_multiome.rds"),
  STARmap = file.path(DIR_OBJS, "starmap_mapped_to_multiome.rds")
)

for (nm in names(paths)) {
  if (!file.exists(paths[[nm]])) {
    log_msg(paste0("Missing object for ", nm, "; skipping."), logfile)
    next
  }
  obj <- readRDS(paths[[nm]])
  assay_use <- select_best_assay_for_scoring(obj, prefer_imputed = TRUE)
  DefaultAssay(obj) <- assay_use
  obj <- safe_join_layers(obj, assay_use)
  if (!has_data_layer(obj, assay_use)) obj <- NormalizeData(obj, verbose = FALSE)

  # Select features: use VariableFeatures if set; else pick most variable genes.
  if (length(VariableFeatures(obj)) < 200) {
    obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  }
  feats <- VariableFeatures(obj)
  feats <- feats[feats %in% rownames(obj[[assay_use]])]
  feats <- head(feats, 2000)

  # Extract matrix (non-negative required for NMF).
  mat <- get_assay_matrix(obj, assay = assay_use, layer = "data")
  mat <- mat[feats, , drop = FALSE]
  mat <- as.matrix(mat)
  mat[mat < 0] <- 0

  # Choose rank (k). You can tune this; we use 8 as a reasonable default.
  k <- 8
  log_msg(paste0("Running NMF for ", nm, " (k=", k, ", genes=", nrow(mat), ", cells=", ncol(mat), ")"), logfile)

  set.seed(SEED)
  nmf_res <- NMF::nmf(mat, rank = k, nrun = 10, .opt = "v")

  W <- NMF::basis(nmf_res)      # genes x k
  H <- NMF::coef(nmf_res)       # k x cells

  # Save metagene definitions (top genes per factor)
  top_genes <- lapply(seq_len(k), function(i) {
    g <- names(sort(W[, i], decreasing = TRUE))
    head(g, 50)
  })
  names(top_genes) <- paste0("meta", seq_len(k))

  out_def <- dplyr::bind_rows(lapply(names(top_genes), function(mg) {
    data.frame(metagene = mg, rank = seq_along(top_genes[[mg]]), gene = top_genes[[mg]])
  }))
  write.csv(out_def, file.path(DIR_TABLES, paste0("nmf_metagenes_", nm, ".csv")), row.names = FALSE)

  # Add metagene scores to metadata (H)
  for (i in seq_len(k)) {
    obj@meta.data[[paste0("score_nmf_meta", i)]] <- as.numeric(H[i, colnames(obj)])
  }

  # Quick visualization (UMAP if exists)
  if ("umap" %in% names(obj@reductions)) {
    p <- FeaturePlot(obj, features = paste0("score_nmf_meta", seq_len(min(4, k))), reduction = "umap", ncol = 2)
    save_plot(p, file.path(DIR_FIGURES, paste0("nmf_metagene_scores_", nm, ".png")), w = 12, h = 10)
  }

  # Persist object (optional)
  saveRDS(obj, file.path(DIR_OBJS, paste0("nmf_scored_", nm, ".rds")))
}

log_msg("08B complete.", logfile)
