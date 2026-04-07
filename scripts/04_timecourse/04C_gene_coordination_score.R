# ======================================================================
# scripts/04_timecourse/04C_gene_coordination_score.R
# OPTIONAL (but informative): "Gene coordination" score inspired by the
# Greenbaum timeline framework.
#
# Concept:
# - For each gene, compute a "center-of-mass" (COM) along gestational age,
#   weighting by expression.
# - For a pathway gene set, measure dispersion of COMs across its genes.
# - Compare to random gene sets of the same size -> coordination z-score.
#
# Output:
# - output/tables/gene_coordination_scores.csv
# - output/figures/gene_coordination_heatmap.png
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("config/config.R")
source("scripts/R/utils.R")

ensure_dir(DIR_TABLES); ensure_dir(DIR_FIGURES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "04C_gene_coordination_score.log")

ref <- readRDS(file.path(DIR_OBJS, "multiome_reference_processed.rds"))
ref <- ensure_week_column(ref, COL_WEEK_CANDIDATES)

assay_ref <- if ((ref@misc$norm_method_use %||% "LogNormalize") == "SCT") "SCT" else "RNA"
DefaultAssay(ref) <- assay_ref
ref <- safe_join_layers(ref, assay = assay_ref)
if (!has_data_layer(ref, assay = assay_ref)) ref <- NormalizeData(ref, verbose = FALSE)

# Gene universe: expressed genes (nonzero mean)
expr_mat <- get_assay_matrix(ref, assay = assay_ref, layer = "data")
gene_means <- Matrix::rowMeans(expr_mat)
genes_universe <- names(gene_means)[gene_means > 0]
log_msg(paste0("Gene universe (mean>0): ", length(genes_universe)), logfile)

# Gini coefficient (dispersion)
gini <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  x <- sort(x)
  n <- length(x)
  if (sum(x) == 0) return(0)
  G <- sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x))
  G
}

# Center-of-mass across week
compute_com <- function(expr_vec, week_vec) {
  # expr_vec: numeric expression for a gene across cells
  # week_vec: numeric week across cells
  w <- week_vec
  e <- expr_vec
  ok <- is.finite(w) & is.finite(e)
  w <- w[ok]; e <- e[ok]
  if (length(e) == 0 || sum(e) == 0) return(NA_real_)
  sum(w * e) / sum(e)
}

# Prepare cell-level week and (optionally) celltype grouping
md <- ref@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  mutate(week_num = as.numeric(.data[[COL_WEEK]]),
         celltype = as.character(celltype_ref))

# Optional: restrict to a subset of cell types for stability (edit as needed)
celltypes_keep <- unique(md$celltype)
# Example: prioritize trophoblast + immune + stromal categories if present
# celltypes_keep <- celltypes_keep[grepl("EVT|CTB|STB|fib|endo|mac|NK", celltypes_keep, ignore.case = TRUE)]

# Settings
N_NULL <- 200  # random gene sets per pathway/celltype (increase for final)
set.seed(SEED)

# Compute coordination per celltype (loop)
results <- list()

for (ct in celltypes_keep) {
  cells_ct <- md$cell[md$celltype == ct & is.finite(md$week_num)]
  if (length(cells_ct) < 200) next
  
  log_msg(paste0("Celltype ", ct, ": n_cells=", length(cells_ct)), logfile)
  
  weeks <- md$week_num[match(cells_ct, md$cell)]
  mat_ct <- expr_mat[, cells_ct, drop = FALSE]
  
  # Precompute COM for all genes (in universe) to avoid repeated work
  com_all <- sapply(genes_universe, function(g) compute_com(as.numeric(mat_ct[g, ]), weeks))
  names(com_all) <- genes_universe
  
  for (gs_name in names(GENESETS_CORE)) {
    geneset <- intersect(GENESETS_CORE[[gs_name]], genes_universe)
    if (length(geneset) < 5) next
    
    com_gs <- com_all[geneset]
    disp_gs <- gini(com_gs)
    
    # Null dispersion distribution
    disp_null <- numeric(N_NULL)
    for (b in seq_len(N_NULL)) {
      samp <- sample(genes_universe, size = length(geneset), replace = FALSE)
      disp_null[b] <- gini(com_all[samp])
    }
    
    z <- (mean(disp_null, na.rm = TRUE) - disp_gs) / (sd(disp_null, na.rm = TRUE) + 1e-9)
    # Interpretation: higher z => more coordinated (tighter COM dispersion) than random
    
    results[[length(results) + 1]] <- tibble::tibble(
      celltype = ct,
      geneset = gs_name,
      n_genes = length(geneset),
      disp_gini = disp_gs,
      null_mean = mean(disp_null, na.rm = TRUE),
      null_sd = sd(disp_null, na.rm = TRUE),
      coordination_z = z
    )
  }
}

res <- bind_rows(results)
write.csv(res, file.path(DIR_TABLES, "gene_coordination_scores.csv"), row.names = FALSE)

# Plot heatmap (celltype x geneset)
p <- res %>%
  ggplot(aes(x = geneset, y = celltype, fill = coordination_z)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene coordination z-score (COM dispersion vs random)",
       x = "Gene set", y = "Cell type", fill = "coordination z")
save_plot(p, file.path(DIR_FIGURES, "gene_coordination_heatmap.png"), w = 10, h = 10)

log_msg("Done gene coordination scoring.", logfile)
