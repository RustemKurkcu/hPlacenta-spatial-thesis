# ======================================================================
# scripts/08_metagenes/08B_metagene_module_discovery.R
# Metagene discovery: correlation -> clustering -> module score
#
# PURPOSE:
#   Discover gene modules (co-expression programs) from the reference data
#   using correlation-based clustering. This is a "WGCNA-lite" approach
#   that's faster and more interpretable for thesis/grant work.
#
# METHOD:
#   1. Calculate gene-gene correlation on variable features
#   2. Cluster genes into modules using hierarchical clustering
#   3. Generate module gene lists
#   4. Score modules in the reference object
#   5. Identify high-variance modules (biologically interesting)
#
# BIOLOGICAL RATIONALE:
#   Gene modules represent coordinated biological programs:
#   - ECM remodeling programs
#   - Immune tolerance programs
#   - Metabolic programs (hypoxia, ethanolamine)
#   - Cell cycle programs
#   These can be mapped spatially and temporally to identify
#   "permissiveness windows" for infection.
#
# OUTPUT:
#   - Module gene lists (CSV)
#   - Module scores in reference object
#   - Module variance plot (identifies important modules)
#   - Updated reference object with metagene scores
#
# USAGE:
#   source("scripts/08_metagenes/08B_metagene_module_discovery.R")
#
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
})

source("config/config.R")
source("scripts/R/utils.R")

ensure_dir(DIR_TABLES); ensure_dir(DIR_FIGURES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "08B_metagene_module_discovery.log")

# Create metagenes subdirectories
dir.create(file.path(DIR_TABLES, "metagenes"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(DIR_FIGURES, "metagenes"), showWarnings = FALSE, recursive = TRUE)

log_msg("Starting metagene module discovery...", logfile)

obj_path <- file.path(DIR_OBJS, "multiome_reference_processed.rds")
if (!file.exists(obj_path)) stop("Missing: ", obj_path, " (run 02A first)")

log_msg(paste0("Loading reference: ", obj_path), logfile)
obj <- readRDS(obj_path)

# Choose normalized layer:
# - If SCT exists, use SCT "data"
# - Else use RNA "data" after NormalizeData
assay_use <- if ("SCT" %in% names(obj@assays)) "SCT" else "RNA"
DefaultAssay(obj) <- assay_use

log_msg(paste0("Using assay: ", assay_use), logfile)

if (assay_use == "RNA" && !("data" %in% Layers(obj[["RNA"]]))) {
  log_msg("RNA data slot empty; running NormalizeData + FindVariableFeatures.", logfile)
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, verbose = FALSE)
}

# Use variable features
genes <- VariableFeatures(obj)
if (length(genes) < 500) {
  log_msg("Too few variable genes; increasing variable features.", logfile)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  genes <- VariableFeatures(obj)
}

# Keep it fast and stable (top 1500 variable genes)
genes <- head(genes, 1500)
log_msg(sprintf("Using %d variable genes for module discovery", length(genes)), logfile)

# Extract expression matrix
mat <- GetAssayData(obj, slot = "data")[genes, , drop=FALSE]
mat <- as.matrix(mat)

log_msg(sprintf("Matrix dimensions: %d genes x %d cells", nrow(mat), ncol(mat)), logfile)

# ---- Correlation between genes ----
log_msg(paste0("Computing gene-gene correlation for ", length(genes), " genes..."), logfile)
log_msg("  (This may take a few minutes)", logfile)

cormat <- cor(t(mat), method="pearson")
cormat[is.na(cormat)] <- 0

log_msg("  Correlation matrix computed", logfile)

# ---- Cluster genes into modules ----
k <- 20  # number of modules (tunable)
log_msg(sprintf("Clustering genes into %d modules...", k), logfile)

hc <- hclust(as.dist(1 - cormat), method="average")
module_id <- cutree(hc, k = k)

# Create module lists
modules <- split(names(module_id), module_id)

# Keep only non-tiny modules (>=10 genes)
modules <- modules[sapply(modules, length) >= 10]

log_msg(sprintf("Created %d modules with >=10 genes", length(modules)), logfile)

# ---- Save module gene lists ----
module_tbl <- bind_rows(lapply(names(modules), function(m) {
  tibble(module = paste0("MG", m), gene = modules[[m]])
}))

output_file <- file.path(DIR_TABLES, "metagenes", "metagene_modules_gene_lists.csv")
write.csv(module_tbl, output_file, row.names = FALSE)
log_msg(paste0("Saved module gene lists: ", output_file), logfile)

# ---- Score modules in the object ----
log_msg("Scoring modules in reference object...", logfile)

# Use safe module scoring if available
if (exists("add_module_score_safe")) {
  for (m in names(modules)) {
    score_name <- paste0("MG_", m)
    obj <- add_module_score_safe(
      obj, 
      genes = modules[[m]], 
      score_name = score_name,
      assay = assay_use,
      seed = SEED
    )
  }
} else {
  # Fallback to standard AddModuleScore
  obj <- AddModuleScore(
    obj,
    features = lapply(modules, unique),
    name = "MG_",
    assay = assay_use,
    search = FALSE
  )
}

# Get module score columns
mg_cols <- grep("^MG_", colnames(obj@meta.data), value = TRUE)
log_msg(sprintf("Added %d module scores to metadata", length(mg_cols)), logfile)

# Save module scores
write.csv(obj@meta.data[, mg_cols, drop=FALSE],
          file.path(DIR_TABLES, "metagenes", "metagene_module_scores_meta.csv"))

# ---- Identify high-variance modules ----
log_msg("Calculating module score variances...", logfile)

var_df <- tibble(
  module = mg_cols,
  variance = sapply(mg_cols, function(cn) var(obj@meta.data[[cn]], na.rm=TRUE))
) %>% 
  arrange(desc(variance))

write.csv(var_df, 
          file.path(DIR_TABLES, "metagenes", "metagene_module_score_variances.csv"), 
          row.names=FALSE)

log_msg("Top 5 high-variance modules:", logfile)
for (i in 1:min(5, nrow(var_df))) {
  log_msg(sprintf("  %d. %s (variance = %.4f)", 
                  i, var_df$module[i], var_df$variance[i]), logfile)
}

# ---- Plot: module score variance ----
p1 <- ggplot(var_df, aes(x=reorder(module, variance), y=variance)) +
  geom_col(fill="steelblue", alpha=0.8) +
  coord_flip() +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y = element_text(size=8),
    plot.title = element_text(face="bold")
  ) +
  labs(
    title="Metagene Modules: Score Variance (Reference)",
    subtitle="Higher variance = more biologically interesting",
    x="Module", 
    y="Variance"
  )

save_plot(p1, 
          file.path(DIR_FIGURES, "metagenes", "metagenes_reference_module_variance.png"), 
          w=8, h=10)

# ---- Save updated reference object ----
output_obj <- file.path(DIR_OBJS, "multiome_reference_processed_with_metagenes.rds")
saveRDS(obj, output_obj)
log_msg(paste0("Saved updated reference: ", output_obj), logfile)

log_msg("08B complete.", logfile)

cat("\n")
cat(strrep("=", 70), "\n")
cat("METAGENE MODULE DISCOVERY COMPLETE\n")
cat(strrep("=", 70), "\n")
cat(sprintf("Modules discovered: %d\n", length(modules)))
cat(sprintf("Module gene lists: %s\n", output_file))
cat(sprintf("Module scores: %s/metagenes/metagene_module_scores_meta.csv\n", DIR_TABLES))
cat(sprintf("Variance plot: %s/metagenes/metagenes_reference_module_variance.png\n", DIR_FIGURES))
cat(sprintf("Updated object: %s\n", output_obj))
cat(strrep("=", 70), "\n\n")

cat("NEXT STEPS:\n")
cat("1. Review high-variance modules in variance plot\n")
cat("2. Examine gene lists for top modules\n")
cat("3. Run 08C to map modules spatially and temporally\n")
cat("4. Interpret modules biologically (ECM, immune, metabolic programs)\n\n")

cat("INTERPRETATION GUIDE:\n")
cat("- High variance modules = biologically interesting programs\n")
cat("- Look for modules enriched in:\n")
cat("  * ECM remodeling genes (MMP2, MMP9, collagens)\n")
cat("  * Immune tolerance genes (HLA-G, IDO1, TGFB1)\n")
cat("  * Metabolic genes (PLD1, ethanolamine pathway)\n")
cat("  * Hypoxia response genes (HIF1A, VEGFA)\n\n")

# ======================================================================
# END OF SCRIPT
# ======================================================================