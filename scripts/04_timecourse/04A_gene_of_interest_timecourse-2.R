# ======================================================================
# scripts/04_timecourse/04B_immune_subsets_refinement.R
# Immune subtyping (lightweight, score-based) for macrophage and NK cells.
#
# Goal:
# - Add immune-focused subtype labels WITHOUT deleting any author labels.
# - Use module-score heuristics so this runs even on STARmap panels.
#
# Outputs:
# - Updated RDS objects in output/objects/*_with_immune_subtypes.rds
# - Summary tables + quick figures
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("config/config.R")
source("scripts/R/utils.R")@@ -24,51 +24,
52 @@ source("scripts/R/utils.R")
ensure_dir(DIR_OBJS); ensure_dir(DIR_FIGURES); ensure_dir(DIR_TABLES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "04B_immune_subsets_refinement.log")

ref <- readRDS(file.path(DIR_OBJS, "multiome_reference_processed.rds"))
slide <- readRDS(if (file.exists(file.path(DIR_OBJS, "slidetags_harmonized.rds"))) file.path(DIR_OBJS, "slidetags_harmonized.rds") else file.path(DIR_OBJS, "slidetags_mapped_to_multiome.rds"))
star <- readRDS(if (file.exists(file.path(DIR_OBJS, "starmap_harmonized.rds"))) file.path(DIR_OBJS, "starmap_harmonized.rds") else file.path(DIR_OBJS, "starmap_mapped_to_multiome.rds"))

# Add broad celltype column used for filtering
ref$celltype_use <- ref$celltype_ref
slide$celltype_use <- slide@meta.data[[COL_PRED_CELLTYPE]] %||% slide$celltype_author
star$celltype_use <- star@meta.data[[COL_PRED_CELLTYPE]]

# Define immune-focused gene sets (reuse from config/gene_sets.R)
IMMUNE_SETS <- list(
  Cytotoxic_NK = GENESETS$Cytotoxic_NK,
  Immune_Tolerance = GENESETS$Immune_Tolerance,
  Myeloid_Inflammation = GENESETS$Myeloid_Inflammation,
  Interferon_Stimulated = GENESETS$Interferon_Stimulated
)

# Helper: add immune scores and assign subtype labels
refine_immune <- function(obj, dataset_name, assay_use) {
  DefaultAssay(obj) <- assay_use
  
  # Ensure normalized data exists
  if (!"data" %in% Layers(obj[[assay_use]])) obj <- NormalizeData(obj, verbose = FALSE)
  obj <- safe_join_layers(obj, assay = assay_use)
  if (!has_data_layer(obj, assay = assay_use)) obj <- NormalizeData(obj, verbose = FALSE)
  
  obj <- add_modules_from_list(obj, IMMUNE_SETS, assay = assay_use, prefix = "imm_", seed = SEED)
  
  md <- obj@meta.data
  
  # Broad immune detection via label text (works across datasets)
  ct <- as.character(md$celltype_use)
  is_nk  <- grepl("NK|dNK", ct, ignore.case = TRUE)
  is_mac <- grepl("mac|mono|myelo", ct, ignore.case = TRUE)
  
  md$immune_class <- ifelse(is_nk, "NK",
                            ifelse(is_mac, "Macrophage", NA_character_))
  
  # Score-based subtyping
  md$immune_subtype <- NA_character_
  
  # NK
  # (very simple: cytotoxic-high vs other)
  if ("imm_Cytotoxic_NK" %in% colnames(md)) {
    thr <- quantile(md$imm_Cytotoxic_NK[is_nk], probs = 0.75, na.rm = TRUE)
    md$immune_subtype[is_nk & md$imm_Cytotoxic_NK >= thr] <- "NK_cytotoxic_high"
    md$immune_subtype[is_nk & is.na(md$immune_subtype)] <- "NK_other"
  }
  
  # Macrophage
  @@ -82,36 +83,36 @@ refine_immune <- function(obj, dataset_name, assay_use) {
  }
  
  obj@meta.data <- md
  
  # Save summary
  tab <- md %>%
    filter(!is.na(immune_class)) %>%
    count(week, immune_class, immune_subtype, name = "n_cells") %>%
    arrange(week, immune_class, desc(n_cells))
  
  write.csv(tab, file.path(DIR_TABLES, paste0(dataset_name, "_immune_subtype_counts.csv")), row.names = FALSE)
  
  # Quick UMAP plot if available
  if ("umap" %in% names(obj@reductions)) {
    p <- DimPlot(obj, reduction = "umap", group.by = "immune_subtype", label = TRUE) +
      ggtitle(paste0(dataset_name, ": immune_subtype (UMAP)"))
    save_plot(p, file.path(DIR_FIGURES, paste0(dataset_name, "_immune_subtype_umap.png")), w = 10, h = 7)
  }
  
  obj
}

# Choose assays
assay_ref <- if ((ref@misc$norm_method_use %||% "LogNormalize") == "SCT") "SCT" else "RNA"
assay_slide <- if ("SCT" %in% names(slide@assays)) "SCT" else "RNA"
assay_star <- "RNA_raw"
assay_star <- if (exists("select_starmap_assay")) select_starmap_assay(star, prefer_imputed = TRUE) else "RNA_raw"

ref2 <- refine_immune(ref, "Multiome", assay_ref)
slide2 <- refine_immune(slide, "SlideTags", assay_slide)
star2 <- refine_immune(star, "STARmap", assay_star)

saveRDS(ref2, file.path(DIR_OBJS, "multiome_reference_with_immune_subtypes.rds"))
saveRDS(slide2, file.path(DIR_OBJS, "slidetags_with_immune_subtypes.rds"))
saveRDS(star2, file.path(DIR_OBJS, "starmap_with_immune_subtypes.rds"))

log_msg("Done immune subtype refinement.", logfile)