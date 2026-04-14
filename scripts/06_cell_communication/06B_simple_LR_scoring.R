<<<<<<< HEAD
# ======================================================================
# scripts/06_cell_communication/06B_simple_LR_scoring.R
# Lightweight ligand-receptor scoring (no heavy dependencies).
#
# For each ligand-receptor pair:
#   score(sender,receiver) = mean_expr(ligand in sender) * mean_expr(receptor in receiver)
#
# This is NOT a full probabilistic communication model, but it is useful
# for grant figures / hypothesis generation and runs fast.
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

source("config/config.R")
source("scripts/R/utils.R")

ensure_dir(DIR_FIGURES); ensure_dir(DIR_TABLES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "06B_simple_LR_scoring.log")

pairs <- readr::read_csv(PATH_LR_PAIRS, show_col_types = FALSE)

load_dataset <- function(name) {
  if (name == "Multiome") return(readRDS(file.path(DIR_OBJS, "multiome_reference_processed.rds")))
  if (name == "SlideTags") return(readRDS(file.path(DIR_OBJS, "slidetags_mapped_to_multiome.rds")))
  if (name == "STARmap") return(readRDS(file.path(DIR_OBJS, "starmap_mapped_to_multiome.rds")))
  stop("Unknown dataset: ", name)
}

choose_assay <- function(obj, dataset) {
  if (dataset == "STARmap") return("RNA_raw")
  if ("SCT" %in% names(obj@assays)) return("SCT")
  "RNA"
}

prep_obj <- function(obj, dataset) {
  obj <- ensure_week_column(obj, COL_WEEK_CANDIDATES)

  assay <- choose_assay(obj, dataset)
  DefaultAssay(obj) <- assay

  # Ensure normalized data exists
  if (!"data" %in% Layers(obj[[assay]])) {
    log_msg(paste0(dataset, ": running NormalizeData on assay ", assay), logfile)
    obj <- NormalizeData(obj, verbose = FALSE)
  }
  obj
}


pick_label_col <- function(md) {
  if ("celltype_final_refined" %in% colnames(md)) return("celltype_final_refined")
  if ("celltype_final_conservative" %in% colnames(md)) return("celltype_final_conservative")
  if (COL_PRED_CELLTYPE %in% colnames(md)) return(COL_PRED_CELLTYPE)
  if ("celltype_author" %in% colnames(md)) return("celltype_author")
  NULL
}

mean_by_celltype <- function(obj, assay, genes, celltype_col) {
  genes <- intersect(genes, rownames(obj))
  if (length(genes) == 0) return(NULL)

  DefaultAssay(obj) <- assay
  mat <- get_assay_matrix(obj, assay = assay, layer = "data")[genes, , drop = FALSE]
  ct <- as.character(obj@meta.data[[celltype_col]])

  # mean expression per cell type
  df <- as.data.frame(t(as.matrix(mat)))
  df$celltype <- ct
  out <- df %>%
    group_by(celltype) %>%
    summarise(across(all_of(genes), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  out
}

score_lr <- function(obj, dataset, celltype_col) {
  obj <- prep_obj(obj, dataset)
  assay <- choose_assay(obj, dataset)

  # pick top-level labels
  if (!(celltype_col %in% colnames(obj@meta.data))) stop("Missing celltype column: ", celltype_col)

  genes_needed <- unique(c(pairs$ligand, pairs$receptor))
  avg <- mean_by_celltype(obj, assay, genes_needed, celltype_col)
  if (is.null(avg)) return(NULL)

  avg_long <- avg %>%
    pivot_longer(-celltype, names_to = "gene", values_to = "mean_expr")

  # Precompute maps
  ligand_map <- avg_long %>% filter(gene %in% pairs$ligand) %>%
    rename(sender = celltype, ligand = gene, mean_ligand = mean_expr)
  receptor_map <- avg_long %>% filter(gene %in% pairs$receptor) %>%
    rename(receiver = celltype, receptor = gene, mean_receptor = mean_expr)

  # Cross-join for each pair (sender x receiver)
  out_list <- list()

  for (i in seq_len(nrow(pairs))) {
    p <- pairs[i,]
    lig <- p$ligand
    rec <- p$receptor

    L <- ligand_map %>% filter(ligand == lig)
    R <- receptor_map %>% filter(receptor == rec)

    if (nrow(L) == 0 || nrow(R) == 0) next

    edges <- tidyr::crossing(L, R) %>%
      mutate(
        pathway = p$pathway,
        note = p$note,
        dataset = dataset,
        score = mean_ligand * mean_receptor
      ) %>%
      select(dataset, pathway, ligand, receptor, sender, receiver, mean_ligand, mean_receptor, score, note)

    out_list[[length(out_list) + 1]] <- edges
  }

  bind_rows(out_list)
}

# Run all datasets
ref <- load_dataset("Multiome")
slide <- load_dataset("SlideTags")
star <- load_dataset("STARmap")

edges_ref <- score_lr(ref, "Multiome", celltype_col = "celltype_ref")
edges_slide <- score_lr(slide, "SlideTags", celltype_col = pick_label_col(slide@meta.data) %||% COL_PRED_CELLTYPE)
edges_star <- score_lr(star, "STARmap", celltype_col = pick_label_col(star@meta.data) %||% COL_PRED_CELLTYPE)

edges_all <- bind_rows(edges_ref, edges_slide, edges_star)
write.csv(edges_all, file.path(DIR_TABLES, "simple_LR_edges_all_datasets.csv"), row.names = FALSE)
log_msg(paste0("Saved LR edges: n=", nrow(edges_all)), logfile)

# Plot: top interactions per dataset/pathway
top_edges <- edges_all %>%
  group_by(dataset, pathway, ligand, receptor, sender, receiver) %>%
  summarise(score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  group_by(dataset, pathway) %>%
  slice_max(order_by = score, n = 25, with_ties = FALSE) %>%
  ungroup()

p <- top_edges %>%
  mutate(pair = paste0(ligand, "->", receptor),
         edge = paste0(sender, "->", receiver)) %>%
  ggplot(aes(x = edge, y = pair, fill = score)) +
  geom_tile() +
  facet_wrap(~dataset + pathway, scales = "free", ncol = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "Top ligand-receptor interaction scores (simple product metric)",
       x = "Sender->Receiver (cell types)", y = "Ligand->Receptor", fill = "score")
save_plot(p, file.path(DIR_FIGURES, "simple_LR_top_edges_heatmap.png"), w = 14, h = 12)

log_msg("Done simple LR scoring.", logfile)
=======
# ======================================================================
# scripts/06_cell_communication/06B_simple_LR_scoring.R
# Lightweight ligand-receptor scoring (no heavy dependencies).
#
# For each ligand-receptor pair:
#   score(sender,receiver) = mean_expr(ligand in sender) * mean_expr(receptor in receiver)
#
# This is NOT a full probabilistic communication model, but it is useful
# for grant figures / hypothesis generation and runs fast.
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

source("config/config.R")
source("scripts/R/utils.R")

ensure_dir(DIR_FIGURES); ensure_dir(DIR_TABLES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "06B_simple_LR_scoring.log")

pairs <- readr::read_csv(PATH_LR_PAIRS, show_col_types = FALSE)

load_dataset <- function(name) {
  if (name == "Multiome") return(readRDS(file.path(DIR_OBJS, "multiome_reference_processed.rds")))
  if (name == "SlideTags") return(readRDS(file.path(DIR_OBJS, "slidetags_mapped_to_multiome.rds")))
  if (name == "STARmap") return(readRDS(file.path(DIR_OBJS, "starmap_mapped_to_multiome.rds")))
  stop("Unknown dataset: ", name)
}

choose_assay <- function(obj, dataset) {
  if (dataset == "STARmap") return("RNA_raw")
  if ("SCT" %in% names(obj@assays)) return("SCT")
  "RNA"
}

prep_obj <- function(obj, dataset) {
  obj <- ensure_week_column(obj, COL_WEEK_CANDIDATES)

  assay <- choose_assay(obj, dataset)
  DefaultAssay(obj) <- assay

  # Ensure normalized data exists
  if (!"data" %in% Layers(obj[[assay]])) {
    log_msg(paste0(dataset, ": running NormalizeData on assay ", assay), logfile)
    obj <- NormalizeData(obj, verbose = FALSE)
  }
  obj
}


pick_label_col <- function(md) {
  if ("celltype_final_refined" %in% colnames(md)) return("celltype_final_refined")
  if ("celltype_final_conservative" %in% colnames(md)) return("celltype_final_conservative")
  if (COL_PRED_CELLTYPE %in% colnames(md)) return(COL_PRED_CELLTYPE)
  if ("celltype_author" %in% colnames(md)) return("celltype_author")
  NULL
}

mean_by_celltype <- function(obj, assay, genes, celltype_col) {
  genes <- intersect(genes, rownames(obj))
  if (length(genes) == 0) return(NULL)

  DefaultAssay(obj) <- assay
  mat <- get_assay_matrix(obj, assay = assay, layer = "data")[genes, , drop = FALSE]
  ct <- as.character(obj@meta.data[[celltype_col]])

  # mean expression per cell type
  df <- as.data.frame(t(as.matrix(mat)))
  df$celltype <- ct
  out <- df %>%
    group_by(celltype) %>%
    summarise(across(all_of(genes), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  out
}

score_lr <- function(obj, dataset, celltype_col) {
  obj <- prep_obj(obj, dataset)
  assay <- choose_assay(obj, dataset)

  # pick top-level labels
  if (!(celltype_col %in% colnames(obj@meta.data))) stop("Missing celltype column: ", celltype_col)

  genes_needed <- unique(c(pairs$ligand, pairs$receptor))
  avg <- mean_by_celltype(obj, assay, genes_needed, celltype_col)
  if (is.null(avg)) return(NULL)

  avg_long <- avg %>%
    pivot_longer(-celltype, names_to = "gene", values_to = "mean_expr")

  # Precompute maps
  ligand_map <- avg_long %>% filter(gene %in% pairs$ligand) %>%
    rename(sender = celltype, ligand = gene, mean_ligand = mean_expr)
  receptor_map <- avg_long %>% filter(gene %in% pairs$receptor) %>%
    rename(receiver = celltype, receptor = gene, mean_receptor = mean_expr)

  # Cross-join for each pair (sender x receiver)
  out_list <- list()

  for (i in seq_len(nrow(pairs))) {
    p <- pairs[i,]
    lig <- p$ligand
    rec <- p$receptor

    L <- ligand_map %>% filter(ligand == lig)
    R <- receptor_map %>% filter(receptor == rec)

    if (nrow(L) == 0 || nrow(R) == 0) next

    edges <- tidyr::crossing(L, R) %>%
      mutate(
        pathway = p$pathway,
        note = p$note,
        dataset = dataset,
        score = mean_ligand * mean_receptor
      ) %>%
      select(dataset, pathway, ligand, receptor, sender, receiver, mean_ligand, mean_receptor, score, note)

    out_list[[length(out_list) + 1]] <- edges
  }

  bind_rows(out_list)
}

# Run all datasets
ref <- load_dataset("Multiome")
slide <- load_dataset("SlideTags")
star <- load_dataset("STARmap")

edges_ref <- score_lr(ref, "Multiome", celltype_col = "celltype_ref")
edges_slide <- score_lr(slide, "SlideTags", celltype_col = pick_label_col(slide@meta.data) %||% COL_PRED_CELLTYPE)
edges_star <- score_lr(star, "STARmap", celltype_col = pick_label_col(star@meta.data) %||% COL_PRED_CELLTYPE)

edges_all <- bind_rows(edges_ref, edges_slide, edges_star)
write.csv(edges_all, file.path(DIR_TABLES, "simple_LR_edges_all_datasets.csv"), row.names = FALSE)
log_msg(paste0("Saved LR edges: n=", nrow(edges_all)), logfile)

# Plot: top interactions per dataset/pathway
top_edges <- edges_all %>%
  group_by(dataset, pathway, ligand, receptor, sender, receiver) %>%
  summarise(score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  group_by(dataset, pathway) %>%
  slice_max(order_by = score, n = 25, with_ties = FALSE) %>%
  ungroup()

p <- top_edges %>%
  mutate(pair = paste0(ligand, "->", receptor),
         edge = paste0(sender, "->", receiver)) %>%
  ggplot(aes(x = edge, y = pair, fill = score)) +
  geom_tile() +
  facet_wrap(~dataset + pathway, scales = "free", ncol = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "Top ligand-receptor interaction scores (simple product metric)",
       x = "Sender->Receiver (cell types)", y = "Ligand->Receptor", fill = "score")
save_plot(p, file.path(DIR_FIGURES, "simple_LR_top_edges_heatmap.png"), w = 14, h = 12)

log_msg("Done simple LR scoring.", logfile)
>>>>>>> 54c1a5c675ca223d183141f8a75de090237ce902
