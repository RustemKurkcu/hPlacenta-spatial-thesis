<<<<<<< HEAD
# ======================================================================
# scripts/05_spatial/05D_neighborhood_DE.R
# ADVANCED MODULE: spatial niches + week-specific rewiring + DE.
#
# This extends 05B/05C:
#   * 05B gives pairwise adjacency enrichment (Z-scores)
#   * 05C gives a scalar “permissiveness” score
#
# Here we cluster cells by local neighborhood composition, inspired by the
# niche logic in spatial placenta atlases.
# Enhancements:
# - Niche clustering from neighborhood composition
# - Week-specific niche rewiring statistics
# - Marker-style DE and pseudo-bulk niche DE by donor/week
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

source("config/config.R")
source("scripts/R/utils.R")

check_required_packages(c("Seurat", "dplyr", "tibble", "ggplot2", "RANN", "Matrix"),
                        context = "05D_neighborhood_DE")

ensure_dir(DIR_FIGURES); ensure_dir(DIR_TABLES); ensure_dir(DIR_LOGS); ensure_dir(DIR_OBJS)
logfile <- file.path(DIR_LOGS, "05D_neighborhood_DE.log")

log_msg("[05D] Loading objects...", logfile)
slide_path <- if (file.exists(file.path(DIR_OBJS, "slidetags_harmonized.rds"))) file.path(DIR_OBJS, "slidetags_harmonized.rds") else file.path(DIR_OBJS, "slidetags_mapped_to_multiome.rds")
star_path <- if (file.exists(file.path(DIR_OBJS, "starmap_harmonized.rds"))) file.path(DIR_OBJS, "starmap_harmonized.rds") else file.path(DIR_OBJS, "starmap_mapped_to_multiome.rds")

objs <- list()
if (file.exists(slide_path)) objs$`Slide-tags` <- readRDS(slide_path)
if (file.exists(star_path)) objs$STARmap <- readRDS(star_path)
if (length(objs) == 0) stop("No spatial objects found. Run mapping first.")

K_NEIGHBORS <- 25
N_NEIGH_CLUSTERS <- 8
PSEUDOBULK_MAX_GENES <- 5000

pick_donor_col <- function(md) {
  cand <- c("donor_id", "sample_id", "orig.ident", "biosample_id", "donor", "patient_id")
  hit <- intersect(cand, colnames(md))
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

summarize_niche_celltype_context <- function(md, dataset, label_col) {
  if (!("niche_cluster" %in% colnames(md)) || !(label_col %in% colnames(md))) return(invisible(NULL))
  
  md2 <- md %>%
    mutate(celltype = as.character(.data[[label_col]]))
  
  comp_ct <- md2 %>%
    filter(!is.na(celltype), nzchar(celltype)) %>%
    count(niche_cluster, celltype, name = "n_cells") %>%
    group_by(niche_cluster) %>%
    mutate(frac = n_cells / sum(n_cells)) %>%
    ungroup() %>%
    arrange(niche_cluster, desc(frac), celltype)
  write.csv(comp_ct,
            file.path(DIR_TABLES, paste0(dataset, "_niche_center_celltype_composition.csv")),
            row.names = FALSE)
  
  top_ct <- comp_ct %>%
    group_by(niche_cluster) %>%
    slice_head(n = 5) %>%
    ungroup()
  write.csv(top_ct,
            file.path(DIR_TABLES, paste0(dataset, "_niche_center_celltype_top5.csv")),
            row.names = FALSE)
  
  p_ct <- comp_ct %>%
    group_by(niche_cluster) %>%
    slice_head(n = 8) %>%
    ungroup() %>%
    ggplot(aes(x = niche_cluster, y = frac, fill = celltype)) +
    geom_col(position = "stack") +
    theme_classic() +
    labs(title = paste0(dataset, ": niche cell-type composition (top 8 per niche)"),
         x = "Niche cluster", y = "Fraction")
  save_plot(p_ct, file.path(DIR_FIGURES, paste0(dataset, "_niche_center_celltype_composition.png")), w = 11, h = 6)
  
  if ("week" %in% colnames(md2)) {
    comp_ct_week <- md2 %>%
      filter(!is.na(celltype), nzchar(celltype), !is.na(week)) %>%
      count(week, niche_cluster, celltype, name = "n_cells") %>%
      group_by(week, niche_cluster) %>%
      mutate(frac = n_cells / sum(n_cells)) %>%
      ungroup()
    write.csv(comp_ct_week,
              file.path(DIR_TABLES, paste0(dataset, "_niche_week_center_celltype_composition.csv")),
              row.names = FALSE)
  }
  
  invisible(comp_ct)
}

summarize_niche_pathways <- function(obj, md, dataset, assay_use) {
  if (!("niche_cluster" %in% colnames(md))) return(invisible(NULL))
  if (!exists("GENESETS_CORE")) return(invisible(NULL))
  
  mat <- get_assay_matrix(obj, assay = assay_use, layer = "data")
  if (is.null(mat)) return(invisible(NULL))
  cells <- intersect(colnames(mat), md$cell)
  if (length(cells) < 50) return(invisible(NULL))
  mat <- mat[, cells, drop = FALSE]
  
  md_use <- md %>%
    filter(cell %in% cells)
  if (!("week" %in% colnames(md_use))) md_use$week <- NA
  
  genesets <- GENESETS_CORE
  gs_scores <- lapply(names(genesets), function(gs) {
    g <- intersect(genesets[[gs]], rownames(mat))
    if (length(g) == 0) {
      return(tibble::tibble(cell = colnames(mat), geneset = gs, score = NA_real_, n_genes = 0L))
    }
    sc <- Matrix::colMeans(mat[g, , drop = FALSE])
    tibble::tibble(cell = names(sc), geneset = gs, score = as.numeric(sc), n_genes = length(g))
  }) %>% bind_rows()
  
  gs_niche <- gs_scores %>%
    left_join(md_use %>% select(cell, niche_cluster, week), by = "cell") %>%
    group_by(niche_cluster, geneset) %>%
    summarize(mean_score = mean(score, na.rm = TRUE),
              median_score = median(score, na.rm = TRUE),
              n_cells = dplyr::n(),
              n_genes = max(n_genes, na.rm = TRUE),
              .groups = "drop") %>%
    group_by(geneset) %>%
    mutate(z_across_niches = (mean_score - mean(mean_score, na.rm = TRUE)) /
             (sd(mean_score, na.rm = TRUE) + 1e-9)) %>%
    ungroup()
  
  write.csv(gs_niche,
            file.path(DIR_TABLES, paste0(dataset, "_niche_geneset_scores.csv")),
            row.names = FALSE)
  
  p_gs <- gs_niche %>%
    ggplot(aes(x = geneset, y = niche_cluster, fill = z_across_niches)) +
    geom_tile() +
    scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick4", midpoint = 0) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0(dataset, ": niche pathway/module profile"),
         subtitle = "z-score of mean module score across niches",
         x = "Gene set", y = "Niche cluster", fill = "z")
  save_plot(p_gs, file.path(DIR_FIGURES, paste0(dataset, "_niche_geneset_heatmap.png")), w = 10, h = 6)
  
  if ("week" %in% colnames(md_use) && any(!is.na(md_use$week))) {
    gs_week <- gs_scores %>%
      left_join(md_use %>% select(cell, niche_cluster, week), by = "cell") %>%
      filter(!is.na(week)) %>%
      group_by(week, niche_cluster, geneset) %>%
      summarize(mean_score = mean(score, na.rm = TRUE), n_cells = dplyr::n(), .groups = "drop")
    write.csv(gs_week,
              file.path(DIR_TABLES, paste0(dataset, "_niche_week_geneset_scores.csv")),
              row.names = FALSE)
  }
  
  invisible(gs_niche)
}

compute_rewiring_stats <- function(md, dataset) {
  if (!("week" %in% colnames(md)) || !("niche_cluster" %in% colnames(md))) return(NULL)
  tab <- table(md$week, md$niche_cluster)
  tab <- tab[rowSums(tab) > 0, , drop = FALSE]
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NULL)
  
  chi <- suppressWarnings(chisq.test(tab))
  out_chi <- as.data.frame(as.table(chi$stdres))
  colnames(out_chi) <- c("week", "niche_cluster", "std_resid")
  out_chi$dataset <- dataset
  write.csv(out_chi, file.path(DIR_TABLES, paste0(dataset, "_niche_week_chisq_stdres.csv")), row.names = FALSE)
  
  comp <- prop.table(tab, margin = 1)
  comp_df <- as.data.frame(as.table(comp))
  colnames(comp_df) <- c("week", "niche_cluster", "fraction")
  comp_df$dataset <- dataset
  write.csv(comp_df, file.path(DIR_TABLES, paste0(dataset, "_niche_week_composition.csv")), row.names = FALSE)
  
  niche_trends <- comp_df %>%
    group_by(niche_cluster) %>%
    summarize(
      n_weeks = n_distinct(week),
      rho_week = if (n_distinct(week) >= 3) suppressWarnings(cor(as.numeric(as.character(week)), fraction, method = "spearman", use = "complete.obs")) else NA_real_,
      p_value = if (n_distinct(week) >= 3) suppressWarnings(cor.test(as.numeric(as.character(week)), fraction, method = "spearman", exact = FALSE)$p.value) else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"),
           trend = case_when(
             !is.na(p_adj) & p_adj < 0.05 & rho_week > 0 ~ "increasing",
             !is.na(p_adj) & p_adj < 0.05 & rho_week < 0 ~ "decreasing",
             TRUE ~ "not_significant"
           ))
  write.csv(niche_trends, file.path(DIR_TABLES, paste0(dataset, "_niche_fraction_trends.csv")), row.names = FALSE)
  
  p <- ggplot(comp_df, aes(x = as.factor(week), y = fraction, fill = niche_cluster)) +
    geom_col(position = "stack") +
    theme_classic() +
    labs(title = paste0(dataset, ": niche composition by week"), x = "Week", y = "Fraction")
  save_plot(p, file.path(DIR_FIGURES, paste0(dataset, "_niche_week_composition.png")), w = 10, h = 6)
  
  # Week-to-week rewiring (Jensen-Shannon divergence)
  jsd <- function(p, q) {
    p <- p / (sum(p) + 1e-12)
    q <- q / (sum(q) + 1e-12)
    m <- 0.5 * (p + q)
    kl <- function(a, b) sum(ifelse(a > 0, a * log2((a + 1e-12) / (b + 1e-12)), 0))
    0.5 * kl(p, m) + 0.5 * kl(q, m)
  }
  
  wk <- as.numeric(rownames(comp))
  ord <- order(wk)
  comp <- comp[ord, , drop = FALSE]
  wk <- wk[ord]
  if (length(wk) >= 2) {
    rw <- lapply(seq_len(length(wk) - 1), function(i) {
      tibble::tibble(dataset = dataset,
                     week_from = wk[i],
                     week_to = wk[i + 1],
                     jsd = jsd(comp[i, ], comp[i + 1, ]))
    }) %>% bind_rows()
    write.csv(rw, file.path(DIR_TABLES, paste0(dataset, "_niche_week_rewiring_jsd.csv")), row.names = FALSE)
  }
  
  list(chisq_p = chi$p.value)
}

run_pseudobulk_niche_de <- function(obj, dataset, assay_use) {
  md <- obj@meta.data %>% rownames_to_column("cell")
  donor_col <- pick_donor_col(md)
  if (is.null(donor_col) || !"week" %in% colnames(md)) {
    log_msg(paste0("[05D] pseudo-bulk skipped for ", dataset, ": donor/week columns unavailable"), logfile)
    return(NULL)
  }
  
  md <- md %>% mutate(donor = as.character(.data[[donor_col]]),
                      week = as.character(week),
                      replicate = paste(donor, week, sep = "__"))
  keep <- !is.na(md$replicate) & !is.na(md$niche_cluster)
  md <- md[keep, , drop = FALSE]
  if (nrow(md) < 200) return(NULL)
  
  mat <- get_assay_matrix(obj, assay = assay_use, layer = "data")
  if (is.null(mat)) return(NULL)
  
  common_cells <- intersect(colnames(mat), md$cell)
  if (length(common_cells) < 50) {
    log_msg(paste0("[05D] pseudo-bulk skipped for ", dataset, ": insufficient overlapping cells in assay matrix"), logfile)
    return(NULL)
  }
  md <- md %>% filter(cell %in% common_cells)
  md <- md[match(common_cells, md$cell), , drop = FALSE]
  mat <- mat[, common_cells, drop = FALSE]
  
  means <- Matrix::rowMeans(mat)
  genes_use <- names(sort(means, decreasing = TRUE))[seq_len(min(PSEUDOBULK_MAX_GENES, length(means)))]
  mat <- mat[genes_use, , drop = FALSE]
  
  grp <- paste(md$replicate, md$niche_cluster, sep = "@@")
  M <- Matrix::sparse.model.matrix(~ 0 + grp)
  colnames(M) <- sub("^grp", "", colnames(M))
  pb <- mat %*% M
  sizes <- Matrix::colSums(M)
  pb <- sweep(pb, 2, sizes, "/")
  
  pb_meta <- tibble::tibble(group = colnames(pb)) %>%
    tidyr::separate(group, into = c("replicate", "niche_cluster"), sep = "@@", remove = FALSE)
  
  out <- list()
  for (niche in sort(unique(pb_meta$niche_cluster))) {
    idx_in <- which(pb_meta$niche_cluster == niche)
    idx_out <- which(pb_meta$niche_cluster != niche)
    if (length(idx_in) < 2 || length(idx_out) < 2) next
    
    pvals <- apply(as.matrix(pb), 1, function(v) {
      suppressWarnings(stats::wilcox.test(v[idx_in], v[idx_out])$p.value)
    })
    lfc <- rowMeans(as.matrix(pb)[, idx_in, drop = FALSE]) - rowMeans(as.matrix(pb)[, idx_out, drop = FALSE])
    
    out[[length(out) + 1]] <- tibble::tibble(
      dataset = dataset,
      niche_cluster = niche,
      gene = rownames(pb),
      avg_diff = lfc,
      p_val = pvals,
      p_adj = p.adjust(pvals, method = "BH")
    )
  }
  
  de_pb <- bind_rows(out) %>% arrange(niche_cluster, p_adj, desc(avg_diff))
  if (nrow(de_pb) > 0) {
    write.csv(de_pb, file.path(DIR_TABLES, paste0(dataset, "_niche_pseudobulk_de.csv")), row.names = FALSE)
  }
  de_pb
}

compute_niche_labels <- function(obj, dataset) {
  obj <- ensure_week_column(obj, COL_WEEK_CANDIDATES)
  obj <- ensure_spatial_coords(obj, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)
  
  md <- obj@meta.data %>% rownames_to_column("cell")
  lab_col <- if ("celltype_final_refined" %in% colnames(md)) "celltype_final_refined" else COL_PRED_CELLTYPE
  if (!lab_col %in% colnames(md)) stop("Missing label column for neighborhood building.")
  if (!all(c("spatial_x_use", "spatial_y_use") %in% colnames(md))) stop("Missing spatial coords")
  
  g <- compute_spatial_knn_graph(md$spatial_x_use, md$spatial_y_use, k = K_NEIGHBORS)
  idx <- g$idx
  
  labs <- as.character(md[[lab_col]])
  uniq <- sort(unique(labs))
  comp <- matrix(0, nrow = length(labs), ncol = length(uniq), dimnames = list(md$cell, uniq))
  for (i in seq_len(nrow(idx))) {
    tab <- table(labs[idx[i, ]])
    comp[i, names(tab)] <- as.numeric(tab) / sum(tab)
  }
  
  # Cluster neighborhoods in composition space
  set.seed(SEED)
  km <- kmeans(comp, centers = min(N_NEIGH_CLUSTERS, nrow(comp)), nstart = 20)
  md$niche_cluster <- paste0("N", km$cluster)
  obj@meta.data <- md %>% column_to_rownames("cell")
  
  # Save table: niche composition
  comp_df <- as.data.frame(comp)
  comp_df$cell <- rownames(comp_df)
  comp_df$niche_cluster <- md$niche_cluster
  niche_comp <- comp_df %>% group_by(niche_cluster) %>% summarize(across(all_of(uniq), mean), n_cells = dplyr::n(), .groups = "drop")
  write.csv(niche_comp, file.path(DIR_TABLES, paste0(dataset, "_niche_composition.csv")), row.names = FALSE)
  
  # Plot: niche clusters in space
  summarize_niche_celltype_context(md, dataset, lab_col)
  
  p <- ggplot(md, aes(x = spatial_x_use, y = spatial_y_use, color = niche_cluster)) +
    geom_point(size = 0.5) + coord_equal() + theme_bw() +
    labs(title = paste0(dataset, " niche clusters (K=", K_NEIGHBORS, ")"), x = "x", y = "y")
  save_plot(p, file.path(DIR_FIGURES, paste0(dataset, "_niche_clusters_spatial.png")), w = 9, h = 7)
  
  assay_use <- if (dataset == "STARmap" && exists("select_starmap_assay")) select_starmap_assay(obj, prefer_imputed = TRUE) else Seurat::DefaultAssay(obj)
  Seurat::DefaultAssay(obj) <- assay_use
  obj <- safe_join_layers(obj, assay = assay_use)
  if (!has_data_layer(obj, assay = assay_use)) obj <- NormalizeData(obj, verbose = FALSE)
  
  summarize_niche_pathways(obj, obj@meta.data %>% rownames_to_column("cell"), dataset, assay_use)
  
  Idents(obj) <- obj$niche_cluster
  de <- tryCatch(FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE),
                 error = function(e) { log_msg(paste0("[05D] marker DE skipped for ", dataset, ": ", e$message), logfile); NULL })
  if (!is.null(de) && nrow(de) > 0) write.csv(de, file.path(DIR_TABLES, paste0(dataset, "_niche_markers.csv")), row.names = FALSE)
  
  compute_rewiring_stats(obj@meta.data %>% rownames_to_column("cell"), dataset)
  run_pseudobulk_niche_de(obj, dataset, assay_use)
  
  obj
}

if (!requireNamespace("RANN", quietly = TRUE)) {
  log_msg("[05D] Package 'RANN' not installed; skipping neighborhood module.", logfile)
} else {
  for (nm in names(objs)) {
    log_msg(paste0("[05D] Computing niches: ", nm), logfile)
    objs[[nm]] <- compute_niche_labels(objs[[nm]], nm)
  }
  
  if (!is.null(objs$`Slide-tags`)) saveRDS(objs$`Slide-tags`, slide_path)
  if (!is.null(objs$STARmap)) saveRDS(objs$STARmap, star_path)
  
  log_msg("[05D] Done.", logfile)
=======
# ======================================================================
# scripts/05_spatial/05D_neighborhood_DE.R
# ADVANCED MODULE: spatial niches + week-specific rewiring + DE.
#
# This extends 05B/05C:
#   * 05B gives pairwise adjacency enrichment (Z-scores)
#   * 05C gives a scalar “permissiveness” score
#
# Here we cluster cells by local neighborhood composition, inspired by the
# niche logic in spatial placenta atlases.
# Enhancements:
# - Niche clustering from neighborhood composition
# - Week-specific niche rewiring statistics
# - Marker-style DE and pseudo-bulk niche DE by donor/week
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

source("config/config.R")
source("scripts/R/utils.R")

check_required_packages(c("Seurat", "dplyr", "tibble", "ggplot2", "RANN", "Matrix"),
                        context = "05D_neighborhood_DE")

ensure_dir(DIR_FIGURES); ensure_dir(DIR_TABLES); ensure_dir(DIR_LOGS); ensure_dir(DIR_OBJS)
logfile <- file.path(DIR_LOGS, "05D_neighborhood_DE.log")

log_msg("[05D] Loading objects...", logfile)
slide_path <- if (file.exists(file.path(DIR_OBJS, "slidetags_harmonized.rds"))) file.path(DIR_OBJS, "slidetags_harmonized.rds") else file.path(DIR_OBJS, "slidetags_mapped_to_multiome.rds")
star_path <- if (file.exists(file.path(DIR_OBJS, "starmap_harmonized.rds"))) file.path(DIR_OBJS, "starmap_harmonized.rds") else file.path(DIR_OBJS, "starmap_mapped_to_multiome.rds")

objs <- list()
if (file.exists(slide_path)) objs$`Slide-tags` <- readRDS(slide_path)
if (file.exists(star_path)) objs$STARmap <- readRDS(star_path)
if (length(objs) == 0) stop("No spatial objects found. Run mapping first.")

K_NEIGHBORS <- 25
N_NEIGH_CLUSTERS <- 8
PSEUDOBULK_MAX_GENES <- 5000

pick_donor_col <- function(md) {
  cand <- c("donor_id", "sample_id", "orig.ident", "biosample_id", "donor", "patient_id")
  hit <- intersect(cand, colnames(md))
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

summarize_niche_celltype_context <- function(md, dataset, label_col) {
  if (!("niche_cluster" %in% colnames(md)) || !(label_col %in% colnames(md))) return(invisible(NULL))
  
  md2 <- md %>%
    mutate(celltype = as.character(.data[[label_col]]))
  
  comp_ct <- md2 %>%
    filter(!is.na(celltype), nzchar(celltype)) %>%
    count(niche_cluster, celltype, name = "n_cells") %>%
    group_by(niche_cluster) %>%
    mutate(frac = n_cells / sum(n_cells)) %>%
    ungroup() %>%
    arrange(niche_cluster, desc(frac), celltype)
  write.csv(comp_ct,
            file.path(DIR_TABLES, paste0(dataset, "_niche_center_celltype_composition.csv")),
            row.names = FALSE)
  
  top_ct <- comp_ct %>%
    group_by(niche_cluster) %>%
    slice_head(n = 5) %>%
    ungroup()
  write.csv(top_ct,
            file.path(DIR_TABLES, paste0(dataset, "_niche_center_celltype_top5.csv")),
            row.names = FALSE)
  
  p_ct <- comp_ct %>%
    group_by(niche_cluster) %>%
    slice_head(n = 8) %>%
    ungroup() %>%
    ggplot(aes(x = niche_cluster, y = frac, fill = celltype)) +
    geom_col(position = "stack") +
    theme_classic() +
    labs(title = paste0(dataset, ": niche cell-type composition (top 8 per niche)"),
         x = "Niche cluster", y = "Fraction")
  save_plot(p_ct, file.path(DIR_FIGURES, paste0(dataset, "_niche_center_celltype_composition.png")), w = 11, h = 6)
  
  if ("week" %in% colnames(md2)) {
    comp_ct_week <- md2 %>%
      filter(!is.na(celltype), nzchar(celltype), !is.na(week)) %>%
      count(week, niche_cluster, celltype, name = "n_cells") %>%
      group_by(week, niche_cluster) %>%
      mutate(frac = n_cells / sum(n_cells)) %>%
      ungroup()
    write.csv(comp_ct_week,
              file.path(DIR_TABLES, paste0(dataset, "_niche_week_center_celltype_composition.csv")),
              row.names = FALSE)
  }
  
  invisible(comp_ct)
}

summarize_niche_pathways <- function(obj, md, dataset, assay_use) {
  if (!("niche_cluster" %in% colnames(md))) return(invisible(NULL))
  if (!exists("GENESETS_CORE")) return(invisible(NULL))
  
  mat <- get_assay_matrix(obj, assay = assay_use, layer = "data")
  if (is.null(mat)) return(invisible(NULL))
  cells <- intersect(colnames(mat), md$cell)
  if (length(cells) < 50) return(invisible(NULL))
  mat <- mat[, cells, drop = FALSE]
  
  md_use <- md %>%
    filter(cell %in% cells)
  if (!("week" %in% colnames(md_use))) md_use$week <- NA
  
  genesets <- GENESETS_CORE
  gs_scores <- lapply(names(genesets), function(gs) {
    g <- intersect(genesets[[gs]], rownames(mat))
    if (length(g) == 0) {
      return(tibble::tibble(cell = colnames(mat), geneset = gs, score = NA_real_, n_genes = 0L))
    }
    sc <- Matrix::colMeans(mat[g, , drop = FALSE])
    tibble::tibble(cell = names(sc), geneset = gs, score = as.numeric(sc), n_genes = length(g))
  }) %>% bind_rows()
  
  gs_niche <- gs_scores %>%
    left_join(md_use %>% select(cell, niche_cluster, week), by = "cell") %>%
    group_by(niche_cluster, geneset) %>%
    summarize(mean_score = mean(score, na.rm = TRUE),
              median_score = median(score, na.rm = TRUE),
              n_cells = dplyr::n(),
              n_genes = max(n_genes, na.rm = TRUE),
              .groups = "drop") %>%
    group_by(geneset) %>%
    mutate(z_across_niches = (mean_score - mean(mean_score, na.rm = TRUE)) /
             (sd(mean_score, na.rm = TRUE) + 1e-9)) %>%
    ungroup()
  
  write.csv(gs_niche,
            file.path(DIR_TABLES, paste0(dataset, "_niche_geneset_scores.csv")),
            row.names = FALSE)
  
  p_gs <- gs_niche %>%
    ggplot(aes(x = geneset, y = niche_cluster, fill = z_across_niches)) +
    geom_tile() +
    scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick4", midpoint = 0) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0(dataset, ": niche pathway/module profile"),
         subtitle = "z-score of mean module score across niches",
         x = "Gene set", y = "Niche cluster", fill = "z")
  save_plot(p_gs, file.path(DIR_FIGURES, paste0(dataset, "_niche_geneset_heatmap.png")), w = 10, h = 6)
  
  if ("week" %in% colnames(md_use) && any(!is.na(md_use$week))) {
    gs_week <- gs_scores %>%
      left_join(md_use %>% select(cell, niche_cluster, week), by = "cell") %>%
      filter(!is.na(week)) %>%
      group_by(week, niche_cluster, geneset) %>%
      summarize(mean_score = mean(score, na.rm = TRUE), n_cells = dplyr::n(), .groups = "drop")
    write.csv(gs_week,
              file.path(DIR_TABLES, paste0(dataset, "_niche_week_geneset_scores.csv")),
              row.names = FALSE)
  }
  
  invisible(gs_niche)
}

compute_rewiring_stats <- function(md, dataset) {
  if (!("week" %in% colnames(md)) || !("niche_cluster" %in% colnames(md))) return(NULL)
  tab <- table(md$week, md$niche_cluster)
  tab <- tab[rowSums(tab) > 0, , drop = FALSE]
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NULL)
  
  chi <- suppressWarnings(chisq.test(tab))
  out_chi <- as.data.frame(as.table(chi$stdres))
  colnames(out_chi) <- c("week", "niche_cluster", "std_resid")
  out_chi$dataset <- dataset
  write.csv(out_chi, file.path(DIR_TABLES, paste0(dataset, "_niche_week_chisq_stdres.csv")), row.names = FALSE)
  
  comp <- prop.table(tab, margin = 1)
  comp_df <- as.data.frame(as.table(comp))
  colnames(comp_df) <- c("week", "niche_cluster", "fraction")
  comp_df$dataset <- dataset
  write.csv(comp_df, file.path(DIR_TABLES, paste0(dataset, "_niche_week_composition.csv")), row.names = FALSE)
  
  niche_trends <- comp_df %>%
    group_by(niche_cluster) %>%
    summarize(
      n_weeks = n_distinct(week),
      rho_week = if (n_distinct(week) >= 3) suppressWarnings(cor(as.numeric(as.character(week)), fraction, method = "spearman", use = "complete.obs")) else NA_real_,
      p_value = if (n_distinct(week) >= 3) suppressWarnings(cor.test(as.numeric(as.character(week)), fraction, method = "spearman", exact = FALSE)$p.value) else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"),
           trend = case_when(
             !is.na(p_adj) & p_adj < 0.05 & rho_week > 0 ~ "increasing",
             !is.na(p_adj) & p_adj < 0.05 & rho_week < 0 ~ "decreasing",
             TRUE ~ "not_significant"
           ))
  write.csv(niche_trends, file.path(DIR_TABLES, paste0(dataset, "_niche_fraction_trends.csv")), row.names = FALSE)
  
  p <- ggplot(comp_df, aes(x = as.factor(week), y = fraction, fill = niche_cluster)) +
    geom_col(position = "stack") +
    theme_classic() +
    labs(title = paste0(dataset, ": niche composition by week"), x = "Week", y = "Fraction")
  save_plot(p, file.path(DIR_FIGURES, paste0(dataset, "_niche_week_composition.png")), w = 10, h = 6)
  
  # Week-to-week rewiring (Jensen-Shannon divergence)
  jsd <- function(p, q) {
    p <- p / (sum(p) + 1e-12)
    q <- q / (sum(q) + 1e-12)
    m <- 0.5 * (p + q)
    kl <- function(a, b) sum(ifelse(a > 0, a * log2((a + 1e-12) / (b + 1e-12)), 0))
    0.5 * kl(p, m) + 0.5 * kl(q, m)
  }
  
  wk <- as.numeric(rownames(comp))
  ord <- order(wk)
  comp <- comp[ord, , drop = FALSE]
  wk <- wk[ord]
  if (length(wk) >= 2) {
    rw <- lapply(seq_len(length(wk) - 1), function(i) {
      tibble::tibble(dataset = dataset,
                     week_from = wk[i],
                     week_to = wk[i + 1],
                     jsd = jsd(comp[i, ], comp[i + 1, ]))
    }) %>% bind_rows()
    write.csv(rw, file.path(DIR_TABLES, paste0(dataset, "_niche_week_rewiring_jsd.csv")), row.names = FALSE)
  }
  
  list(chisq_p = chi$p.value)
}

run_pseudobulk_niche_de <- function(obj, dataset, assay_use) {
  md <- obj@meta.data %>% rownames_to_column("cell")
  donor_col <- pick_donor_col(md)
  if (is.null(donor_col) || !"week" %in% colnames(md)) {
    log_msg(paste0("[05D] pseudo-bulk skipped for ", dataset, ": donor/week columns unavailable"), logfile)
    return(NULL)
  }
  
  md <- md %>% mutate(donor = as.character(.data[[donor_col]]),
                      week = as.character(week),
                      replicate = paste(donor, week, sep = "__"))
  keep <- !is.na(md$replicate) & !is.na(md$niche_cluster)
  md <- md[keep, , drop = FALSE]
  if (nrow(md) < 200) return(NULL)
  
  mat <- get_assay_matrix(obj, assay = assay_use, layer = "data")
  if (is.null(mat)) return(NULL)
  
  common_cells <- intersect(colnames(mat), md$cell)
  if (length(common_cells) < 50) {
    log_msg(paste0("[05D] pseudo-bulk skipped for ", dataset, ": insufficient overlapping cells in assay matrix"), logfile)
    return(NULL)
  }
  md <- md %>% filter(cell %in% common_cells)
  md <- md[match(common_cells, md$cell), , drop = FALSE]
  mat <- mat[, common_cells, drop = FALSE]
  
  means <- Matrix::rowMeans(mat)
  genes_use <- names(sort(means, decreasing = TRUE))[seq_len(min(PSEUDOBULK_MAX_GENES, length(means)))]
  mat <- mat[genes_use, , drop = FALSE]
  
  grp <- paste(md$replicate, md$niche_cluster, sep = "@@")
  M <- Matrix::sparse.model.matrix(~ 0 + grp)
  colnames(M) <- sub("^grp", "", colnames(M))
  pb <- mat %*% M
  sizes <- Matrix::colSums(M)
  pb <- sweep(pb, 2, sizes, "/")
  
  pb_meta <- tibble::tibble(group = colnames(pb)) %>%
    tidyr::separate(group, into = c("replicate", "niche_cluster"), sep = "@@", remove = FALSE)
  
  out <- list()
  for (niche in sort(unique(pb_meta$niche_cluster))) {
    idx_in <- which(pb_meta$niche_cluster == niche)
    idx_out <- which(pb_meta$niche_cluster != niche)
    if (length(idx_in) < 2 || length(idx_out) < 2) next
    
    pvals <- apply(as.matrix(pb), 1, function(v) {
      suppressWarnings(stats::wilcox.test(v[idx_in], v[idx_out])$p.value)
    })
    lfc <- rowMeans(as.matrix(pb)[, idx_in, drop = FALSE]) - rowMeans(as.matrix(pb)[, idx_out, drop = FALSE])
    
    out[[length(out) + 1]] <- tibble::tibble(
      dataset = dataset,
      niche_cluster = niche,
      gene = rownames(pb),
      avg_diff = lfc,
      p_val = pvals,
      p_adj = p.adjust(pvals, method = "BH")
    )
  }
  
  de_pb <- bind_rows(out) %>% arrange(niche_cluster, p_adj, desc(avg_diff))
  if (nrow(de_pb) > 0) {
    write.csv(de_pb, file.path(DIR_TABLES, paste0(dataset, "_niche_pseudobulk_de.csv")), row.names = FALSE)
  }
  de_pb
}

compute_niche_labels <- function(obj, dataset) {
  obj <- ensure_week_column(obj, COL_WEEK_CANDIDATES)
  obj <- ensure_spatial_coords(obj, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)
  
  md <- obj@meta.data %>% rownames_to_column("cell")
  lab_col <- if ("celltype_final_refined" %in% colnames(md)) "celltype_final_refined" else COL_PRED_CELLTYPE
  if (!lab_col %in% colnames(md)) stop("Missing label column for neighborhood building.")
  if (!all(c("spatial_x_use", "spatial_y_use") %in% colnames(md))) stop("Missing spatial coords")
  
  g <- compute_spatial_knn_graph(md$spatial_x_use, md$spatial_y_use, k = K_NEIGHBORS)
  idx <- g$idx
  
  labs <- as.character(md[[lab_col]])
  uniq <- sort(unique(labs))
  comp <- matrix(0, nrow = length(labs), ncol = length(uniq), dimnames = list(md$cell, uniq))
  for (i in seq_len(nrow(idx))) {
    tab <- table(labs[idx[i, ]])
    comp[i, names(tab)] <- as.numeric(tab) / sum(tab)
  }
  
  # Cluster neighborhoods in composition space
  set.seed(SEED)
  km <- kmeans(comp, centers = min(N_NEIGH_CLUSTERS, nrow(comp)), nstart = 20)
  md$niche_cluster <- paste0("N", km$cluster)
  obj@meta.data <- md %>% column_to_rownames("cell")
  
  # Save table: niche composition
  comp_df <- as.data.frame(comp)
  comp_df$cell <- rownames(comp_df)
  comp_df$niche_cluster <- md$niche_cluster
  niche_comp <- comp_df %>% group_by(niche_cluster) %>% summarize(across(all_of(uniq), mean), n_cells = dplyr::n(), .groups = "drop")
  write.csv(niche_comp, file.path(DIR_TABLES, paste0(dataset, "_niche_composition.csv")), row.names = FALSE)
  
  # Plot: niche clusters in space
  summarize_niche_celltype_context(md, dataset, lab_col)
  
  p <- ggplot(md, aes(x = spatial_x_use, y = spatial_y_use, color = niche_cluster)) +
    geom_point(size = 0.5) + coord_equal() + theme_bw() +
    labs(title = paste0(dataset, " niche clusters (K=", K_NEIGHBORS, ")"), x = "x", y = "y")
  save_plot(p, file.path(DIR_FIGURES, paste0(dataset, "_niche_clusters_spatial.png")), w = 9, h = 7)
  
  assay_use <- if (dataset == "STARmap" && exists("select_starmap_assay")) select_starmap_assay(obj, prefer_imputed = TRUE) else Seurat::DefaultAssay(obj)
  Seurat::DefaultAssay(obj) <- assay_use
  obj <- safe_join_layers(obj, assay = assay_use)
  if (!has_data_layer(obj, assay = assay_use)) obj <- NormalizeData(obj, verbose = FALSE)
  
  summarize_niche_pathways(obj, obj@meta.data %>% rownames_to_column("cell"), dataset, assay_use)
  
  Idents(obj) <- obj$niche_cluster
  de <- tryCatch(FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE),
                 error = function(e) { log_msg(paste0("[05D] marker DE skipped for ", dataset, ": ", e$message), logfile); NULL })
  if (!is.null(de) && nrow(de) > 0) write.csv(de, file.path(DIR_TABLES, paste0(dataset, "_niche_markers.csv")), row.names = FALSE)
  
  compute_rewiring_stats(obj@meta.data %>% rownames_to_column("cell"), dataset)
  run_pseudobulk_niche_de(obj, dataset, assay_use)
  
  obj
}

if (!requireNamespace("RANN", quietly = TRUE)) {
  log_msg("[05D] Package 'RANN' not installed; skipping neighborhood module.", logfile)
} else {
  for (nm in names(objs)) {
    log_msg(paste0("[05D] Computing niches: ", nm), logfile)
    objs[[nm]] <- compute_niche_labels(objs[[nm]], nm)
  }
  
  if (!is.null(objs$`Slide-tags`)) saveRDS(objs$`Slide-tags`, slide_path)
  if (!is.null(objs$STARmap)) saveRDS(objs$STARmap, star_path)
  
  log_msg("[05D] Done.", logfile)
>>>>>>> 54c1a5c675ca223d183141f8a75de090237ce902
}