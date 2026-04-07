======================================================================
  # scripts/05_spatial/05E_spatial_lr_proximity.R
  # Spatial ligand-receptor proximity analysis (distance-constrained edges).
  #
  # Enhanced methods:
  # - K sensitivity panel (k = 8, 15, 25, 40)
  # - Week-stratified LR summaries
  # - Spatial-density constrained null model (receiver label shuffle)
  # ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

source("config/config.R")
source("scripts/R/utils.R")

check_required_packages(c("Seurat", "dplyr", "ggplot2", "tidyr", "tibble", "RANN"),
                        context = "05E_spatial_lr_proximity")

ensure_dir(DIR_TABLES); ensure_dir(DIR_FIGURES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "05E_spatial_lr_proximity.log")

K_SENS <- c(8, 15, 25, 40)
K_SENS <- unique(pmax(2, as.integer(K_SENS)))
LR_NULL_N_PERM <- 100
LR_MAX_PAIRS <- 150

slidetags_path <- if (file.exists(file.path(DIR_OBJS, "slidetags_harmonized.rds"))) {
  file.path(DIR_OBJS, "slidetags_harmonized.rds")
} else {
  file.path(DIR_OBJS, "slidetags_mapped_to_multiome.rds")
}
starmap_path <- if (file.exists(file.path(DIR_OBJS, "starmap_harmonized.rds"))) {
  file.path(DIR_OBJS, "starmap_harmonized.rds")
} else {
  file.path(DIR_OBJS, "starmap_mapped_to_multiome.rds")
}

pick_label_col <- function(md) {
  if ("celltype_final_refined" %in% colnames(md)) return("celltype_final_refined")
  if ("celltype_final_conservative" %in% colnames(md)) return("celltype_final_conservative")
  if (COL_PRED_CELLTYPE %in% colnames(md)) return(COL_PRED_CELLTYPE)
  if ("celltype_author" %in% colnames(md)) return("celltype_author")
  NULL
}

permute_receiver_labels <- function(df, n_perm = 100, seed = 1) {
  set.seed(seed)
  grp <- interaction(df$week, df$density_bin, drop = TRUE)
  splits <- split(seq_len(nrow(df)), grp)
  
  replicate(n_perm, {
    lab <- df$receiver_label
    for (ix in splits) {
      if (length(ix) > 1) lab[ix] <- sample(lab[ix], replace = FALSE)
    }
    lab
  }, simplify = FALSE)
}

score_lr_spatial <- function(obj, dataset_name, assay_use) {
  DefaultAssay(obj) <- assay_use
  obj <- safe_join_layers(obj, assay = assay_use)
  if (!has_data_layer(obj, assay = assay_use)) obj <- NormalizeData(obj, verbose = FALSE)
  
  md <- obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    mutate(week = .data[[COL_WEEK]],
           x = .data[["spatial_x_use"]],
           y = .data[["spatial_y_use"]])
  
  label_col <- pick_label_col(md)
  if (is.null(label_col)) stop("No label column found for ", dataset_name)
  md$label <- as.character(md[[label_col]])
  
  keep <- is.finite(md$x) & is.finite(md$y) & !is.na(md$label)
  md <- md[keep, , drop = FALSE]
  if (nrow(md) < 100) {
    log_msg(paste0("[05E] ", dataset_name, ": too few spatial cells; skipping"), logfile)
    return(invisible(NULL))
  }
  
  lr_pairs <- read.csv(PATH_LR_PAIRS, stringsAsFactors = FALSE)
  colnames(lr_pairs) <- tolower(colnames(lr_pairs))
  if (!all(c("ligand", "receptor") %in% colnames(lr_pairs))) stop("PATH_LR_PAIRS must include ligand,receptor")
  
  feats <- rownames(obj[[assay_use]])
  lr_pairs <- lr_pairs %>% filter(ligand %in% feats, receptor %in% feats)
  if (nrow(lr_pairs) == 0) {
    log_msg(paste0("[05E] ", dataset_name, ": no LR pairs in assay ", assay_use), logfile)
    return(invisible(NULL))
  }
  
  lr_pairs <- lr_pairs %>%
    mutate(pair = paste0(ligand, "->", receptor)) %>%
    distinct(pair, .keep_all = TRUE) %>%
    slice_head(n = LR_MAX_PAIRS)
  
  expr_lig <- FetchData(obj, vars = unique(lr_pairs$ligand)) %>% tibble::rownames_to_column("cell")
  expr_rec <- FetchData(obj, vars = unique(lr_pairs$receptor)) %>% tibble::rownames_to_column("cell")
  
  all_out <- list()
  
  for (k in K_SENS) {
    if (nrow(md) <= k + 1) next
    log_msg(paste0("[05E] ", dataset_name, ": k=", k), logfile)
    
    g <- compute_spatial_knn_graph(md$x, md$y, k = k)
    idx <- g$idx
    local_density <- 1 / (rowMeans(g$dist, na.rm = TRUE) + 1e-9)
    md$density_bin <- dplyr::ntile(local_density, 4)
    
    edge_df <- tibble::tibble(
      sender = rep(md$cell, each = ncol(idx)),
      receiver = as.vector(matrix(md$cell[idx], nrow = nrow(idx), ncol = ncol(idx)))
    ) %>%
      left_join(md %>% select(cell, week, density_bin, sender_label = label), by = c("sender" = "cell")) %>%
      left_join(md %>% select(cell, receiver_label = label), by = c("receiver" = "cell"))
    
    perm_receiver_labels <- permute_receiver_labels(edge_df, n_perm = LR_NULL_N_PERM, seed = SEED)
    
    per_pair <- lapply(seq_len(nrow(lr_pairs)), function(i) {
      L <- lr_pairs$ligand[i]
      R <- lr_pairs$receptor[i]
      lr_name <- lr_pairs$pair[i]
      
      d <- edge_df %>%
        left_join(expr_lig %>% select(cell, ligand_expr = .data[[L]]), by = c("sender" = "cell")) %>%
        left_join(expr_rec %>% select(cell, receptor_expr = .data[[R]]), by = c("receiver" = "cell")) %>%
        mutate(lr = lr_name,
               lr_score = log1p(pmax(0, ligand_expr)) * log1p(pmax(0, receptor_expr)))
      
      obs <- d %>%
        group_by(week, sender_label, receiver_label, lr) %>%
        summarize(mean_lr_score = mean(lr_score, na.rm = TRUE),
                  median_lr_score = median(lr_score, na.rm = TRUE),
                  n_edges = n(), .groups = "drop")
      
      if (nrow(obs) == 0) return(NULL)
      
      null_means <- lapply(perm_receiver_labels, function(lab_perm) {
        d$receiver_label_perm <- lab_perm
        d %>% group_by(week, sender_label, receiver_label_perm, lr) %>%
          summarize(null_mean = mean(lr_score, na.rm = TRUE), .groups = "drop") %>%
          rename(receiver_label = receiver_label_perm)
      })
      null_df <- bind_rows(null_means, .id = "perm_id")
      
      stats_df <- null_df %>%
        group_by(week, sender_label, receiver_label, lr) %>%
        summarize(null_mean = mean(null_mean, na.rm = TRUE),
                  null_sd = sd(null_mean, na.rm = TRUE),
                  .groups = "drop")
      
      out <- obs %>%
        left_join(stats_df, by = c("week", "sender_label", "receiver_label", "lr")) %>%
        mutate(z_null = (mean_lr_score - null_mean) / (null_sd + 1e-9),
               p_empirical = pnorm(-abs(z_null)) * 2,
               k = k)
      out
    })
    
    all_out[[length(all_out) + 1]] <- bind_rows(per_pair)
  }
  
  out <- bind_rows(all_out)
  if (nrow(out) == 0) return(invisible(NULL))
  out <- out %>% mutate(p_adj = p.adjust(p_empirical, method = "BH")) %>% arrange(p_adj, desc(mean_lr_score))
  
  write.csv(out, file.path(DIR_TABLES, paste0(dataset_name, "_spatial_lr_edges.csv")), row.names = FALSE)
  
  top <- out %>% group_by(k, week, sender_label, receiver_label) %>%
    slice_max(mean_lr_score, n = 3, with_ties = FALSE) %>% ungroup()
  write.csv(top, file.path(DIR_TABLES, paste0(dataset_name, "_spatial_lr_summary.csv")), row.names = FALSE)
  
  robust <- out %>%
    group_by(week, sender_label, receiver_label, lr) %>%
    summarize(n_k = n_distinct(k),
              n_k_sig = sum(p_adj < 0.05, na.rm = TRUE),
              median_z = median(z_null, na.rm = TRUE),
              robust_fraction = n_k_sig / pmax(1, n_k), .groups = "drop")
  write.csv(robust, file.path(DIR_TABLES, paste0(dataset_name, "_spatial_lr_k_robustness.csv")), row.names = FALSE)
  
  p <- robust %>%
    arrange(desc(robust_fraction), desc(median_z)) %>%
    slice_head(n = 20) %>%
    mutate(pair = paste0(sender_label, "->", receiver_label, " :: ", lr)) %>%
    ggplot(aes(x = reorder(pair, robust_fraction), y = robust_fraction, fill = median_z)) +
    geom_col() + coord_flip() + theme_classic() +
    labs(title = paste0(dataset_name, ": robust spatial LR (across k)"), x = "sender->receiver :: LR", y = "Robust fraction")
  save_plot(p, file.path(DIR_FIGURES, paste0(dataset_name, "_spatial_lr_top20.png")), w = 11, h = 8)
  
  log_msg(paste0("[05E] Done spatial LR for ", dataset_name, " using assay ", assay_use), logfile)
  invisible(out)
}

slide <- readRDS(slidetags_path)
star <- readRDS(starmap_path)
slide <- ensure_week_column(slide, COL_WEEK_CANDIDATES)
star <- ensure_week_column(star, COL_WEEK_CANDIDATES)
slide <- ensure_spatial_coords(slide, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)
star <- ensure_spatial_coords(star, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)

assay_slide <- if ("SCT" %in% names(slide@assays)) "SCT" else "RNA"
assay_star <- if (exists("select_starmap_assay")) select_starmap_assay(star, prefer_imputed = TRUE) else "RNA_raw"

score_lr_spatial(slide, "SlideTags", assay_slide)
score_lr_spatial(star, "STARmap", assay_star)

log_msg("[05E] Complete.", logfile)