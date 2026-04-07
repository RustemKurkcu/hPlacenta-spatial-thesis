# ======================================================================
# scripts/05_spatial/05B_neighborhood_enrichment.R
# Spatial neighborhood enrichment (cell-type adjacency) for Slide-tags and STARmap.
#
# Method:
# - Build a KNN graph in (x,y) space.
# - Count directed edges between labels.
# - Compare to a label-permutation null to compute z-score enrichment.
# Enhanced methods:
# - KNN adjacency + permutation null (z-score)
# - K sensitivity panel (default k = 8, 15, 25, 40)
# - Week stratification
# - Spatial density stratification (quartiles)
# - Dual heatmap palettes (white->red and blue-white-red)
# - Effect-size and trend tables (observed-vs-expected and increasing/decreasing over time)
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("config/config.R")
source("scripts/R/utils.R")

check_required_packages(c("Seurat", "dplyr", "tidyr", "ggplot2", "tibble", "RANN"),
                        context = "05B_neighborhood_enrichment")

ensure_dir(DIR_FIGURES); ensure_dir(DIR_TABLES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "05B_neighborhood_enrichment.log")

# Prefer harmonized objects (if available)
K_SENS <- c(8, 15, 25, 40)
K_SENS <- unique(pmax(2, as.integer(K_SENS)))


# Backward-compatibility shim: some environments may have older utils.R
get_knn_graph <- function(x, y, k) {
  if (exists("compute_spatial_knn_graph")) {
    return(compute_spatial_knn_graph(x, y, k = k))
  }
  idx <- compute_spatial_knn(x, y, k = k)
  coords <- cbind(as.numeric(x), as.numeric(y))
  keep <- is.finite(coords[, 1]) & is.finite(coords[, 2])
  coords <- coords[keep, , drop = FALSE]
  k_use <- max(1L, min(as.integer(k), nrow(coords) - 1L))
  nn <- RANN::nn2(coords, k = k_use + 1L)
  list(idx = idx, dist = nn$nn.dists[, -1, drop = FALSE], k = k_use)
}

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

log_msg(paste0("Loading SlideTags object: ", slidetags_path), logfile)
slidetags <- readRDS(slidetags_path)
starmap   <- readRDS(starmap_path)
log_msg(paste0("Loading STARmap object: ", starmap_path), logfile)
starmap <- readRDS(starmap_path)

slidetags <- ensure_week_column(slidetags, COL_WEEK_CANDIDATES)
starmap   <- ensure_week_column(starmap, COL_WEEK_CANDIDATES)

starmap <- ensure_week_column(starmap, COL_WEEK_CANDIDATES)
slidetags <- ensure_spatial_coords(slidetags, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)
starmap   <- ensure_spatial_coords(starmap, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)
starmap <- ensure_spatial_coords(starmap, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)

pick_label_col <- function(md) {
  if ("celltype_final_refined" %in% colnames(md)) return("celltype_final_refined")
  if ("celltype_final_conservative" %in% colnames(md)) return("celltype_final_conservative")
  if (COL_PRED_CELLTYPE %in% colnames(md)) return(COL_PRED_CELLTYPE)
  if ("celltype_author" %in% colnames(md)) return("celltype_author")
  NULL
}

order_labels_biological <- function(labels) {
  labels <- unique(as.character(labels))
  groups <- dplyr::case_when(
    grepl("EVT|CTB|STB|Troph", labels, ignore.case = TRUE) ~ "01_trophoblast",
    grepl("NK|Tcell|T cell|Bcell|B cell|mac|mono|myelo|immune", labels, ignore.case = TRUE) ~ "02_immune",
    grepl("Endo|vascular|vessel", labels, ignore.case = TRUE) ~ "03_endothelial",
    grepl("Fib|strom|mesench", labels, ignore.case = TRUE) ~ "04_stromal",
    grepl("gland|epith", labels, ignore.case = TRUE) ~ "05_epithelial",
    TRUE ~ "99_other"
  )
  tibble::tibble(label = labels, group = groups) %>%
    arrange(group, label) %>%
    pull(label)
}

plot_neighbor_heatmap_versions <- function(zmat, dataset_name, week, k, out_prefix) {
  zdf <- as.data.frame(as.table(zmat))
  colnames(zdf) <- c("from", "to", "z")
  
  # Version 1: biologically ordered labels
  fixed_order <- order_labels_biological(c(as.character(zdf$from), as.character(zdf$to)))
  zdf_fixed <- zdf %>%
    mutate(from = factor(from, levels = fixed_order),
           to = factor(to, levels = fixed_order))
  
  p_fixed <- ggplot(zdf_fixed, aes(x = to, y = from, fill = z)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "firebrick4") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0(dataset_name, " neighbor enrichment z (week ", week, ", k=", k, ")"),
         subtitle = "Biologically ordered labels",
         x = "Neighbor label", y = "Center label", fill = "z")
  save_plot(p_fixed, paste0(out_prefix, "_neighbor_z_heatmap_fixed.png"), w = 12, h = 10)
  
  p_fixed_bwr <- ggplot(zdf_fixed, aes(x = to, y = from, fill = z)) +
    geom_tile() +
    scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick4", midpoint = 0) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0(dataset_name, " neighbor enrichment z (week ", week, ", k=", k, ")"),
         subtitle = "Biologically ordered labels (diverging palette)",
         x = "Neighbor label", y = "Center label", fill = "z")
  save_plot(p_fixed_bwr, paste0(out_prefix, "_neighbor_z_heatmap_fixed_bwr.png"), w = 12, h = 10)
  
  # Version 2: clustered labels
  row_ord <- rownames(zmat)
  col_ord <- colnames(zmat)
  if (nrow(zmat) >= 2) row_ord <- rownames(zmat)[hclust(dist(zmat))$order]
  if (ncol(zmat) >= 2) col_ord <- colnames(zmat)[hclust(dist(t(zmat)))$order]
  
  zdf_clustered <- zdf %>%
    mutate(from = factor(from, levels = row_ord),
           to = factor(to, levels = col_ord))
  
  p_clustered <- ggplot(zdf_clustered, aes(x = to, y = from, fill = z)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "firebrick4") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0(dataset_name, " neighbor enrichment z (week ", week, ", k=", k, ")"),
         subtitle = "Hierarchically clustered labels",
         x = "Neighbor label", y = "Center label", fill = "z")
  save_plot(p_clustered, paste0(out_prefix, "_neighbor_z_heatmap_clustered.png"), w = 12, h = 10)
  
  p_clustered_bwr <- ggplot(zdf_clustered, aes(x = to, y = from, fill = z)) +
    geom_tile() +
    scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick4", midpoint = 0) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0(dataset_name, " neighbor enrichment z (week ", week, ", k=", k, ")"),
         subtitle = "Hierarchically clustered labels (diverging palette)",
         x = "Neighbor label", y = "Center label", fill = "z")
  save_plot(p_clustered_bwr, paste0(out_prefix, "_neighbor_z_heatmap_clustered_bwr.png"), w = 12, h = 10)
}

enrichment_to_long <- function(enr, week, k, scope = "global", density_bin = NA_integer_) {
  zdf <- as.data.frame(as.table(enr$z))
  colnames(zdf) <- c("from", "to", "z")
  
  obsdf <- as.data.frame(as.table(enr$observed))
  expdf <- as.data.frame(as.table(enr$expected))
  sddf <- as.data.frame(as.table(enr$sd_null))
  colnames(obsdf) <- c("from", "to", "observed")
  colnames(expdf) <- c("from", "to", "expected")
  colnames(sddf) <- c("from", "to", "sd_null")
  
  zdf %>%
    left_join(obsdf, by = c("from", "to")) %>%
    left_join(expdf, by = c("from", "to")) %>%
    left_join(sddf, by = c("from", "to")) %>%
    mutate(
      week = week,
      k = k,
      scope = scope,
      density_bin = density_bin,
      obs_minus_exp = observed - expected,
      obs_over_exp = observed / (expected + 1e-9),
      log2_obs_over_exp = log2(obs_over_exp + 1e-9)
    )
}

run_one <- function(obj, dataset_name) {
  md0 <- obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    mutate(week = .data[[COL_WEEK]],
           x = .data[["spatial_x_use"]],
           y = .data[["spatial_y_use"]])
  
  label_col <- pick_label_col(md0)
  if (is.null(label_col)) stop("No label column found for ", dataset_name)
  
  md <- md0 %>%
    mutate(label = as.character(.data[[label_col]])) %>%
    filter(is.finite(x), is.finite(y), !is.na(label))
  
  weeks <- sort(unique(md$week[!is.na(md$week)]))
  if (length(weeks) == 0) weeks <- unique(md$week)
  
  all_rows <- list()
  for (w in weeks) {
    d <- md %>% filter(week == w)
    if (nrow(d) < 50) next
    
    log_msg(paste0(dataset_name, " week ", w, ": KNN k=", SPATIAL_KNN_K, " n=", nrow(d)), logfile)
    idx <- compute_spatial_knn(d$x, d$y, k = SPATIAL_KNN_K)
    
    enr <- neighbor_enrichment(d$label, idx, n_perm = SPATIAL_N_PERM, seed = SEED)
    
    # Save matrices
    write.csv(enr$observed, file.path(DIR_TABLES, paste0(dataset_name, "_week_", w, "_neighbor_observed.csv")))
    write.csv(enr$expected, file.path(DIR_TABLES, paste0(dataset_name, "_week_", w, "_neighbor_expected.csv")))
    write.csv(enr$z,        file.path(DIR_TABLES, paste0(dataset_name, "_week_", w, "_neighbor_z.csv")))
    
    # Plot z-score heatmap
    zdf <- as.data.frame(as.table(enr$z))
    colnames(zdf) <- c("from","to","z")
    p <- ggplot(zdf, aes(x = to, y = from, fill = z)) +
      geom_tile() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste0(dataset_name, " neighbor enrichment z (week ", w, ")"),
           x = "Neighbor label", y = "Center label", fill = "z")
    save_plot(p, file.path(DIR_FIGURES, paste0(dataset_name, "_week_", w, "_neighbor_z_heatmap.png")), w = 12, h = 10)
    log_msg(sprintf("%s week %s: n=%d", dataset_name, as.character(w), nrow(d)), logfile)
    
    for (k in K_SENS) {
      if (nrow(d) <= k + 1) next
      log_msg(sprintf("  k=%d", k), logfile)
      
      g <- get_knn_graph(d$x, d$y, k = k)
      idx <- g$idx
      local_density <- 1 / (rowMeans(g$dist, na.rm = TRUE) + 1e-9)
      d$density_bin <- dplyr::ntile(local_density, 4)
      
      # Global enrichment (all density bins)
      enr <- neighbor_enrichment(d$label, idx, n_perm = SPATIAL_N_PERM, seed = SEED)
      
      write.csv(enr$observed, file.path(DIR_TABLES, paste0(dataset_name, "_week_", w, "_k", k, "_neighbor_observed.csv")))
      write.csv(enr$expected, file.path(DIR_TABLES, paste0(dataset_name, "_week_", w, "_k", k, "_neighbor_expected.csv")))
      write.csv(enr$z, file.path(DIR_TABLES, paste0(dataset_name, "_week_", w, "_k", k, "_neighbor_z.csv")))
      
      all_rows[[length(all_rows) + 1]] <- enrichment_to_long(
        enr = enr, week = w, k = k, scope = "global", density_bin = NA_integer_
      )
      
      out_prefix <- file.path(DIR_FIGURES, paste0(dataset_name, "_week_", w, "_k", k))
      plot_neighbor_heatmap_versions(
        zmat = enr$z,
        dataset_name = dataset_name,
        week = w,
        k = k,
        out_prefix = out_prefix
      )
      
      # Density-stratified enrichment (recompute KNN within each density bin)
      for (b in sort(unique(d$density_bin))) {
        d_b <- d %>% filter(density_bin == b)
        if (nrow(d_b) < max(60, k + 2)) next
        
        idx_b <- compute_spatial_knn(d_b$x, d_b$y, k = min(k, nrow(d_b) - 1))
        labels_b <- d_b$label
        enr_b <- tryCatch(neighbor_enrichment(labels_b, idx_b,
                                              n_perm = max(50, floor(SPATIAL_N_PERM / 2)), seed = SEED),
                          error = function(e) NULL)
        if (is.null(enr_b)) next
        
        all_rows[[length(all_rows) + 1]] <- enrichment_to_long(
          enr = enr_b, week = w, k = k, scope = "density_stratified", density_bin = b
        )
      }
    }
  }
  
  if (length(all_rows) > 0) {
    edges <- dplyr::bind_rows(all_rows) %>%
      mutate(p_adj = stats::p.adjust(2 * pnorm(-abs(z)), method = "BH"),
             significant = p_adj < 0.05)
    write.csv(edges, file.path(DIR_TABLES, paste0(dataset_name, "_neighbor_enrichment_long.csv")), row.names = FALSE)
    
    robust <- edges %>%
      filter(scope == "global") %>%
      group_by(week, from, to) %>%
      summarize(n_k_tested = n_distinct(k),
                n_k_sig = sum(significant, na.rm = TRUE),
                median_z = median(z, na.rm = TRUE),
                robust_fraction = n_k_sig / pmax(1, n_k_tested),
                .groups = "drop")
    write.csv(robust, file.path(DIR_TABLES, paste0(dataset_name, "_neighbor_enrichment_k_robustness.csv")), row.names = FALSE)
    
    # Human-interpretable effect summary
    effect_summary <- edges %>%
      filter(scope == "global") %>%
      group_by(from, to) %>%
      summarize(
        n_obs = n(),
        median_z = median(z, na.rm = TRUE),
        median_obs_over_exp = median(obs_over_exp, na.rm = TRUE),
        median_obs_minus_exp = median(obs_minus_exp, na.rm = TRUE),
        frac_significant = mean(significant, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(frac_significant), desc(median_z))
    write.csv(effect_summary,
              file.path(DIR_TABLES, paste0(dataset_name, "_neighbor_enrichment_effect_summary.csv")),
              row.names = FALSE)
    
    # Week-trend analysis for interactions: increasing/decreasing over time
    trend_tbl <- edges %>%
      filter(scope == "global") %>%
      group_by(k, from, to) %>%
      summarize(
        n_weeks = n_distinct(week),
        rho_z = if (n_distinct(week) >= 3) suppressWarnings(cor(week, z, method = "spearman", use = "complete.obs")) else NA_real_,
        p_z = if (n_distinct(week) >= 3) suppressWarnings(cor.test(week, z, method = "spearman", exact = FALSE)$p.value) else NA_real_,
        rho_l2fc = if (n_distinct(week) >= 3) suppressWarnings(cor(week, log2_obs_over_exp, method = "spearman", use = "complete.obs")) else NA_real_,
        p_l2fc = if (n_distinct(week) >= 3) suppressWarnings(cor.test(week, log2_obs_over_exp, method = "spearman", exact = FALSE)$p.value) else NA_real_,
        median_obs_over_exp = median(obs_over_exp, na.rm = TRUE),
        median_z = median(z, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        p_adj_z = p.adjust(p_z, method = "BH"),
        p_adj_l2fc = p.adjust(p_l2fc, method = "BH"),
        trend = case_when(
          !is.na(p_adj_z) & p_adj_z < 0.05 & rho_z > 0 ~ "increasing",
          !is.na(p_adj_z) & p_adj_z < 0.05 & rho_z < 0 ~ "decreasing",
          TRUE ~ "not_significant"
        )
      )
    
    write.csv(trend_tbl,
              file.path(DIR_TABLES, paste0(dataset_name, "_neighbor_enrichment_trends.csv")),
              row.names = FALSE)
    
    shortlist <- trend_tbl %>%
      filter(trend %in% c("increasing", "decreasing"), median_obs_over_exp > 1) %>%
      arrange(trend, p_adj_z, desc(abs(rho_z)), desc(median_obs_over_exp))
    write.csv(shortlist,
              file.path(DIR_TABLES, paste0(dataset_name, "_neighbor_enrichment_trend_shortlist.csv")),
              row.names = FALSE)
    
    # Human-friendly trend figure (top changing interactions)
    n_top <- min(20L, sum(trend_tbl$trend %in% c("increasing", "decreasing"), na.rm = TRUE))
    top_plot_tbl <- trend_tbl %>%
      filter(trend %in% c("increasing", "decreasing")) %>%
      arrange(p_adj_z, desc(abs(rho_z))) %>%
      slice_head(n = n_top) %>%
      mutate(pair = paste(from, to, sep = " -> "))
    if (nrow(top_plot_tbl) > 0) {
      p_trend <- ggplot(top_plot_tbl, aes(x = reorder(pair, rho_z), y = rho_z, fill = trend)) +
        geom_col() +
        coord_flip() +
        theme_classic() +
        labs(title = paste0(dataset_name, ": interactions changing over time"),
             subtitle = "Spearman trend of neighborhood z-score across weeks",
             x = "Center -> Neighbor", y = "rho(week, z)")
      save_plot(p_trend,
                file.path(DIR_FIGURES, paste0(dataset_name, "_neighbor_enrichment_trend_top.png")),
                w = 12, h = 8)
    }
  }
}

run_one(slidetags, "SlideTags")
run_one(starmap, "STARmap")

log_msg("Done neighborhood enrichment.", logfile)
log_msg("Done neighborhood enrichment (with k-sensitivity and density stratification).", logfile)
