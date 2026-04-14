<<<<<<< HEAD
# ======================================================================
# scripts/05_spatial/05C_permissiveness_score_maps.R
# Compute and visualize a composite "permissiveness" score:
#   permissiveness = z(MMP/ECM remodeling) + z(Immune tolerance) - z(NK cytotoxic)
#                   + z(Ethanolamine metabolism)
#
# Rationale:
# - Mirrors the Greenbaum concept of a temporally gated, locally regulated
#   permissive milieu (immune tolerance + remodeling around EVT/arteries),
#   and operationalizes the project's "Spatiotemporal Window" hypothesis.
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

source("config/config.R")
source("scripts/R/utils.R")

check_required_packages(c("Seurat", "dplyr", "ggplot2", "tibble"), context = "05C_permissiveness_score_maps")

ensure_dir(DIR_FIGURES); ensure_dir(DIR_TABLES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "05C_permissiveness_score_maps.log")

# Prefer harmonized objects
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

slide <- readRDS(slidetags_path)
star  <- readRDS(starmap_path)

slide <- ensure_week_column(slide, COL_WEEK_CANDIDATES)
star  <- ensure_week_column(star, COL_WEEK_CANDIDATES)

slide <- ensure_spatial_coords(slide, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)
star  <- ensure_spatial_coords(star, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)

pick_label_col <- function(md) {
  if ("celltype_final_refined" %in% colnames(md)) return("celltype_final_refined")
  if ("celltype_final_conservative" %in% colnames(md)) return("celltype_final_conservative")
  if (COL_PRED_CELLTYPE %in% colnames(md)) return(COL_PRED_CELLTYPE)
  if ("celltype_author" %in% colnames(md)) return("celltype_author")
  NULL
}

compute_perm <- function(obj, dataset_name, assay_use) {
  DefaultAssay(obj) <- assay_use
  obj <- safe_join_layers(obj, assay = assay_use)
  if (!has_data_layer(obj, assay = assay_use)) obj <- NormalizeData(obj, verbose = FALSE)
  
  # Ensure required module scores exist
  obj <- add_module_score_safe(obj, GENESETS$MMP_ECM_Remodeling, "score_MMP_ECM_Remodeling", assay = assay_use, seed = SEED)
  obj <- add_module_score_safe(obj, GENESETS$Immune_Tolerance, "score_Immune_Tolerance", assay = assay_use, seed = SEED)
  obj <- add_module_score_safe(obj, GENESETS$Cytotoxic_NK, "score_Cytotoxic_NK", assay = assay_use, seed = SEED)
  
  nk_candidates <- list(
    Cytotoxic_NK = GENESETS$Cytotoxic_NK,
    NK_Expanded = if ("NK_Expanded" %in% names(GENESETS)) GENESETS$NK_Expanded else NULL,
    dNK_Residency = if ("dNK_Residency" %in% names(GENESETS)) GENESETS$dNK_Residency else NULL
  )
  nk_candidates <- nk_candidates[!vapply(nk_candidates, is.null, logical(1))]
  
  nk_info <- NULL
  nk_selected <- NULL
  for (nk_nm in names(nk_candidates)) {
    obj <- add_module_score_safe(
      obj,
      nk_candidates[[nk_nm]],
      "score_Cytotoxic_NK",
      assay = assay_use,
      seed = SEED,
      min_genes = 2
    )
    info <- attr(obj, "last_module_score_info")
    if (!is.null(info) && identical(info$score_name, "score_Cytotoxic_NK")) {
      nk_info <- info
      if (isTRUE(info$found_n >= 2)) {
        nk_selected <- nk_nm
        break
      }
    }
  }
  if (is.null(nk_selected)) {
    nk_selected <- "none"
    warning(sprintf("%s: insufficient NK coverage across all NK modules; NK contribution will default to 0 in z_safe.", dataset_name))
  }
  log_msg(sprintf("[05C] %s NK module used: %s (%d genes found)",
                  dataset_name, nk_selected, nk_info$found_n %||% 0L), logfile)
  
  obj <- add_module_score_safe(obj, GENESETS$Ethanolamine_Metabolism, "score_Ethanolamine_Metabolism", assay = assay_use, seed = SEED)
  
  md <- obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    mutate(
      week = .data[[COL_WEEK]],
      x = .data[["spatial_x_use"]],
      y = .data[["spatial_y_use"]]
    )
  
  label_col <- pick_label_col(md) %||% "celltype_use"
  md$label <- md[[label_col]]
  
  # z-score within dataset (robust to sparse/all-NA modules)
  z_safe <- function(v, nm = "module") {
    v <- suppressWarnings(as.numeric(v))
    ok <- is.finite(v)
    if (sum(ok) < 2) {
      warning(sprintf("%s: <2 finite values for '%s'; using 0 contribution for permissiveness.", dataset_name, nm))
      return(rep(0, length(v)))
    }
    mu <- mean(v[ok], na.rm = TRUE)
    sdv <- stats::sd(v[ok], na.rm = TRUE)
    if (!is.finite(sdv) || sdv < 1e-12) {
      warning(sprintf("%s: near-zero SD for '%s'; using 0 contribution for permissiveness.", dataset_name, nm))
      return(rep(0, length(v)))
    }
    z <- (v - mu) / sdv
    z[!is.finite(z)] <- 0
    z
  }
  
  md$z_tol <- z_safe(md$score_Immune_Tolerance, "score_Immune_Tolerance")
  md$z_nk  <- z_safe(md$score_Cytotoxic_NK, "score_Cytotoxic_NK")
  md$z_ea  <- z_safe(md$score_Ethanolamine_Metabolism, "score_Ethanolamine_Metabolism")
  
  md$permissiveness <- md$z_mmp + md$z_tol - md$z_nk + md$z_ea
  
  nk_qc <- data.frame(
    dataset = dataset_name,
    assay = assay_use,
    nk_gene_set_used = nk_selected,
    nk_requested_genes = nk_info$requested_n %||% NA_integer_,
    nk_found_genes = nk_info$found_n %||% NA_integer_,
    nk_used_aliases = nk_info$used_aliases %||% FALSE,
    nk_found_features = paste(nk_info$found_features %||% character(0), collapse = ";"),
    nk_missing_genes = paste(nk_info$missing_genes %||% character(0), collapse = ";"),
    stringsAsFactors = FALSE
  )
  write.csv(nk_qc, file.path(DIR_TABLES, paste0(dataset_name, "_NK_module_coverage_qc.csv")), row.names = FALSE)
  
  # Save cell-level table (lightweight)
  out_csv <- file.path(DIR_TABLES, paste0(dataset_name, "_permissiveness_cell_level.csv"))
  write.csv(md %>% select(cell, week, label, permissiveness, z_mmp, z_tol, z_nk, z_ea, x, y),
            out_csv, row.names = FALSE)
  
  # Plot spatial map per week (continuous)
  weeks <- sort(unique(md$week[!is.na(md$week)]))
  if (length(weeks) == 0) weeks <- unique(md$week)
  
  for (w in weeks) {
    d <- md %>% filter(week == w) %>% filter(is.finite(x), is.finite(y))
    if (nrow(d) == 0) next
    p <- ggplot(d, aes(x = x, y = y, color = permissiveness)) +
      geom_point(size = 0.6, alpha = 0.9) +
      coord_fixed() +
      theme_classic() +
      labs(title = paste0(dataset_name, " permissiveness score (week ", w, ")"),
           x = "spatial x", y = "spatial y", color = "permissiveness")
    save_plot(p, file.path(DIR_FIGURES, paste0(dataset_name, "_permissiveness_week_", w, ".png")), w = 10, h = 8)
  }
  
  # Composition of top permissive cells
  thr <- quantile(md$permissiveness, probs = 0.90, na.rm = TRUE)
  comp <- md %>%
    mutate(is_top = permissiveness >= thr) %>%
    count(is_top, label, name = "n_cells") %>%
    group_by(is_top) %>%
    mutate(frac = n_cells / sum(n_cells)) %>%
    ungroup()
  write.csv(comp, file.path(DIR_TABLES, paste0(dataset_name, "_permissiveness_top10pct_composition.csv")), row.names = FALSE)
  
  p2 <- comp %>% filter(is_top) %>%
    ggplot(aes(x = reorder(label, -frac), y = frac)) +
    geom_col() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0(dataset_name, ": composition of top 10% permissive cells"),
         x = "Label", y = "Fraction of top cells")
  save_plot(p2, file.path(DIR_FIGURES, paste0(dataset_name, "_permissiveness_top10pct_composition.png")), w = 10, h = 6)
  
  obj@meta.data$permissiveness <- md$permissiveness[match(rownames(obj@meta.data), md$cell)]
  obj
}

assay_slide <- if ("SCT" %in% names(slide@assays)) "SCT" else "RNA"
assay_star  <- if (exists("select_starmap_assay")) select_starmap_assay(star, prefer_imputed = TRUE) else "RNA_raw"
log_msg(paste0("05C assays -> SlideTags: ", assay_slide, ", STARmap: ", assay_star), logfile)

slide2 <- compute_perm(slide, "SlideTags", assay_slide)
star2  <- compute_perm(star,  "STARmap", assay_star)

saveRDS(slide2, file.path(DIR_OBJS, "slidetags_with_permissiveness.rds"))
saveRDS(star2,  file.path(DIR_OBJS, "starmap_with_permissiveness.rds"))

=======
# ======================================================================
# scripts/05_spatial/05C_permissiveness_score_maps.R
# Compute and visualize a composite "permissiveness" score:
#   permissiveness = z(MMP/ECM remodeling) + z(Immune tolerance) - z(NK cytotoxic)
#                   + z(Ethanolamine metabolism)
#
# Rationale:
# - Mirrors the Greenbaum concept of a temporally gated, locally regulated
#   permissive milieu (immune tolerance + remodeling around EVT/arteries),
#   and operationalizes the project's "Spatiotemporal Window" hypothesis.
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

source("config/config.R")
source("scripts/R/utils.R")

check_required_packages(c("Seurat", "dplyr", "ggplot2", "tibble"), context = "05C_permissiveness_score_maps")

ensure_dir(DIR_FIGURES); ensure_dir(DIR_TABLES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "05C_permissiveness_score_maps.log")

# Prefer harmonized objects
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

slide <- readRDS(slidetags_path)
star  <- readRDS(starmap_path)

slide <- ensure_week_column(slide, COL_WEEK_CANDIDATES)
star  <- ensure_week_column(star, COL_WEEK_CANDIDATES)

slide <- ensure_spatial_coords(slide, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)
star  <- ensure_spatial_coords(star, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)

pick_label_col <- function(md) {
  if ("celltype_final_refined" %in% colnames(md)) return("celltype_final_refined")
  if ("celltype_final_conservative" %in% colnames(md)) return("celltype_final_conservative")
  if (COL_PRED_CELLTYPE %in% colnames(md)) return(COL_PRED_CELLTYPE)
  if ("celltype_author" %in% colnames(md)) return("celltype_author")
  NULL
}

compute_perm <- function(obj, dataset_name, assay_use) {
  DefaultAssay(obj) <- assay_use
  obj <- safe_join_layers(obj, assay = assay_use)
  if (!has_data_layer(obj, assay = assay_use)) obj <- NormalizeData(obj, verbose = FALSE)
  
  # Ensure required module scores exist
  obj <- add_module_score_safe(obj, GENESETS$MMP_ECM_Remodeling, "score_MMP_ECM_Remodeling", assay = assay_use, seed = SEED)
  obj <- add_module_score_safe(obj, GENESETS$Immune_Tolerance, "score_Immune_Tolerance", assay = assay_use, seed = SEED)
  obj <- add_module_score_safe(obj, GENESETS$Cytotoxic_NK, "score_Cytotoxic_NK", assay = assay_use, seed = SEED)
  
  nk_candidates <- list(
    Cytotoxic_NK = GENESETS$Cytotoxic_NK,
    NK_Expanded = if ("NK_Expanded" %in% names(GENESETS)) GENESETS$NK_Expanded else NULL,
    dNK_Residency = if ("dNK_Residency" %in% names(GENESETS)) GENESETS$dNK_Residency else NULL
  )
  nk_candidates <- nk_candidates[!vapply(nk_candidates, is.null, logical(1))]
  
  nk_info <- NULL
  nk_selected <- NULL
  for (nk_nm in names(nk_candidates)) {
    obj <- add_module_score_safe(
      obj,
      nk_candidates[[nk_nm]],
      "score_Cytotoxic_NK",
      assay = assay_use,
      seed = SEED,
      min_genes = 2
    )
    info <- attr(obj, "last_module_score_info")
    if (!is.null(info) && identical(info$score_name, "score_Cytotoxic_NK")) {
      nk_info <- info
      if (isTRUE(info$found_n >= 2)) {
        nk_selected <- nk_nm
        break
      }
    }
  }
  if (is.null(nk_selected)) {
    nk_selected <- "none"
    warning(sprintf("%s: insufficient NK coverage across all NK modules; NK contribution will default to 0 in z_safe.", dataset_name))
  }
  log_msg(sprintf("[05C] %s NK module used: %s (%d genes found)",
                  dataset_name, nk_selected, nk_info$found_n %||% 0L), logfile)
  
  obj <- add_module_score_safe(obj, GENESETS$Ethanolamine_Metabolism, "score_Ethanolamine_Metabolism", assay = assay_use, seed = SEED)
  
  md <- obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    mutate(
      week = .data[[COL_WEEK]],
      x = .data[["spatial_x_use"]],
      y = .data[["spatial_y_use"]]
    )
  
  label_col <- pick_label_col(md) %||% "celltype_use"
  md$label <- md[[label_col]]
  
  # z-score within dataset (robust to sparse/all-NA modules)
  z_safe <- function(v, nm = "module") {
    v <- suppressWarnings(as.numeric(v))
    ok <- is.finite(v)
    if (sum(ok) < 2) {
      warning(sprintf("%s: <2 finite values for '%s'; using 0 contribution for permissiveness.", dataset_name, nm))
      return(rep(0, length(v)))
    }
    mu <- mean(v[ok], na.rm = TRUE)
    sdv <- stats::sd(v[ok], na.rm = TRUE)
    if (!is.finite(sdv) || sdv < 1e-12) {
      warning(sprintf("%s: near-zero SD for '%s'; using 0 contribution for permissiveness.", dataset_name, nm))
      return(rep(0, length(v)))
    }
    z <- (v - mu) / sdv
    z[!is.finite(z)] <- 0
    z
  }
  
  md$z_tol <- z_safe(md$score_Immune_Tolerance, "score_Immune_Tolerance")
  md$z_nk  <- z_safe(md$score_Cytotoxic_NK, "score_Cytotoxic_NK")
  md$z_ea  <- z_safe(md$score_Ethanolamine_Metabolism, "score_Ethanolamine_Metabolism")
  
  md$permissiveness <- md$z_mmp + md$z_tol - md$z_nk + md$z_ea
  
  nk_qc <- data.frame(
    dataset = dataset_name,
    assay = assay_use,
    nk_gene_set_used = nk_selected,
    nk_requested_genes = nk_info$requested_n %||% NA_integer_,
    nk_found_genes = nk_info$found_n %||% NA_integer_,
    nk_used_aliases = nk_info$used_aliases %||% FALSE,
    nk_found_features = paste(nk_info$found_features %||% character(0), collapse = ";"),
    nk_missing_genes = paste(nk_info$missing_genes %||% character(0), collapse = ";"),
    stringsAsFactors = FALSE
  )
  write.csv(nk_qc, file.path(DIR_TABLES, paste0(dataset_name, "_NK_module_coverage_qc.csv")), row.names = FALSE)
  
  # Save cell-level table (lightweight)
  out_csv <- file.path(DIR_TABLES, paste0(dataset_name, "_permissiveness_cell_level.csv"))
  write.csv(md %>% select(cell, week, label, permissiveness, z_mmp, z_tol, z_nk, z_ea, x, y),
            out_csv, row.names = FALSE)
  
  # Plot spatial map per week (continuous)
  weeks <- sort(unique(md$week[!is.na(md$week)]))
  if (length(weeks) == 0) weeks <- unique(md$week)
  
  for (w in weeks) {
    d <- md %>% filter(week == w) %>% filter(is.finite(x), is.finite(y))
    if (nrow(d) == 0) next
    p <- ggplot(d, aes(x = x, y = y, color = permissiveness)) +
      geom_point(size = 0.6, alpha = 0.9) +
      coord_fixed() +
      theme_classic() +
      labs(title = paste0(dataset_name, " permissiveness score (week ", w, ")"),
           x = "spatial x", y = "spatial y", color = "permissiveness")
    save_plot(p, file.path(DIR_FIGURES, paste0(dataset_name, "_permissiveness_week_", w, ".png")), w = 10, h = 8)
  }
  
  # Composition of top permissive cells
  thr <- quantile(md$permissiveness, probs = 0.90, na.rm = TRUE)
  comp <- md %>%
    mutate(is_top = permissiveness >= thr) %>%
    count(is_top, label, name = "n_cells") %>%
    group_by(is_top) %>%
    mutate(frac = n_cells / sum(n_cells)) %>%
    ungroup()
  write.csv(comp, file.path(DIR_TABLES, paste0(dataset_name, "_permissiveness_top10pct_composition.csv")), row.names = FALSE)
  
  p2 <- comp %>% filter(is_top) %>%
    ggplot(aes(x = reorder(label, -frac), y = frac)) +
    geom_col() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0(dataset_name, ": composition of top 10% permissive cells"),
         x = "Label", y = "Fraction of top cells")
  save_plot(p2, file.path(DIR_FIGURES, paste0(dataset_name, "_permissiveness_top10pct_composition.png")), w = 10, h = 6)
  
  obj@meta.data$permissiveness <- md$permissiveness[match(rownames(obj@meta.data), md$cell)]
  obj
}

assay_slide <- if ("SCT" %in% names(slide@assays)) "SCT" else "RNA"
assay_star  <- if (exists("select_starmap_assay")) select_starmap_assay(star, prefer_imputed = TRUE) else "RNA_raw"
log_msg(paste0("05C assays -> SlideTags: ", assay_slide, ", STARmap: ", assay_star), logfile)

slide2 <- compute_perm(slide, "SlideTags", assay_slide)
star2  <- compute_perm(star,  "STARmap", assay_star)

saveRDS(slide2, file.path(DIR_OBJS, "slidetags_with_permissiveness.rds"))
saveRDS(star2,  file.path(DIR_OBJS, "starmap_with_permissiveness.rds"))

>>>>>>> 54c1a5c675ca223d183141f8a75de090237ce902
log_msg("Done permissiveness scoring.", logfile)