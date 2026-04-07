======================================================================
  # scripts/05_spatial/05H_spatial_permissiveness_panels_global.R
  # Build per-sample spatial permissiveness panel figures from 05C outputs,
  # including global cross-sample calibration and protected-region maps.
  #
  # Inputs (expected from prior steps):
  # - output/objects/*_with_permissiveness.rds (preferred) or *_harmonized.rds
  # - output/tables/*_permissiveness_cell_level.csv (optional; copied if present)
  #
  # Outputs:
  # - output/figures/05_spatial/<dataset>/<sample>/
  # - output/tables/05_spatial/<dataset>/<sample>/
# - output/tables/05_spatial/summary_table_with_highlights.csv
# - output/tables/05_spatial/permissiveness_global_allcells.csv
# - output/tables/05_spatial/PROCESS_LOG.txt
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(MASS)
  library(png)
  library(jsonlite)
})

source("config/config.R")
source("scripts/R/utils.R")

check_required_packages(
  c("Seurat", "dplyr", "tibble", "ggplot2", "patchwork", "viridis", "MASS", "png", "jsonlite"),
  context = "05H_spatial_permissiveness_panels_global"
)

# Config knobs ---------------------------------------------------------------
TOP_FRAC <- 0.10
BOTTOM_FRAC <- 0.10
MIN_POINTS_FOR_KDE <- 30L
KDE_GRID_N <- 250L
BIO_IMPORTANT <- c("EVT", "Endothelial_cells", "Macrophage", "NK", "Fibroblast", "vCTB", "Hofbauer cells")
TOP_N_DYNAMIC_LABELS <- 3L

DIR_FIG <- file.path(DIR_FIGURES, "05_spatial")
DIR_TAB <- file.path(DIR_TABLES, "05_spatial")
ensure_dir(DIR_FIG); ensure_dir(DIR_TAB); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "05H_spatial_permissiveness_panels_global.log")

plot_bg_points <- function(df, col = "grey80", size = 0.35, alpha = 0.45) {
  ggplot2::geom_point(data = df, aes(x = x, y = y), color = col, size = size, alpha = alpha)
}

z_global <- function(v) {
  v <- suppressWarnings(as.numeric(v))
  ok <- is.finite(v)
  if (sum(ok) < 2) return(rep(0, length(v)))
  mu <- mean(v[ok], na.rm = TRUE)
  sdv <- stats::sd(v[ok], na.rm = TRUE)
  if (!is.finite(sdv) || sdv < 1e-12) return(rep(0, length(v)))
  z <- (v - mu) / sdv
  z[!is.finite(z)] <- 0
  z
}

safe_ecdf <- function(v) {
  v <- suppressWarnings(as.numeric(v))
  out <- rep(NA_real_, length(v))
  ok <- is.finite(v)
  if (sum(ok) == 0) return(out)
  f <- stats::ecdf(v[ok])
  out[ok] <- f(v[ok])
  out
}

safe_kde_df <- function(df, frac = 0.10, mode = c("top", "bottom"), n_grid = 250L) {
  mode <- match.arg(mode)
  d <- df %>% filter(is.finite(x), is.finite(y), is.finite(permissiveness_global))
  if (nrow(d) < MIN_POINTS_FOR_KDE) return(NULL)
  
  k <- max(20L, as.integer(round(frac * nrow(d))))
  sel <- if (mode == "top") {
    d %>% arrange(desc(permissiveness_global)) %>% slice(seq_len(min(k, nrow(d))))
  } else {
    d %>% arrange(permissiveness_global) %>% slice(seq_len(min(k, nrow(d))))
  }
  
  if (nrow(sel) < MIN_POINTS_FOR_KDE) return(NULL)
  
  xr <- range(d$x, na.rm = TRUE)
  yr <- range(d$y, na.rm = TRUE)
  kd <- MASS::kde2d(sel$x, sel$y, n = n_grid, lims = c(xr[1], xr[2], yr[1], yr[2]))
  kd_df <- expand.grid(x = kd$x, y = kd$y)
  kd_df$z <- as.vector(kd$z)
  kd_df
}

pick_dataset_name <- function(path) {
  b <- basename(path)
  b <- sub("_with_permissiveness\\.rds$", "", b)
  b <- sub("_harmonized\\.rds$", "", b)
  b
}

pick_sample_col <- function(md) {
  cands <- c("sample_id", "sample", "orig.ident", "slide", "slide_id", "section")
  hit <- intersect(cands, colnames(md))
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

# Discover objects -----------------------------------------------------------
obj_paths <- c(
  list.files(DIR_OBJS, pattern = "_with_permissiveness\\.rds$", full.names = TRUE),
  list.files(DIR_OBJS, pattern = "_harmonized\\.rds$", full.names = TRUE)
)
obj_paths <- unique(obj_paths)

if (length(obj_paths) == 0) {
  stop("[05H] No *_with_permissiveness.rds or *_harmonized.rds files found in ", DIR_OBJS)
}

log_msg(sprintf("[05H] Found %d object candidates.", length(obj_paths)), logfile)

# Pass 1: collect all cells for global comparability -------------------------
all_cells <- list()
objects <- list()

for (pth in obj_paths) {
  obj <- tryCatch(readRDS(pth), error = function(e) NULL)
  if (is.null(obj) || !inherits(obj, "Seurat")) next
  
  ds <- pick_dataset_name(pth)
  md <- obj@meta.data %>% rownames_to_column("cell")
  md$dataset <- ds
  
  # Ensure week + spatial fields are available consistently.
  obj <- ensure_week_column(obj, COL_WEEK_CANDIDATES)
  obj <- ensure_spatial_coords(obj, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)
  md <- obj@meta.data %>% rownames_to_column("cell") %>% mutate(dataset = ds)
  
  if (!("permissiveness" %in% colnames(md))) next
  
  all_cells[[length(all_cells) + 1L]] <- md %>%
    transmute(
      dataset,
      cell,
      permissiveness = suppressWarnings(as.numeric(permissiveness)),
      score_MMP_ECM_Remodeling = suppressWarnings(as.numeric(.data$score_MMP_ECM_Remodeling)),
      score_Immune_Tolerance = suppressWarnings(as.numeric(.data$score_Immune_Tolerance)),
      score_Cytotoxic_NK = suppressWarnings(as.numeric(.data$score_Cytotoxic_NK)),
      score_Ethanolamine_Metabolism = suppressWarnings(as.numeric(.data$score_Ethanolamine_Metabolism))
    )
  
  objects[[ds]] <- obj
}

if (length(all_cells) == 0) {
  stop("[05H] No permissiveness metadata found in candidate objects. Run 05C first.")
}

all_md <- bind_rows(all_cells)

all_md <- all_md %>% mutate(
  z_mmp_global = z_global(score_MMP_ECM_Remodeling),
  z_tol_global = z_global(score_Immune_Tolerance),
  z_nk_global = z_global(score_Cytotoxic_NK),
  z_ea_global = z_global(score_Ethanolamine_Metabolism),
  permissiveness_global = z_mmp_global + z_tol_global - z_nk_global + z_ea_global,
  permissiveness_global_pct = safe_ecdf(permissiveness_global)
)

write.csv(all_md, file.path(DIR_TAB, "permissiveness_global_allcells.csv"), row.names = FALSE)
log_msg(sprintf("[05H] Wrote global calibration table with %d cells.", nrow(all_md)), logfile)

# Pass 2: sample-level figures -----------------------------------------------
summary_rows <- list()
process_log <- c(sprintf("Run at: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

for (ds in names(objects)) {
  obj <- objects[[ds]]
  md <- obj@meta.data %>% rownames_to_column("cell")
  
  sample_col <- pick_sample_col(md)
  samples <- if (is.null(sample_col)) "all" else unique(as.character(md[[sample_col]]))
  
  for (sp in samples) {
    md_s <- if (is.null(sample_col)) {
      md
    } else {
      md %>% filter(as.character(.data[[sample_col]]) == sp)
    }
    if (nrow(md_s) == 0) next
    
    md_s <- md_s %>%
      left_join(all_md %>% select(dataset, cell, permissiveness_global, permissiveness_global_pct),
                by = c("dataset" = "dataset", "cell" = "cell"))
    
    md_s <- md_s %>%
      mutate(
        x = suppressWarnings(as.numeric(spatial_x_use)),
        y = suppressWarnings(as.numeric(spatial_y_use))
      ) %>%
      filter(is.finite(x), is.finite(y), is.finite(permissiveness_global))
    
    if (nrow(md_s) == 0) {
      process_log <- c(process_log, sprintf("[%s/%s] skipped: no finite x/y/permissiveness_global", ds, sp))
      next
    }
    
    label_col <- pick_first_present(md_s, c("celltype_final_refined", "celltype_final_conservative", "celltype_author", "predicted.id"))
    md_s$label <- if (!is.null(label_col)) as.character(md_s[[label_col]]) else "unknown"
    
    dyn_labels <- md_s %>% group_by(label) %>% summarize(mean_perm = mean(permissiveness_global, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(mean_perm)) %>% slice_head(n = TOP_N_DYNAMIC_LABELS) %>% pull(label)
    hl <- unique(c(intersect(BIO_IMPORTANT, unique(md_s$label)), dyn_labels))
    
    fig_dir <- file.path(DIR_FIG, ds, as.character(sp))
    tab_dir <- file.path(DIR_TAB, ds, as.character(sp))
    ensure_dir(fig_dir); ensure_dir(tab_dir)
    
    # Copy 05C per-dataset tables into sample folder for provenance.
    copy_tbls <- list.files(DIR_TABLES, pattern = paste0("^", ds, "_.*(permissiveness_cell_level|NK_module_coverage_qc)\\.csv$"),
                            full.names = TRUE)
    if (length(copy_tbls) > 0) file.copy(copy_tbls, tab_dir, overwrite = TRUE)
    
    p_base <- ggplot(md_s, aes(x = x, y = y)) +
      theme_classic() +
      coord_fixed() +
      labs(x = "spatial x", y = "spatial y")
    
    p_overlay_global <- p_base +
      geom_point(aes(color = permissiveness_global), size = 0.55, alpha = 0.9) +
      scale_color_viridis_c(option = "viridis", trans = "sqrt") +
      ggtitle(sprintf("%s / %s: Global permissiveness", ds, sp))
    
    p_overlay_white_red <- p_base +
      geom_point(aes(color = permissiveness_global), size = 0.55, alpha = 0.9) +
      scale_color_gradient2(low = "white", mid = "#ffd6c7", high = "red", midpoint = 0) +
      ggtitle(sprintf("%s / %s: Global permissiveness (white->red)", ds, sp))
    
    p_overlay_hl <- p_base +
      geom_point(aes(color = permissiveness_global), size = 0.45, alpha = 0.65) +
      scale_color_viridis_c(option = "viridis", trans = "sqrt") +
      ggtitle(sprintf("%s / %s: Highlights + global permissiveness", ds, sp))
    
    if (length(hl) > 0) {
      p_overlay_hl <- p_overlay_hl +
        geom_point(data = md_s %>% filter(label %in% hl),
                   aes(fill = label), shape = 21, size = 2.4, color = "black", stroke = 0.5, alpha = 0.9)
    }
    
    kd_top <- safe_kde_df(md_s, frac = TOP_FRAC, mode = "top", n_grid = KDE_GRID_N)
    kd_bottom <- safe_kde_df(md_s, frac = BOTTOM_FRAC, mode = "bottom", n_grid = KDE_GRID_N)
    
    p_top <- p_base + ggtitle(sprintf("%s / %s: Top %.0f%% permissive density", ds, sp, TOP_FRAC * 100))
    if (!is.null(kd_top)) {
      p_top <- p_top + geom_raster(data = kd_top, aes(x = x, y = y, fill = z), interpolate = TRUE) +
        scale_fill_viridis_c(option = "magma")
    }
    
    p_bottom <- p_base + ggtitle(sprintf("%s / %s: Bottom %.0f%% (protected) density", ds, sp, BOTTOM_FRAC * 100))
    if (!is.null(kd_bottom)) {
      p_bottom <- p_bottom + geom_raster(data = kd_bottom, aes(x = x, y = y, fill = z), interpolate = TRUE) +
        scale_fill_viridis_c(option = "plasma")
    }
    
    p_diff <- p_base + ggtitle(sprintf("%s / %s: Hotspot - Protected density", ds, sp))
    if (!is.null(kd_top) && !is.null(kd_bottom)) {
      kd_diff <- kd_top %>% select(x, y, z_top = z) %>%
        left_join(kd_bottom %>% select(x, y, z_bottom = z), by = c("x", "y")) %>%
        mutate(z_diff = z_top - z_bottom)
      lim <- max(abs(kd_diff$z_diff), na.rm = TRUE)
      if (!is.finite(lim) || lim <= 0) lim <- 1
      p_diff <- p_diff +
        geom_raster(data = kd_diff, aes(x = x, y = y, fill = z_diff), interpolate = TRUE) +
        scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                             limits = c(-lim, lim))
    }
    
    ggsave(file.path(fig_dir, paste0(ds, "_", sp, "_3panel_global.png")),
           (p_overlay_global + p_top + p_overlay_hl) + plot_layout(nrow = 1), width = 18, height = 6, dpi = 300)
    ggsave(file.path(fig_dir, paste0(ds, "_", sp, "_3panel_global_white_red.png")),
           (p_overlay_white_red + p_top + p_overlay_hl) + plot_layout(nrow = 1), width = 18, height = 6, dpi = 300)
    ggsave(file.path(fig_dir, paste0(ds, "_", sp, "_protected_and_diff.png")),
           (p_top + p_bottom + p_diff) + plot_layout(nrow = 1), width = 18, height = 6, dpi = 300)
    
    # Per-week global maps.
    if ("week" %in% colnames(md_s)) {
      weeks <- sort(unique(md_s$week[!is.na(md_s$week)]))
      for (w in weeks) {
        d <- md_s %>% filter(week == w)
        if (nrow(d) == 0) next
        p_w <- ggplot(d, aes(x = x, y = y, color = permissiveness_global)) +
          geom_point(size = 0.55, alpha = 0.9) +
          scale_color_viridis_c(option = "viridis", trans = "sqrt") +
          coord_fixed() + theme_classic() +
          labs(title = sprintf("%s / %s week %s: global permissiveness", ds, sp, w), x = "spatial x", y = "spatial y")
        ggsave(file.path(fig_dir, paste0(ds, "_", sp, "_permissiveness_global_week_", w, ".png")),
               p_w, width = 10, height = 8, dpi = 300)
      }
    }
    
    # Export hotspot/protected cell lists for downstream DE/CellChat.
    top_thr <- stats::quantile(md_s$permissiveness_global, probs = 1 - TOP_FRAC, na.rm = TRUE)
    bot_thr <- stats::quantile(md_s$permissiveness_global, probs = BOTTOM_FRAC, na.rm = TRUE)
    hotspots <- md_s %>% filter(permissiveness_global >= top_thr)
    protected <- md_s %>% filter(permissiveness_global <= bot_thr)
    write.csv(hotspots, file.path(tab_dir, paste0(ds, "_", sp, "_hotspot_cells_top", as.integer(TOP_FRAC * 100), "pct.csv")), row.names = FALSE)
    write.csv(protected, file.path(tab_dir, paste0(ds, "_", sp, "_protected_cells_bottom", as.integer(BOTTOM_FRAC * 100), "pct.csv")), row.names = FALSE)
    
    write.csv(md_s, file.path(tab_dir, paste0(ds, "_", sp, "_mapped_coords_global.csv")), row.names = FALSE)
    
    meta <- list(
      dataset = ds,
      sample = as.character(sp),
      created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      n_cells = nrow(md_s),
      highlights = hl,
      global_score_definition = "z_global(MMP) + z_global(Tolerance) - z_global(NK) + z_global(EA)",
      top_fraction = TOP_FRAC,
      bottom_fraction = BOTTOM_FRAC,
      comparability_note = "permissiveness_global is calibrated across pooled cells from all datasets loaded in this run"
    )
    write_json(meta, file.path(tab_dir, "meta.json"), pretty = TRUE, auto_unbox = TRUE)
    
    readme <- c(
      paste0("Dataset: ", ds),
      paste0("Sample: ", sp),
      "This folder contains global-calibrated spatial permissiveness outputs.",
      "Global score = z_global(MMP) + z_global(Immune_Tolerance) - z_global(Cytotoxic_NK) + z_global(Ethanolamine).",
      "protected_and_diff figure includes bottom-percentile protected density and hotspot-protected difference.",
      paste0("Top fraction: ", TOP_FRAC, "; Bottom fraction: ", BOTTOM_FRAC)
    )
    writeLines(readme, file.path(tab_dir, "README_sample.txt"))
    
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      dataset = ds,
      sample = as.character(sp),
      n_cells = nrow(md_s),
      highlight_labels = paste(hl, collapse = ";"),
      fig_dir = fig_dir,
      tab_dir = tab_dir,
      stringsAsFactors = FALSE
    )
    
    process_log <- c(process_log, sprintf("[%s/%s] done (cells=%d)", ds, sp, nrow(md_s)))
  }
}

if (length(summary_rows) > 0) {
  summary_df <- bind_rows(summary_rows)
  write.csv(summary_df, file.path(DIR_TAB, "summary_table_with_highlights.csv"), row.names = FALSE)
  log_msg(sprintf("[05H] Summary written: %d rows.", nrow(summary_df)), logfile)
} else {
  log_msg("[05H] No sample-level outputs produced.", logfile)
}

writeLines(process_log, file.path(DIR_TAB, "PROCESS_LOG.txt"))
log_msg("[05H] Done.", logfile)