<<<<<<< HEAD
# ======================================================================
# scripts/05_spatial/05H_spatial_permissiveness_panels_global.R
# Build per-sample spatial permissiveness panel figures from 05C outputs,
# including global cross-sample calibration and protected-region maps.
#
# Run from repo root:
# source("scripts/05_spatial/05H_spatial_permissiveness_panels_global.R")
# ======================================================================
# scripts/05_spatial/05H_spatial_permissiveness_panels_global.R
# Build per-sample spatial permissiveness panel figures from 05C outputs,
# including global cross-sample calibration and protected-region maps.
#
# Run from repo root:
# source("scripts/05_spatial/05H_spatial_permissiveness_panels_global.R")

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
TOP_FRAC <- 0.10          # top X% = hotspot
BOTTOM_FRAC <- 0.10       # bottom X% = protected
MIN_POINTS_FOR_KDE <- 30L
KDE_GRID_N <- 250L
BIO_IMPORTANT <- c("EVT", "Endothelial_cells", "Macrophage", "NK", "Fibroblast", "vCTB", "Hofbauer cells")
TOP_N_DYNAMIC_LABELS <- 3L

DIR_FIG <- file.path(DIR_FIGURES, "05_spatial")
DIR_TAB <- file.path(DIR_TABLES, "05_spatial")
ensure_dir(DIR_FIG); ensure_dir(DIR_TAB); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "05H_spatial_permissiveness_panels_global.log")

# small logger
log_msg2 <- function(msg, logfile = logfile) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  tryCatch(cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg), file = logfile, append = TRUE),
           error = function(e) NULL)
}

plot_bg_points <- function(df, col = "grey80", size = 0.35, alpha = 0.45) {
  ggplot2::geom_point(data = df, aes(x = x, y = y), color = col, size = size, alpha = alpha)
}

# robust/global z
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
  xr <- range(d$x, na.rm = TRUE); yr <- range(d$y, na.rm = TRUE)
  kd <- tryCatch(MASS::kde2d(sel$x, sel$y, n = n_grid, lims = c(xr[1], xr[2], yr[1], yr[2])),
                 error = function(e) NULL)
  if (is.null(kd)) return(NULL)
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

# New helper used in script
pick_first_present <- function(df_or_cols, cols) {
  if (is.data.frame(df_or_cols)) {
    avail <- intersect(cols, colnames(df_or_cols))
  } else {
    avail <- intersect(cols, df_or_cols)
  }
  if (length(avail) == 0) return(NULL)
  avail[[1]]
}

# Helper: write per-figure caption text
write_caption <- function(fig_path, caption_text) {
  cap_path <- sub("\\.(png|pdf)$", "_caption.txt", fig_path, ignore.case = TRUE)
  tryCatch(writeLines(caption_text, cap_path), error = function(e) NULL)
}

# Create a repo-level description (05F + 05H) to place in scripts/05_spatial
write_repo_description <- function(out_path) {
  txt <- c(
    "# 05F and 05H description",
    "",
    "This file describes the outputs and interpretation guidelines for:",
    "- scripts/05_spatial/05C_permissiveness_score_maps.R (05C)",
    "- scripts/05_spatial/05H_spatial_permissiveness_panels_global.R (05H)",
    "",
    "### Purpose",
    "05C: compute dataset-relative permissiveness per cell: permissiveness = z_MMP + z_Tolerance - z_NK + z_Ethanolamine (z computed within each dataset).",
    "05H: pooled global calibration across datasets to allow cross-sample/time comparisons: z_global computed across pooled cells.",
    "",
    "### Figure types created by 05H (per-sample)",
    "1. 3-panel global: (A) image/background, (B) KDE of top 10% permissive cells, (C) global-permissiveness overlay (viridis).",
    "2. 3-panel white->red: (A) image, (B) same KDE, (C) white->red overlay emphasizing positive permissiveness.",
    "3. Cell-type highlighted variant: (A) image, (B) same KDE, (C) overlay with highlighted biologically-important cell types (filled circles) on top of permissiveness.",
    "4. Protected/diff panels: (A) top-density, (B) bottom-density (protected), (C) hotspot - protected difference map.",
    "",
    "### Interpretation guidance (for figure legends / methods)",
    "- Hotspot defined as KDE of the top 10% permissive cells (dataset-pooled percentile when using permissiveness_global).",
    "- Protected regions: KDE of the bottom 10% permissive cells.",
    "- Difference map = hotspot_density - protected_density (positive = permissive, negative = protected).",
    "- Include a note whether permissiveness is dataset-relative (05C) or globally-calibrated (05H).",
    "",
    "### Files produced",
    "- output/figures/05_spatial/<dataset>/<sample>/ : PNG/PDF figures and caption .txt files",
    "- output/tables/05_spatial/ : permissiveness_global_allcells.csv, summary_table_with_highlights.csv, PROCESS_LOG.txt, per-sample hotspot/protected CSVs and meta.json files",
    "",
    "### Methods notes",
    "See the Methods bullets saved to scripts/05_spatial/METHODS_05C_05H.md for full reproducible details."
  )
  tryCatch(writeLines(txt, out_path), error = function(e) NULL)
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

log_msg2(sprintf("[05H] Found %d object candidates.", length(obj_paths)), logfile)

# Pass 1: collect all cells for global comparability -------------------------
all_cells <- list()
objects <- list()
process_log <- list()
process_log[[length(process_log) + 1L]] <- paste0("Run at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

for (pth in obj_paths) {
  obj <- tryCatch(readRDS(pth), error = function(e) { log_msg2(paste0("  ERROR reading ", pth, ": ", e$message)); NULL })
  if (is.null(obj) || !inherits(obj, "Seurat")) {
    process_log[[length(process_log) + 1L]] <- paste0("Skipping non-Seurat or unreadable: ", pth); next
  }
  ds <- pick_dataset_name(pth)
  
  # ensure week + spatial columns present (try/catch)
  obj <- tryCatch({
    obj2 <- ensure_week_column(obj, COL_WEEK_CANDIDATES)
    obj2 <- ensure_spatial_coords(obj2, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)
    obj2
  }, error = function(e) { log_msg2(paste0("  WARNING: ensure_week/coords failed for ", pth, ": ", e$message)); obj })
  
  # md from object; ensure dataset column
  md <- obj@meta.data %>% rownames_to_column("cell") %>% mutate(dataset = ds)
  
  if (!("permissiveness" %in% colnames(md))) {
    process_log[[length(process_log) + 1L]] <- paste0("[", ds, "] skipping: no 'permissiveness' column in metadata")
    objects[[ds]] <- obj  # keep object for plotting if needed
    next
  }
  
  # Safely collect score columns; warn if missing
  required_scores <- c("score_MMP_ECM_Remodeling", "score_Immune_Tolerance", "score_Cytotoxic_NK", "score_Ethanolamine_Metabolism")
  missing_scores <- setdiff(required_scores, colnames(md))
  if (length(missing_scores) > 0) {
    process_log[[length(process_log) + 1L]] <- paste0("[", ds, "] WARNING: missing score columns: ", paste(missing_scores, collapse = "; "))
    for (mc in missing_scores) md[[mc]] <- NA_real_
  }
  
  # push to all_cells
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
  process_log[[length(process_log) + 1L]] <- paste0("[", ds, "] queued (n_cells=", nrow(md), ")")
}

if (length(all_cells) == 0) {
  log_msg2("[05H] No permissiveness metadata found in candidate objects. Run 05C first.", logfile)
  writeLines(unlist(process_log), file.path(DIR_TAB, "PROCESS_LOG.txt"))
  stop("[05H] No permissiveness metadata found in candidate objects. Run 05C first.")
}

all_md <- bind_rows(all_cells)

# Check required score columns exist in pooled table before z-scoring
required_scores <- c("score_MMP_ECM_Remodeling", "score_Immune_Tolerance", "score_Cytotoxic_NK", "score_Ethanolamine_Metabolism")
if (!all(required_scores %in% colnames(all_md))) {
  missing_scores <- setdiff(required_scores, colnames(all_md))
  log_msg2(paste0("[05H] WARNING: pooled table missing score columns: ", paste(missing_scores, collapse = "; "), " - those components will be NA"), logfile)
}

# compute global z and permissiveness (NA-safe)
all_md <- all_md %>% mutate(
  z_mmp_global = if ("score_MMP_ECM_Remodeling" %in% colnames(all_md)) z_global(score_MMP_ECM_Remodeling) else NA_real_,
  z_tol_global = if ("score_Immune_Tolerance" %in% colnames(all_md)) z_global(score_Immune_Tolerance) else NA_real_,
  z_nk_global  = if ("score_Cytotoxic_NK" %in% colnames(all_md)) z_global(score_Cytotoxic_NK) else NA_real_,
  z_ea_global  = if ("score_Ethanolamine_Metabolism" %in% colnames(all_md)) z_global(score_Ethanolamine_Metabolism) else NA_real_,
  permissiveness_global = rowSums(cbind(z_mmp_global, z_tol_global, -z_nk_global, z_ea_global), na.rm = FALSE),
  permissiveness_global_pct = safe_ecdf(permissiveness_global)
)

# If permissiveness_global is NA for many rows, warn
if (sum(is.finite(all_md$permissiveness_global)) < 2) {
  log_msg2("[05H] WARNING: permissiveness_global could not be computed (not enough finite components).", logfile)
}

# Write pooled table
tryCatch({
  write.csv(all_md, file.path(DIR_TAB, "permissiveness_global_allcells.csv"), row.names = FALSE)
  log_msg2(sprintf("[05H] Wrote global calibration table with %d cells.", nrow(all_md)), logfile)
}, error = function(e) log_msg2(paste0("[05H] Warning: cannot write permissiveness_global_allcells.csv: ", e$message), logfile))

# Create repo-level description file in scripts/05_spatial
tryCatch({
  ensure_dir(file.path("scripts", "05_spatial"))
  write_repo_description(file.path("scripts", "05_spatial", "05F_and_05H_description.md"))
}, error = function(e) NULL)

# Pass 2: sample-level figures -----------------------------------------------
summary_rows <- list()
process_log[[length(process_log) + 1L]] <- paste0("Pass2 start: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

for (ds in names(objects)) {
  obj <- objects[[ds]]
  if (!inherits(obj, "Seurat")) {
    process_log[[length(process_log) + 1L]] <- paste0("Skipping non-Seurat object for dataset: ", ds); next
  }
  
  # Ensure md carries dataset name for robust joins later
  md <- obj@meta.data %>% rownames_to_column("cell") %>% mutate(dataset = ds)
  
  sample_col <- pick_sample_col(md)
  samples <- if (is.null(sample_col)) "all" else unique(as.character(md[[sample_col]]))
  
  for (sp in samples) {
    md_s <- if (is.null(sample_col)) { md } else { md %>% filter(as.character(.data[[sample_col]]) == sp) }
    if (nrow(md_s) == 0) {
      process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] no cells -> skip")
      next
    }
    
    # left_join with pooled all_md; ensure join keys present
    if (!("dataset" %in% colnames(md_s)) || !("cell" %in% colnames(md_s))) {
      process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] missing dataset/cell columns; skipping join")
      next
    }
    
    joined <- tryCatch({
      md_s %>%
        left_join(all_md %>% select(dataset, cell, permissiveness_global, permissiveness_global_pct),
                  by = c("dataset" = "dataset", "cell" = "cell"))
    }, error = function(e) {
      process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] left_join failed: ", e$message)
      NULL
    })
    if (is.null(joined)) {
      process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] join returned NULL -> skipping sample")
      next
    }
    md_s <- joined
    
    if (!("permissiveness_global" %in% colnames(md_s)) || sum(is.finite(md_s$permissiveness_global)) < 2) {
      process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] insufficient permissiveness_global after join -> skipping sample")
      next
    }
    
    md_s <- md_s %>%
      mutate(x = suppressWarnings(as.numeric(spatial_x_use)),
             y = suppressWarnings(as.numeric(spatial_y_use))) %>%
      filter(is.finite(x), is.finite(y), is.finite(permissiveness_global))
    
    if (nrow(md_s) == 0) {
      process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] skipped: no finite x/y/permissiveness_global")
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
      theme_classic() + coord_fixed() + labs(x = "spatial x", y = "spatial y")
    
    # overlay variants
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
      # color palette for labels (Set2)
      n_hl <- length(hl)
      pal <- RColorBrewer::brewer.pal(max(3, min(8, n_hl)), "Set2")
      fill_vals <- setNames(rep(pal, length.out = n_hl), hl)
      p_overlay_hl <- p_overlay_hl +
        geom_point(data = md_s %>% filter(label %in% hl),
                   aes(fill = label), shape = 21, size = 2.4, color = "black", stroke = 0.5, alpha = 0.9) +
        scale_fill_manual(name = "Highlighted cell types", values = fill_vals)
    }
    
    kd_top <- safe_kde_df(md_s, frac = TOP_FRAC, mode = "top", n_grid = KDE_GRID_N)
    kd_bottom <- safe_kde_df(md_s, frac = BOTTOM_FRAC, mode = "bottom", n_grid = KDE_GRID_N)
    
    p_top <- p_base + ggtitle(sprintf("%s / %s: Top %.0f%% permissive density", ds, sp, TOP_FRAC * 100))
    if (!is.null(kd_top)) {
      p_top <- p_top + geom_raster(data = kd_top, aes(x = x, y = y, fill = z), interpolate = TRUE) +
        scale_fill_viridis_c(option = "magma")
    } else {
      p_top <- p_top + ggtitle(paste0("Too few points for top density (need >= ", MIN_POINTS_FOR_KDE, ")"))
    }
    
    p_bottom <- p_base + ggtitle(sprintf("%s / %s: Bottom %.0f%% (protected) density", ds, sp, BOTTOM_FRAC * 100))
    if (!is.null(kd_bottom)) {
      p_bottom <- p_bottom + geom_raster(data = kd_bottom, aes(x = x, y = y, fill = z), interpolate = TRUE) +
        scale_fill_viridis_c(option = "plasma")
    } else {
      p_bottom <- p_bottom + ggtitle(paste0("Too few points for bottom density (need >= ", MIN_POINTS_FOR_KDE, ")"))
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
    } else {
      p_diff <- p_diff + ggtitle("Diff not available (top or bottom density missing)")
    }
    
    # figure paths
    f3_global <- file.path(fig_dir, paste0(ds, "_", sp, "_3panel_global.png"))
    f3_white_red <- file.path(fig_dir, paste0(ds, "_", sp, "_3panel_global_white_red.png"))
    f_prot_diff <- file.path(fig_dir, paste0(ds, "_", sp, "_protected_and_diff.png"))
    
    # Save panels (with try-catch)
    tryCatch({
      ggsave(f3_global, (p_overlay_global + p_top + p_overlay_hl) + plot_layout(nrow = 1),
             width = 18, height = 6, dpi = 300)
      # caption
      write_caption(f3_global,
                    c(sprintf("3-panel global (dataset=%s sample=%s)", ds, sp),
                      "- Panel A: original sample coordinates (background); Panel B: KDE of top 10% permissive cells (magma); Panel C: per-cell global permissiveness (viridis).",
                      "- Interpretation: hotspots in Panel B correspond to clusters where permissiveness is concentrated; Panel C shows per-cell global-level permissiveness."))
    }, error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing _3panel_global.png: ", e$message))
    
    tryCatch({
      ggsave(f3_white_red, (p_overlay_white_red + p_top + p_overlay_hl) + plot_layout(nrow = 1),
             width = 18, height = 6, dpi = 300)
      write_caption(f3_white_red,
                    c(sprintf("3-panel white->red (dataset=%s sample=%s)", ds, sp),
                      "- Panel C uses white->red to emphasize positive permissiveness (centered at zero).",
                      "- Use this panel to emphasize regions with strongly positive permissiveness."))
    }, error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing _3panel_global_white_red.png: ", e$message))
    
    tryCatch({
      ggsave(f_prot_diff, (p_top + p_bottom + p_diff) + plot_layout(nrow = 1),
             width = 18, height = 6, dpi = 300)
      write_caption(f_prot_diff,
                    c(sprintf("Protected and difference maps (dataset=%s sample=%s)", ds, sp),
                      "- Left: KDE of top 10% (hotspot). Middle: KDE of bottom 10% (protected). Right: hotspot - protected difference (blue negative = protected; red positive = permissive)."))
    }, error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing _protected_and_diff.png: ", e$message))
    
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
        out_w <- file.path(fig_dir, paste0(ds, "_", sp, "_permissiveness_global_week_", w, ".png"))
        tryCatch({
          ggsave(out_w, p_w, width = 10, height = 8, dpi = 300)
          write_caption(out_w, c(sprintf("Global permissiveness (dataset=%s sample=%s week=%s)", ds, sp, w),
                                 "- Interpret relative to other weeks/samples only when using permissiveness_global (pooled calibration)."))
        }, error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing weekly png: ", e$message))
      }
    }
    
    # Export hotspot/protected cell lists
    top_thr <- tryCatch(stats::quantile(md_s$permissiveness_global, probs = 1 - TOP_FRAC, na.rm = TRUE), error = function(e) NA_real_)
    bot_thr <- tryCatch(stats::quantile(md_s$permissiveness_global, probs = BOTTOM_FRAC, na.rm = TRUE), error = function(e) NA_real_)
    hotspots <- if (!is.na(top_thr)) md_s %>% filter(permissiveness_global >= top_thr) else md_s[0,]
    protected <- if (!is.na(bot_thr)) md_s %>% filter(permissiveness_global <= bot_thr) else md_s[0,]
    tryCatch(write.csv(hotspots, file.path(tab_dir, paste0(ds, "_", sp, "_hotspot_cells_top", as.integer(TOP_FRAC * 100), "pct.csv")), row.names = FALSE),
             error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing hotspots CSV: ", e$message))
    tryCatch(write.csv(protected, file.path(tab_dir, paste0(ds, "_", sp, "_protected_cells_bottom", as.integer(BOTTOM_FRAC * 100), "pct.csv")), row.names = FALSE),
             error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing protected CSV: ", e$message))
    
    tryCatch(write.csv(md_s, file.path(tab_dir, paste0(ds, "_", sp, "_mapped_coords_global.csv")), row.names = FALSE),
             error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing mapped coords CSV: ", e$message))
    
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
    tryCatch(write_json(meta, file.path(tab_dir, "meta.json"), pretty = TRUE, auto_unbox = TRUE),
             error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing meta.json: ", e$message))
    
    readme <- c(
      paste0("Dataset: ", ds),
      paste0("Sample: ", sp),
      "This folder contains global-calibrated spatial permissiveness outputs.",
      "Global score = z_global(MMP) + z_global(Immune_Tolerance) - z_global(Cytotoxic_NK) + z_global(Ethanolamine).",
      "protected_and_diff figure includes bottom-percentile protected density and hotspot-protected difference.",
      paste0("Top fraction: ", TOP_FRAC, "; Bottom fraction: ", BOTTOM_FRAC)
    )
    tryCatch(writeLines(readme, file.path(tab_dir, "README_sample.txt")),
             error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing README: ", e$message))
    
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      dataset = ds,
      sample = as.character(sp),
      n_cells = nrow(md_s),
      highlight_labels = paste(hl, collapse = ";"),
      fig_dir = fig_dir,
      tab_dir = tab_dir,
      stringsAsFactors = FALSE
    )
    
    process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] done (cells=", nrow(md_s), ")")
  }
}

# summary + logs
if (length(summary_rows) > 0) {
  summary_df <- bind_rows(summary_rows)
  tryCatch(write.csv(summary_df, file.path(DIR_TAB, "summary_table_with_highlights.csv"), row.names = FALSE),
           error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("Error writing summary_table_with_highlights.csv: ", e$message))
  log_msg2(sprintf("[05H] Summary written: %d rows.", nrow(summary_df)), logfile)
} else {
  log_msg2("[05H] No sample-level outputs produced.", logfile)
}

# PROCESS LOG
proc_file <- file.path(DIR_TAB, "PROCESS_LOG.txt")
tryCatch(writeLines(unlist(c(process_log)), proc_file), error = function(e) log_msg2(paste0("Failed to write PROCESS_LOG: ", e$message), logfile))
log_msg2("[05H] Done.", logfile)

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

# small logger (prints + appends to logfile)
log_msg2 <- function(msg, logfile = logfile) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  tryCatch(cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg), file = logfile, append = TRUE),
           error = function(e) NULL)
}

# ---- Helpers ----
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
  kd <- tryCatch(MASS::kde2d(sel$x, sel$y, n = n_grid, lims = c(xr[1], xr[2], yr[1], yr[2])),
                 error = function(e) NULL)
  if (is.null(kd)) return(NULL)
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

# New helper used in script
pick_first_present <- function(df_or_cols, cols) {
  if (is.data.frame(df_or_cols)) {
    avail <- intersect(cols, colnames(df_or_cols))
  } else {
    avail <- intersect(cols, df_or_cols)
  }
  if (length(avail) == 0) return(NULL)
  avail[[1]]
}

# ---- Discover objects ----
obj_paths <- c(
  list.files(DIR_OBJS, pattern = "_with_permissiveness\\.rds$", full.names = TRUE),
  list.files(DIR_OBJS, pattern = "_harmonized\\.rds$", full.names = TRUE)
)
obj_paths <- unique(obj_paths)

if (length(obj_paths) == 0) {
  stop("[05H] No *_with_permissiveness.rds or *_harmonized.rds files found in ", DIR_OBJS)
}

log_msg2(sprintf("[05H] Found %d object candidates.", length(obj_paths)), logfile)

# ---- Pass 1: collect all cells for global comparability ----
all_cells <- list()
objects <- list()
process_log <- list()
process_log[[length(process_log) + 1L]] <- paste0("Run at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

for (pth in obj_paths) {
  obj <- tryCatch(readRDS(pth), error = function(e) { log_msg2(paste0("  ERROR reading ", pth, ": ", e$message)); NULL })
  if (is.null(obj) || !inherits(obj, "Seurat")) {
    process_log[[length(process_log) + 1L]] <- paste0("Skipping non-Seurat or unreadable: ", pth); next
  }
  ds <- pick_dataset_name(pth)
  # ensure week + spatial columns present
  obj <- tryCatch({
    obj2 <- ensure_week_column(obj, COL_WEEK_CANDIDATES)
    obj2 <- ensure_spatial_coords(obj2, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)
    obj2
  }, error = function(e) { log_msg2(paste0("  WARNING: ensure_week/coords failed for ", pth, ": ", e$message)); obj })
  
  # md from object; ensure dataset column
  md <- obj@meta.data %>% rownames_to_column("cell") %>% mutate(dataset = ds)
  
  if (!("permissiveness" %in% colnames(md))) {
    process_log[[length(process_log) + 1L]] <- paste0("[", ds, "] skipping: no 'permissiveness' column in metadata")
    objects[[ds]] <- obj  # still keep for plotting (maybe permissiveness propagated later), but don't add cells
    next
  }
  
  # Safely collect columns; if missing, fill NA and warn
  required_scores <- c("score_MMP_ECM_Remodeling", "score_Immune_Tolerance", "score_Cytotoxic_NK", "score_Ethanolamine_Metabolism")
  missing_scores <- setdiff(required_scores, colnames(md))
  if (length(missing_scores) > 0) {
    process_log[[length(process_log) + 1L]] <- paste0("[", ds, "] WARNING: missing score columns: ", paste(missing_scores, collapse = "; "))
    # create NA columns to keep table shape
    for (mc in missing_scores) md[[mc]] <- NA_real_
  }
  
  # push to all_cells
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
  process_log[[length(process_log) + 1L]] <- paste0("[", ds, "] queued (n_cells=", nrow(md), ")")
}

if (length(all_cells) == 0) {
  log_msg2("[05H] No permissiveness metadata found in candidate objects. Run 05C first.", logfile)
  writeLines(unlist(process_log), file.path(DIR_TAB, "PROCESS_LOG.txt"))
  stop("[05H] No permissiveness metadata found in candidate objects. Run 05C first.")
}

all_md <- bind_rows(all_cells)

# Check required score columns exist in pooled table before z-scoring
required_scores <- c("score_MMP_ECM_Remodeling", "score_Immune_Tolerance", "score_Cytotoxic_NK", "score_Ethanolamine_Metabolism")
if (!all(required_scores %in% colnames(all_md))) {
  missing_scores <- setdiff(required_scores, colnames(all_md))
  log_msg2(paste0("[05H] WARNING: pooled table missing score columns: ", paste(missing_scores, collapse = "; "), " - those components will be NA"), logfile)
}

# compute global z and permissiveness (NA-safe)
all_md <- all_md %>% mutate(
  z_mmp_global = if ("score_MMP_ECM_Remodeling" %in% colnames(all_md)) z_global(score_MMP_ECM_Remodeling) else NA_real_,
  z_tol_global = if ("score_Immune_Tolerance" %in% colnames(all_md)) z_global(score_Immune_Tolerance) else NA_real_,
  z_nk_global  = if ("score_Cytotoxic_NK" %in% colnames(all_md)) z_global(score_Cytotoxic_NK) else NA_real_,
  z_ea_global  = if ("score_Ethanolamine_Metabolism" %in% colnames(all_md)) z_global(score_Ethanolamine_Metabolism) else NA_real_,
  permissiveness_global = rowSums(cbind(z_mmp_global, z_tol_global, -z_nk_global, z_ea_global), na.rm = FALSE),
  permissiveness_global_pct = safe_ecdf(permissiveness_global)
)

# If permissiveness_global is NA for many rows, warn
if (sum(is.finite(all_md$permissiveness_global)) < 2) {
  log_msg2("[05H] WARNING: permissiveness_global could not be computed (not enough finite components).", logfile)
}

# === global stats & symmetric color limit for white->red ===
perm_maxabs <- NA_real_
if (exists("all_md") && "permissiveness_global" %in% colnames(all_md)) {
  perm_maxabs <- suppressWarnings(max(abs(all_md$permissiveness_global), na.rm = TRUE))
  if (!is.finite(perm_maxabs)) perm_maxabs <- NA_real_
}
# Save global stats (provenance)
global_stats <- list(
  n_cells_pooled = ifelse(exists("all_md"), nrow(all_md), NA_integer_),
  perm_maxabs = perm_maxabs,
  TOP_FRAC = TOP_FRAC,
  BOTTOM_FRAC = BOTTOM_FRAC
)
# create dir if missing
tryCatch({
  if (!dir.exists(DIR_TAB)) dir.create(DIR_TAB, recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(global_stats, file.path(DIR_TAB, "permissiveness_global_stats.json"), pretty = TRUE, auto_unbox = TRUE)
  log_msg2(sprintf("[05H] Saved permissiveness_global_stats.json (perm_maxabs=%s)", as.character(perm_maxabs)), logfile)
}, error = function(e) log_msg2(paste0("[05H] Warning: failed to write global stats: ", e$message), logfile))

# Write pooled table
tryCatch({
  write.csv(all_md, file.path(DIR_TAB, "permissiveness_global_allcells.csv"), row.names = FALSE)
  log_msg2(sprintf("[05H] Wrote global calibration table with %d cells.", nrow(all_md)), logfile)
}, error = function(e) log_msg2(paste0("[05H] Warning: cannot write permissiveness_global_allcells.csv: ", e$message), logfile))

# ---- Pass 2: sample-level figures ----
summary_rows <- list()
process_log[[length(process_log) + 1L]] <- paste0("Pass2 start: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

for (ds in names(objects)) {
  obj <- objects[[ds]]
  if (!inherits(obj, "Seurat")) {
    process_log[[length(process_log) + 1L]] <- paste0("Skipping non-Seurat object for dataset: ", ds); next
  }
  
  # Ensure md carries dataset name for robust joins later
  md <- obj@meta.data %>% rownames_to_column("cell") %>% mutate(dataset = ds)
  
  sample_col <- pick_sample_col(md)
  samples <- if (is.null(sample_col)) "all" else unique(as.character(md[[sample_col]]))
  
  for (sp in samples) {
    md_s <- if (is.null(sample_col)) {
      md
    } else {
      md %>% filter(as.character(.data[[sample_col]]) == sp)
    }
    if (nrow(md_s) == 0) {
      process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] no cells -> skip")
      next
    }
    
    # left_join with pooled all_md; but ensure join keys present
    if (!("dataset" %in% colnames(md_s)) || !("cell" %in% colnames(md_s))) {
      process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] missing dataset/cell columns; skipping join")
      next
    }
    
    # perform join; if no matching rows, warn but continue
    joined <- tryCatch({
      md_s %>%
        left_join(all_md %>% select(dataset, cell, permissiveness_global, permissiveness_global_pct),
                  by = c("dataset" = "dataset", "cell" = "cell"))
    }, error = function(e) {
      process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] left_join failed: ", e$message)
      NULL
    })
    if (is.null(joined)) {
      process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] join returned NULL -> skipping sample")
      next
    }
    md_s <- joined
    
    # If permissiveness_global didn't join (all NA), warn and skip plotting
    if (!("permissiveness_global" %in% colnames(md_s)) || sum(is.finite(md_s$permissiveness_global)) < 2) {
      process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] insufficient permissiveness_global after join -> skipping sample")
      next
    }
    
    md_s <- md_s %>%
      mutate(
        x = suppressWarnings(as.numeric(spatial_x_use)),
        y = suppressWarnings(as.numeric(spatial_y_use))
      ) %>%
      filter(is.finite(x), is.finite(y), is.finite(permissiveness_global))
    
    if (nrow(md_s) == 0) {
      process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] skipped: no finite x/y/permissiveness_global")
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
    
    # white->red overlay using symmetric limits (perm_maxabs) for cross-sample comparability
    if (is.finite(perm_maxabs) && perm_maxabs > 0) {
      p_overlay_white_red <- p_base +
        geom_point(aes(color = permissiveness_global), size = 0.55, alpha = 0.9) +
        scale_color_gradient2(low = "white", mid = "#ffd6c7", high = "red", midpoint = 0,
                              limits = c(-perm_maxabs, perm_maxabs), oob = scales::squish) +
        ggtitle(sprintf("%s / %s: Global permissiveness (white->red; symmetric limits)", ds, sp))
    } else {
      p_overlay_white_red <- p_base +
        geom_point(aes(color = permissiveness_global), size = 0.55, alpha = 0.9) +
        scale_color_gradient2(low = "white", mid = "#ffd6c7", high = "red", midpoint = 0) +
        ggtitle(sprintf("%s / %s: Global permissiveness (white->red)", ds, sp))
    }
    
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
    } else {
      p_top <- p_top + ggtitle(paste0("Too few points for top density (need >= ", MIN_POINTS_FOR_KDE, ")"))
    }
    
    p_bottom <- p_base + ggtitle(sprintf("%s / %s: Bottom %.0f%% (protected) density", ds, sp, BOTTOM_FRAC * 100))
    if (!is.null(kd_bottom)) {
      p_bottom <- p_bottom + geom_raster(data = kd_bottom, aes(x = x, y = y, fill = z), interpolate = TRUE) +
        scale_fill_viridis_c(option = "plasma")
    } else {
      p_bottom <- p_bottom + ggtitle(paste0("Too few points for bottom density (need >= ", MIN_POINTS_FOR_KDE, ")"))
    }
    
    p_diff <- p_base + ggtitle(sprintf("%s / %s: Hotspot - Protected density", ds, sp))
    if (!is.null(kd_top) && !is.null(kd_bottom)) {
      kd_diff <- kd_top %>% select(x, y, z_top = z) %>%
        left_join(kd_bottom %>% select(x, y, z_bottom = z), by = c("x", "y")) %>%
        mutate(z_diff = z_top - z_bottom)
      lim <- max(abs(kd_diff$z_diff), na.rm = TRUE)
      if (!is.finite(lim) || lim <= 0) lim <- 1
      # Use kd_diff-specific lim if present, otherwise fall back to perm_maxabs
      if (is.finite(perm_maxabs) && perm_maxabs > 0) {
        lim_use <- max(lim, perm_maxabs/10, na.rm = TRUE)
      } else {
        lim_use <- lim
      }
      p_diff <- p_diff +
        geom_raster(data = kd_diff, aes(x = x, y = y, fill = z_diff), interpolate = TRUE) +
        scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                             limits = c(-lim_use, lim_use), oob = scales::squish)
    } else {
      p_diff <- p_diff + ggtitle("Diff not available (top or bottom density missing)")
    }
    
    # Save panels with try-catch so script continues on intermittent file errors
    tryCatch({
      ggsave(file.path(fig_dir, paste0(ds, "_", sp, "_3panel_global.png")),
             (p_overlay_global + p_top + p_overlay_hl) + plot_layout(nrow = 1), width = 18, height = 6, dpi = 300)
    }, error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing _3panel_global.png: ", e$message))
    
    tryCatch({
      ggsave(file.path(fig_dir, paste0(ds, "_", sp, "_3panel_global_white_red.png")),
             (p_overlay_white_red + p_top + p_overlay_hl) + plot_layout(nrow = 1), width = 18, height = 6, dpi = 300)
    }, error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing _3panel_global_white_red.png: ", e$message))
    
    tryCatch({
      ggsave(file.path(fig_dir, paste0(ds, "_", sp, "_protected_and_diff.png")),
             (p_top + p_bottom + p_diff) + plot_layout(nrow = 1), width = 18, height = 6, dpi = 300)
    }, error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing _protected_and_diff.png: ", e$message))
    
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
        tryCatch({
          ggsave(file.path(fig_dir, paste0(ds, "_", sp, "_permissiveness_global_week_", w, ".png")), p_w, width = 10, height = 8, dpi = 300)
        }, error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing weekly png: ", e$message))
      }
    }
    
    # Export hotspot/protected cell lists for downstream DE/CellChat.
    top_thr <- tryCatch(stats::quantile(md_s$permissiveness_global, probs = 1 - TOP_FRAC, na.rm = TRUE), error = function(e) NA_real_)
    bot_thr <- tryCatch(stats::quantile(md_s$permissiveness_global, probs = BOTTOM_FRAC, na.rm = TRUE), error = function(e) NA_real_)
    hotspots <- if (!is.na(top_thr)) md_s %>% filter(permissiveness_global >= top_thr) else md_s[0,]
    protected <- if (!is.na(bot_thr)) md_s %>% filter(permissiveness_global <= bot_thr) else md_s[0,]
    tryCatch(write.csv(hotspots, file.path(tab_dir, paste0(ds, "_", sp, "_hotspot_cells_top", as.integer(TOP_FRAC * 100), "pct.csv")), row.names = FALSE),
             error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing hotspots CSV: ", e$message))
    tryCatch(write.csv(protected, file.path(tab_dir, paste0(ds, "_", sp, "_protected_cells_bottom", as.integer(BOTTOM_FRAC * 100), "pct.csv")), row.names = FALSE),
             error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing protected CSV: ", e$message))
    
    tryCatch(write.csv(md_s, file.path(tab_dir, paste0(ds, "_", sp, "_mapped_coords_global.csv")), row.names = FALSE),
             error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing mapped coords CSV: ", e$message))
    
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
    tryCatch(write_json(meta, file.path(tab_dir, "meta.json"), pretty = TRUE, auto_unbox = TRUE),
             error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing meta.json: ", e$message))
    
    readme <- c(
      paste0("Dataset: ", ds),
      paste0("Sample: ", sp),
      "This folder contains global-calibrated spatial permissiveness outputs.",
      "Global score = z_global(MMP) + z_global(Immune_Tolerance) - z_global(Cytotoxic_NK) + z_global(Ethanolamine).",
      "protected_and_diff figure includes bottom-percentile protected density and hotspot-protected difference.",
      paste0("Top fraction: ", TOP_FRAC, "; Bottom fraction: ", BOTTOM_FRAC)
    )
    tryCatch(writeLines(readme, file.path(tab_dir, "README_sample.txt")),
             error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing README: ", e$message))
    
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      dataset = ds,
      sample = as.character(sp),
      n_cells = nrow(md_s),
      highlight_labels = paste(hl, collapse = ";"),
      fig_dir = fig_dir,
      tab_dir = tab_dir,
      stringsAsFactors = FALSE
    )
    
    process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] done (cells=", nrow(md_s), ")")
  }
}

# summary + logs
if (length(summary_rows) > 0) {
  summary_df <- bind_rows(summary_rows)
  tryCatch(write.csv(summary_df, file.path(DIR_TAB, "summary_table_with_highlights.csv"), row.names = FALSE),
           error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("Error writing summary_table_with_highlights.csv: ", e$message))
  log_msg2(sprintf("[05H] Summary written: %d rows.", nrow(summary_df)), logfile)
} else {
  log_msg2("[05H] No sample-level outputs produced.", logfile)
}

# PROCESS LOG
proc_file <- file.path(DIR_TAB, "PROCESS_LOG.txt")
tryCatch(writeLines(unlist(c(process_log)), proc_file), error = function(e) log_msg2(paste0("Failed to write PROCESS_LOG: ", e$message), logfile))
log_msg2("[05H] Done.", logfile)

# ---- Write description files for methods & interpretation ----
# Scripts descriptions (methods notes)
desc_dir <- "scripts/05_spatial"
tryCatch({
  if (!dir.exists(desc_dir)) dir.create(desc_dir, recursive = TRUE, showWarnings = FALSE)
  f05f <- file.path(desc_dir, "05F_description.txt")
  f05h <- file.path(desc_dir, "05H_description.txt")
  fig_readme <- file.path(DIR_FIG, "README_plot_interpretation.txt")
  
  txt_05f <- c(
    "05F — Spatial Permissiveness: purpose and usage",
    "-----------------------------------------------",
    "",
    "Purpose",
    "-------",
    "Step 05F (companion to 05C/05H) documents the conceptual and practical",
    "purpose of the spatial permissiveness workflow.",
    "",
    "05C computes per-dataset 'permissiveness' by z-scoring 4 modules within each",
    "dataset and combining: permissiveness = z(MMP/ECM remodeling) + z(Immune tolerance)",
    "- z(NK cytotoxic) + z(Ethanolamine metabolism).",
    "",
    "05H performs cross-sample global calibration (pooled z-scores across datasets),",
    "derives 'permissiveness_global', computes hotspot (top N%) and protected",
    "(bottom N%) kernel densities, and produces per-sample figures (vanilla,",
    "white→red symmetric, and highlight variants) for interpretation and downstream",
    "analyses.",
    "",
    "Inputs",
    "------",
    "- output/objects/*_with_permissiveness.rds (preferred), or *_harmonized.rds",
    "- output/tables/*_permissiveness_cell_level.csv (optional but helpful)",
    "- raw images in data/raw/ (image_map.csv optionally maps samples->images)",
    "",
    "Outputs",
    "-------",
    "- output/figures/05_spatial/<dataset>/<sample>/*  (3-panel, white->red, protected/diff, weeks)",
    "- output/tables/05_spatial/<dataset>/<sample>/*  (mapped coords, hotspots, protected lists, meta.json)",
    "- output/tables/05_spatial/permissiveness_global_allcells.csv",
    "- output/tables/05_spatial/permissiveness_global_stats.json",
    "- output/tables/05_spatial/PROCESS_LOG.txt",
    "",
    "Interpretation notes",
    "--------------------",
    "- 05C scores are dataset-relative (z within dataset). 05H provides a globally",
    "  calibrated view (permissiveness_global) that is comparable across samples.",
    "- White→red panels use a symmetric color scale centered at zero for comparability:",
    "  +/- perm_maxabs (stored in permissiveness_global_stats.json).",
    "- Hotspots = kernel-density of top N% permissive cells (default 10%); Protected =",
    "  kernel-density of bottom N% (default 10%). Hotspot - Protected = difference map.",
    "- Use hotspot cell lists for DE, pathway analysis, and CellChat.",
    "",
    "Caveats",
    "-------",
    "- If a dataset lacks enough genes for some modules (STARmap or targeted assays),",
    "  NK or MMP contributions may be weak — inspect NK_module_coverage_qc.csv.",
    "- Global calibration assumes module scores are computed on compatible assays or",
    "  on an integrated assay. If module scores are heterogeneous (different assays),",
    "  interpret global comparisons with caution."
  )
  writeLines(txt_05f, f05f)
  
  txt_05h <- c(
    "05H — Global-calibrated spatial permissiveness panels (README)",
    "---------------------------------------------------------------",
    "",
    "What 05H does",
    "-------------",
    "1. Reads harmonized Seurat objects (preferably *_with_permissiveness.rds).",
    "2. Pools per-cell module scores across datasets and computes global z-scores:",
    "   z_mmp_global, z_tol_global, z_nk_global, z_ea_global.",
    "3. Calculates permissiveness_global = z_mmp_global + z_tol_global - z_nk_global + z_ea_global.",
    "4. Exports a pooled table (permissiveness_global_allcells.csv) and stats file with",
    "   perm_maxabs (permissiveness_global_stats.json).",
    "5. For each dataset and sample:",
    "   - Maps coordinates to image pixels using heuristics,",
    "   - Draws three-panel figures:",
    "       A) Vanilla: per-cell permissiveness_global (viridis)",
    "       B) Density: top-permissive kernel density (top N%)",
    "       C) Highlights: permissiveness + highlighted cell types",
    "     and a white->red variant with symmetric color limits.",
    "   - Creates protected (bottom N%) density maps and a Hotspot-Protected difference map.",
    "6. Writes per-sample files: mapped_coords_global.csv, hotspot/protected lists,",
    "   meta.json, README_sample.txt, and figure PNG/PDFs.",
    "",
    "Why symmetric color limits?",
    "---------------------------",
    "We compute a pooled `perm_maxabs = max(abs(permissiveness_global))`. Using",
    "symmetric limits (`[-perm_maxabs, +perm_maxabs]`) centered at 0 allows direct",
    "visual comparison between samples: a red pixel of the same intensity means",
    "the same relative permissiveness across samples.",
    "",
    "How to interpret each panel",
    "---------------------------",
    "- Panel A (vanilla): per-cell score (viridis). Use for local cell-level patterns.",
    "- Panel B (density): spatial concentration of top permissive cells (hotspots).",
    "  Useful for locating regions to analyze (DE genes, neighborhood composition).",
    "- Panel C (highlights): overlay of highlighted biologically-important cell types",
    "  (e.g. EVT, NK, Endothelial) on permissiveness — helps connect biology to score.",
    "- White->Red: same cell-level data but with symmetric red ramp — good for cross-sample comparisons.",
    "- Protected and Diff: regions enriched for bottom-permissiveness (protected) and Hotspot - Protected difference.",
    "",
    "Suggested next analyses",
    "-----------------------",
    "- Differential expression between hotspot cells and background cells (Seurat FindMarkers).",
    "- CellChat / ligand-receptor analysis on hotspot cells vs background cells.",
    "- Spatially-aware tests (SPARK, spatialDE) focusing on hotspot regions."
  )
  writeLines(txt_05h, f05h)
  
  txt_figreadme <- c(
    "Plot interpretation guide (short)",
    "--------------------------------",
    "",
    "Files produced per sample:",
    "- <dataset>_<sample>_3panel_global.png : (Viridis) Global permissiveness, hotspot density, highlights.",
    "- <dataset>_<sample>_3panel_global_white_red.png : same but white->red symmetric color scale for cross-sample comparability.",
    "- <dataset>_<sample>_protected_and_diff.png : Top density, bottom(protected) density, and difference (hotspot - protected).",
    "- <dataset>_<sample>_permissiveness_global_week_<w>.png : per-week global permissiveness maps.",
    "",
    "Key points:",
    "- Permissiveness_global is pooled across datasets and comparable across samples.",
    "- Hotspots = highest-permissive cells (top 10% by default). Protected = lowest 10%.",
    "- Difference map (hotspot - protected) highlights spatially differential permissiveness.",
    "- Use hotspot CSVs for downstream DE/CellChat.",
    "- Check meta.json and permissiveness_global_stats.json for provenance and the perm_maxabs used for symmetric color scales.",
    "",
    "Caveats:",
    "- If module scores are missing (e.g., due to limited gene panels), permissiveness contributions may be NA; check NK_module_coverage_qc files.",
    "- Interpret differences across platforms (STARmap vs SlideTags) carefully even after global calibration; module computation assumptions matter."
  )
  writeLines(txt_figreadme, fig_readme)
  
  log_msg2("[05H] Wrote descriptive README files in scripts/05_spatial and output/figures/05_spatial", logfile)
=======
# ======================================================================
# scripts/05_spatial/05H_spatial_permissiveness_panels_global.R
# Build per-sample spatial permissiveness panel figures from 05C outputs,
# including global cross-sample calibration and protected-region maps.
#
# Run from repo root:
# source("scripts/05_spatial/05H_spatial_permissiveness_panels_global.R")
# ======================================================================
# scripts/05_spatial/05H_spatial_permissiveness_panels_global.R
# Build per-sample spatial permissiveness panel figures from 05C outputs,
# including global cross-sample calibration and protected-region maps.
#
# Run from repo root:
# source("scripts/05_spatial/05H_spatial_permissiveness_panels_global.R")

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
TOP_FRAC <- 0.10          # top X% = hotspot
BOTTOM_FRAC <- 0.10       # bottom X% = protected
MIN_POINTS_FOR_KDE <- 30L
KDE_GRID_N <- 250L
BIO_IMPORTANT <- c("EVT", "Endothelial_cells", "Macrophage", "NK", "Fibroblast", "vCTB", "Hofbauer cells")
TOP_N_DYNAMIC_LABELS <- 3L

DIR_FIG <- file.path(DIR_FIGURES, "05_spatial")
DIR_TAB <- file.path(DIR_TABLES, "05_spatial")
ensure_dir(DIR_FIG); ensure_dir(DIR_TAB); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "05H_spatial_permissiveness_panels_global.log")

# small logger
log_msg2 <- function(msg, logfile = logfile) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  tryCatch(cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg), file = logfile, append = TRUE),
           error = function(e) NULL)
}

plot_bg_points <- function(df, col = "grey80", size = 0.35, alpha = 0.45) {
  ggplot2::geom_point(data = df, aes(x = x, y = y), color = col, size = size, alpha = alpha)
}

# robust/global z
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
  xr <- range(d$x, na.rm = TRUE); yr <- range(d$y, na.rm = TRUE)
  kd <- tryCatch(MASS::kde2d(sel$x, sel$y, n = n_grid, lims = c(xr[1], xr[2], yr[1], yr[2])),
                 error = function(e) NULL)
  if (is.null(kd)) return(NULL)
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

# New helper used in script
pick_first_present <- function(df_or_cols, cols) {
  if (is.data.frame(df_or_cols)) {
    avail <- intersect(cols, colnames(df_or_cols))
  } else {
    avail <- intersect(cols, df_or_cols)
  }
  if (length(avail) == 0) return(NULL)
  avail[[1]]
}

# Helper: write per-figure caption text
write_caption <- function(fig_path, caption_text) {
  cap_path <- sub("\\.(png|pdf)$", "_caption.txt", fig_path, ignore.case = TRUE)
  tryCatch(writeLines(caption_text, cap_path), error = function(e) NULL)
}

# Create a repo-level description (05F + 05H) to place in scripts/05_spatial
write_repo_description <- function(out_path) {
  txt <- c(
    "# 05F and 05H description",
    "",
    "This file describes the outputs and interpretation guidelines for:",
    "- scripts/05_spatial/05C_permissiveness_score_maps.R (05C)",
    "- scripts/05_spatial/05H_spatial_permissiveness_panels_global.R (05H)",
    "",
    "### Purpose",
    "05C: compute dataset-relative permissiveness per cell: permissiveness = z_MMP + z_Tolerance - z_NK + z_Ethanolamine (z computed within each dataset).",
    "05H: pooled global calibration across datasets to allow cross-sample/time comparisons: z_global computed across pooled cells.",
    "",
    "### Figure types created by 05H (per-sample)",
    "1. 3-panel global: (A) image/background, (B) KDE of top 10% permissive cells, (C) global-permissiveness overlay (viridis).",
    "2. 3-panel white->red: (A) image, (B) same KDE, (C) white->red overlay emphasizing positive permissiveness.",
    "3. Cell-type highlighted variant: (A) image, (B) same KDE, (C) overlay with highlighted biologically-important cell types (filled circles) on top of permissiveness.",
    "4. Protected/diff panels: (A) top-density, (B) bottom-density (protected), (C) hotspot - protected difference map.",
    "",
    "### Interpretation guidance (for figure legends / methods)",
    "- Hotspot defined as KDE of the top 10% permissive cells (dataset-pooled percentile when using permissiveness_global).",
    "- Protected regions: KDE of the bottom 10% permissive cells.",
    "- Difference map = hotspot_density - protected_density (positive = permissive, negative = protected).",
    "- Include a note whether permissiveness is dataset-relative (05C) or globally-calibrated (05H).",
    "",
    "### Files produced",
    "- output/figures/05_spatial/<dataset>/<sample>/ : PNG/PDF figures and caption .txt files",
    "- output/tables/05_spatial/ : permissiveness_global_allcells.csv, summary_table_with_highlights.csv, PROCESS_LOG.txt, per-sample hotspot/protected CSVs and meta.json files",
    "",
    "### Methods notes",
    "See the Methods bullets saved to scripts/05_spatial/METHODS_05C_05H.md for full reproducible details."
  )
  tryCatch(writeLines(txt, out_path), error = function(e) NULL)
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

log_msg2(sprintf("[05H] Found %d object candidates.", length(obj_paths)), logfile)

# Pass 1: collect all cells for global comparability -------------------------
all_cells <- list()
objects <- list()
process_log <- list()
process_log[[length(process_log) + 1L]] <- paste0("Run at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

for (pth in obj_paths) {
  obj <- tryCatch(readRDS(pth), error = function(e) { log_msg2(paste0("  ERROR reading ", pth, ": ", e$message)); NULL })
  if (is.null(obj) || !inherits(obj, "Seurat")) {
    process_log[[length(process_log) + 1L]] <- paste0("Skipping non-Seurat or unreadable: ", pth); next
  }
  ds <- pick_dataset_name(pth)
  
  # ensure week + spatial columns present (try/catch)
  obj <- tryCatch({
    obj2 <- ensure_week_column(obj, COL_WEEK_CANDIDATES)
    obj2 <- ensure_spatial_coords(obj2, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)
    obj2
  }, error = function(e) { log_msg2(paste0("  WARNING: ensure_week/coords failed for ", pth, ": ", e$message)); obj })
  
  # md from object; ensure dataset column
  md <- obj@meta.data %>% rownames_to_column("cell") %>% mutate(dataset = ds)
  
  if (!("permissiveness" %in% colnames(md))) {
    process_log[[length(process_log) + 1L]] <- paste0("[", ds, "] skipping: no 'permissiveness' column in metadata")
    objects[[ds]] <- obj  # keep object for plotting if needed
    next
  }
  
  # Safely collect score columns; warn if missing
  required_scores <- c("score_MMP_ECM_Remodeling", "score_Immune_Tolerance", "score_Cytotoxic_NK", "score_Ethanolamine_Metabolism")
  missing_scores <- setdiff(required_scores, colnames(md))
  if (length(missing_scores) > 0) {
    process_log[[length(process_log) + 1L]] <- paste0("[", ds, "] WARNING: missing score columns: ", paste(missing_scores, collapse = "; "))
    for (mc in missing_scores) md[[mc]] <- NA_real_
  }
  
  # push to all_cells
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
  process_log[[length(process_log) + 1L]] <- paste0("[", ds, "] queued (n_cells=", nrow(md), ")")
}

if (length(all_cells) == 0) {
  log_msg2("[05H] No permissiveness metadata found in candidate objects. Run 05C first.", logfile)
  writeLines(unlist(process_log), file.path(DIR_TAB, "PROCESS_LOG.txt"))
  stop("[05H] No permissiveness metadata found in candidate objects. Run 05C first.")
}

all_md <- bind_rows(all_cells)

# Check required score columns exist in pooled table before z-scoring
required_scores <- c("score_MMP_ECM_Remodeling", "score_Immune_Tolerance", "score_Cytotoxic_NK", "score_Ethanolamine_Metabolism")
if (!all(required_scores %in% colnames(all_md))) {
  missing_scores <- setdiff(required_scores, colnames(all_md))
  log_msg2(paste0("[05H] WARNING: pooled table missing score columns: ", paste(missing_scores, collapse = "; "), " - those components will be NA"), logfile)
}

# compute global z and permissiveness (NA-safe)
all_md <- all_md %>% mutate(
  z_mmp_global = if ("score_MMP_ECM_Remodeling" %in% colnames(all_md)) z_global(score_MMP_ECM_Remodeling) else NA_real_,
  z_tol_global = if ("score_Immune_Tolerance" %in% colnames(all_md)) z_global(score_Immune_Tolerance) else NA_real_,
  z_nk_global  = if ("score_Cytotoxic_NK" %in% colnames(all_md)) z_global(score_Cytotoxic_NK) else NA_real_,
  z_ea_global  = if ("score_Ethanolamine_Metabolism" %in% colnames(all_md)) z_global(score_Ethanolamine_Metabolism) else NA_real_,
  permissiveness_global = rowSums(cbind(z_mmp_global, z_tol_global, -z_nk_global, z_ea_global), na.rm = FALSE),
  permissiveness_global_pct = safe_ecdf(permissiveness_global)
)

# If permissiveness_global is NA for many rows, warn
if (sum(is.finite(all_md$permissiveness_global)) < 2) {
  log_msg2("[05H] WARNING: permissiveness_global could not be computed (not enough finite components).", logfile)
}

# Write pooled table
tryCatch({
  write.csv(all_md, file.path(DIR_TAB, "permissiveness_global_allcells.csv"), row.names = FALSE)
  log_msg2(sprintf("[05H] Wrote global calibration table with %d cells.", nrow(all_md)), logfile)
}, error = function(e) log_msg2(paste0("[05H] Warning: cannot write permissiveness_global_allcells.csv: ", e$message), logfile))

# Create repo-level description file in scripts/05_spatial
tryCatch({
  ensure_dir(file.path("scripts", "05_spatial"))
  write_repo_description(file.path("scripts", "05_spatial", "05F_and_05H_description.md"))
}, error = function(e) NULL)

# Pass 2: sample-level figures -----------------------------------------------
summary_rows <- list()
process_log[[length(process_log) + 1L]] <- paste0("Pass2 start: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

for (ds in names(objects)) {
  obj <- objects[[ds]]
  if (!inherits(obj, "Seurat")) {
    process_log[[length(process_log) + 1L]] <- paste0("Skipping non-Seurat object for dataset: ", ds); next
  }
  
  # Ensure md carries dataset name for robust joins later
  md <- obj@meta.data %>% rownames_to_column("cell") %>% mutate(dataset = ds)
  
  sample_col <- pick_sample_col(md)
  samples <- if (is.null(sample_col)) "all" else unique(as.character(md[[sample_col]]))
  
  for (sp in samples) {
    md_s <- if (is.null(sample_col)) { md } else { md %>% filter(as.character(.data[[sample_col]]) == sp) }
    if (nrow(md_s) == 0) {
      process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] no cells -> skip")
      next
    }
    
    # left_join with pooled all_md; ensure join keys present
    if (!("dataset" %in% colnames(md_s)) || !("cell" %in% colnames(md_s))) {
      process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] missing dataset/cell columns; skipping join")
      next
    }
    
    joined <- tryCatch({
      md_s %>%
        left_join(all_md %>% select(dataset, cell, permissiveness_global, permissiveness_global_pct),
                  by = c("dataset" = "dataset", "cell" = "cell"))
    }, error = function(e) {
      process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] left_join failed: ", e$message)
      NULL
    })
    if (is.null(joined)) {
      process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] join returned NULL -> skipping sample")
      next
    }
    md_s <- joined
    
    if (!("permissiveness_global" %in% colnames(md_s)) || sum(is.finite(md_s$permissiveness_global)) < 2) {
      process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] insufficient permissiveness_global after join -> skipping sample")
      next
    }
    
    md_s <- md_s %>%
      mutate(x = suppressWarnings(as.numeric(spatial_x_use)),
             y = suppressWarnings(as.numeric(spatial_y_use))) %>%
      filter(is.finite(x), is.finite(y), is.finite(permissiveness_global))
    
    if (nrow(md_s) == 0) {
      process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] skipped: no finite x/y/permissiveness_global")
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
      theme_classic() + coord_fixed() + labs(x = "spatial x", y = "spatial y")
    
    # overlay variants
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
      # color palette for labels (Set2)
      n_hl <- length(hl)
      pal <- RColorBrewer::brewer.pal(max(3, min(8, n_hl)), "Set2")
      fill_vals <- setNames(rep(pal, length.out = n_hl), hl)
      p_overlay_hl <- p_overlay_hl +
        geom_point(data = md_s %>% filter(label %in% hl),
                   aes(fill = label), shape = 21, size = 2.4, color = "black", stroke = 0.5, alpha = 0.9) +
        scale_fill_manual(name = "Highlighted cell types", values = fill_vals)
    }
    
    kd_top <- safe_kde_df(md_s, frac = TOP_FRAC, mode = "top", n_grid = KDE_GRID_N)
    kd_bottom <- safe_kde_df(md_s, frac = BOTTOM_FRAC, mode = "bottom", n_grid = KDE_GRID_N)
    
    p_top <- p_base + ggtitle(sprintf("%s / %s: Top %.0f%% permissive density", ds, sp, TOP_FRAC * 100))
    if (!is.null(kd_top)) {
      p_top <- p_top + geom_raster(data = kd_top, aes(x = x, y = y, fill = z), interpolate = TRUE) +
        scale_fill_viridis_c(option = "magma")
    } else {
      p_top <- p_top + ggtitle(paste0("Too few points for top density (need >= ", MIN_POINTS_FOR_KDE, ")"))
    }
    
    p_bottom <- p_base + ggtitle(sprintf("%s / %s: Bottom %.0f%% (protected) density", ds, sp, BOTTOM_FRAC * 100))
    if (!is.null(kd_bottom)) {
      p_bottom <- p_bottom + geom_raster(data = kd_bottom, aes(x = x, y = y, fill = z), interpolate = TRUE) +
        scale_fill_viridis_c(option = "plasma")
    } else {
      p_bottom <- p_bottom + ggtitle(paste0("Too few points for bottom density (need >= ", MIN_POINTS_FOR_KDE, ")"))
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
    } else {
      p_diff <- p_diff + ggtitle("Diff not available (top or bottom density missing)")
    }
    
    # figure paths
    f3_global <- file.path(fig_dir, paste0(ds, "_", sp, "_3panel_global.png"))
    f3_white_red <- file.path(fig_dir, paste0(ds, "_", sp, "_3panel_global_white_red.png"))
    f_prot_diff <- file.path(fig_dir, paste0(ds, "_", sp, "_protected_and_diff.png"))
    
    # Save panels (with try-catch)
    tryCatch({
      ggsave(f3_global, (p_overlay_global + p_top + p_overlay_hl) + plot_layout(nrow = 1),
             width = 18, height = 6, dpi = 300)
      # caption
      write_caption(f3_global,
                    c(sprintf("3-panel global (dataset=%s sample=%s)", ds, sp),
                      "- Panel A: original sample coordinates (background); Panel B: KDE of top 10% permissive cells (magma); Panel C: per-cell global permissiveness (viridis).",
                      "- Interpretation: hotspots in Panel B correspond to clusters where permissiveness is concentrated; Panel C shows per-cell global-level permissiveness."))
    }, error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing _3panel_global.png: ", e$message))
    
    tryCatch({
      ggsave(f3_white_red, (p_overlay_white_red + p_top + p_overlay_hl) + plot_layout(nrow = 1),
             width = 18, height = 6, dpi = 300)
      write_caption(f3_white_red,
                    c(sprintf("3-panel white->red (dataset=%s sample=%s)", ds, sp),
                      "- Panel C uses white->red to emphasize positive permissiveness (centered at zero).",
                      "- Use this panel to emphasize regions with strongly positive permissiveness."))
    }, error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing _3panel_global_white_red.png: ", e$message))
    
    tryCatch({
      ggsave(f_prot_diff, (p_top + p_bottom + p_diff) + plot_layout(nrow = 1),
             width = 18, height = 6, dpi = 300)
      write_caption(f_prot_diff,
                    c(sprintf("Protected and difference maps (dataset=%s sample=%s)", ds, sp),
                      "- Left: KDE of top 10% (hotspot). Middle: KDE of bottom 10% (protected). Right: hotspot - protected difference (blue negative = protected; red positive = permissive)."))
    }, error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing _protected_and_diff.png: ", e$message))
    
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
        out_w <- file.path(fig_dir, paste0(ds, "_", sp, "_permissiveness_global_week_", w, ".png"))
        tryCatch({
          ggsave(out_w, p_w, width = 10, height = 8, dpi = 300)
          write_caption(out_w, c(sprintf("Global permissiveness (dataset=%s sample=%s week=%s)", ds, sp, w),
                                 "- Interpret relative to other weeks/samples only when using permissiveness_global (pooled calibration)."))
        }, error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing weekly png: ", e$message))
      }
    }
    
    # Export hotspot/protected cell lists
    top_thr <- tryCatch(stats::quantile(md_s$permissiveness_global, probs = 1 - TOP_FRAC, na.rm = TRUE), error = function(e) NA_real_)
    bot_thr <- tryCatch(stats::quantile(md_s$permissiveness_global, probs = BOTTOM_FRAC, na.rm = TRUE), error = function(e) NA_real_)
    hotspots <- if (!is.na(top_thr)) md_s %>% filter(permissiveness_global >= top_thr) else md_s[0,]
    protected <- if (!is.na(bot_thr)) md_s %>% filter(permissiveness_global <= bot_thr) else md_s[0,]
    tryCatch(write.csv(hotspots, file.path(tab_dir, paste0(ds, "_", sp, "_hotspot_cells_top", as.integer(TOP_FRAC * 100), "pct.csv")), row.names = FALSE),
             error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing hotspots CSV: ", e$message))
    tryCatch(write.csv(protected, file.path(tab_dir, paste0(ds, "_", sp, "_protected_cells_bottom", as.integer(BOTTOM_FRAC * 100), "pct.csv")), row.names = FALSE),
             error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing protected CSV: ", e$message))
    
    tryCatch(write.csv(md_s, file.path(tab_dir, paste0(ds, "_", sp, "_mapped_coords_global.csv")), row.names = FALSE),
             error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing mapped coords CSV: ", e$message))
    
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
    tryCatch(write_json(meta, file.path(tab_dir, "meta.json"), pretty = TRUE, auto_unbox = TRUE),
             error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing meta.json: ", e$message))
    
    readme <- c(
      paste0("Dataset: ", ds),
      paste0("Sample: ", sp),
      "This folder contains global-calibrated spatial permissiveness outputs.",
      "Global score = z_global(MMP) + z_global(Immune_Tolerance) - z_global(Cytotoxic_NK) + z_global(Ethanolamine).",
      "protected_and_diff figure includes bottom-percentile protected density and hotspot-protected difference.",
      paste0("Top fraction: ", TOP_FRAC, "; Bottom fraction: ", BOTTOM_FRAC)
    )
    tryCatch(writeLines(readme, file.path(tab_dir, "README_sample.txt")),
             error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing README: ", e$message))
    
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      dataset = ds,
      sample = as.character(sp),
      n_cells = nrow(md_s),
      highlight_labels = paste(hl, collapse = ";"),
      fig_dir = fig_dir,
      tab_dir = tab_dir,
      stringsAsFactors = FALSE
    )
    
    process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] done (cells=", nrow(md_s), ")")
  }
}

# summary + logs
if (length(summary_rows) > 0) {
  summary_df <- bind_rows(summary_rows)
  tryCatch(write.csv(summary_df, file.path(DIR_TAB, "summary_table_with_highlights.csv"), row.names = FALSE),
           error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("Error writing summary_table_with_highlights.csv: ", e$message))
  log_msg2(sprintf("[05H] Summary written: %d rows.", nrow(summary_df)), logfile)
} else {
  log_msg2("[05H] No sample-level outputs produced.", logfile)
}

# PROCESS LOG
proc_file <- file.path(DIR_TAB, "PROCESS_LOG.txt")
tryCatch(writeLines(unlist(c(process_log)), proc_file), error = function(e) log_msg2(paste0("Failed to write PROCESS_LOG: ", e$message), logfile))
log_msg2("[05H] Done.", logfile)

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

# small logger (prints + appends to logfile)
log_msg2 <- function(msg, logfile = logfile) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  tryCatch(cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg), file = logfile, append = TRUE),
           error = function(e) NULL)
}

# ---- Helpers ----
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
  kd <- tryCatch(MASS::kde2d(sel$x, sel$y, n = n_grid, lims = c(xr[1], xr[2], yr[1], yr[2])),
                 error = function(e) NULL)
  if (is.null(kd)) return(NULL)
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

# New helper used in script
pick_first_present <- function(df_or_cols, cols) {
  if (is.data.frame(df_or_cols)) {
    avail <- intersect(cols, colnames(df_or_cols))
  } else {
    avail <- intersect(cols, df_or_cols)
  }
  if (length(avail) == 0) return(NULL)
  avail[[1]]
}

# ---- Discover objects ----
obj_paths <- c(
  list.files(DIR_OBJS, pattern = "_with_permissiveness\\.rds$", full.names = TRUE),
  list.files(DIR_OBJS, pattern = "_harmonized\\.rds$", full.names = TRUE)
)
obj_paths <- unique(obj_paths)

if (length(obj_paths) == 0) {
  stop("[05H] No *_with_permissiveness.rds or *_harmonized.rds files found in ", DIR_OBJS)
}

log_msg2(sprintf("[05H] Found %d object candidates.", length(obj_paths)), logfile)

# ---- Pass 1: collect all cells for global comparability ----
all_cells <- list()
objects <- list()
process_log <- list()
process_log[[length(process_log) + 1L]] <- paste0("Run at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

for (pth in obj_paths) {
  obj <- tryCatch(readRDS(pth), error = function(e) { log_msg2(paste0("  ERROR reading ", pth, ": ", e$message)); NULL })
  if (is.null(obj) || !inherits(obj, "Seurat")) {
    process_log[[length(process_log) + 1L]] <- paste0("Skipping non-Seurat or unreadable: ", pth); next
  }
  ds <- pick_dataset_name(pth)
  # ensure week + spatial columns present
  obj <- tryCatch({
    obj2 <- ensure_week_column(obj, COL_WEEK_CANDIDATES)
    obj2 <- ensure_spatial_coords(obj2, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)
    obj2
  }, error = function(e) { log_msg2(paste0("  WARNING: ensure_week/coords failed for ", pth, ": ", e$message)); obj })
  
  # md from object; ensure dataset column
  md <- obj@meta.data %>% rownames_to_column("cell") %>% mutate(dataset = ds)
  
  if (!("permissiveness" %in% colnames(md))) {
    process_log[[length(process_log) + 1L]] <- paste0("[", ds, "] skipping: no 'permissiveness' column in metadata")
    objects[[ds]] <- obj  # still keep for plotting (maybe permissiveness propagated later), but don't add cells
    next
  }
  
  # Safely collect columns; if missing, fill NA and warn
  required_scores <- c("score_MMP_ECM_Remodeling", "score_Immune_Tolerance", "score_Cytotoxic_NK", "score_Ethanolamine_Metabolism")
  missing_scores <- setdiff(required_scores, colnames(md))
  if (length(missing_scores) > 0) {
    process_log[[length(process_log) + 1L]] <- paste0("[", ds, "] WARNING: missing score columns: ", paste(missing_scores, collapse = "; "))
    # create NA columns to keep table shape
    for (mc in missing_scores) md[[mc]] <- NA_real_
  }
  
  # push to all_cells
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
  process_log[[length(process_log) + 1L]] <- paste0("[", ds, "] queued (n_cells=", nrow(md), ")")
}

if (length(all_cells) == 0) {
  log_msg2("[05H] No permissiveness metadata found in candidate objects. Run 05C first.", logfile)
  writeLines(unlist(process_log), file.path(DIR_TAB, "PROCESS_LOG.txt"))
  stop("[05H] No permissiveness metadata found in candidate objects. Run 05C first.")
}

all_md <- bind_rows(all_cells)

# Check required score columns exist in pooled table before z-scoring
required_scores <- c("score_MMP_ECM_Remodeling", "score_Immune_Tolerance", "score_Cytotoxic_NK", "score_Ethanolamine_Metabolism")
if (!all(required_scores %in% colnames(all_md))) {
  missing_scores <- setdiff(required_scores, colnames(all_md))
  log_msg2(paste0("[05H] WARNING: pooled table missing score columns: ", paste(missing_scores, collapse = "; "), " - those components will be NA"), logfile)
}

# compute global z and permissiveness (NA-safe)
all_md <- all_md %>% mutate(
  z_mmp_global = if ("score_MMP_ECM_Remodeling" %in% colnames(all_md)) z_global(score_MMP_ECM_Remodeling) else NA_real_,
  z_tol_global = if ("score_Immune_Tolerance" %in% colnames(all_md)) z_global(score_Immune_Tolerance) else NA_real_,
  z_nk_global  = if ("score_Cytotoxic_NK" %in% colnames(all_md)) z_global(score_Cytotoxic_NK) else NA_real_,
  z_ea_global  = if ("score_Ethanolamine_Metabolism" %in% colnames(all_md)) z_global(score_Ethanolamine_Metabolism) else NA_real_,
  permissiveness_global = rowSums(cbind(z_mmp_global, z_tol_global, -z_nk_global, z_ea_global), na.rm = FALSE),
  permissiveness_global_pct = safe_ecdf(permissiveness_global)
)

# If permissiveness_global is NA for many rows, warn
if (sum(is.finite(all_md$permissiveness_global)) < 2) {
  log_msg2("[05H] WARNING: permissiveness_global could not be computed (not enough finite components).", logfile)
}

# === global stats & symmetric color limit for white->red ===
perm_maxabs <- NA_real_
if (exists("all_md") && "permissiveness_global" %in% colnames(all_md)) {
  perm_maxabs <- suppressWarnings(max(abs(all_md$permissiveness_global), na.rm = TRUE))
  if (!is.finite(perm_maxabs)) perm_maxabs <- NA_real_
}
# Save global stats (provenance)
global_stats <- list(
  n_cells_pooled = ifelse(exists("all_md"), nrow(all_md), NA_integer_),
  perm_maxabs = perm_maxabs,
  TOP_FRAC = TOP_FRAC,
  BOTTOM_FRAC = BOTTOM_FRAC
)
# create dir if missing
tryCatch({
  if (!dir.exists(DIR_TAB)) dir.create(DIR_TAB, recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(global_stats, file.path(DIR_TAB, "permissiveness_global_stats.json"), pretty = TRUE, auto_unbox = TRUE)
  log_msg2(sprintf("[05H] Saved permissiveness_global_stats.json (perm_maxabs=%s)", as.character(perm_maxabs)), logfile)
}, error = function(e) log_msg2(paste0("[05H] Warning: failed to write global stats: ", e$message), logfile))

# Write pooled table
tryCatch({
  write.csv(all_md, file.path(DIR_TAB, "permissiveness_global_allcells.csv"), row.names = FALSE)
  log_msg2(sprintf("[05H] Wrote global calibration table with %d cells.", nrow(all_md)), logfile)
}, error = function(e) log_msg2(paste0("[05H] Warning: cannot write permissiveness_global_allcells.csv: ", e$message), logfile))

# ---- Pass 2: sample-level figures ----
summary_rows <- list()
process_log[[length(process_log) + 1L]] <- paste0("Pass2 start: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

for (ds in names(objects)) {
  obj <- objects[[ds]]
  if (!inherits(obj, "Seurat")) {
    process_log[[length(process_log) + 1L]] <- paste0("Skipping non-Seurat object for dataset: ", ds); next
  }
  
  # Ensure md carries dataset name for robust joins later
  md <- obj@meta.data %>% rownames_to_column("cell") %>% mutate(dataset = ds)
  
  sample_col <- pick_sample_col(md)
  samples <- if (is.null(sample_col)) "all" else unique(as.character(md[[sample_col]]))
  
  for (sp in samples) {
    md_s <- if (is.null(sample_col)) {
      md
    } else {
      md %>% filter(as.character(.data[[sample_col]]) == sp)
    }
    if (nrow(md_s) == 0) {
      process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] no cells -> skip")
      next
    }
    
    # left_join with pooled all_md; but ensure join keys present
    if (!("dataset" %in% colnames(md_s)) || !("cell" %in% colnames(md_s))) {
      process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] missing dataset/cell columns; skipping join")
      next
    }
    
    # perform join; if no matching rows, warn but continue
    joined <- tryCatch({
      md_s %>%
        left_join(all_md %>% select(dataset, cell, permissiveness_global, permissiveness_global_pct),
                  by = c("dataset" = "dataset", "cell" = "cell"))
    }, error = function(e) {
      process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] left_join failed: ", e$message)
      NULL
    })
    if (is.null(joined)) {
      process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] join returned NULL -> skipping sample")
      next
    }
    md_s <- joined
    
    # If permissiveness_global didn't join (all NA), warn and skip plotting
    if (!("permissiveness_global" %in% colnames(md_s)) || sum(is.finite(md_s$permissiveness_global)) < 2) {
      process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] insufficient permissiveness_global after join -> skipping sample")
      next
    }
    
    md_s <- md_s %>%
      mutate(
        x = suppressWarnings(as.numeric(spatial_x_use)),
        y = suppressWarnings(as.numeric(spatial_y_use))
      ) %>%
      filter(is.finite(x), is.finite(y), is.finite(permissiveness_global))
    
    if (nrow(md_s) == 0) {
      process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] skipped: no finite x/y/permissiveness_global")
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
    
    # white->red overlay using symmetric limits (perm_maxabs) for cross-sample comparability
    if (is.finite(perm_maxabs) && perm_maxabs > 0) {
      p_overlay_white_red <- p_base +
        geom_point(aes(color = permissiveness_global), size = 0.55, alpha = 0.9) +
        scale_color_gradient2(low = "white", mid = "#ffd6c7", high = "red", midpoint = 0,
                              limits = c(-perm_maxabs, perm_maxabs), oob = scales::squish) +
        ggtitle(sprintf("%s / %s: Global permissiveness (white->red; symmetric limits)", ds, sp))
    } else {
      p_overlay_white_red <- p_base +
        geom_point(aes(color = permissiveness_global), size = 0.55, alpha = 0.9) +
        scale_color_gradient2(low = "white", mid = "#ffd6c7", high = "red", midpoint = 0) +
        ggtitle(sprintf("%s / %s: Global permissiveness (white->red)", ds, sp))
    }
    
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
    } else {
      p_top <- p_top + ggtitle(paste0("Too few points for top density (need >= ", MIN_POINTS_FOR_KDE, ")"))
    }
    
    p_bottom <- p_base + ggtitle(sprintf("%s / %s: Bottom %.0f%% (protected) density", ds, sp, BOTTOM_FRAC * 100))
    if (!is.null(kd_bottom)) {
      p_bottom <- p_bottom + geom_raster(data = kd_bottom, aes(x = x, y = y, fill = z), interpolate = TRUE) +
        scale_fill_viridis_c(option = "plasma")
    } else {
      p_bottom <- p_bottom + ggtitle(paste0("Too few points for bottom density (need >= ", MIN_POINTS_FOR_KDE, ")"))
    }
    
    p_diff <- p_base + ggtitle(sprintf("%s / %s: Hotspot - Protected density", ds, sp))
    if (!is.null(kd_top) && !is.null(kd_bottom)) {
      kd_diff <- kd_top %>% select(x, y, z_top = z) %>%
        left_join(kd_bottom %>% select(x, y, z_bottom = z), by = c("x", "y")) %>%
        mutate(z_diff = z_top - z_bottom)
      lim <- max(abs(kd_diff$z_diff), na.rm = TRUE)
      if (!is.finite(lim) || lim <= 0) lim <- 1
      # Use kd_diff-specific lim if present, otherwise fall back to perm_maxabs
      if (is.finite(perm_maxabs) && perm_maxabs > 0) {
        lim_use <- max(lim, perm_maxabs/10, na.rm = TRUE)
      } else {
        lim_use <- lim
      }
      p_diff <- p_diff +
        geom_raster(data = kd_diff, aes(x = x, y = y, fill = z_diff), interpolate = TRUE) +
        scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                             limits = c(-lim_use, lim_use), oob = scales::squish)
    } else {
      p_diff <- p_diff + ggtitle("Diff not available (top or bottom density missing)")
    }
    
    # Save panels with try-catch so script continues on intermittent file errors
    tryCatch({
      ggsave(file.path(fig_dir, paste0(ds, "_", sp, "_3panel_global.png")),
             (p_overlay_global + p_top + p_overlay_hl) + plot_layout(nrow = 1), width = 18, height = 6, dpi = 300)
    }, error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing _3panel_global.png: ", e$message))
    
    tryCatch({
      ggsave(file.path(fig_dir, paste0(ds, "_", sp, "_3panel_global_white_red.png")),
             (p_overlay_white_red + p_top + p_overlay_hl) + plot_layout(nrow = 1), width = 18, height = 6, dpi = 300)
    }, error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing _3panel_global_white_red.png: ", e$message))
    
    tryCatch({
      ggsave(file.path(fig_dir, paste0(ds, "_", sp, "_protected_and_diff.png")),
             (p_top + p_bottom + p_diff) + plot_layout(nrow = 1), width = 18, height = 6, dpi = 300)
    }, error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing _protected_and_diff.png: ", e$message))
    
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
        tryCatch({
          ggsave(file.path(fig_dir, paste0(ds, "_", sp, "_permissiveness_global_week_", w, ".png")), p_w, width = 10, height = 8, dpi = 300)
        }, error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing weekly png: ", e$message))
      }
    }
    
    # Export hotspot/protected cell lists for downstream DE/CellChat.
    top_thr <- tryCatch(stats::quantile(md_s$permissiveness_global, probs = 1 - TOP_FRAC, na.rm = TRUE), error = function(e) NA_real_)
    bot_thr <- tryCatch(stats::quantile(md_s$permissiveness_global, probs = BOTTOM_FRAC, na.rm = TRUE), error = function(e) NA_real_)
    hotspots <- if (!is.na(top_thr)) md_s %>% filter(permissiveness_global >= top_thr) else md_s[0,]
    protected <- if (!is.na(bot_thr)) md_s %>% filter(permissiveness_global <= bot_thr) else md_s[0,]
    tryCatch(write.csv(hotspots, file.path(tab_dir, paste0(ds, "_", sp, "_hotspot_cells_top", as.integer(TOP_FRAC * 100), "pct.csv")), row.names = FALSE),
             error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing hotspots CSV: ", e$message))
    tryCatch(write.csv(protected, file.path(tab_dir, paste0(ds, "_", sp, "_protected_cells_bottom", as.integer(BOTTOM_FRAC * 100), "pct.csv")), row.names = FALSE),
             error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing protected CSV: ", e$message))
    
    tryCatch(write.csv(md_s, file.path(tab_dir, paste0(ds, "_", sp, "_mapped_coords_global.csv")), row.names = FALSE),
             error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing mapped coords CSV: ", e$message))
    
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
    tryCatch(write_json(meta, file.path(tab_dir, "meta.json"), pretty = TRUE, auto_unbox = TRUE),
             error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing meta.json: ", e$message))
    
    readme <- c(
      paste0("Dataset: ", ds),
      paste0("Sample: ", sp),
      "This folder contains global-calibrated spatial permissiveness outputs.",
      "Global score = z_global(MMP) + z_global(Immune_Tolerance) - z_global(Cytotoxic_NK) + z_global(Ethanolamine).",
      "protected_and_diff figure includes bottom-percentile protected density and hotspot-protected difference.",
      paste0("Top fraction: ", TOP_FRAC, "; Bottom fraction: ", BOTTOM_FRAC)
    )
    tryCatch(writeLines(readme, file.path(tab_dir, "README_sample.txt")),
             error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("[", ds, "/", sp, "] error writing README: ", e$message))
    
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      dataset = ds,
      sample = as.character(sp),
      n_cells = nrow(md_s),
      highlight_labels = paste(hl, collapse = ";"),
      fig_dir = fig_dir,
      tab_dir = tab_dir,
      stringsAsFactors = FALSE
    )
    
    process_log[[length(process_log) + 1L]] <- paste0("[", ds, "/", sp, "] done (cells=", nrow(md_s), ")")
  }
}

# summary + logs
if (length(summary_rows) > 0) {
  summary_df <- bind_rows(summary_rows)
  tryCatch(write.csv(summary_df, file.path(DIR_TAB, "summary_table_with_highlights.csv"), row.names = FALSE),
           error = function(e) process_log[[length(process_log) + 1L]] <<- paste0("Error writing summary_table_with_highlights.csv: ", e$message))
  log_msg2(sprintf("[05H] Summary written: %d rows.", nrow(summary_df)), logfile)
} else {
  log_msg2("[05H] No sample-level outputs produced.", logfile)
}

# PROCESS LOG
proc_file <- file.path(DIR_TAB, "PROCESS_LOG.txt")
tryCatch(writeLines(unlist(c(process_log)), proc_file), error = function(e) log_msg2(paste0("Failed to write PROCESS_LOG: ", e$message), logfile))
log_msg2("[05H] Done.", logfile)

# ---- Write description files for methods & interpretation ----
# Scripts descriptions (methods notes)
desc_dir <- "scripts/05_spatial"
tryCatch({
  if (!dir.exists(desc_dir)) dir.create(desc_dir, recursive = TRUE, showWarnings = FALSE)
  f05f <- file.path(desc_dir, "05F_description.txt")
  f05h <- file.path(desc_dir, "05H_description.txt")
  fig_readme <- file.path(DIR_FIG, "README_plot_interpretation.txt")
  
  txt_05f <- c(
    "05F — Spatial Permissiveness: purpose and usage",
    "-----------------------------------------------",
    "",
    "Purpose",
    "-------",
    "Step 05F (companion to 05C/05H) documents the conceptual and practical",
    "purpose of the spatial permissiveness workflow.",
    "",
    "05C computes per-dataset 'permissiveness' by z-scoring 4 modules within each",
    "dataset and combining: permissiveness = z(MMP/ECM remodeling) + z(Immune tolerance)",
    "- z(NK cytotoxic) + z(Ethanolamine metabolism).",
    "",
    "05H performs cross-sample global calibration (pooled z-scores across datasets),",
    "derives 'permissiveness_global', computes hotspot (top N%) and protected",
    "(bottom N%) kernel densities, and produces per-sample figures (vanilla,",
    "white→red symmetric, and highlight variants) for interpretation and downstream",
    "analyses.",
    "",
    "Inputs",
    "------",
    "- output/objects/*_with_permissiveness.rds (preferred), or *_harmonized.rds",
    "- output/tables/*_permissiveness_cell_level.csv (optional but helpful)",
    "- raw images in data/raw/ (image_map.csv optionally maps samples->images)",
    "",
    "Outputs",
    "-------",
    "- output/figures/05_spatial/<dataset>/<sample>/*  (3-panel, white->red, protected/diff, weeks)",
    "- output/tables/05_spatial/<dataset>/<sample>/*  (mapped coords, hotspots, protected lists, meta.json)",
    "- output/tables/05_spatial/permissiveness_global_allcells.csv",
    "- output/tables/05_spatial/permissiveness_global_stats.json",
    "- output/tables/05_spatial/PROCESS_LOG.txt",
    "",
    "Interpretation notes",
    "--------------------",
    "- 05C scores are dataset-relative (z within dataset). 05H provides a globally",
    "  calibrated view (permissiveness_global) that is comparable across samples.",
    "- White→red panels use a symmetric color scale centered at zero for comparability:",
    "  +/- perm_maxabs (stored in permissiveness_global_stats.json).",
    "- Hotspots = kernel-density of top N% permissive cells (default 10%); Protected =",
    "  kernel-density of bottom N% (default 10%). Hotspot - Protected = difference map.",
    "- Use hotspot cell lists for DE, pathway analysis, and CellChat.",
    "",
    "Caveats",
    "-------",
    "- If a dataset lacks enough genes for some modules (STARmap or targeted assays),",
    "  NK or MMP contributions may be weak — inspect NK_module_coverage_qc.csv.",
    "- Global calibration assumes module scores are computed on compatible assays or",
    "  on an integrated assay. If module scores are heterogeneous (different assays),",
    "  interpret global comparisons with caution."
  )
  writeLines(txt_05f, f05f)
  
  txt_05h <- c(
    "05H — Global-calibrated spatial permissiveness panels (README)",
    "---------------------------------------------------------------",
    "",
    "What 05H does",
    "-------------",
    "1. Reads harmonized Seurat objects (preferably *_with_permissiveness.rds).",
    "2. Pools per-cell module scores across datasets and computes global z-scores:",
    "   z_mmp_global, z_tol_global, z_nk_global, z_ea_global.",
    "3. Calculates permissiveness_global = z_mmp_global + z_tol_global - z_nk_global + z_ea_global.",
    "4. Exports a pooled table (permissiveness_global_allcells.csv) and stats file with",
    "   perm_maxabs (permissiveness_global_stats.json).",
    "5. For each dataset and sample:",
    "   - Maps coordinates to image pixels using heuristics,",
    "   - Draws three-panel figures:",
    "       A) Vanilla: per-cell permissiveness_global (viridis)",
    "       B) Density: top-permissive kernel density (top N%)",
    "       C) Highlights: permissiveness + highlighted cell types",
    "     and a white->red variant with symmetric color limits.",
    "   - Creates protected (bottom N%) density maps and a Hotspot-Protected difference map.",
    "6. Writes per-sample files: mapped_coords_global.csv, hotspot/protected lists,",
    "   meta.json, README_sample.txt, and figure PNG/PDFs.",
    "",
    "Why symmetric color limits?",
    "---------------------------",
    "We compute a pooled `perm_maxabs = max(abs(permissiveness_global))`. Using",
    "symmetric limits (`[-perm_maxabs, +perm_maxabs]`) centered at 0 allows direct",
    "visual comparison between samples: a red pixel of the same intensity means",
    "the same relative permissiveness across samples.",
    "",
    "How to interpret each panel",
    "---------------------------",
    "- Panel A (vanilla): per-cell score (viridis). Use for local cell-level patterns.",
    "- Panel B (density): spatial concentration of top permissive cells (hotspots).",
    "  Useful for locating regions to analyze (DE genes, neighborhood composition).",
    "- Panel C (highlights): overlay of highlighted biologically-important cell types",
    "  (e.g. EVT, NK, Endothelial) on permissiveness — helps connect biology to score.",
    "- White->Red: same cell-level data but with symmetric red ramp — good for cross-sample comparisons.",
    "- Protected and Diff: regions enriched for bottom-permissiveness (protected) and Hotspot - Protected difference.",
    "",
    "Suggested next analyses",
    "-----------------------",
    "- Differential expression between hotspot cells and background cells (Seurat FindMarkers).",
    "- CellChat / ligand-receptor analysis on hotspot cells vs background cells.",
    "- Spatially-aware tests (SPARK, spatialDE) focusing on hotspot regions."
  )
  writeLines(txt_05h, f05h)
  
  txt_figreadme <- c(
    "Plot interpretation guide (short)",
    "--------------------------------",
    "",
    "Files produced per sample:",
    "- <dataset>_<sample>_3panel_global.png : (Viridis) Global permissiveness, hotspot density, highlights.",
    "- <dataset>_<sample>_3panel_global_white_red.png : same but white->red symmetric color scale for cross-sample comparability.",
    "- <dataset>_<sample>_protected_and_diff.png : Top density, bottom(protected) density, and difference (hotspot - protected).",
    "- <dataset>_<sample>_permissiveness_global_week_<w>.png : per-week global permissiveness maps.",
    "",
    "Key points:",
    "- Permissiveness_global is pooled across datasets and comparable across samples.",
    "- Hotspots = highest-permissive cells (top 10% by default). Protected = lowest 10%.",
    "- Difference map (hotspot - protected) highlights spatially differential permissiveness.",
    "- Use hotspot CSVs for downstream DE/CellChat.",
    "- Check meta.json and permissiveness_global_stats.json for provenance and the perm_maxabs used for symmetric color scales.",
    "",
    "Caveats:",
    "- If module scores are missing (e.g., due to limited gene panels), permissiveness contributions may be NA; check NK_module_coverage_qc files.",
    "- Interpret differences across platforms (STARmap vs SlideTags) carefully even after global calibration; module computation assumptions matter."
  )
  writeLines(txt_figreadme, fig_readme)
  
  log_msg2("[05H] Wrote descriptive README files in scripts/05_spatial and output/figures/05_spatial", logfile)
>>>>>>> 54c1a5c675ca223d183141f8a75de090237ce902
}, error = function(e) log_msg2(paste0("[05H] Warning: failed to write description files: ", e$message), logfile))