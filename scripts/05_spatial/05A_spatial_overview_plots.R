# ======================================================================
# scripts/05_spatial/05A_spatial_overview_plots.R
# Spatial maps for Slide-tags and STARmap (colored by best available label).
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

source("config/config.R")
source("scripts/R/utils.R")

check_required_packages(c("Seurat", "dplyr", "ggplot2", "tibble"), context = "05A_spatial_overview_plots")

ensure_dir(DIR_FIGURES); ensure_dir(DIR_TABLES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "05A_spatial_overview_plots.log")

# Prefer harmonized objects (if you ran 03C)
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
log_msg(paste0("Loading STARmap object: ", starmap_path), logfile)
starmap   <- readRDS(starmap_path)

# Ensure coords
slidetags <- ensure_week_column(slidetags, COL_WEEK_CANDIDATES)
starmap   <- ensure_week_column(starmap, COL_WEEK_CANDIDATES)

slidetags <- ensure_spatial_coords(slidetags, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)
starmap   <- ensure_spatial_coords(starmap, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)

pick_label_col <- function(md) {
  if ("celltype_final_refined" %in% colnames(md)) return("celltype_final_refined")
  if ("celltype_final_conservative" %in% colnames(md)) return("celltype_final_conservative")
  if (COL_PRED_CELLTYPE %in% colnames(md)) return(COL_PRED_CELLTYPE)
  if ("celltype_author" %in% colnames(md)) return("celltype_author")
  NULL
}

plot_spatial_by_week <- function(obj, dataset_name) {
  md <- obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    mutate(week = .data[[COL_WEEK]])
  
  label_col <- pick_label_col(md)
  if (is.null(label_col)) stop("No label column found for ", dataset_name)
  
  md <- md %>%
    mutate(
      celltype = as.character(.data[[label_col]]),
      x = suppressWarnings(as.numeric(spatial_x_use)),
      y = suppressWarnings(as.numeric(spatial_y_use))
    )
  
  n_finite <- sum(is.finite(md$x) & is.finite(md$y))
  log_msg(sprintf("%s: %d/%d cells have finite spatial coordinates",
                  dataset_name, n_finite, nrow(md)), logfile)
  if (n_finite == 0) {
    log_msg(sprintf("%s: no finite coordinates found; skipping spatial plots", dataset_name), logfile)
    return(invisible(0L))
  }
  
  weeks <- sort(unique(md$week[!is.na(md$week)]))
  if (length(weeks) == 0) weeks <- unique(md$week)
  
  n_saved <- 0L
  for (w in weeks) {
    d <- md %>%
      filter(week == w, is.finite(x), is.finite(y), !is.na(celltype), nzchar(celltype))
    if (nrow(d) == 0) {
      log_msg(sprintf("%s week %s: 0 plottable cells; skipped", dataset_name, as.character(w)), logfile)
      next
    }
    p <- scatter_spatial(d, x = "x", y = "y", color = "celltype",
                         title = paste0(dataset_name, " spatial (week ", w, ")"))
    path <- file.path(DIR_FIGURES, paste0(dataset_name, "_spatial_week_", w, ".png"))
    save_plot(p, path, w = 10, h = 8)
    n_saved <- n_saved + 1L
    log_msg(sprintf("Saved: %s (n=%d cells)", path, nrow(d)), logfile)
  }
  
  if (n_saved == 0) {
    log_msg(sprintf("%s: no figures were saved (check week/label/coordinate columns)", dataset_name), logfile)
  }
  invisible(n_saved)
}


plot_spatial_by_week(slidetags, "SlideTags")
plot_spatial_by_week(starmap, "STARmap")
n_slide <- plot_spatial_by_week(slidetags, "SlideTags")
n_star <- plot_spatial_by_week(starmap, "STARmap")

log_msg("Saved spatial overview plots.", logfile)
log_msg(sprintf("Saved spatial overview plots. SlideTags=%d, STARmap=%d", n_slide, n_star), logfile)

