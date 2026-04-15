#!/usr/bin/env Rscript

source("R/config.R")
source("R/helpers_io.R")
source("R/celltype_dictionary.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(qs)
  library(readr)
})

log_file <- file.path(CFG$dirs$logs, "pipeline.log")
log_msg("16c_plot_misi_embeddings_V2 starting.", log_file = log_file)

# ---------------------------------------------------------------------------
# Output directory structure:
#   misi_embeddings/
#   ├── magma/                  ← viridis magma (sequential, perceptually uniform)
#   ├── white_to_darkred/       ← white → dark red (sequential, publication classic)
#   └── diverging_blue_red/     ← navy → white → dark red (diverging, centered on median)
# ---------------------------------------------------------------------------
base_out_dir <- file.path(CFG$dirs$figures, "misi_embeddings")
scale_dirs <- list(
  A_magma              = file.path(base_out_dir, "magma"),
  B_white_to_darkred   = file.path(base_out_dir, "white_to_darkred"),
  C_diverging_blue_red = file.path(base_out_dir, "diverging_blue_red")
)
for (d in scale_dirs) ensure_dir(d)

canonical_levels <- c("UI", "Lm", "Pf", "Tg")
misi_col_priority <- c("MISIv2_MISI", "Target_MISI_Score", "MISI", "MISI_weighted_lambda1se", "MISI_naive")

# ---------------------------------------------------------------------------
# Three color scales
# ---------------------------------------------------------------------------
color_versions <- list(
  A_magma = function(limits, midpoint = NULL) {
    scale_color_viridis_c(option = "magma", limits = limits,
                          oob = scales::squish, name = "MISI")
  },
  B_white_to_darkred = function(limits, midpoint = NULL) {
    scale_color_gradientn(
      colours = c("#FFFFFF", "#FEE5D9", "#FC9272", "#CB181D", "#8B0000"),
      limits = limits, oob = scales::squish, name = "MISI")
  },
  C_diverging_blue_red = function(limits, midpoint = NULL) {
    mid <- if (!is.null(midpoint)) midpoint else mean(limits)
    scale_color_gradientn(
      colours = c("#08306B", "#2171B5", "#6BAED6", "#DEEBF7",
                  "#FFFFFF",
                  "#FEE0D2", "#FC9272", "#CB181D", "#67000D"),
      limits = limits, values = scales::rescale(
        c(limits[1], limits[1] + 0.25 * (mid - limits[1]),
          limits[1] + 0.75 * (mid - limits[1]), mid - 0.001 * diff(limits),
          mid,
          mid + 0.001 * diff(limits), mid + 0.25 * (limits[2] - mid),
          mid + 0.75 * (limits[2] - mid), limits[2]),
        to = c(0, 1)),
      oob = scales::squish, name = "MISI")
  }
)

# ---------------------------------------------------------------------------
# Helper: load Seurat — PRIORITIZE seu_with_misi.qs
# ---------------------------------------------------------------------------
load_any_seurat <- function() {
  primary <- file.path(CFG$dirs$objects, "seu_with_misi_v2.qs")
  if (file.exists(primary)) {
    obj <- tryCatch(qs::qread(primary, nthreads = 1), error = function(e) NULL)
    if (inherits(obj, "Seurat")) {
      log_msg("Loaded Seurat object from ", primary, log_file = log_file)
      return(obj)
    }
    stop("Primary MISI object exists but could not be loaded as Seurat: ", primary)
  }

  candidates <- c(
    file.path(CFG$dirs$objects, "seu_with_misi_v2.qs")
  )
  if (exists("CFG") && !is.null(CFG$data$seu_qs_path))  candidates <- c(candidates, CFG$data$seu_qs_path)
  if (exists("CFG") && !is.null(CFG$data$seu_rds_path)) candidates <- c(candidates, CFG$data$seu_rds_path)
  candidates <- unique(candidates[file.exists(candidates)])
  for (p in candidates) {
    obj <- tryCatch(
      if (grepl("\\.qs$", p, ignore.case = TRUE)) qs::qread(p, nthreads = 1) else readRDS(p),
      error = function(e) NULL)
    if (inherits(obj, "Seurat")) {
      log_msg("Loaded fallback Seurat object from ", p, log_file = log_file)
      return(obj)
    }
  }
  stop("No loadable Seurat object found. Expected primary: ", primary)
}

resolve_reduction <- function(obj, target) {
  aliases <- list(
    umap = c("umap", "harmony_umap", "UMAP", "ref.umap"),
    tsne = c("tsne", "harmony_tsne", "tSNE", "ref.tsne")
  )
  for (nm in aliases[[target]]) if (nm %in% names(obj@reductions)) return(nm)
  NA_character_
}

# ---------------------------------------------------------------------------
# Helper: harmonize condition — handles BOTH "infection" and "infection_hpi"
# ---------------------------------------------------------------------------
canon_condition <- function(x) {
  raw <- trimws(as.character(x))
  # If values look like "Lm_24h" or "UI_0h", extract just the infection part
  # by splitting on underscore and taking the first token
  first_token <- sub("^([^_]+)_.*$", "\\1", raw)
  # Use first_token for matching if it looks like an infection_hpi format
  is_compound <- grepl("^[A-Za-z]+_\\d+", raw)
  to_match <- ifelse(is_compound, first_token, raw)
  low <- tolower(to_match)
  out <- dplyr::case_when(
    grepl("^ui$|uninfect|control|mock", low) ~ "UI",
    grepl("^lm$|lister", low)                ~ "Lm",
    grepl("^pf$|falcip|plasmod", low)        ~ "Pf",
    grepl("^tg$|toxop", low)                 ~ "Tg",
    TRUE ~ to_match
  )
  out[out == "" | is.na(out)] <- "Unknown"
  out
}

# ---------------------------------------------------------------------------
# Load and prepare data
# ---------------------------------------------------------------------------
seu <- load_any_seurat()

# Pull MISI from object or from table and join by cell barcode
misi_col <- misi_col_priority[misi_col_priority %in% colnames(seu@meta.data)][1]

if (is.na(misi_col)) {
  misi_tables <- c(
    file.path(CFG$dirs$tables, "misi_v2_component_scores_per_cell.csv"),
    file.path(CFG$dirs$tables, "misi_per_cell_with_weights.csv"),
    file.path(CFG$dirs$tables, "misi_per_cell.csv")
  )
  misi_tables <- misi_tables[file.exists(misi_tables)]
  if (length(misi_tables) > 0) {
    mt <- readr::read_csv(misi_tables[1], show_col_types = FALSE)
    score_hit <- c("MISIv2_MISI", "Target_MISI_Score", "MISI_weighted_lambda1se", "MISI", "MISI_naive")
    score_hit <- score_hit[score_hit %in% colnames(mt)][1]
    if (!is.na(score_hit) && "cell" %in% colnames(mt)) {
      md <- seu@meta.data %>% tibble::rownames_to_column("cell")
      md <- md %>% left_join(mt[, c("cell", score_hit)], by = "cell")
      seu[["MISI_for_plot"]] <- as.numeric(md[[score_hit]])
      misi_col <- "MISI_for_plot"
    }
  }
}
if (is.na(misi_col) || !misi_col %in% colnames(seu@meta.data)) stop("MISI score not found in object or MISI tables.")
log_msg("16c: using MISI column ", misi_col, log_file = log_file)

cell_col <- CFG$cols$cell_type
if (!cell_col %in% colnames(seu@meta.data)) {
  for (alt in c("celltype_refined", "celltype_author", "cell_type", "celltype"))
    if (alt %in% colnames(seu@meta.data)) seu[[cell_col]] <- seu@meta.data[[alt]]
}
if (!cell_col %in% colnames(seu@meta.data)) stop("No cell-type column found.")

# Determine best source for condition:
# 1) If infection column exists, use it directly (cleanest)
# 2) Else if condition column exists, parse it
# 3) Else fallback
if (CFG$cols$infection %in% colnames(seu@meta.data)) {
  cond_source <- as.character(seu@meta.data[[CFG$cols$infection]])
} else if ("condition" %in% colnames(seu@meta.data)) {
  cond_source <- as.character(seu@meta.data$condition)
} else {
  cond_source <- rep("Unknown", ncol(seu))
}
seu$condition_harmonized <- canon_condition(cond_source)
seu$condition_harmonized <- factor(
  seu$condition_harmonized,
  levels = c(canonical_levels, setdiff(sort(unique(seu$condition_harmonized)), canonical_levels)))

# Also keep the detailed condition (infection_hpi) for optional per-timepoint facets
if ("condition" %in% colnames(seu@meta.data)) {
  seu$condition_detailed <- as.character(seu@meta.data$condition)
} else if (all(c(CFG$cols$infection, CFG$cols$hpi) %in% colnames(seu@meta.data))) {
  seu$condition_detailed <- paste0(
    as.character(seu@meta.data[[CFG$cols$infection]]), "_",
    as.character(seu@meta.data[[CFG$cols$hpi]]))
} else {
  seu$condition_detailed <- seu$condition_harmonized
}

seu$celltype_plot <- rename_with_true_names(as.character(seu@meta.data[[cell_col]]))
cell_levels <- sort(unique(as.character(seu$celltype_plot)))
cell_cols <- universal_colors[cell_levels]
cell_cols[is.na(cell_cols)] <- "#BDBDBD"

if (!any(c("umap", "harmony_umap", "UMAP", "ref.umap") %in% names(seu@reductions)) &&
    !any(c("tsne", "harmony_tsne", "tSNE", "ref.tsne") %in% names(seu@reductions))) {
  stop("No UMAP/tSNE reductions found. Run script 03 core architecture first.")
}

build_df <- function(obj, red) {
  emb <- Embeddings(obj, red)
  tibble::tibble(
    cell = rownames(emb), x = emb[, 1], y = emb[, 2],
    condition = obj$condition_harmonized,
    condition_detailed = obj$condition_detailed,
    celltype_plot = obj$celltype_plot,
    Target_MISI_Score = as.numeric(obj@meta.data[[misi_col]]))
}

save_both <- function(p, base, w, h) {
  ggsave(paste0(base, ".pdf"), p, width = w, height = h, bg = "white")
  ggsave(paste0(base, ".png"), p, width = w, height = h, dpi = 320, bg = "white")
}

# --- Aggregated (all cells, one panel) ---
plot_aggregated <- function(df, red_label, scale_fun, midpoint = NULL) {
  lim <- range(df$Target_MISI_Score, na.rm = TRUE)

  p1 <- ggplot(df, aes(x, y, color = celltype_plot)) +
    geom_point(size = 0.2, alpha = 0.9) +
    scale_color_manual(values = cell_cols, drop = FALSE) +
    coord_equal() + theme_bw(11) +
    labs(title = paste0(red_label, " \u2014 Cell Types"))

  p2 <- ggplot(df, aes(x, y, color = Target_MISI_Score)) +
    geom_point(size = 0.2, alpha = 0.9) +
    scale_fun(lim, midpoint = midpoint) +
    coord_equal() + theme_bw(11) +
    labs(title = paste0(red_label, " \u2014 MISI"))

  (p1 | p2) & theme(axis.text = element_blank(), axis.ticks = element_blank())
}

# --- Faceted by 4 canonical conditions (UI, Lm, Pf, Tg) ---
plot_faceted_condition <- function(df, red_label, scale_fun, midpoint = NULL) {
  lim <- range(df$Target_MISI_Score, na.rm = TRUE)

  p3 <- ggplot(df, aes(x, y, color = celltype_plot)) +
    geom_point(size = 0.2, alpha = 0.9) +
    scale_color_manual(values = cell_cols, drop = FALSE) +
    facet_wrap(~condition, nrow = 1, drop = FALSE) +
    coord_equal() + theme_bw(10) +
    labs(title = paste0(red_label, " \u2014 Cell Types by Condition"))

  p4 <- ggplot(df, aes(x, y, color = Target_MISI_Score)) +
    geom_point(size = 0.2, alpha = 0.9) +
    scale_fun(lim, midpoint = midpoint) +
    facet_wrap(~condition, nrow = 1, drop = FALSE) +
    coord_equal() + theme_bw(10) +
    labs(title = paste0(red_label, " \u2014 MISI by Condition"))

  (p3 / p4) & theme(axis.text = element_blank(), axis.ticks = element_blank(),
                     strip.text = element_text(face = "bold"))
}

# --- Faceted by detailed condition (infection × hpi) ---
plot_faceted_detailed <- function(df, red_label, scale_fun, midpoint = NULL) {
  lim <- range(df$Target_MISI_Score, na.rm = TRUE)
  n_det <- length(unique(df$condition_detailed))

  p5 <- ggplot(df, aes(x, y, color = celltype_plot)) +
    geom_point(size = 0.15, alpha = 0.8) +
    scale_color_manual(values = cell_cols, drop = FALSE) +
    facet_wrap(~condition_detailed, nrow = 2) +
    coord_equal() + theme_bw(8) +
    labs(title = paste0(red_label, " \u2014 Cell Types by Condition \u00d7 Timepoint"))

  p6 <- ggplot(df, aes(x, y, color = Target_MISI_Score)) +
    geom_point(size = 0.15, alpha = 0.8) +
    scale_fun(lim, midpoint = midpoint) +
    facet_wrap(~condition_detailed, nrow = 2) +
    coord_equal() + theme_bw(8) +
    labs(title = paste0(red_label, " \u2014 MISI by Condition \u00d7 Timepoint"))

  (p5 / p6) & theme(axis.text = element_blank(), axis.ticks = element_blank(),
                     strip.text = element_text(face = "bold", size = 7))
}

# ---------------------------------------------------------------------------
# Generate all plots
# For each scale × each reduction:
#   1. Aggregated (all cells, side-by-side celltype + MISI)
#   2. Faceted by 4 canonical conditions (UI, Lm, Pf, Tg)
#   3. Faceted by detailed condition (infection × hpi) — only if >4 unique
# ---------------------------------------------------------------------------
for (target in c("umap", "tsne")) {
  red <- resolve_reduction(seu, target)
  if (is.na(red)) {
    log_msg("Skipping ", target, ": missing reduction.", log_file = log_file)
    next
  }

  df <- build_df(seu, red) %>%
    filter(!is.na(Target_MISI_Score), !is.na(condition), !is.na(x), !is.na(y))
  if (nrow(df) == 0) {
    log_msg("Skipping ", target, ": no rows after filtering.", log_file = log_file)
    next
  }

  misi_median <- median(df$Target_MISI_Score, na.rm = TRUE)
  n_cond <- nlevels(df$condition)
  n_detailed <- length(unique(df$condition_detailed))

  for (ver in names(color_versions)) {
    out_dir <- scale_dirs[[ver]]
    midpt <- if (ver == "C_diverging_blue_red") misi_median else NULL

    # 1. Aggregated
    p_agg <- plot_aggregated(df, toupper(target), color_versions[[ver]], midpoint = midpt)
    save_both(p_agg,
              file.path(out_dir, paste0("Fig16c_", toupper(target), "_Aggregated")),
              15, 6.8)

    # 2. Faceted by canonical condition (4 panels)
    p_fac <- plot_faceted_condition(df, toupper(target), color_versions[[ver]], midpoint = midpt)
    save_both(p_fac,
              file.path(out_dir, paste0("Fig16c_", toupper(target), "_By_Condition")),
              max(14, 3.4 * n_cond), 10)

    # 3. Faceted by detailed condition (only if more than canonical 4)
    if (n_detailed > n_cond) {
      p_det <- plot_faceted_detailed(df, toupper(target), color_versions[[ver]], midpoint = midpt)
      w_det <- max(16, 3.0 * ceiling(n_detailed / 2))
      save_both(p_det,
                file.path(out_dir, paste0("Fig16c_", toupper(target), "_By_Condition_Timepoint")),
                w_det, 14)
    }

    log_msg("  Saved ", ver, " / ", toupper(target), " plots.", log_file = log_file)
  }
}

log_msg("16c_plot_misi_embeddings_V2 done. Outputs in: ", base_out_dir, log_file = log_file)