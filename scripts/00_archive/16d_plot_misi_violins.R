#!/usr/bin/env Rscript

source("R/config.R")
source("R/helpers_io.R")
source("R/celltype_dictionary.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(qs)
})

log_file <- file.path(CFG$dirs$logs, "pipeline.log")
log_msg("16d_plot_misi_violins starting.", log_file = log_file)

# ---------------------------------------------------------------------------
# Output directory structure:
#   misi_violins/
#   ├── by_condition/              ← condition-colored violins (categorical palette)
#   ├── white_to_darkred/          ← MISI-mapped fill: white → dark red
#   └── diverging_blue_red/        ← MISI-mapped fill: navy → white → dark red
# ---------------------------------------------------------------------------
base_out_dir <- file.path(CFG$dirs$figures, "misi_violins")
dir_condition  <- file.path(base_out_dir, "by_condition")
dir_wtr        <- file.path(base_out_dir, "white_to_darkred")
dir_div        <- file.path(base_out_dir, "diverging_blue_red")
for (d in c(dir_condition, dir_wtr, dir_div)) ensure_dir(d)

condition_palette <- c(UI = "#CCCCCC", Lm = "#E31A1C", Pf = "#FF7F00", Tg = "#33A02C")
condition_levels <- names(condition_palette)
misi_col_priority <- c("MISIv2_MISI", "Target_MISI_Score", "MISI", "MISI_weighted_lambda1se", "MISI_naive")

# ---------------------------------------------------------------------------
# Helper: load Seurat — PRIORITIZE seu_with_misi.qs
# ---------------------------------------------------------------------------
load_any_seurat <- function() {
  candidates <- c(
    file.path(CFG$dirs$objects, "seu_with_misi_v2.qs"),
    file.path(CFG$dirs$objects, "seu_with_misi.qs"),
    file.path(CFG$dirs$objects, "seu_with_misi_v2.qs"),
    file.path(CFG$dirs$objects, "seu_core_architecture.qs"),
    file.path(CFG$dirs$objects, "seu_with_scores.qs"),
    file.path(CFG$dirs$objects, "seu_clean.qs"),
    file.path(CFG$dirs$objects, "seu_clean.rds")
  )
  if (exists("CFG") && !is.null(CFG$data$seu_qs_path))  candidates <- c(candidates, CFG$data$seu_qs_path)
  if (exists("CFG") && !is.null(CFG$data$seu_rds_path)) candidates <- c(candidates, CFG$data$seu_rds_path)
  for (p in unique(candidates[file.exists(candidates)])) {
    obj <- tryCatch(if (grepl("\\.qs$", p, ignore.case = TRUE)) qs::qread(p, nthreads = 1) else readRDS(p), error = function(e) NULL)
    if (inherits(obj, "Seurat")) {
      log_msg(sprintf("16d: loaded Seurat object from %s", p), log_file = log_file)
      return(obj)
    }
  }
  stop("No loadable Seurat object found.")
}

# ---------------------------------------------------------------------------
# Helper: harmonize condition — handles BOTH "infection" and "infection_hpi"
# ---------------------------------------------------------------------------
canon_condition <- function(x) {
  raw <- trimws(as.character(x))
  first_token <- sub("^([^_]+)_.*$", "\\1", raw)
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

seu <- load_any_seurat()
misi_col <- misi_col_priority[misi_col_priority %in% colnames(seu@meta.data)][1]
if (is.na(misi_col)) stop("MISI metadata missing.")
log_msg(sprintf("16d: using MISI column '%s'", misi_col), log_file = log_file)

cell_col <- CFG$cols$cell_type
if (!cell_col %in% colnames(seu@meta.data)) {
  for (alt in c("celltype_refined", "celltype_author", "cell_type", "celltype"))
    if (alt %in% colnames(seu@meta.data)) seu[[cell_col]] <- seu@meta.data[[alt]]
}
if (!cell_col %in% colnames(seu@meta.data)) stop("No cell-type column found.")

# Use infection column directly when available (cleanest 4-level condition)
if (CFG$cols$infection %in% colnames(seu@meta.data)) {
  cond_source <- as.character(seu@meta.data[[CFG$cols$infection]])
} else if ("condition" %in% colnames(seu@meta.data)) {
  cond_source <- as.character(seu@meta.data$condition)
} else {
  cond_source <- rep("Unknown", ncol(seu))
}

df <- seu@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  transmute(
    cell,
    condition = canon_condition(cond_source),
    celltype_plot = rename_with_true_names(as.character(.data[[cell_col]])),
    Target_MISI_Score = as.numeric(.data[[misi_col]])
  ) %>%
  filter(!is.na(Target_MISI_Score), !is.na(celltype_plot))

df$condition <- factor(
  df$condition,
  levels = c(condition_levels, setdiff(sort(unique(df$condition)), condition_levels)))

log_msg("  Conditions: ", paste(levels(df$condition), collapse = ", "), log_file = log_file)
log_msg("  Cell types: ", length(unique(df$celltype_plot)), log_file = log_file)

# Filter to cell types with at least 2 conditions having >= 20 cells
valid_ct <- df %>%
  count(celltype_plot, condition, name = "n") %>%
  group_by(celltype_plot) %>%
  summarise(n_cond_ge20 = sum(n >= 20), .groups = "drop") %>%
  filter(n_cond_ge20 >= 2) %>%
  pull(celltype_plot)
df <- df %>% filter(celltype_plot %in% valid_ct)

present_conds <- unique(as.character(df$condition))
pairwise_comparisons <- if (length(present_conds) >= 2) {
  combn(present_conds, 2, simplify = FALSE)
} else {
  list()
}

# Extend condition palette for unexpected conditions
fill_pal <- condition_palette
extra_conds <- setdiff(unique(as.character(df$condition)), names(fill_pal))
if (length(extra_conds) > 0) {
  extra_cols <- grDevices::hcl.colors(length(extra_conds), palette = "Set 2")
  names(extra_cols) <- extra_conds
  fill_pal <- c(fill_pal, extra_cols)
}

# Compute per-celltype × condition median MISI for heatmap-fill violins
ct_medians <- df %>%
  group_by(celltype_plot, condition) %>%
  summarise(median_misi = median(Target_MISI_Score, na.rm = TRUE), .groups = "drop")

df <- df %>%
  left_join(ct_medians, by = c("celltype_plot", "condition"))

misi_range <- range(df$Target_MISI_Score, na.rm = TRUE)
misi_median_global <- median(df$Target_MISI_Score, na.rm = TRUE)

# ---------------------------------------------------------------------------
# stat_compare_means layers (reused across plots)
# ---------------------------------------------------------------------------
add_stats <- function(p) {
  if (length(pairwise_comparisons) >= 1) {
    p <- p +
      ggpubr::stat_compare_means(
        comparisons = pairwise_comparisons, method = "wilcox.test",
        label = "p.signif", hide.ns = FALSE, size = 2.7,
        bracket.size = 0.22, step.increase = 0.055) +
      ggpubr::stat_compare_means(
        comparisons = pairwise_comparisons, method = "wilcox.test",
        label = "p.format", hide.ns = TRUE, size = 2.1,
        vjust = -0.35, step.increase = 0.055)
  } else {
    log_msg("16d: <2 conditions; skipping stat_compare_means.", log_file = log_file)
  }
  p
}

base_theme <- theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = "bold", size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ===========================================================================
# PLOT 1: Condition-colored violins (original style)
# ===========================================================================
p_cond <- ggplot(df, aes(x = condition, y = Target_MISI_Score, fill = condition)) +
  geom_violin(trim = FALSE, scale = "width", color = "black",
              linewidth = 0.25, alpha = 0.95) +
  geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.22,
               color = "black", linewidth = 0.25) +
  scale_fill_manual(values = fill_pal, drop = FALSE) +
  facet_wrap(~celltype_plot, scales = "free_y") +
  labs(title = "MISI vulnerability by condition and cell type",
       subtitle = "Violin + boxplot | per-cell-type pairwise Wilcoxon tests",
       x = "Condition", y = "MISI Score") +
  base_theme

p_cond <- add_stats(p_cond)

save_both(p_cond, file.path(dir_condition, "Fig16d_MISI_Violins_By_CellType"), 16, 11)
log_msg("  Saved condition-colored violins.", log_file = log_file)

# ===========================================================================
# PLOT 2: White-to-dark-red fill mapped to median MISI per violin
# ===========================================================================
p_wtr <- ggplot(df, aes(x = condition, y = Target_MISI_Score, fill = median_misi)) +
  geom_violin(trim = FALSE, scale = "width", color = "black",
              linewidth = 0.25, alpha = 0.95) +
  geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.22,
               color = "black", linewidth = 0.25, fill = NA) +
  scale_fill_gradientn(
    colours = c("#FFFFFF", "#FEE5D9", "#FC9272", "#CB181D", "#8B0000"),
    limits = misi_range, oob = scales::squish, name = "Median\nMISI") +
  facet_wrap(~celltype_plot, scales = "free_y") +
  labs(title = "MISI vulnerability \u2014 white to dark red heatmap fill",
       subtitle = "Violin fill = median MISI per cell-type \u00d7 condition | Wilcoxon tests",
       x = "Condition", y = "MISI Score") +
  base_theme +
  theme(legend.position = "right")

p_wtr <- add_stats(p_wtr)

save_both(p_wtr, file.path(dir_wtr, "Fig16d_MISI_Violins_WhiteToDarkRed"), 17, 11)
log_msg("  Saved white-to-darkred heatmap violins.", log_file = log_file)

# ===========================================================================
# PLOT 3: Diverging navy-blue → white → dark-red fill (centered on median)
# ===========================================================================
p_div <- ggplot(df, aes(x = condition, y = Target_MISI_Score, fill = median_misi)) +
  geom_violin(trim = FALSE, scale = "width", color = "black",
              linewidth = 0.25, alpha = 0.95) +
  geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.22,
               color = "black", linewidth = 0.25, fill = NA) +
  scale_fill_gradientn(
    colours = c("#08306B", "#2171B5", "#6BAED6", "#DEEBF7",
                "#FFFFFF",
                "#FEE0D2", "#FC9272", "#CB181D", "#67000D"),
    limits = misi_range,
    values = scales::rescale(
      c(misi_range[1],
        misi_range[1] + 0.25 * (misi_median_global - misi_range[1]),
        misi_range[1] + 0.75 * (misi_median_global - misi_range[1]),
        misi_median_global - 0.001 * diff(misi_range),
        misi_median_global,
        misi_median_global + 0.001 * diff(misi_range),
        misi_median_global + 0.25 * (misi_range[2] - misi_median_global),
        misi_median_global + 0.75 * (misi_range[2] - misi_median_global),
        misi_range[2]),
      to = c(0, 1)),
    oob = scales::squish, name = "Median\nMISI") +
  facet_wrap(~celltype_plot, scales = "free_y") +
  labs(title = "MISI vulnerability \u2014 diverging blue-white-red heatmap fill",
       subtitle = "Centered on global median | blue = below, red = above | Wilcoxon tests",
       x = "Condition", y = "MISI Score") +
  base_theme +
  theme(legend.position = "right")

p_div <- add_stats(p_div)

save_both(p_div, file.path(dir_div, "Fig16d_MISI_Violins_DivergingBlueRed"), 17, 11)
log_msg("  Saved diverging blue-red heatmap violins.", log_file = log_file)

log_msg("16d_plot_misi_violins done. Outputs in: ", base_out_dir, log_file = log_file)