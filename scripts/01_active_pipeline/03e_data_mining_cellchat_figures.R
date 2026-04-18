#!/usr/bin/env Rscript

# =============================================================================
# Script: 03e_data_mining_cellchat_figures.R
# Purpose: Exploratory/"fun" comparative figures for thesis-pathway CellChat outputs.
#          Saves figures under output/figures/data_mining_cellchat/.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(jsonlite)
})

source("R/celltype_dictionary.R")

PIPELINE_NAME <- "03e_data_mining_cellchat_figures"
PIPELINE_VERSION <- "1.0.0"
RUN_TIMESTAMP <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
RUN_SEED <- 42L
set.seed(RUN_SEED)

OUT_ROOT <- "output"
DIR_LOGS <- file.path(OUT_ROOT, "logs")
DIR_OBJECTS <- file.path(OUT_ROOT, "objects")
DIR_FIGURES <- file.path(OUT_ROOT, "figures", "data_mining_cellchat")
DIR_REPORTS <- file.path(OUT_ROOT, "reports")
for (d in c(DIR_LOGS, DIR_OBJECTS, DIR_FIGURES, DIR_REPORTS)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
LOG_FILE <- file.path(DIR_LOGS, "03e_data_mining_cellchat_figures.log")

log_msg <- function(..., .level = "INFO") {
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  txt <- paste0(..., collapse = "")
  line <- paste0(stamp, " [", .level, "] ", txt)
  cat(line, "\n")
  cat(line, "\n", file = LOG_FILE, append = TRUE)
}

resolve_object_path <- function(path_rds) {
  path_qs <- sub("\\.rds$", ".qs", path_rds)
  if (file.exists(path_qs)) return(path_qs)
  if (file.exists(path_rds)) return(path_rds)
  NA_character_
}

read_object <- function(path) {
  if (is.na(path)) return(NULL)
  if (grepl("\\.qs$", path)) {
    if (!requireNamespace("qs", quietly = TRUE)) {
      log_msg("qs path found but package 'qs' is not installed: ", path, .level = "WARN")
      return(NULL)
    }
    return(qs::qread(path))
  }
  readRDS(path)
}

record_artifact_manifest <- function(manifest_path, source_data, figures_output, notes = NULL) {
  manifest <- list(
    pipeline = PIPELINE_NAME,
    version = PIPELINE_VERSION,
    run_timestamp = RUN_TIMESTAMP,
    seed = RUN_SEED,
    source_data = source_data,
    figures_output = figures_output,
    script_path = "scripts/01_active_pipeline/03e_data_mining_cellchat_figures.R",
    notes = notes
  )
  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
}

if (requireNamespace("SpatialCellChat", quietly = TRUE)) {
  subset_comm_fn <- SpatialCellChat::subsetCommunication
  rank_net_fn <- SpatialCellChat::rankNet
  log_msg("Using SpatialCellChat::subsetCommunication backend.")
} else if (requireNamespace("CellChat", quietly = TRUE)) {
  subset_comm_fn <- CellChat::subsetCommunication
  rank_net_fn <- CellChat::rankNet
  log_msg("Using CellChat::subsetCommunication backend.")
} else {
  stop("Neither SpatialCellChat nor CellChat is installed.")
}

thesis_pathways <- c("MMP", "WNT", "TGFb", "CXCL", "CCL", "SPP1", "NOTCH", "FN1", "MHC-I", "MHC-II", "CD45", "TIGIT", "PD-L1")
week_ids <- c("W7", "W8-2", "W9", "W11")

week_paths <- setNames(
  vapply(file.path(DIR_OBJECTS, paste0("03_spatial_cellchat_", week_ids, ".rds")), resolve_object_path, character(1)),
  week_ids
)
merged_path <- resolve_object_path(file.path(DIR_OBJECTS, "03_spatial_cellchat_merged.rds"))

if (all(is.na(week_paths))) stop("No week-level Script03 CellChat outputs found.")

cellchat_weeks <- lapply(week_paths, read_object)
cellchat_weeks <- cellchat_weeks[!vapply(cellchat_weeks, is.null, logical(1))]
if (length(cellchat_weeks) == 0) stop("Failed to load any week-level CellChat objects.")

cellchat_merged <- read_object(merged_path)
if (is.null(cellchat_merged)) log_msg("Merged object not found/loaded; merged-based plots will be skipped.", .level = "WARN")

# -----------------------------------------------------------------------------
# Summaries
# -----------------------------------------------------------------------------

summarize_week <- function(wk, obj) {
  weight_mat <- tryCatch(obj@net$weight, error = function(e) NULL)
  total_weight <- if (!is.null(weight_mat)) sum(weight_mat, na.rm = TRUE) else NA_real_

  total_links <- if (!is.null(weight_mat)) sum(weight_mat > 0, na.rm = TRUE) else NA_real_
  n_groups <- length(levels(obj@idents))
  n_cells <- length(obj@idents)

  tibble(
    week = wk,
    n_cells = n_cells,
    n_groups = n_groups,
    total_weight = total_weight,
    total_links = total_links
  )
}

summary_tbl <- bind_rows(lapply(names(cellchat_weeks), function(wk) summarize_week(wk, cellchat_weeks[[wk]])))
summary_csv <- file.path(DIR_REPORTS, "03e_cellchat_week_summary.csv")
write.csv(summary_tbl, summary_csv, row.names = FALSE)

pathway_stats <- bind_rows(lapply(names(cellchat_weeks), function(wk) {
  obj <- cellchat_weeks[[wk]]
  bind_rows(lapply(thesis_pathways, function(pw) {
    df <- tryCatch(
      subset_comm_fn(object = obj, signaling = pw),
      error = function(e) NULL
    )
    if (is.null(df) || nrow(df) == 0) {
      tibble(week = wk, pathway = pw, n_edges = 0L, total_prob = 0)
    } else {
      prob_col <- if ("prob" %in% colnames(df)) "prob" else if ("prob.combined" %in% colnames(df)) "prob.combined" else NULL
      total_prob <- if (!is.null(prob_col)) sum(as.numeric(df[[prob_col]]), na.rm = TRUE) else 0
      tibble(week = wk, pathway = pw, n_edges = nrow(df), total_prob = total_prob)
    }
  }))
}))
pathway_csv <- file.path(DIR_REPORTS, "03e_cellchat_thesis_pathway_stats_by_week.csv")
write.csv(pathway_stats, pathway_csv, row.names = FALSE)

edge_tbl <- bind_rows(lapply(names(cellchat_weeks), function(wk) {
  obj <- cellchat_weeks[[wk]]
  bind_rows(lapply(thesis_pathways, function(pw) {
    df <- tryCatch(subset_comm_fn(object = obj, signaling = pw), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) return(tibble())

    prob_col <- if ("prob" %in% colnames(df)) "prob" else if ("prob.combined" %in% colnames(df)) "prob.combined" else NULL
    src_col <- colnames(df)[grepl("^(source|sender)", colnames(df), ignore.case = TRUE)][1]
    tgt_col <- colnames(df)[grepl("^(target|receiver)", colnames(df), ignore.case = TRUE)][1]
    if (is.na(src_col) || is.na(tgt_col) || is.null(prob_col)) return(tibble())

    tibble(
      week = wk,
      pathway = pw,
      source = as.character(df[[src_col]]),
      target = as.character(df[[tgt_col]]),
      prob = as.numeric(df[[prob_col]])
    )
  }))
}))

edge_summary <- if (nrow(edge_tbl) > 0 && all(c("source", "target", "prob") %in% colnames(edge_tbl))) {
  edge_tbl %>%
    mutate(edge = paste(source, "→", target)) %>%
    group_by(week, pathway, edge) %>%
    summarise(total_prob = sum(prob, na.rm = TRUE), .groups = "drop")
} else {
  log_msg("No valid source/target/prob columns found for edge summary; creating empty edge table.", .level = "WARN")
  tibble(week = character(0), pathway = character(0), edge = character(0), total_prob = numeric(0))
}
edge_csv <- file.path(DIR_REPORTS, "03e_cellchat_top_edges_by_week.csv")
write.csv(edge_summary, edge_csv, row.names = FALSE)

# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------

fig_paths <- character(0)

# 1) Total communication weight and links by week
p1 <- summary_tbl %>%
  pivot_longer(cols = c(total_weight, total_links), names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = week, y = value, fill = week)) +
  geom_col(width = 0.7, color = "black", alpha = 0.9) +
  facet_wrap(~metric, scales = "free_y") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  labs(
    title = "CellChat Communication Burden by Week",
    subtitle = "Total network weight and link counts",
    x = "Week", y = "Value"
  )
p1_pdf <- file.path(DIR_FIGURES, "03e_week_total_weight_and_links.pdf")
p1_png <- file.path(DIR_FIGURES, "03e_week_total_weight_and_links.png")
ggsave(p1_pdf, p1, width = 12, height = 6)
ggsave(p1_png, p1, width = 12, height = 6, dpi = 300)
fig_paths <- c(fig_paths, p1_pdf, p1_png)

# 2) Thesis pathway heatmap by week
p2 <- pathway_stats %>%
  mutate(pathway = factor(pathway, levels = thesis_pathways)) %>%
  ggplot(aes(x = week, y = pathway, fill = log1p(total_prob))) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "C", name = "log1p(total prob)") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Thesis Pathway Activity Heatmap by Week",
    subtitle = "Computed from subsetCommunication signaling summaries",
    x = "Week", y = "Pathway"
  )
p2_pdf <- file.path(DIR_FIGURES, "03e_thesis_pathway_heatmap_by_week.pdf")
p2_png <- file.path(DIR_FIGURES, "03e_thesis_pathway_heatmap_by_week.png")
ggsave(p2_pdf, p2, width = 10, height = 7)
ggsave(p2_png, p2, width = 10, height = 7, dpi = 300)
fig_paths <- c(fig_paths, p2_pdf, p2_png)

# 3) Bubble plot: top edges across thesis pathways
top_edges <- edge_summary %>%
  group_by(pathway, edge) %>%
  summarise(global_prob = sum(total_prob, na.rm = TRUE), .groups = "drop") %>%
  group_by(pathway) %>%
  slice_max(order_by = global_prob, n = 6, with_ties = FALSE) %>%
  ungroup() %>%
  select(pathway, edge)

if (nrow(top_edges) > 0) {
  p3_data <- edge_summary %>%
    inner_join(top_edges, by = c("pathway", "edge")) %>%
    mutate(pathway = factor(pathway, levels = thesis_pathways))

  p3 <- ggplot(p3_data, aes(x = week, y = fct_reorder(edge, total_prob, .fun = max), size = total_prob, color = pathway)) +
    geom_point(alpha = 0.85) +
    theme_minimal(base_size = 12) +
    labs(
      title = "Top Sender→Receiver Edges Across Thesis Pathways",
      subtitle = "Largest communication probabilities across weeks",
      x = "Week", y = "Edge", size = "Total prob"
    )
  p3_pdf <- file.path(DIR_FIGURES, "03e_top_edges_bubble_by_week.pdf")
  p3_png <- file.path(DIR_FIGURES, "03e_top_edges_bubble_by_week.png")
  ggsave(p3_pdf, p3, width = 13, height = 9)
  ggsave(p3_png, p3, width = 13, height = 9, dpi = 300)
  fig_paths <- c(fig_paths, p3_pdf, p3_png)
}

# 4) Merged object pathway rank barplot (if available)
if (!is.null(cellchat_merged)) {
  merged_pathway_tbl <- bind_rows(lapply(thesis_pathways, function(pw) {
    df <- tryCatch(subset_comm_fn(object = cellchat_merged, signaling = pw), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) {
      tibble(pathway = pw, total_prob = 0)
    } else {
      prob_col <- if ("prob" %in% colnames(df)) "prob" else if ("prob.combined" %in% colnames(df)) "prob.combined" else NULL
      total_prob <- if (!is.null(prob_col)) sum(as.numeric(df[[prob_col]]), na.rm = TRUE) else 0
      tibble(pathway = pw, total_prob = total_prob)
    }
  }))

  p4 <- merged_pathway_tbl %>%
    mutate(pathway = fct_reorder(pathway, total_prob)) %>%
    ggplot(aes(x = pathway, y = total_prob, fill = total_prob)) +
    geom_col() +
    coord_flip() +
    scale_fill_viridis_c(option = "B") +
    theme_minimal(base_size = 13) +
    labs(
      title = "Merged CellChat: Thesis Pathway Ranking",
      subtitle = "Total communication probability across merged datasets",
      x = "Pathway", y = "Total prob"
    )

  p4_pdf <- file.path(DIR_FIGURES, "03e_merged_thesis_pathway_ranking.pdf")
  p4_png <- file.path(DIR_FIGURES, "03e_merged_thesis_pathway_ranking.png")
  ggsave(p4_pdf, p4, width = 9, height = 6)
  ggsave(p4_png, p4, width = 9, height = 6, dpi = 300)
  fig_paths <- c(fig_paths, p4_pdf, p4_png)

  # CellChat-native pathway information flow ranking (thesis-relevant summary)
  p5_pdf <- file.path(DIR_FIGURES, "03e_merged_rankNet_information_flow.pdf")
  p5_png <- file.path(DIR_FIGURES, "03e_merged_rankNet_information_flow.png")
  pdf(p5_pdf, width = 12, height = 8)
  tryCatch(
    rank_net_fn(cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = FALSE),
    error = function(e) {
      log_msg("rankNet(comparison) failed on merged object; trying mode='single'. ", conditionMessage(e), .level = "WARN")
      rank_net_fn(cellchat_merged, mode = "single")
    }
  )
  dev.off()
  png(p5_png, width = 3600, height = 2400, res = 300)
  tryCatch(
    rank_net_fn(cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = FALSE),
    error = function(e) rank_net_fn(cellchat_merged, mode = "single")
  )
  dev.off()
  fig_paths <- c(fig_paths, p5_pdf, p5_png)
}

manifest_path <- file.path(DIR_REPORTS, "03e_data_mining_cellchat_figures_manifest.json")
record_artifact_manifest(
  manifest_path = manifest_path,
  source_data = c(unname(week_paths[!is.na(week_paths)]), if (!is.na(merged_path)) merged_path else character(0)),
  figures_output = fig_paths,
  notes = c(
    "Exploratory comparative CellChat figures for thesis pathways",
    "Outputs saved under output/figures/data_mining_cellchat",
    "Includes week-level and merged-level summaries"
  )
)

log_msg("Saved data-mining CellChat figures + manifest: ", manifest_path)
