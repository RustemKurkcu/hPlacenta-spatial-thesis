# ======================================================================
# scripts/04_timecourse/04A_timecourse_posthoc_significance_plots.R
# Lightweight post-hoc plotting from 04A enhanced summary table.
#
# Purpose:
# - Build "all" and "significant-only" timecourse figures without rerunning 04A.
#
# Inputs:
# - output/tables/timecourse_gene_module_summaries_enhanced.csv
#
# Outputs:
# - output/tables/timecourse_trend_significance_posthoc.csv
# - output/tables/timecourse_trend_significance_posthoc_significant_only.csv
# - output/figures/posthoc_timecourse_<metric>_all.png
# - output/figures/posthoc_timecourse_<metric>_significant_only.png
# ======================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

source("config/config.R")
source("scripts/R/utils.R")

ensure_dir(DIR_TABLES); ensure_dir(DIR_FIGURES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "04A_timecourse_posthoc_significance_plots.log")

infile <- file.path(DIR_TABLES, "timecourse_gene_module_summaries_enhanced.csv")
if (!file.exists(infile)) stop("Missing input: ", infile, ". Run 04A enhanced first.")

sumdf <- read.csv(infile, stringsAsFactors = FALSE)
if (!all(c("dataset", "week", "celltype_use") %in% colnames(sumdf))) {
  stop("Required columns missing from summary table.")
}

metric_cols <- c(
  "PLD1_mean", "MMP2_mean", "MMP9_mean", "HLA.G_mean", "NKG7_mean",
  "score_MMP_ECM_Remodeling_mean", "score_Immune_Tolerance_mean",
  "score_Cytotoxic_NK_mean", "score_Ethanolamine_Metabolism_mean"
)
metric_cols <- intersect(metric_cols, colnames(sumdf))
if (length(metric_cols) == 0) stop("No target *_mean metric columns found.")

trend_stats <- function(metric_col) {
  sumdf %>%
    filter(!is.na(.data[[metric_col]]), is.finite(week)) %>%
    group_by(dataset, celltype_use) %>%
    summarize(
      n_weeks = n_distinct(week),
      rho = if (n_weeks >= 3) suppressWarnings(cor(week, .data[[metric_col]], method = "spearman", use = "complete.obs")) else NA_real_,
      p_value = if (n_weeks >= 3) suppressWarnings(cor.test(week, .data[[metric_col]], method = "spearman", exact = FALSE)$p.value) else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      trend_sig = ifelse(!is.na(p_adj) & p_adj < 0.05, "FDR<0.05", "NS"),
      metric = metric_col
    )
}

trend_all <- bind_rows(lapply(metric_cols, trend_stats))
write.csv(trend_all, file.path(DIR_TABLES, "timecourse_trend_significance_posthoc.csv"), row.names = FALSE)
write.csv(trend_all %>% filter(trend_sig == "FDR<0.05"),
          file.path(DIR_TABLES, "timecourse_trend_significance_posthoc_significant_only.csv"), row.names = FALSE)

for (metric_col in metric_cols) {
  df_metric <- sumdf %>%
    filter(!is.na(.data[[metric_col]])) %>%
    left_join(trend_all %>% filter(metric == metric_col), by = c("dataset", "celltype_use"))
  
  nice_name <- gsub("_mean$", "", metric_col)
  n_sig <- sum(unique(df_metric[c("dataset", "celltype_use", "trend_sig")])$trend_sig == "FDR<0.05", na.rm = TRUE)
  
  p_all <- ggplot(df_metric, aes(x = week, y = .data[[metric_col]], color = celltype_use, group = celltype_use, linetype = trend_sig)) +
    geom_line(alpha = 0.7) +
    geom_point(size = 1.8) +
    facet_wrap(~dataset, scales = "free_y", ncol = 1) +
    scale_linetype_manual(values = c("FDR<0.05" = "solid", "NS" = "22")) +
    theme_classic(base_size = 12) +
    labs(title = paste0(nice_name, " Over Gestational Time (all)"),
         subtitle = "Line type = Spearman trend significance",
         caption = paste0("Significant trends (BH FDR<0.05): ", n_sig),
         x = "Gestational week", y = nice_name, color = "Cell type", linetype = "Trend")
  save_plot(p_all, file.path(DIR_FIGURES, paste0("posthoc_timecourse_", nice_name, "_all.png")), w = 10, h = 8)
  
  sig_ct <- trend_all %>% filter(metric == metric_col, trend_sig == "FDR<0.05") %>% pull(celltype_use) %>% unique()
  if (length(sig_ct) > 0) {
    p_sig <- p_all
    p_sig$data <- p_sig$data %>% filter(celltype_use %in% sig_ct)
    p_sig <- p_sig + labs(title = paste0(nice_name, " Over Gestational Time (significant only)"))
    save_plot(p_sig, file.path(DIR_FIGURES, paste0("posthoc_timecourse_", nice_name, "_significant_only.png")), w = 10, h = 8)
  }
}

log_msg("Done post-hoc timecourse significance plotting.", logfile)