<<<<<<< HEAD
# ======================================================================
# scripts/05_spatial/05G_niche_susceptibility_scoring.R
#
# Build niche-level susceptibility from 05D niche clusters + 05C permissiveness.
#
# Outputs per dataset (SlideTags, STARmap):
# - *_niche_susceptibility_scores.csv
# - *_niche_susceptibility_classification.csv
# - *_niche_susceptibility_week_composition.csv
# - *_niche_susceptibility_week_trends.csv
# - *_with_niche_susceptibility.rds
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
})

source("config/config.R")
source("scripts/R/utils.R")
if (file.exists("config/vulnerability_config.R")) {
  source("config/vulnerability_config.R")
}

check_required_packages(c("Seurat", "dplyr", "tidyr", "ggplot2", "tibble"),
                        context = "05G_niche_susceptibility_scoring")

ensure_dir(DIR_OBJS); ensure_dir(DIR_TABLES); ensure_dir(DIR_FIGURES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "05G_niche_susceptibility_scoring.log")

SUSCEPT_MODE <- get0("SUSCEPT_MODE", ifnotfound = "quantile")
SUSCEPT_Q <- as.numeric(get0("SUSCEPT_Q", ifnotfound = 0.80))
SUSCEPT_Z_THRESHOLD <- as.numeric(get0("SUSCEPT_Z_THRESHOLD", ifnotfound = 0.0))
VULN_MIN_WEEKS_FOR_TREND <- as.integer(get0("VULN_MIN_WEEKS_FOR_TREND", ifnotfound = 3L))

safe_read_rds_with_fallback <- function(candidates, label) {
  for (pth in unique(candidates)) {
    if (is.null(pth) || is.na(pth) || !nzchar(pth) || !file.exists(pth)) next
    finfo <- file.info(pth)
    if (is.na(finfo$size) || finfo$size <= 0) {
      log_msg(paste0("[05G] ", label, ": empty file skipped: ", pth), logfile)
      next
    }
    obj <- tryCatch(readRDS(pth), error = function(e) e)
    if (!inherits(obj, "error")) {
      log_msg(paste0("[05G] Loaded ", label, ": ", pth), logfile)
      return(obj)
    }
    log_msg(paste0("[05G] Failed reading ", label, " at ", pth, " (", obj$message, ")"), logfile)
  }
  stop("[05G] Could not load ", label, " from candidate paths")
}

safe_read_csv_with_fallback <- function(candidates, label) {
  for (pth in unique(candidates)) {
    if (is.null(pth) || is.na(pth) || !nzchar(pth) || !file.exists(pth)) next
    tab <- tryCatch(read.csv(pth, stringsAsFactors = FALSE), error = function(e) e)
    if (!inherits(tab, "error")) {
      log_msg(paste0("[05G] Loaded ", label, " table: ", pth), logfile)
      return(tab)
    }
  }
  NULL
}

parse_week_numeric <- function(x) {
  x <- as.character(x)
  out <- suppressWarnings(as.numeric(x))
  miss <- !is.finite(out)
  if (any(miss)) {
    y <- x[miss]
    y <- gsub("[Ww]", "", y)
    y <- gsub("week[_-]?", "", y)
    y <- sub("^(\\d+)-(\\d+)$", "\\1.\\2", y)
    y <- sub(".*?(\\d+(?:\\.\\d+)?).*$", "\\1", y)
    out[miss] <- suppressWarnings(as.numeric(y))
  }
  out
}

resolve_dataset_paths <- function(dataset) {
  ds <- tolower(dataset)
  if (ds == "slidetags") {
    list(
      niche_obj = c(file.path(DIR_OBJS, "slidetags_harmonized.rds"),
                    file.path(DIR_OBJS, "slidetags_mapped.rds"),
                    file.path(DIR_OBJS, "slidetags_mapped_to_multiome.rds"),
                    PATH_SLIDETAGS_RDS, PATH_SLIDETAGS_RAW),
      permiss_obj = c(file.path(DIR_OBJS, "slidetags_with_permissiveness.rds"),
                      file.path(DIR_OBJS, "slidetags_harmonized.rds"),
                      PATH_SLIDETAGS_RDS, PATH_SLIDETAGS_RAW),
      permiss_table = c(file.path(DIR_TABLES, "SlideTags_permissiveness_cell_level.csv"),
                        file.path(DIR_TABLES, "Slide-tags_permissiveness_cell_level.csv"))
    )
  } else {
    list(
      niche_obj = c(file.path(DIR_OBJS, "starmap_harmonized.rds"),
                    file.path(DIR_OBJS, "starmap_mapped.rds"),
                    file.path(DIR_OBJS, "starmap_mapped_to_multiome.rds"),
                    PATH_STARMAP_RDS),
      permiss_obj = c(file.path(DIR_OBJS, "starmap_with_permissiveness.rds"),
                      file.path(DIR_OBJS, "starmap_harmonized.rds"),
                      PATH_STARMAP_RDS),
      permiss_table = c(file.path(DIR_TABLES, "STARmap_permissiveness_cell_level.csv"))
    )
  }
}

compute_threshold <- function(scores) {
  scores <- scores[is.finite(scores)]
  if (length(scores) == 0) return(Inf)
  if (tolower(SUSCEPT_MODE) == "fixed") {
    return(SUSCEPT_Z_THRESHOLD)
  }
  stats::quantile(scores, probs = SUSCEPT_Q, na.rm = TRUE)
}

score_one_dataset <- function(dataset_name) {
  paths <- resolve_dataset_paths(dataset_name)
  
  obj_niche <- safe_read_rds_with_fallback(paths$niche_obj, paste0(dataset_name, " niche object"))
  obj_niche <- ensure_week_column(obj_niche, COL_WEEK_CANDIDATES)
  
  md <- obj_niche@meta.data %>%
    rownames_to_column("cell") %>%
    mutate(week = as.character(.data[[COL_WEEK]]),
           week_num = parse_week_numeric(week))
  
  if (!("niche_cluster" %in% colnames(md))) {
    stop("[05G] ", dataset_name, ": niche_cluster not found. Run 05D first.")
  }
  
  if (!("permissiveness" %in% colnames(md))) {
    obj_perm <- safe_read_rds_with_fallback(paths$permiss_obj, paste0(dataset_name, " permissiveness object"))
    md_perm <- obj_perm@meta.data %>%
      rownames_to_column("cell")
    
    if ("permissiveness" %in% colnames(md_perm)) {
      md <- md %>% left_join(md_perm %>% select(cell, permissiveness), by = "cell")
    }
  }
  
  if (!("permissiveness" %in% colnames(md)) || all(is.na(md$permissiveness))) {
    tbl_perm <- safe_read_csv_with_fallback(paths$permiss_table, paste0(dataset_name, " permissiveness"))
    if (!is.null(tbl_perm) && all(c("cell", "permissiveness") %in% colnames(tbl_perm))) {
      md <- md %>%
        select(-any_of("permissiveness")) %>%
        left_join(tbl_perm %>% select(cell, permissiveness), by = "cell")
    }
  }
  
  if (!("permissiveness" %in% colnames(md)) || all(!is.finite(md$permissiveness))) {
    stop("[05G] ", dataset_name, ": permissiveness unavailable after fallbacks. Run 05C first.")
  }
  
  md <- md %>% mutate(permissiveness = as.numeric(permissiveness))
  
  niche_scores <- md %>%
    filter(!is.na(niche_cluster), is.finite(permissiveness)) %>%
    group_by(niche_cluster) %>%
    summarize(
      dataset = dataset_name,
      n_cells = dplyr::n(),
      mean_perm = mean(permissiveness, na.rm = TRUE),
      median_perm = median(permissiveness, na.rm = TRUE),
      q75_perm = quantile(permissiveness, 0.75, na.rm = TRUE),
      q90_perm = quantile(permissiveness, 0.90, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      susceptibility_score = mean_perm,
      susceptibility_rank = dplyr::dense_rank(dplyr::desc(susceptibility_score))
    ) %>%
    arrange(susceptibility_rank, niche_cluster)
  
  if (nrow(niche_scores) == 0) {
    stop("[05G] ", dataset_name, ": no finite niche permissiveness scores.")
  }
  
  thr <- compute_threshold(niche_scores$susceptibility_score)
  
  niche_class <- niche_scores %>%
    mutate(
      threshold_used = thr,
      threshold_mode = ifelse(tolower(SUSCEPT_MODE) == "fixed", "fixed", paste0("quantile_", SUSCEPT_Q)),
      susceptibility_class = ifelse(susceptibility_score >= thr, "susceptible", "non_susceptible")
    ) %>%
    select(dataset, niche_cluster, susceptibility_score, threshold_mode, threshold_used, susceptibility_class)
  
  md2 <- md %>%
    left_join(niche_class %>%
                select(niche_cluster, susceptibility_score, susceptibility_class), by = "niche_cluster") %>%
    mutate(
      susceptibility_score = as.numeric(susceptibility_score),
      susceptibility_class = ifelse(is.na(susceptibility_class), "non_susceptible", susceptibility_class)
    )
  
  week_comp <- md2 %>%
    filter(!is.na(week), nzchar(week), !is.na(niche_cluster)) %>%
    count(dataset = dataset_name, week, week_num, niche_cluster, susceptibility_class, name = "n_cells") %>%
    group_by(dataset, week, week_num) %>%
    mutate(frac_of_week = n_cells / sum(n_cells)) %>%
    ungroup() %>%
    arrange(week_num, niche_cluster)
  
  trend_tbl <- week_comp %>%
    group_by(dataset, niche_cluster, susceptibility_class) %>%
    summarize(
      n_weeks = n_distinct(week),
      rho_week = if (sum(is.finite(week_num) & is.finite(frac_of_week)) >= VULN_MIN_WEEKS_FOR_TREND) {
        suppressWarnings(cor(week_num, frac_of_week, method = "spearman", use = "complete.obs"))
      } else {
        NA_real_
      },
      slope_week = if (sum(is.finite(week_num) & is.finite(frac_of_week)) >= 2) {
        suppressWarnings(unname(stats::coef(stats::lm(frac_of_week ~ week_num))[2]))
      } else {
        NA_real_
      },
      p_value = if (sum(is.finite(week_num) & is.finite(frac_of_week)) >= VULN_MIN_WEEKS_FOR_TREND) {
        suppressWarnings(stats::cor.test(week_num, frac_of_week, method = "spearman", exact = FALSE)$p.value)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      trend = dplyr::case_when(
        !is.na(p_adj) & p_adj < 0.05 & rho_week > 0 ~ "increasing",
        !is.na(p_adj) & p_adj < 0.05 & rho_week < 0 ~ "decreasing",
        TRUE ~ "not_significant"
      )
    )
  
  # Persist tables
  write.csv(niche_scores,
            file.path(DIR_TABLES, paste0(dataset_name, "_niche_susceptibility_scores.csv")),
            row.names = FALSE)
  write.csv(niche_class,
            file.path(DIR_TABLES, paste0(dataset_name, "_niche_susceptibility_classification.csv")),
            row.names = FALSE)
  write.csv(week_comp,
            file.path(DIR_TABLES, paste0(dataset_name, "_niche_susceptibility_week_composition.csv")),
            row.names = FALSE)
  write.csv(trend_tbl,
            file.path(DIR_TABLES, paste0(dataset_name, "_niche_susceptibility_week_trends.csv")),
            row.names = FALSE)
  
  # Add metadata back to object and save
  obj_niche@meta.data$niche_susceptibility_score <- md2$susceptibility_score[match(rownames(obj_niche@meta.data), md2$cell)]
  obj_niche@meta.data$niche_susceptibility_class <- md2$susceptibility_class[match(rownames(obj_niche@meta.data), md2$cell)]
  
  out_rds <- file.path(DIR_OBJS, paste0(tolower(dataset_name), "_with_niche_susceptibility.rds"))
  saveRDS(obj_niche, out_rds)
  
  log_msg(sprintf("[05G] %s: niches=%d susceptible=%d saved=%s",
                  dataset_name,
                  nrow(niche_scores),
                  sum(niche_class$susceptibility_class == "susceptible", na.rm = TRUE),
                  out_rds), logfile)
  
  list(scores = niche_scores, class = niche_class, week = week_comp, trend = trend_tbl)
}

log_msg("[05G] Starting niche susceptibility scoring...", logfile)
res_slide <- score_one_dataset("SlideTags")
res_star <- score_one_dataset("STARmap")

# Combined convenience outputs
scores_all <- bind_rows(res_slide$scores, res_star$scores) %>%
  mutate(timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
class_all <- bind_rows(res_slide$class, res_star$class) %>%
  mutate(timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
trend_all <- bind_rows(res_slide$trend, res_star$trend) %>%
  mutate(timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

write.csv(scores_all, file.path(DIR_TABLES, "combined_niche_susceptibility_scores.csv"), row.names = FALSE)
write.csv(class_all, file.path(DIR_TABLES, "combined_niche_susceptibility_classification.csv"), row.names = FALSE)
write.csv(trend_all, file.path(DIR_TABLES, "combined_niche_susceptibility_week_trends.csv"), row.names = FALSE)

log_msg("[05G] Done niche susceptibility scoring.", logfile)
=======
# ======================================================================
# scripts/05_spatial/05G_niche_susceptibility_scoring.R
#
# Build niche-level susceptibility from 05D niche clusters + 05C permissiveness.
#
# Outputs per dataset (SlideTags, STARmap):
# - *_niche_susceptibility_scores.csv
# - *_niche_susceptibility_classification.csv
# - *_niche_susceptibility_week_composition.csv
# - *_niche_susceptibility_week_trends.csv
# - *_with_niche_susceptibility.rds
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
})

source("config/config.R")
source("scripts/R/utils.R")
if (file.exists("config/vulnerability_config.R")) {
  source("config/vulnerability_config.R")
}

check_required_packages(c("Seurat", "dplyr", "tidyr", "ggplot2", "tibble"),
                        context = "05G_niche_susceptibility_scoring")

ensure_dir(DIR_OBJS); ensure_dir(DIR_TABLES); ensure_dir(DIR_FIGURES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "05G_niche_susceptibility_scoring.log")

SUSCEPT_MODE <- get0("SUSCEPT_MODE", ifnotfound = "quantile")
SUSCEPT_Q <- as.numeric(get0("SUSCEPT_Q", ifnotfound = 0.80))
SUSCEPT_Z_THRESHOLD <- as.numeric(get0("SUSCEPT_Z_THRESHOLD", ifnotfound = 0.0))
VULN_MIN_WEEKS_FOR_TREND <- as.integer(get0("VULN_MIN_WEEKS_FOR_TREND", ifnotfound = 3L))

safe_read_rds_with_fallback <- function(candidates, label) {
  for (pth in unique(candidates)) {
    if (is.null(pth) || is.na(pth) || !nzchar(pth) || !file.exists(pth)) next
    finfo <- file.info(pth)
    if (is.na(finfo$size) || finfo$size <= 0) {
      log_msg(paste0("[05G] ", label, ": empty file skipped: ", pth), logfile)
      next
    }
    obj <- tryCatch(readRDS(pth), error = function(e) e)
    if (!inherits(obj, "error")) {
      log_msg(paste0("[05G] Loaded ", label, ": ", pth), logfile)
      return(obj)
    }
    log_msg(paste0("[05G] Failed reading ", label, " at ", pth, " (", obj$message, ")"), logfile)
  }
  stop("[05G] Could not load ", label, " from candidate paths")
}

safe_read_csv_with_fallback <- function(candidates, label) {
  for (pth in unique(candidates)) {
    if (is.null(pth) || is.na(pth) || !nzchar(pth) || !file.exists(pth)) next
    tab <- tryCatch(read.csv(pth, stringsAsFactors = FALSE), error = function(e) e)
    if (!inherits(tab, "error")) {
      log_msg(paste0("[05G] Loaded ", label, " table: ", pth), logfile)
      return(tab)
    }
  }
  NULL
}

parse_week_numeric <- function(x) {
  x <- as.character(x)
  out <- suppressWarnings(as.numeric(x))
  miss <- !is.finite(out)
  if (any(miss)) {
    y <- x[miss]
    y <- gsub("[Ww]", "", y)
    y <- gsub("week[_-]?", "", y)
    y <- sub("^(\\d+)-(\\d+)$", "\\1.\\2", y)
    y <- sub(".*?(\\d+(?:\\.\\d+)?).*$", "\\1", y)
    out[miss] <- suppressWarnings(as.numeric(y))
  }
  out
}

resolve_dataset_paths <- function(dataset) {
  ds <- tolower(dataset)
  if (ds == "slidetags") {
    list(
      niche_obj = c(file.path(DIR_OBJS, "slidetags_harmonized.rds"),
                    file.path(DIR_OBJS, "slidetags_mapped.rds"),
                    file.path(DIR_OBJS, "slidetags_mapped_to_multiome.rds"),
                    PATH_SLIDETAGS_RDS, PATH_SLIDETAGS_RAW),
      permiss_obj = c(file.path(DIR_OBJS, "slidetags_with_permissiveness.rds"),
                      file.path(DIR_OBJS, "slidetags_harmonized.rds"),
                      PATH_SLIDETAGS_RDS, PATH_SLIDETAGS_RAW),
      permiss_table = c(file.path(DIR_TABLES, "SlideTags_permissiveness_cell_level.csv"),
                        file.path(DIR_TABLES, "Slide-tags_permissiveness_cell_level.csv"))
    )
  } else {
    list(
      niche_obj = c(file.path(DIR_OBJS, "starmap_harmonized.rds"),
                    file.path(DIR_OBJS, "starmap_mapped.rds"),
                    file.path(DIR_OBJS, "starmap_mapped_to_multiome.rds"),
                    PATH_STARMAP_RDS),
      permiss_obj = c(file.path(DIR_OBJS, "starmap_with_permissiveness.rds"),
                      file.path(DIR_OBJS, "starmap_harmonized.rds"),
                      PATH_STARMAP_RDS),
      permiss_table = c(file.path(DIR_TABLES, "STARmap_permissiveness_cell_level.csv"))
    )
  }
}

compute_threshold <- function(scores) {
  scores <- scores[is.finite(scores)]
  if (length(scores) == 0) return(Inf)
  if (tolower(SUSCEPT_MODE) == "fixed") {
    return(SUSCEPT_Z_THRESHOLD)
  }
  stats::quantile(scores, probs = SUSCEPT_Q, na.rm = TRUE)
}

score_one_dataset <- function(dataset_name) {
  paths <- resolve_dataset_paths(dataset_name)
  
  obj_niche <- safe_read_rds_with_fallback(paths$niche_obj, paste0(dataset_name, " niche object"))
  obj_niche <- ensure_week_column(obj_niche, COL_WEEK_CANDIDATES)
  
  md <- obj_niche@meta.data %>%
    rownames_to_column("cell") %>%
    mutate(week = as.character(.data[[COL_WEEK]]),
           week_num = parse_week_numeric(week))
  
  if (!("niche_cluster" %in% colnames(md))) {
    stop("[05G] ", dataset_name, ": niche_cluster not found. Run 05D first.")
  }
  
  if (!("permissiveness" %in% colnames(md))) {
    obj_perm <- safe_read_rds_with_fallback(paths$permiss_obj, paste0(dataset_name, " permissiveness object"))
    md_perm <- obj_perm@meta.data %>%
      rownames_to_column("cell")
    
    if ("permissiveness" %in% colnames(md_perm)) {
      md <- md %>% left_join(md_perm %>% select(cell, permissiveness), by = "cell")
    }
  }
  
  if (!("permissiveness" %in% colnames(md)) || all(is.na(md$permissiveness))) {
    tbl_perm <- safe_read_csv_with_fallback(paths$permiss_table, paste0(dataset_name, " permissiveness"))
    if (!is.null(tbl_perm) && all(c("cell", "permissiveness") %in% colnames(tbl_perm))) {
      md <- md %>%
        select(-any_of("permissiveness")) %>%
        left_join(tbl_perm %>% select(cell, permissiveness), by = "cell")
    }
  }
  
  if (!("permissiveness" %in% colnames(md)) || all(!is.finite(md$permissiveness))) {
    stop("[05G] ", dataset_name, ": permissiveness unavailable after fallbacks. Run 05C first.")
  }
  
  md <- md %>% mutate(permissiveness = as.numeric(permissiveness))
  
  niche_scores <- md %>%
    filter(!is.na(niche_cluster), is.finite(permissiveness)) %>%
    group_by(niche_cluster) %>%
    summarize(
      dataset = dataset_name,
      n_cells = dplyr::n(),
      mean_perm = mean(permissiveness, na.rm = TRUE),
      median_perm = median(permissiveness, na.rm = TRUE),
      q75_perm = quantile(permissiveness, 0.75, na.rm = TRUE),
      q90_perm = quantile(permissiveness, 0.90, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      susceptibility_score = mean_perm,
      susceptibility_rank = dplyr::dense_rank(dplyr::desc(susceptibility_score))
    ) %>%
    arrange(susceptibility_rank, niche_cluster)
  
  if (nrow(niche_scores) == 0) {
    stop("[05G] ", dataset_name, ": no finite niche permissiveness scores.")
  }
  
  thr <- compute_threshold(niche_scores$susceptibility_score)
  
  niche_class <- niche_scores %>%
    mutate(
      threshold_used = thr,
      threshold_mode = ifelse(tolower(SUSCEPT_MODE) == "fixed", "fixed", paste0("quantile_", SUSCEPT_Q)),
      susceptibility_class = ifelse(susceptibility_score >= thr, "susceptible", "non_susceptible")
    ) %>%
    select(dataset, niche_cluster, susceptibility_score, threshold_mode, threshold_used, susceptibility_class)
  
  md2 <- md %>%
    left_join(niche_class %>%
                select(niche_cluster, susceptibility_score, susceptibility_class), by = "niche_cluster") %>%
    mutate(
      susceptibility_score = as.numeric(susceptibility_score),
      susceptibility_class = ifelse(is.na(susceptibility_class), "non_susceptible", susceptibility_class)
    )
  
  week_comp <- md2 %>%
    filter(!is.na(week), nzchar(week), !is.na(niche_cluster)) %>%
    count(dataset = dataset_name, week, week_num, niche_cluster, susceptibility_class, name = "n_cells") %>%
    group_by(dataset, week, week_num) %>%
    mutate(frac_of_week = n_cells / sum(n_cells)) %>%
    ungroup() %>%
    arrange(week_num, niche_cluster)
  
  trend_tbl <- week_comp %>%
    group_by(dataset, niche_cluster, susceptibility_class) %>%
    summarize(
      n_weeks = n_distinct(week),
      rho_week = if (sum(is.finite(week_num) & is.finite(frac_of_week)) >= VULN_MIN_WEEKS_FOR_TREND) {
        suppressWarnings(cor(week_num, frac_of_week, method = "spearman", use = "complete.obs"))
      } else {
        NA_real_
      },
      slope_week = if (sum(is.finite(week_num) & is.finite(frac_of_week)) >= 2) {
        suppressWarnings(unname(stats::coef(stats::lm(frac_of_week ~ week_num))[2]))
      } else {
        NA_real_
      },
      p_value = if (sum(is.finite(week_num) & is.finite(frac_of_week)) >= VULN_MIN_WEEKS_FOR_TREND) {
        suppressWarnings(stats::cor.test(week_num, frac_of_week, method = "spearman", exact = FALSE)$p.value)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      trend = dplyr::case_when(
        !is.na(p_adj) & p_adj < 0.05 & rho_week > 0 ~ "increasing",
        !is.na(p_adj) & p_adj < 0.05 & rho_week < 0 ~ "decreasing",
        TRUE ~ "not_significant"
      )
    )
  
  # Persist tables
  write.csv(niche_scores,
            file.path(DIR_TABLES, paste0(dataset_name, "_niche_susceptibility_scores.csv")),
            row.names = FALSE)
  write.csv(niche_class,
            file.path(DIR_TABLES, paste0(dataset_name, "_niche_susceptibility_classification.csv")),
            row.names = FALSE)
  write.csv(week_comp,
            file.path(DIR_TABLES, paste0(dataset_name, "_niche_susceptibility_week_composition.csv")),
            row.names = FALSE)
  write.csv(trend_tbl,
            file.path(DIR_TABLES, paste0(dataset_name, "_niche_susceptibility_week_trends.csv")),
            row.names = FALSE)
  
  # Add metadata back to object and save
  obj_niche@meta.data$niche_susceptibility_score <- md2$susceptibility_score[match(rownames(obj_niche@meta.data), md2$cell)]
  obj_niche@meta.data$niche_susceptibility_class <- md2$susceptibility_class[match(rownames(obj_niche@meta.data), md2$cell)]
  
  out_rds <- file.path(DIR_OBJS, paste0(tolower(dataset_name), "_with_niche_susceptibility.rds"))
  saveRDS(obj_niche, out_rds)
  
  log_msg(sprintf("[05G] %s: niches=%d susceptible=%d saved=%s",
                  dataset_name,
                  nrow(niche_scores),
                  sum(niche_class$susceptibility_class == "susceptible", na.rm = TRUE),
                  out_rds), logfile)
  
  list(scores = niche_scores, class = niche_class, week = week_comp, trend = trend_tbl)
}

log_msg("[05G] Starting niche susceptibility scoring...", logfile)
res_slide <- score_one_dataset("SlideTags")
res_star <- score_one_dataset("STARmap")

# Combined convenience outputs
scores_all <- bind_rows(res_slide$scores, res_star$scores) %>%
  mutate(timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
class_all <- bind_rows(res_slide$class, res_star$class) %>%
  mutate(timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
trend_all <- bind_rows(res_slide$trend, res_star$trend) %>%
  mutate(timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

write.csv(scores_all, file.path(DIR_TABLES, "combined_niche_susceptibility_scores.csv"), row.names = FALSE)
write.csv(class_all, file.path(DIR_TABLES, "combined_niche_susceptibility_classification.csv"), row.names = FALSE)
write.csv(trend_all, file.path(DIR_TABLES, "combined_niche_susceptibility_week_trends.csv"), row.names = FALSE)

log_msg("[05G] Done niche susceptibility scoring.", logfile)
>>>>>>> 54c1a5c675ca223d183141f8a75de090237ce902
