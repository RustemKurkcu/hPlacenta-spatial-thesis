<<<<<<< HEAD
# ======================================================================
# scripts/08_metagene/08A_housekeeping_qc.R
#
# ADVANCED MODULE (optional): Housekeeping/QC gene selection.
#
# Motivation:
#   Housekeeping gene stability is useful for:
#     (1) detecting sample-to-sample technical shifts,
#     (2) confirming that normalization did not introduce artifacts,
#     (3) providing a stable baseline for metagene...

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

source("config/config.R")
source("scripts/R/utils.R")

ensure_dir(DIR_TABLES); ensure_dir(DIR_FIGURES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "08A_housekeeping_qc.log")

log_msg("Loading processed objects (for housekeeping QC).", logfile)

ref_path <- file.path(DIR_OBJS, "multiome_reference_processed.rds")
sl_path  <- file.path(DIR_OBJS, "slidetags_mapped_to_multiome.rds")
st_path  <- file.path(DIR_OBJS, "starmap_mapped_to_multiome.rds")

objs <- list()
if (file.exists(ref_path)) objs[["Multiome"]] <- readRDS(ref_path)
if (file.exists(sl_path))  objs[["Slide-tags"]] <- readRDS(sl_path)
if (file.exists(st_path))  objs[["STARmap"]] <- readRDS(st_path)

if (length(objs) == 0) {
  stop("No processed objects found. Run the minimal pipeline first.")
}

# --- Helper: compute housekeeping candidates within an assay ---
hk_candidates <- function(obj, assay = NULL, layer = "data") {
  assay <- assay %||% DefaultAssay(obj)
  safe_join_layers(obj, assay = assay)
  DefaultAssay(obj) <- assay
  if (!has_data_layer(obj, assay, layer = layer)) {
    obj <- NormalizeData(obj, verbose = FALSE)
  }
  mat <- Seurat::GetAssayData(obj, assay = assay, layer = layer)
  # Detection frequency (fraction of cells >0) and mean expression
  det <- Matrix::rowSums(mat > 0) / ncol(mat)
  mu  <- Matrix::rowMeans(mat)
  tibble::tibble(gene = names(mu), mean = as.numeric(mu), detect_frac = as.numeric(det)) %>%
    mutate(rank = dplyr::min_rank(desc(detect_frac)))
}

all_hk <- list()
for (nm in names(objs)) {
  obj <- objs[[nm]]
  assay_use <- if (nm == "STARmap" && "imputed" %in% names(obj@assays)) "imputed" else DefaultAssay(obj)
  log_msg(paste0("Housekeeping scan: ", nm, " (assay=", assay_use, ")"), logfile)
  all_hk[[nm]] <- hk_candidates(obj, assay = assay_use)
}

# Combine and identify robust genes: high detection across datasets
merged <- Reduce(function(a, b) full_join(a, b, by = "gene", suffix = c("", "_y")),
                 lapply(names(all_hk), function(nm) {
                   all_hk[[nm]] %>%
                     select(gene, detect_frac, mean) %>%
                     rename_with(~paste0(., "_", nm), c(detect_frac, mean))
                 }))

# A simple robustness score: average detection across available datasets
detect_cols <- grep("^detect_frac_", colnames(merged), value = TRUE)
merged$detect_mean_across <- rowMeans(merged[, detect_cols], na.rm = TRUE)

hk_top <- merged %>%
  arrange(desc(detect_mean_across)) %>%
  slice(1:200)

out_csv <- file.path(DIR_TABLES, "housekeeping_top200_by_detection.csv")
write.csv(hk_top, out_csv, row.names = FALSE)
log_msg(paste0("Wrote: ", out_csv), logfile)

# Plot: detection distributions (per dataset)
plot_df <- bind_rows(lapply(names(all_hk), function(nm) {
  all_hk[[nm]] %>% mutate(dataset = nm)
}))

p <- ggplot(plot_df, aes(x = detect_frac)) +
  geom_histogram(bins = 50) +
  facet_wrap(~dataset, scales = "free_y") +
  labs(title = "Detection frequency distributions (candidate housekeeping scan)",
       x = "fraction of cells with expression > 0", y = "#genes")
save_plot(p, file.path(DIR_FIGURES, "housekeeping_detection_distributions.png"), w = 10, h = 5)

=======
# ======================================================================
# scripts/08_metagene/08A_housekeeping_qc.R
#
# ADVANCED MODULE (optional): Housekeeping/QC gene selection.
#
# Motivation:
#   Housekeeping gene stability is useful for:
#     (1) detecting sample-to-sample technical shifts,
#     (2) confirming that normalization did not introduce artifacts,
#     (3) providing a stable baseline for metagene...

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

source("config/config.R")
source("scripts/R/utils.R")

ensure_dir(DIR_TABLES); ensure_dir(DIR_FIGURES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "08A_housekeeping_qc.log")

log_msg("Loading processed objects (for housekeeping QC).", logfile)

ref_path <- file.path(DIR_OBJS, "multiome_reference_processed.rds")
sl_path  <- file.path(DIR_OBJS, "slidetags_mapped_to_multiome.rds")
st_path  <- file.path(DIR_OBJS, "starmap_mapped_to_multiome.rds")

objs <- list()
if (file.exists(ref_path)) objs[["Multiome"]] <- readRDS(ref_path)
if (file.exists(sl_path))  objs[["Slide-tags"]] <- readRDS(sl_path)
if (file.exists(st_path))  objs[["STARmap"]] <- readRDS(st_path)

if (length(objs) == 0) {
  stop("No processed objects found. Run the minimal pipeline first.")
}

# --- Helper: compute housekeeping candidates within an assay ---
hk_candidates <- function(obj, assay = NULL, layer = "data") {
  assay <- assay %||% DefaultAssay(obj)
  safe_join_layers(obj, assay = assay)
  DefaultAssay(obj) <- assay
  if (!has_data_layer(obj, assay, layer = layer)) {
    obj <- NormalizeData(obj, verbose = FALSE)
  }
  mat <- Seurat::GetAssayData(obj, assay = assay, layer = layer)
  # Detection frequency (fraction of cells >0) and mean expression
  det <- Matrix::rowSums(mat > 0) / ncol(mat)
  mu  <- Matrix::rowMeans(mat)
  tibble::tibble(gene = names(mu), mean = as.numeric(mu), detect_frac = as.numeric(det)) %>%
    mutate(rank = dplyr::min_rank(desc(detect_frac)))
}

all_hk <- list()
for (nm in names(objs)) {
  obj <- objs[[nm]]
  assay_use <- if (nm == "STARmap" && "imputed" %in% names(obj@assays)) "imputed" else DefaultAssay(obj)
  log_msg(paste0("Housekeeping scan: ", nm, " (assay=", assay_use, ")"), logfile)
  all_hk[[nm]] <- hk_candidates(obj, assay = assay_use)
}

# Combine and identify robust genes: high detection across datasets
merged <- Reduce(function(a, b) full_join(a, b, by = "gene", suffix = c("", "_y")),
                 lapply(names(all_hk), function(nm) {
                   all_hk[[nm]] %>%
                     select(gene, detect_frac, mean) %>%
                     rename_with(~paste0(., "_", nm), c(detect_frac, mean))
                 }))

# A simple robustness score: average detection across available datasets
detect_cols <- grep("^detect_frac_", colnames(merged), value = TRUE)
merged$detect_mean_across <- rowMeans(merged[, detect_cols], na.rm = TRUE)

hk_top <- merged %>%
  arrange(desc(detect_mean_across)) %>%
  slice(1:200)

out_csv <- file.path(DIR_TABLES, "housekeeping_top200_by_detection.csv")
write.csv(hk_top, out_csv, row.names = FALSE)
log_msg(paste0("Wrote: ", out_csv), logfile)

# Plot: detection distributions (per dataset)
plot_df <- bind_rows(lapply(names(all_hk), function(nm) {
  all_hk[[nm]] %>% mutate(dataset = nm)
}))

p <- ggplot(plot_df, aes(x = detect_frac)) +
  geom_histogram(bins = 50) +
  facet_wrap(~dataset, scales = "free_y") +
  labs(title = "Detection frequency distributions (candidate housekeeping scan)",
       x = "fraction of cells with expression > 0", y = "#genes")
save_plot(p, file.path(DIR_FIGURES, "housekeeping_detection_distributions.png"), w = 10, h = 5)

>>>>>>> 54c1a5c675ca223d183141f8a75de090237ce902
log_msg("08A complete.", logfile)