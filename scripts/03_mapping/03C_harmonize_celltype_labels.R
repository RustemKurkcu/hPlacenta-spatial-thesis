<<<<<<< HEAD
# ======================================================================
# scripts/03_mapping/03C_harmonize_celltype_labels.R
# Create "celltype_final" without deleting any original/author annotations.
#
# KEY FIXES:
#   1. Adaptive AddModuleScore (no more "insufficient bins" error)
#   2. Smart assay selection for STARmap (prefer imputed)
#   3. JoinLayers before scoring
#   4. Graceful degradation if genes are missing
#
# Rules (conservative):
#   1) If author label exists and is not unknown -> keep as celltype_final
#   2) Else, if predicted.id exists and score >= PRED_SCORE_HIGH -> use it
#   3) Else -> "unknown"
#
# Additional: "lineage rescue" for unknown cells via marker scoring.
#
# Output:
#   output/objects/slidetags_harmonized.rds
#   output/objects/starmap_harmonized.rds
#   output/tables/*_celltype_harmonization_summary.csv
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

source("config/config.R")
source("scripts/R/utils.R")

ensure_dir(DIR_OBJS); ensure_dir(DIR_TABLES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "03C_harmonize_celltype_labels.log")

log_msg("Starting cell type harmonization...", logfile)

# ---- Load mapped objects ----
slide_path <- file.path(DIR_OBJS, "slidetags_mapped.rds")
star_path  <- file.path(DIR_OBJS, "starmap_mapped.rds")

if (!file.exists(slide_path)) stop("Missing: ", slide_path, " (run 03A first)")
if (!file.exists(star_path))  stop("Missing: ", star_path,  " (run 03B first)")

log_msg(paste0("Loading Slide-tags: ", slide_path), logfile)
slide <- readRDS(slide_path)

log_msg(paste0("Loading STARmap: ", star_path), logfile)
star <- readRDS(star_path)

# ======================================================================
# HARMONIZATION FUNCTION
# ======================================================================
harmonize_one <- function(obj, dataset_name) {
  log_msg(paste0("\n=== Harmonizing ", dataset_name, " ==="), logfile)
  log_msg(sprintf("  Cells: %d", ncol(obj)), logfile)

  md <- obj@meta.data

  # Ensure columns exist
  if (!("celltype_author" %in% colnames(md))) md$celltype_author <- NA_character_
  if (!(COL_PRED_CELLTYPE %in% colnames(md))) md[[COL_PRED_CELLTYPE]] <- NA_character_

  # prediction.score.max
  if (!(COL_PRED_SCORE_MAX %in% colnames(md))) {
    score_cols <- grep("^prediction\\.score\\.", colnames(md), value = TRUE)
    if (length(score_cols) > 0) {
      md[[COL_PRED_SCORE_MAX]] <- apply(md[, score_cols, drop = FALSE], 1, max, na.rm = TRUE)
    } else {
      # Try predicted.celltype.score (STARmap)
      if ("predicted.celltype.score" %in% colnames(md)) {
        md[[COL_PRED_SCORE_MAX]] <- as.numeric(md$predicted.celltype.score)
      } else {
        md[[COL_PRED_SCORE_MAX]] <- NA_real_
      }
    }
  }

  # Unknown author labels
  is_unknown_author <- is.na(md$celltype_author) |
    md$celltype_author %in% c("", "NA", "Unknown", "unknown", "unassigned", "Unassigned")

  # Conservative final
  md$celltype_final_conservative <- ifelse(
    !is_unknown_author, md$celltype_author,
    ifelse(!is.na(md[[COL_PRED_CELLTYPE]]) &
             !is.na(md[[COL_PRED_SCORE_MAX]]) &
             md[[COL_PRED_SCORE_MAX]] >= PRED_SCORE_HIGH,
           as.character(md[[COL_PRED_CELLTYPE]]),
           "unknown"))

  md$celltype_final_source <- ifelse(
    !is_unknown_author, "author",
    ifelse(md$celltype_final_conservative != "unknown", "predicted", "unknown"))

  log_msg(sprintf("  Author labels: %d, Predicted: %d, Unknown: %d",
                  sum(md$celltype_final_source == "author"),
                  sum(md$celltype_final_source == "predicted"),
                  sum(md$celltype_final_source == "unknown")), logfile)

  # Optional refinement using a mapping table
  map_path <- file.path("config", "celltype_refinement_map.csv")
  md$celltype_final_refined <- md$celltype_final_conservative
  md$celltype_refined_source <- md$celltype_final_source

  if (file.exists(map_path)) {
    log_msg("  Loading refinement map...", logfile)
    map <- readr::read_csv(map_path, show_col_types = FALSE)
    if (all(c("author_label", "allow_refine") %in% colnames(map))) {
      allow <- map$author_label[map$allow_refine %in% c(TRUE, "TRUE", "True", 1)]
      can_refine <- md$celltype_author %in% allow &
        !is.na(md[[COL_PRED_SCORE_MAX]]) &
        md[[COL_PRED_SCORE_MAX]] >= PRED_SCORE_REFINE &
        !is.na(md[[COL_PRED_CELLTYPE]])
      md$celltype_final_refined[can_refine] <- as.character(md[[COL_PRED_CELLTYPE]])[can_refine]
      md$celltype_refined_source[can_refine] <- "author_refined_to_predicted"
      log_msg(sprintf("  Refined %d cells using mapping table", sum(can_refine)), logfile)
    }
  }

  obj@meta.data <- md

  # ---- Lineage rescue for unknown cells ----
  log_msg("  Computing lineage scores...", logfile)

  # Choose assay -- prefer imputed for STARmap
  if (dataset_name == "STARmap") {
    assay_use <- select_starmap_assay(obj, prefer_imputed = TRUE)
  } else {
    assay_use <- if ("SCT" %in% names(obj@assays)) "SCT" else "RNA"
  }

  log_msg(paste0("  Using assay: ", assay_use), logfile)
  DefaultAssay(obj) <- assay_use

  # Join layers before scoring

  obj <- safe_join_layers(obj, assay_use)

  # Ensure data layer exists
  if (!has_data_layer(obj, assay_use)) {
    log_msg("  Normalizing data for scoring...", logfile)
    obj <- NormalizeData(obj, verbose = FALSE)
  }

  # Score each lineage
  for (nm in names(LINEAGE_MARKERS)) {
    log_msg(sprintf("  Scoring lineage: %s (%d markers)",
                    nm, length(LINEAGE_MARKERS[[nm]])), logfile)
    obj <- add_module_score_safe(
      obj, LINEAGE_MARKERS[[nm]],
      score_name = paste0("lineage_", nm),
      assay = assay_use, seed = SEED)
  }

  md <- obj@meta.data
  lineage_cols <- grep("^lineage_", colnames(md), value = TRUE)

  if (length(lineage_cols) == 0) {
    warning("No lineage scores computed -- all genes may be missing")
    md$celltype_lineage_rescue <- NA_character_
    md$lineage_score_max <- NA_real_
    md$lineage_label_max <- NA_character_
  } else {
    m <- as.matrix(md[, lineage_cols, drop = FALSE])
    # Replace NA with -Inf for max computation
    m[is.na(m)] <- -Inf
    max_score <- apply(m, 1, max)
    max_score[!is.finite(max_score)] <- NA_real_
    max_idx <- apply(m, 1, which.max)
    max_lab <- gsub("^lineage_", "", lineage_cols[max_idx])

    # Threshold: top quartile among finite scores
    finite_scores <- max_score[is.finite(max_score)]
    thr <- if (length(finite_scores) > 0) quantile(finite_scores, 0.75, na.rm = TRUE) else Inf

    md$celltype_lineage_rescue <- ifelse(
      md$celltype_final_conservative == "unknown" &
        !is.na(max_score) & max_score >= thr,
      paste0("lineage_", max_lab),
      NA_character_)
    md$lineage_score_max <- max_score
    md$lineage_label_max <- max_lab

    n_rescued <- sum(!is.na(md$celltype_lineage_rescue))
    log_msg(sprintf("  Rescued %d unknown cells with lineage labels", n_rescued), logfile)
  }

  obj@meta.data <- md

  # Export summary
  log_msg("  Exporting harmonization summary...", logfile)
  tab <- md %>%
    count(celltype_author, .data[[COL_PRED_CELLTYPE]],
          celltype_final_conservative, celltype_final_refined,
          celltype_lineage_rescue, name = "n_cells") %>%
    arrange(desc(n_cells))

  out_csv <- file.path(DIR_TABLES, paste0(dataset_name, "_celltype_harmonization_summary.csv"))
  write.csv(tab, out_csv, row.names = FALSE)
  log_msg(paste0("  Saved: ", out_csv), logfile)

  log_msg(sprintf("  Final: %d cells, %d conservative types, %d refined types",
                  ncol(obj),
                  length(unique(md$celltype_final_conservative)),
                  length(unique(md$celltype_final_refined))), logfile)
  obj
}

# ---- Process both datasets ----
slide2 <- harmonize_one(slide, "SlideTags")
star2  <- harmonize_one(star,  "STARmap")

# ---- Save ----
log_msg("\nSaving harmonized objects...", logfile)

slide_out <- file.path(DIR_OBJS, "slidetags_harmonized.rds")
saveRDS(slide2, slide_out)
log_msg(paste0("Saved: ", slide_out), logfile)

star_out <- file.path(DIR_OBJS, "starmap_harmonized.rds")
saveRDS(star2, star_out)
log_msg(paste0("Saved: ", star_out), logfile)

log_msg("\n03C complete.", logfile)

cat("\n", strrep("=", 70), "\n")
cat("Cell Type Harmonization Complete\n")
cat(strrep("=", 70), "\n")
cat(sprintf("Slide-tags: %d cells, %d types\n",
            ncol(slide2), length(unique(slide2$celltype_final_refined))))
cat(sprintf("STARmap:    %d cells, %d types\n",
            ncol(star2), length(unique(star2$celltype_final_refined))))
cat(sprintf("\nOutputs:\n  %s\n  %s\n", slide_out, star_out))
=======
# ======================================================================
# scripts/03_mapping/03C_harmonize_celltype_labels.R
# Create "celltype_final" without deleting any original/author annotations.
#
# KEY FIXES:
#   1. Adaptive AddModuleScore (no more "insufficient bins" error)
#   2. Smart assay selection for STARmap (prefer imputed)
#   3. JoinLayers before scoring
#   4. Graceful degradation if genes are missing
#
# Rules (conservative):
#   1) If author label exists and is not unknown -> keep as celltype_final
#   2) Else, if predicted.id exists and score >= PRED_SCORE_HIGH -> use it
#   3) Else -> "unknown"
#
# Additional: "lineage rescue" for unknown cells via marker scoring.
#
# Output:
#   output/objects/slidetags_harmonized.rds
#   output/objects/starmap_harmonized.rds
#   output/tables/*_celltype_harmonization_summary.csv
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

source("config/config.R")
source("scripts/R/utils.R")

ensure_dir(DIR_OBJS); ensure_dir(DIR_TABLES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "03C_harmonize_celltype_labels.log")

log_msg("Starting cell type harmonization...", logfile)

# ---- Load mapped objects ----
slide_path <- file.path(DIR_OBJS, "slidetags_mapped.rds")
star_path  <- file.path(DIR_OBJS, "starmap_mapped.rds")

if (!file.exists(slide_path)) stop("Missing: ", slide_path, " (run 03A first)")
if (!file.exists(star_path))  stop("Missing: ", star_path,  " (run 03B first)")

log_msg(paste0("Loading Slide-tags: ", slide_path), logfile)
slide <- readRDS(slide_path)

log_msg(paste0("Loading STARmap: ", star_path), logfile)
star <- readRDS(star_path)

# ======================================================================
# HARMONIZATION FUNCTION
# ======================================================================
harmonize_one <- function(obj, dataset_name) {
  log_msg(paste0("\n=== Harmonizing ", dataset_name, " ==="), logfile)
  log_msg(sprintf("  Cells: %d", ncol(obj)), logfile)

  md <- obj@meta.data

  # Ensure columns exist
  if (!("celltype_author" %in% colnames(md))) md$celltype_author <- NA_character_
  if (!(COL_PRED_CELLTYPE %in% colnames(md))) md[[COL_PRED_CELLTYPE]] <- NA_character_

  # prediction.score.max
  if (!(COL_PRED_SCORE_MAX %in% colnames(md))) {
    score_cols <- grep("^prediction\\.score\\.", colnames(md), value = TRUE)
    if (length(score_cols) > 0) {
      md[[COL_PRED_SCORE_MAX]] <- apply(md[, score_cols, drop = FALSE], 1, max, na.rm = TRUE)
    } else {
      # Try predicted.celltype.score (STARmap)
      if ("predicted.celltype.score" %in% colnames(md)) {
        md[[COL_PRED_SCORE_MAX]] <- as.numeric(md$predicted.celltype.score)
      } else {
        md[[COL_PRED_SCORE_MAX]] <- NA_real_
      }
    }
  }

  # Unknown author labels
  is_unknown_author <- is.na(md$celltype_author) |
    md$celltype_author %in% c("", "NA", "Unknown", "unknown", "unassigned", "Unassigned")

  # Conservative final
  md$celltype_final_conservative <- ifelse(
    !is_unknown_author, md$celltype_author,
    ifelse(!is.na(md[[COL_PRED_CELLTYPE]]) &
             !is.na(md[[COL_PRED_SCORE_MAX]]) &
             md[[COL_PRED_SCORE_MAX]] >= PRED_SCORE_HIGH,
           as.character(md[[COL_PRED_CELLTYPE]]),
           "unknown"))

  md$celltype_final_source <- ifelse(
    !is_unknown_author, "author",
    ifelse(md$celltype_final_conservative != "unknown", "predicted", "unknown"))

  log_msg(sprintf("  Author labels: %d, Predicted: %d, Unknown: %d",
                  sum(md$celltype_final_source == "author"),
                  sum(md$celltype_final_source == "predicted"),
                  sum(md$celltype_final_source == "unknown")), logfile)

  # Optional refinement using a mapping table
  map_path <- file.path("config", "celltype_refinement_map.csv")
  md$celltype_final_refined <- md$celltype_final_conservative
  md$celltype_refined_source <- md$celltype_final_source

  if (file.exists(map_path)) {
    log_msg("  Loading refinement map...", logfile)
    map <- readr::read_csv(map_path, show_col_types = FALSE)
    if (all(c("author_label", "allow_refine") %in% colnames(map))) {
      allow <- map$author_label[map$allow_refine %in% c(TRUE, "TRUE", "True", 1)]
      can_refine <- md$celltype_author %in% allow &
        !is.na(md[[COL_PRED_SCORE_MAX]]) &
        md[[COL_PRED_SCORE_MAX]] >= PRED_SCORE_REFINE &
        !is.na(md[[COL_PRED_CELLTYPE]])
      md$celltype_final_refined[can_refine] <- as.character(md[[COL_PRED_CELLTYPE]])[can_refine]
      md$celltype_refined_source[can_refine] <- "author_refined_to_predicted"
      log_msg(sprintf("  Refined %d cells using mapping table", sum(can_refine)), logfile)
    }
  }

  obj@meta.data <- md

  # ---- Lineage rescue for unknown cells ----
  log_msg("  Computing lineage scores...", logfile)

  # Choose assay -- prefer imputed for STARmap
  if (dataset_name == "STARmap") {
    assay_use <- select_starmap_assay(obj, prefer_imputed = TRUE)
  } else {
    assay_use <- if ("SCT" %in% names(obj@assays)) "SCT" else "RNA"
  }

  log_msg(paste0("  Using assay: ", assay_use), logfile)
  DefaultAssay(obj) <- assay_use

  # Join layers before scoring

  obj <- safe_join_layers(obj, assay_use)

  # Ensure data layer exists
  if (!has_data_layer(obj, assay_use)) {
    log_msg("  Normalizing data for scoring...", logfile)
    obj <- NormalizeData(obj, verbose = FALSE)
  }

  # Score each lineage
  for (nm in names(LINEAGE_MARKERS)) {
    log_msg(sprintf("  Scoring lineage: %s (%d markers)",
                    nm, length(LINEAGE_MARKERS[[nm]])), logfile)
    obj <- add_module_score_safe(
      obj, LINEAGE_MARKERS[[nm]],
      score_name = paste0("lineage_", nm),
      assay = assay_use, seed = SEED)
  }

  md <- obj@meta.data
  lineage_cols <- grep("^lineage_", colnames(md), value = TRUE)

  if (length(lineage_cols) == 0) {
    warning("No lineage scores computed -- all genes may be missing")
    md$celltype_lineage_rescue <- NA_character_
    md$lineage_score_max <- NA_real_
    md$lineage_label_max <- NA_character_
  } else {
    m <- as.matrix(md[, lineage_cols, drop = FALSE])
    # Replace NA with -Inf for max computation
    m[is.na(m)] <- -Inf
    max_score <- apply(m, 1, max)
    max_score[!is.finite(max_score)] <- NA_real_
    max_idx <- apply(m, 1, which.max)
    max_lab <- gsub("^lineage_", "", lineage_cols[max_idx])

    # Threshold: top quartile among finite scores
    finite_scores <- max_score[is.finite(max_score)]
    thr <- if (length(finite_scores) > 0) quantile(finite_scores, 0.75, na.rm = TRUE) else Inf

    md$celltype_lineage_rescue <- ifelse(
      md$celltype_final_conservative == "unknown" &
        !is.na(max_score) & max_score >= thr,
      paste0("lineage_", max_lab),
      NA_character_)
    md$lineage_score_max <- max_score
    md$lineage_label_max <- max_lab

    n_rescued <- sum(!is.na(md$celltype_lineage_rescue))
    log_msg(sprintf("  Rescued %d unknown cells with lineage labels", n_rescued), logfile)
  }

  obj@meta.data <- md

  # Export summary
  log_msg("  Exporting harmonization summary...", logfile)
  tab <- md %>%
    count(celltype_author, .data[[COL_PRED_CELLTYPE]],
          celltype_final_conservative, celltype_final_refined,
          celltype_lineage_rescue, name = "n_cells") %>%
    arrange(desc(n_cells))

  out_csv <- file.path(DIR_TABLES, paste0(dataset_name, "_celltype_harmonization_summary.csv"))
  write.csv(tab, out_csv, row.names = FALSE)
  log_msg(paste0("  Saved: ", out_csv), logfile)

  log_msg(sprintf("  Final: %d cells, %d conservative types, %d refined types",
                  ncol(obj),
                  length(unique(md$celltype_final_conservative)),
                  length(unique(md$celltype_final_refined))), logfile)
  obj
}

# ---- Process both datasets ----
slide2 <- harmonize_one(slide, "SlideTags")
star2  <- harmonize_one(star,  "STARmap")

# ---- Save ----
log_msg("\nSaving harmonized objects...", logfile)

slide_out <- file.path(DIR_OBJS, "slidetags_harmonized.rds")
saveRDS(slide2, slide_out)
log_msg(paste0("Saved: ", slide_out), logfile)

star_out <- file.path(DIR_OBJS, "starmap_harmonized.rds")
saveRDS(star2, star_out)
log_msg(paste0("Saved: ", star_out), logfile)

log_msg("\n03C complete.", logfile)

cat("\n", strrep("=", 70), "\n")
cat("Cell Type Harmonization Complete\n")
cat(strrep("=", 70), "\n")
cat(sprintf("Slide-tags: %d cells, %d types\n",
            ncol(slide2), length(unique(slide2$celltype_final_refined))))
cat(sprintf("STARmap:    %d cells, %d types\n",
            ncol(star2), length(unique(star2$celltype_final_refined))))
cat(sprintf("\nOutputs:\n  %s\n  %s\n", slide_out, star_out))
>>>>>>> 54c1a5c675ca223d183141f8a75de090237ce902
cat(strrep("=", 70), "\n\n")