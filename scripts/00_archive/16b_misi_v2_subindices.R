# =============================================================================
# scripts/16b_misi_v2_subindices.R
# MISI v2: 4-Subindex Decomposition, Bootstrap CIs, Gene Jackknife, Spatial
# =============================================================================
#
# PURPOSE:
#   Compute the MISI v2 modular subindices for every cell, attach them to the
#   Seurat object, and export comprehensive tables with uncertainty (bootstrap
#   CIs), sensitivity (gene jackknife), and optional spatial calibration
#   (Moran's I radius sweep).
#
# INPUTS:
#   - outputs/objects/seu_with_misi.qs  (from script 16_compute_misi.R)
#
# OUTPUTS:
#   - outputs/objects/seu_with_misi_v2.qs
#   - outputs/tables/misi_v2_gene_presence.csv
#   - outputs/tables/misi_v2_component_scores_per_cell.csv
#   - outputs/tables/misi_v2_summary_by_celltype_condition.csv
#   - outputs/tables/misi_v2_bootstrap_ci.csv
#   - outputs/tables/misi_v2_gene_jackknife.csv
#   - outputs/tables/misi_v2_module_correlations.csv
#   - outputs/tables/misi_v2_weight_sensitivity.csv
#   - outputs/tables/misi_v2_spatial_radius_sweep.csv  (if coords exist)
#
# DEPENDENCIES:
#   - R/config.R, R/helpers_io.R, R/helpers_methods.R (pipeline standard)
#   - R/misi_v2_subindices.R (this PR)
#   - UCell (Bioconductor), dplyr, tibble
#   - Optional: spdep (for spatial sweep), irlba (for PC1 scoring)
#
# EXECUTION:
#   Rscript scripts/16b_misi_v2_subindices.R
#
# AUTHOR:  Shan Kurkcu
# DATE:    2025-04-01
# VERSION: 1.0.0
# =============================================================================

cat("
╔══════════════════════════════════════════════════════════════════╗
║  16b — MISI v2: 4-Subindex Decomposition + Uncertainty + QC    ║
╚══════════════════════════════════════════════════════════════════╝
")

# =============================================================================
# 0. SETUP
# =============================================================================

source("R/config.R")
source("R/helpers_io.R")
source("R/helpers_methods.R")
source("R/misi_v2_subindices.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(qs)
})

log_file <- file.path(CFG$dirs$logs, "pipeline.log")
log_msg("16b_misi_v2_subindices starting.", log_file = log_file)

# --- Configuration ---
SCORING_METHOD  <- "ucell"   # "ucell" (default, recommended) or "pc1"
BOOTSTRAP_B     <- 500L      # Number of bootstrap resamples
BOOTSTRAP_SEED  <- 42L       # Seed for reproducibility
MIN_CELLS_CI    <- 20L       # Minimum cells per stratum for valid CI
NCORES_UCELL    <- 1L        # Cores for UCell (increase on HPC)
RUN_JACKKNIFE   <- TRUE      # Gene leave-one-out jackknife (slow but important)
RUN_SPATIAL     <- TRUE      # Moran's I spatial sweep (only if coords exist)
SPATIAL_RADII   <- c(25, 50, 100, 200, 400)

# --- Directory setup ---
ensure_dir(CFG$dirs$tables)
ensure_dir(CFG$dirs$figures)
ensure_dir(CFG$dirs$objects)
dir_objects <- CFG$dirs$objects

# --- safe_write_csv (matching Script 16 convention) ---
safe_write_csv <- function(df, path) {
  tryCatch(readr::write_csv(df, path), error = function(e) {
    write.csv(df, path, row.names = FALSE)
  })
}

log_msg("  Scoring method: ", SCORING_METHOD, log_file = log_file)
log_msg("  Bootstrap B: ", BOOTSTRAP_B, log_file = log_file)
log_msg("  Jackknife: ", RUN_JACKKNIFE, log_file = log_file)
log_msg("  Spatial sweep: ", RUN_SPATIAL, log_file = log_file)


# =============================================================================
# 1. LOAD SEURAT OBJECT
# =============================================================================

log_msg("Loading Seurat object...", log_file = log_file)

# Try candidates in priority order (Script 16 outputs seu_with_misi.qs)
obj_candidates <- c(
  file.path(dir_objects, "seu_with_misi.qs")
)

seu <- NULL
obj_path <- NULL
for (cand in obj_candidates) {
  if (file.exists(cand)) {
    s <- tryCatch(qs::qread(cand, nthreads = 1), error = function(e) NULL)
    if (!is.null(s) && inherits(s, "Seurat")) {
      seu <- s
      obj_path <- cand
      break
    }
  }
}
if (is.null(seu)) {
  stop("No suitable Seurat object found. Tried:\n  ",
       paste(obj_candidates, collapse = "\n  "))
}

log_msg("  Loaded: ", obj_path, log_file = log_file)
log_msg("  Cells: ", ncol(seu), " | Genes: ", nrow(seu), log_file = log_file)
DefaultAssay(seu) <- "RNA"


# =============================================================================
# 2. RESOLVE METADATA COLUMNS (matching Script 16 conventions)
# =============================================================================

# Use CFG$cols for standard column names (same as Script 16)
ct_col    <- CFG$cols$cell_type   # "cell_type"
donor_col <- CFG$cols$donor_id    # "donor_id"

# cell_type: Script 16 ensures this exists (with fallback to celltype_refined/celltype_author)
if (!ct_col %in% colnames(seu@meta.data)) {
  if ("celltype_refined" %in% colnames(seu@meta.data)) {
    seu[[ct_col]] <- seu$celltype_refined
  } else if ("celltype_author" %in% colnames(seu@meta.data)) {
    seu[[ct_col]] <- seu$celltype_author
  } else if ("seurat_clusters" %in% colnames(seu@meta.data)) {
    seu[[ct_col]] <- seu$seurat_clusters
    log_msg("  WARNING: Using seurat_clusters as cell_type fallback.", log_file = log_file)
  } else {
    stop("No cell-type column found in metadata.")
  }
}

# donor_id: Script 16 ensures this exists
if (!donor_col %in% colnames(seu@meta.data)) {
  seu[[donor_col]] <- "donor_unknown"
  log_msg("  WARNING: No donor_id found; set to 'donor_unknown'.", log_file = log_file)
}

# condition: Script 16 creates this as paste0(infection, "_", hpi)
if (!"condition" %in% colnames(seu@meta.data)) {
  if (all(c(CFG$cols$infection, CFG$cols$hpi) %in% colnames(seu@meta.data))) {
    seu$condition <- paste0(
      seu[[CFG$cols$infection]][, 1], "_",
      seu[[CFG$cols$hpi]][, 1]
    )
    log_msg("  Created 'condition' from infection + hpi.", log_file = log_file)
  } else {
    seu$condition <- "unknown"
    log_msg("  WARNING: Cannot create condition; set to 'unknown'.", log_file = log_file)
  }
}
cond_col <- "condition"

log_msg("  Donor column: ", donor_col, " (", length(unique(seu@meta.data[[donor_col]])),
        " donors)", log_file = log_file)
log_msg("  CellType column: ", ct_col, " (", length(unique(seu@meta.data[[ct_col]])),
        " types)", log_file = log_file)
log_msg("  Condition column: ", cond_col, " (", length(unique(seu@meta.data[[cond_col]])),
        " conditions)", log_file = log_file)


# =============================================================================
# 3. COMPUTE MISI v2 MODULE SCORES
# =============================================================================

log_msg("Computing MISI v2 module scores (", SCORING_METHOD, ")...", log_file = log_file)

misi_result <- run_misi_v2(
  seu,
  scoring_method = SCORING_METHOD,
  assay          = "RNA",
  ncores         = NCORES_UCELL
)

raw_scores    <- misi_result$raw_scores
z_scores      <- misi_result$z_scores
subindices    <- misi_result$subindices
gene_presence <- misi_result$gene_presence

log_msg("  Modules scored: ", paste(colnames(raw_scores), collapse = ", "),
        log_file = log_file)
log_msg("  Gene presence: ", sum(gene_presence$present), "/",
        nrow(gene_presence), " genes detected", log_file = log_file)

# --- Save gene presence report ---
safe_write_csv(gene_presence, file.path(CFG$dirs$tables, "misi_v2_gene_presence.csv"))
log_msg("  Saved: misi_v2_gene_presence.csv", log_file = log_file)


# =============================================================================
# 4. ATTACH SCORES TO SEURAT OBJECT
# =============================================================================

log_msg("Attaching MISI v2 scores to Seurat metadata...", log_file = log_file)

# Module z-scores (prefixed MISIv2_ to avoid collision with Script 16's MISI column)
for (mod in colnames(z_scores)) {
  seu@meta.data[[paste0("MISIv2_", mod, "_z")]] <- z_scores[[mod]]
}

# Subindices + overall MISI v2
for (si in colnames(subindices)) {
  seu@meta.data[[paste0("MISIv2_", si)]] <- subindices[[si]]
}

log_msg("  Added ", ncol(z_scores) + ncol(subindices),
        " MISIv2_ columns to metadata.", log_file = log_file)


# =============================================================================
# 5. BUILD PER-CELL COMPONENT SCORES TABLE
# =============================================================================

log_msg("Building per-cell component scores table...", log_file = log_file)

meta <- seu@meta.data
cell_table <- tibble::tibble(
  cell       = colnames(seu),
  donor      = as.character(meta[[donor_col]]),
  cell_type  = as.character(meta[[ct_col]]),
  condition  = as.character(meta[[cond_col]])
)

# Add raw scores
for (mod in colnames(raw_scores)) {
  cell_table[[paste0(mod, "_raw")]] <- raw_scores[[mod]]
}
# Add z-scores
for (mod in colnames(z_scores)) {
  cell_table[[paste0(mod, "_z")]] <- z_scores[[mod]]
}
# Add subindices
for (si in colnames(subindices)) {
  cell_table[[si]] <- subindices[[si]]
}

safe_write_csv(cell_table, file.path(CFG$dirs$tables, "misi_v2_component_scores_per_cell.csv"))
log_msg("  Saved: misi_v2_component_scores_per_cell.csv (",
        nrow(cell_table), " cells × ", ncol(cell_table), " columns)",
        log_file = log_file)


# =============================================================================
# 6. SUMMARY TABLE BY CELLTYPE × CONDITION
# =============================================================================

log_msg("Computing summary statistics by celltype × condition...", log_file = log_file)

score_cols <- intersect(
  c("MISI", "MISI_Homing", "MISI_Nutrient", "MISI_SwitchPressure",
    "MISI_VascularFragility", "Invasion_MMP_composite"),
  colnames(subindices)
)

summary_input <- data.frame(
  cell_type = as.character(meta[[ct_col]]),
  condition = as.character(meta[[cond_col]]),
  stringsAsFactors = FALSE
)
for (sc in score_cols) summary_input[[sc]] <- subindices[[sc]]

summary_tbl <- summary_input %>%
  dplyr::group_by(cell_type, condition) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    dplyr::across(
      dplyr::all_of(score_cols),
      list(mean = ~mean(.x, na.rm = TRUE),
           sd   = ~sd(.x, na.rm = TRUE),
           median = ~median(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

safe_write_csv(summary_tbl, file.path(CFG$dirs$tables, "misi_v2_summary_by_celltype_condition.csv"))
log_msg("  Saved: misi_v2_summary_by_celltype_condition.csv (",
        nrow(summary_tbl), " strata)", log_file = log_file)


# =============================================================================
# 7. BOOTSTRAP CONFIDENCE INTERVALS
# =============================================================================

log_msg("Computing bootstrap CIs (B=", BOOTSTRAP_B, ")...", log_file = log_file)

boot_value_cols <- intersect(
  c("MISI", "MISI_Homing", "MISI_Nutrient",
    "MISI_SwitchPressure", "MISI_VascularFragility"),
  colnames(subindices)
)

boot_input <- data.frame(
  donor     = as.character(meta[[donor_col]]),
  cell_type = as.character(meta[[ct_col]]),
  condition = as.character(meta[[cond_col]]),
  stringsAsFactors = FALSE
)
for (vc in boot_value_cols) boot_input[[vc]] <- subindices[[vc]]

boot_ci <- bootstrap_group_ci(
  df         = boot_input,
  group_cols = c("donor", "cell_type", "condition"),
  value_cols = boot_value_cols,
  B          = BOOTSTRAP_B,
  min_cells  = MIN_CELLS_CI,
  seed       = BOOTSTRAP_SEED
)

safe_write_csv(boot_ci, file.path(CFG$dirs$tables, "misi_v2_bootstrap_ci.csv"))
log_msg("  Saved: misi_v2_bootstrap_ci.csv (", nrow(boot_ci), " rows)",
        log_file = log_file)

# Report strata with wide CIs
if ("ci_width" %in% colnames(boot_ci)) {
  wide_ci <- boot_ci %>%
    dplyr::filter(!is.na(ci_width)) %>%
    dplyr::arrange(dplyr::desc(ci_width))
  if (nrow(wide_ci) > 0) {
    log_msg("  Top 5 widest CIs (low-confidence strata):", log_file = log_file)
    for (i in seq_len(min(nrow(wide_ci), 5))) {
      log_msg("    ", wide_ci$score[i], " | width=",
              round(wide_ci$ci_width[i], 3), " | n=", wide_ci$n_cells[i],
              log_file = log_file)
    }
  }
}


# =============================================================================
# 8. GENE JACKKNIFE (LEAVE-ONE-OUT SENSITIVITY)
# =============================================================================

if (RUN_JACKKNIFE) {
  log_msg("Running gene jackknife analysis (this may take several minutes)...",
          log_file = log_file)

  jk_results <- jackknife_all_modules(
    seu            = seu,
    raw_scores     = raw_scores,
    scoring_method = SCORING_METHOD,
    assay          = "RNA",
    ncores         = NCORES_UCELL
  )

  safe_write_csv(jk_results, file.path(CFG$dirs$tables, "misi_v2_gene_jackknife.csv"))
  log_msg("  Saved: misi_v2_gene_jackknife.csv (", nrow(jk_results),
          " gene-module pairs)", log_file = log_file)

  # Flag influential genes (correlation with full < 0.9)
  influential <- jk_results %>%
    dplyr::filter(!is.na(cor_with_full) & cor_with_full < 0.9) %>%
    dplyr::arrange(cor_with_full)

  if (nrow(influential) > 0) {
    log_msg("  WARNING: Influential genes (LOO cor < 0.9):", log_file = log_file)
    for (i in seq_len(min(nrow(influential), 10))) {
      log_msg("    ", influential$module[i], " -> ",
              influential$dropped_gene[i], " (r=",
              round(influential$cor_with_full[i], 3), ")",
              log_file = log_file)
    }
  } else {
    log_msg("  All genes show LOO correlation >= 0.9 — no single-gene dominance.",
            log_file = log_file)
  }
} else {
  log_msg("Skipping gene jackknife (RUN_JACKKNIFE = FALSE).", log_file = log_file)
}


# =============================================================================
# 9. MODULE CORRELATION QC
# =============================================================================

log_msg("Computing module correlation matrix (QC)...", log_file = log_file)

cor_qc <- module_correlation_qc(raw_scores, method = "spearman")

safe_write_csv(
  tibble::rownames_to_column(as.data.frame(cor_qc$cor_matrix), var = "module"),
  file.path(CFG$dirs$tables, "misi_v2_module_correlations.csv")
)
log_msg("  Saved: misi_v2_module_correlations.csv", log_file = log_file)

if (nrow(cor_qc$flagged_pairs) > 0) {
  log_msg("  WARNING: Highly correlated module pairs (|rho| > 0.8):",
          log_file = log_file)
  for (i in seq_len(nrow(cor_qc$flagged_pairs))) {
    fp <- cor_qc$flagged_pairs
    log_msg("    ", fp$module_1[i], " <-> ", fp$module_2[i],
            " (rho=", round(fp$correlation[i], 3), ")",
            log_file = log_file)
  }
} else {
  log_msg("  No module pairs with |rho| > 0.8 — acceptable independence.",
          log_file = log_file)
}


# =============================================================================
# 10. WEIGHT SENSITIVITY COMPARISON
# =============================================================================

log_msg("Running weight sensitivity comparison...", log_file = log_file)

# Fresh z-scores (compose_misi mutates the input)
z_scores_fresh <- zscore_cols(raw_scores)

# Literature-informed prior weights (emphasize nutrient + homing for Fn)
prior_weights <- c(
  MISI_Homing            = 1.5,
  MISI_Nutrient          = 1.5,
  MISI_SwitchPressure    = 1.0,
  MISI_VascularFragility = 1.0,
  Invasion_composite     = 0.75
)

wt_sensitivity <- weight_sensitivity_compare(
  z_scores_fresh,
  learned_weights = NULL,
  prior_weights   = prior_weights
)

safe_write_csv(wt_sensitivity$summary,
               file.path(CFG$dirs$tables, "misi_v2_weight_sensitivity.csv"))
log_msg("  Saved: misi_v2_weight_sensitivity.csv", log_file = log_file)

# Report correlations between weighting schemes
cor_schemes <- wt_sensitivity$correlation_matrix
for (i in seq_len(nrow(cor_schemes))) {
  for (j in seq_len(ncol(cor_schemes))) {
    if (i < j) {
      log_msg("  ", rownames(cor_schemes)[i], " <-> ", colnames(cor_schemes)[j],
              " r=", round(cor_schemes[i, j], 4), log_file = log_file)
    }
  }
}


# =============================================================================
# 11. SPATIAL RADIUS SWEEP (OPTIONAL)
# =============================================================================

has_spatial <- FALSE
if (RUN_SPATIAL) {
  log_msg("Checking for spatial coordinates...", log_file = log_file)

  xy <- NULL

  # Check for spatial reduction
  if ("spatial" %in% names(seu@reductions)) {
    xy <- Embeddings(seu, "spatial")[, 1:2]
    has_spatial <- TRUE
    log_msg("  Found 'spatial' reduction.", log_file = log_file)
  }

  # Check metadata coordinate columns (matching Script 16's convention)
  if (!has_spatial) {
    coord_pairs <- list(
      c("x", "y"), c("X", "Y"),
      c("imagerow", "imagecol"),
      c("x_centroid", "y_centroid"),
      c("x_coord", "y_coord"),
      c("spatial_x", "spatial_y")
    )
    for (cp in coord_pairs) {
      if (all(cp %in% colnames(meta))) {
        xy <- as.matrix(meta[, cp])
        has_spatial <- TRUE
        log_msg("  Found spatial coordinates: ", paste(cp, collapse = ", "),
                log_file = log_file)
        break
      }
    }
  }

  if (has_spatial) {
    log_msg("  Running Moran's I sweep across radii: ",
            paste(SPATIAL_RADII, collapse = ", "), log_file = log_file)

    spatial_results <- spatial_sweep_all_scores(
      xy          = xy,
      subindex_df = subindices,
      radii       = SPATIAL_RADII
    )

    safe_write_csv(spatial_results,
                   file.path(CFG$dirs$tables, "misi_v2_spatial_radius_sweep.csv"))
    log_msg("  Saved: misi_v2_spatial_radius_sweep.csv", log_file = log_file)

    # Report peak Moran's I per score
    for (sc in unique(spatial_results$score)) {
      sc_data <- spatial_results %>%
        dplyr::filter(score == sc, !is.na(moran_I))
      if (nrow(sc_data) > 0) {
        peak <- sc_data[which.max(sc_data$moran_I), ]
        log_msg("    ", sc, ": peak I=", round(peak$moran_I, 4),
                " at radius=", peak$radius,
                " (p=", format(peak$p_value, digits = 3), ")",
                log_file = log_file)
      }
    }
  } else {
    log_msg("  No spatial coordinates found. Skipping spatial sweep.",
            log_file = log_file)
    log_msg("  (Expected for dissociated scRNA-seq data.)", log_file = log_file)
  }
} else {
  log_msg("Skipping spatial sweep (RUN_SPATIAL = FALSE).", log_file = log_file)
}


# =============================================================================
# 12. SAVE SEURAT OBJECT
# =============================================================================

log_msg("Saving Seurat object with MISI v2 scores...", log_file = log_file)
out_path <- file.path(CFG$dirs$objects, "seu_with_misi_v2.qs")
qs::qsave(seu, out_path, preset = "high")
log_msg("  Saved: ", out_path, log_file = log_file)


# =============================================================================
# 13. METHODS DRAFT
# =============================================================================

append_methods_draft(
  script_name = "16b_misi_v2_subindices.R",
  hypothesis = "Decomposing MISI into 4 mechanistic subindices (Homing, Nutrient, SwitchPressure, VascularFragility) enables pathway-specific vulnerability attribution and stabilizes the overall index through module-level uncertainty quantification.",
  method_details = list(
    "Score 11 mechanistic gene modules using UCell rank-based signatures (Andreatta & Carmona 2021)",
    "Z-score normalize per module; invert protective modules (Barrier, Defense) for vulnerability orientation",
    "Compose 4 subindices: MISI_Homing (Glyco + inverted Barrier), MISI_Nutrient (Liberation + Sequestration), MISI_SwitchPressure (Immune + Hypoxia), MISI_VascularFragility (Angiogenic + EndoDysfx)",
    "Invasion composite = MMP z-score minus TIMP z-score; overall MISI = mean of 5 components",
    "Cell-bootstrap CIs (B=500) per donor x celltype x condition stratum",
    "Leave-one-gene-out jackknife: recompute module scores dropping each gene, correlate with full",
    "Module correlation QC: Spearman rho matrix to detect redundancy",
    "Weight sensitivity: compare equal vs literature-prior weights",
    "Optional Moran's I spatial autocorrelation sweep at radii 25/50/100/200/400"
  ),
  rationale = "A single composite MISI number obscures which biological axis drives vulnerability. Decomposition into mechanistic subindices enables pathway-level hypothesis testing (e.g., 'Is Fn vulnerability primarily driven by nutrient availability or barrier failure?') and provides interpretable targets for experimental follow-up. Bootstrap CIs and gene jackknife transform MISI from a descriptive visualization into a defensible statistic.",
  citations = c(
    "Andreatta M, Carmona SJ. UCell: Robust and scalable single-cell gene signature scoring. Comput Struct Biotechnol J. 2021;19:3796-3798.",
    "Hao Y, et al. Integrated analysis of multimodal single-cell data. Cell. 2021;184(13):3573-3587.",
    "Palla G, et al. Squidpy: a scalable framework for spatial omics analysis. Nat Methods. 2022;19:171-178."
  )
)


# =============================================================================
# 14. FINAL REPORT
# =============================================================================

cat("
╔══════════════════════════════════════════════════════════════════╗
║  16b — MISI v2 COMPUTATION COMPLETE                            ║
╠══════════════════════════════════════════════════════════════════╣
")

log_msg("=== OUTPUT MANIFEST ===", log_file = log_file)
log_msg("  Object: ", out_path, log_file = log_file)
log_msg("  Tables:", log_file = log_file)

output_tables <- c(
  "misi_v2_gene_presence.csv",
  "misi_v2_component_scores_per_cell.csv",
  "misi_v2_summary_by_celltype_condition.csv",
  "misi_v2_bootstrap_ci.csv"
)
if (RUN_JACKKNIFE) output_tables <- c(output_tables, "misi_v2_gene_jackknife.csv")
output_tables <- c(output_tables,
                   "misi_v2_module_correlations.csv",
                   "misi_v2_weight_sensitivity.csv")
if (has_spatial) output_tables <- c(output_tables, "misi_v2_spatial_radius_sweep.csv")

for (tbl in output_tables) {
  full <- file.path(CFG$dirs$tables, tbl)
  if (file.exists(full)) {
    sz <- file.size(full)
    log_msg("    OK  ", tbl, " (", round(sz / 1024, 1), " KB)", log_file = log_file)
  } else {
    log_msg("    ERR ", tbl, " (NOT FOUND)", log_file = log_file)
  }
}

log_msg("", log_file = log_file)
log_msg("=== MISI v2 SUBINDEX SUMMARY ===", log_file = log_file)
for (si in c("MISI_Homing", "MISI_Nutrient", "MISI_SwitchPressure",
             "MISI_VascularFragility", "MISI")) {
  if (si %in% colnames(subindices)) {
    vals <- subindices[[si]]
    log_msg("  ", si, ": mean=", round(mean(vals, na.rm = TRUE), 4),
            " sd=", round(sd(vals, na.rm = TRUE), 4),
            " range=[", round(min(vals, na.rm = TRUE), 3), ", ",
            round(max(vals, na.rm = TRUE), 3), "]",
            log_file = log_file)
  }
}

cat("
╚══════════════════════════════════════════════════════════════════╝
")

log_msg("16b_misi_v2_subindices done. MISIv2_ columns added: ",
        ncol(z_scores) + ncol(subindices), log_file = log_file)
