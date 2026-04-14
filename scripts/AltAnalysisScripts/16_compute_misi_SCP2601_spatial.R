# =============================================================================
# scripts/16_compute_misi_SCP2601_spatial.R
# MISI v2 + IDO1_Tolerogenic_Shield for SCP2601 Spatial Data
# =============================================================================
#
# PURPOSE:
#   Adapt the canonical 16b_misi_v2_subindices.R pipeline to the SCP2601
#   spatial multiomic dataset (Ounadjela et al. 2024, Nature Medicine).
#   Adds the 12-gene IDO1_Tolerogenic_Shield module to the UCell scoring
#   block, and ensures the RUN_SPATIAL <- TRUE block executes the Moran's I
#   radius sweep over the physical STARmap/Slide-tag coordinates.
#
# KEY ADDITIONS vs. canonical 16b:
#   1. IDO1_Tolerogenic_Shield module (12 genes from Part 3 PDF):
#      IDO1, TGFB1, IL10, CD274 (PD-L1), PDCD1LG2 (PD-L2),
#      HAVCR2 (TIM-3), LGALS9, CD80, CD86, ENTPD1 (CD39),
#      NT5E (CD73), FOXP3
#   2. Physical-coordinate Moran's I (not UMAP) using "spatial" DimReduc
#   3. Spatial feature maps per subindex + IDO1 shield
#   4. Cross-sample (week) comparisons
#
# THESIS HYPOTHESIS (Guardrail 3):
#   EPAS1 (HIF-2α), NOT HIF-1α, drives FLT1/sFLT1 in trophoblasts.
#   The IDO1 Tolerogenic Shield collapses at the maternal-fetal interface
#   when Fn colonizes via Fap2/Gal-GalNAc → this script measures that
#   shield's spatial integrity.
#
# INPUT:
#   - output/objects/SCP2601_combined_spatial.rds  (from Script 05B)
#     OR: provide path via command-line arg
#
# OUTPUTS:
#   - output/objects/SCP2601_with_misi_v2.rds
#   - output/tables/SCP2601_misi_v2_gene_presence.csv
#   - output/tables/SCP2601_misi_v2_per_cell_scores.csv
#   - output/tables/SCP2601_misi_v2_summary_by_celltype_week.csv
#   - output/tables/SCP2601_misi_v2_spatial_radius_sweep.csv
#   - output/tables/SCP2601_ido1_shield_spatial_sweep.csv
#   - output/figures/SCP2601_misi_spatial_featuremaps.pdf
#   - output/figures/SCP2601_ido1_shield_spatial_map.pdf
#   - output/figures/SCP2601_misi_moranI_curves.pdf
#   - output/figures/SCP2601_misi_violins_by_week.pdf
#   - output/figures/SCP2601_misi_vs_ido1_scatter.pdf
#
# EXECUTION:
#   Rscript scripts/16_compute_misi_SCP2601_spatial.R
#
# AUTHOR:  Shan Kurkcu (PhD Thesis Pipeline)
# DATE:    2025-06-19
# VERSION: 1.0.0
# =============================================================================

cat("\n╔══════════════════════════════════════════════════════════════════════╗\n")
cat(  "║  16 — MISI v2 + IDO1 Shield for SCP2601 Spatial Data              ║\n")
cat(  "╚══════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# 0. SETUP
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
})

# --- Pipeline infrastructure ---
tryCatch({
  source("config/config.R")
  source("scripts/R/utils.R")
}, error = function(e) {
  cat("[16] Running in standalone mode.\n")
})

# --- Try to load canonical MISI module (from placenta-infection-pipeline) ---
misi_module_loaded <- FALSE
tryCatch({
  source("R/misi_v2_subindices.R")
  misi_module_loaded <- TRUE
  cat("[16] Canonical misi_v2_subindices.R loaded.\n")
}, error = function(e) {
  cat("[16] misi_v2_subindices.R not found — using embedded definitions.\n")
})

# --- Directories ---
DIR_OBJS    <- "output/objects"
DIR_TABLES  <- "output/tables"
DIR_FIGURES <- "output/figures"
DIR_LOGS    <- "output/logs"
for (d in c(DIR_OBJS, DIR_TABLES, DIR_FIGURES, DIR_LOGS)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

log_file <- file.path(DIR_LOGS, "16_compute_misi_SCP2601.log")
log_msg <- function(...) {
  msg <- paste0(Sys.time(), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}
log_msg("16_compute_misi_SCP2601_spatial starting.")

# --- Configuration ---
SCORING_METHOD  <- "ucell"
RUN_SPATIAL     <- TRUE     # FORCE TRUE — this is a spatial dataset
SPATIAL_RADII   <- c(25, 50, 100, 150, 200, 300, 400, 500)
SEED            <- 42
set.seed(SEED)


# =============================================================================
# 1. MODULE DEFINITIONS (Extended with IDO1_Tolerogenic_Shield)
# =============================================================================

log_msg("Defining MISI v2 modules + IDO1_Tolerogenic_Shield...")

# The canonical 11 MISI modules (from misi_v2_subindices.R)
if (misi_module_loaded && exists("misi_modules_v2", mode = "function")) {
  base_modules <- misi_modules_v2()
} else {
  # Embedded fallback (frozen copy from R/misi_v2_subindices.R)
  base_modules <- list(
    Nutrient_Liberation    = c("PLD1", "PLD2", "GDPD1", "ETNK1", "FAAH"),
    Nutrient_Sequestration = c("PCYT2", "CHPT1", "CEPT1"),
    Homing_Glyco           = c("GALNT1", "GALNT2", "GALNT3", "C1GALT1",
                                "C1GALT1C1", "GCNT1", "MGAT5"),
    Homing_Barrier         = c("CDH1", "CDH5", "OCLN", "TJP1",
                                "CLDN1", "CLDN3", "EPCAM"),
    Switch_Immune          = c("TLR2", "TLR4", "TNF", "IL6", "IL1B",
                                "NFKB1", "NFKB2", "IL18"),
    Switch_Hypoxia         = c("HIF1A", "VEGFA", "LDHA", "PGK1", "ENO1"),
    Invasion_MMP           = c("MMP1", "MMP2", "MMP9", "MMP14"),
    Invasion_TIMP          = c("TIMP1", "TIMP2", "TIMP3"),
    Defense                = c("DEFB1", "DEFB4A", "CAMP", "C3", "C5AR1", "LCN2"),
    Vascular_Angiogenic    = c("FLT1", "PGF", "KDR", "ENG", "VEGFA"),
    Vascular_EndoDysfx     = c("NOS3", "EDN1", "VCAM1", "ICAM1", "SELE")
  )
}

# ╔══════════════════════════════════════════════════════════════════════╗
# ║  NEW: IDO1_Tolerogenic_Shield — 12 genes from Part 3 PDF           ║
# ║                                                                    ║
# ║  The IDO1 tolerogenic shield maintains maternal immune tolerance    ║
# ║  at the decidua. Fn (via MegL → H₂S) collapses this shield,       ║
# ║  leading to loss of decidual immune privilege.                     ║
# ║                                                                    ║
# ║  IMPORTANT: IDO1 shield is PROTECTIVE — higher = more tolerant.    ║
# ║  In MISI composition, it should be INVERTED (like Barrier/Defense) ║
# ║  so that shield collapse → higher MISI vulnerability.              ║
# ╚══════════════════════════════════════════════════════════════════════╝

IDO1_TOLEROGENIC_SHIELD <- c(
  "IDO1",      # Indoleamine 2,3-dioxygenase 1 — master tolerogenic enzyme
  "TGFB1",     # TGF-β1 — immunosuppressive cytokine
  "IL10",      # Interleukin-10 — anti-inflammatory
  "CD274",     # PD-L1 (B7-H1) — T cell exhaustion ligand
  "PDCD1LG2",  # PD-L2 (B7-DC) — secondary PD-1 ligand
  "HAVCR2",    # TIM-3 — checkpoint receptor on decidual immune cells
  "LGALS9",    # Galectin-9 — TIM-3 ligand, immunoregulatory
  "CD80",      # B7-1 — co-stimulatory/co-inhibitory (CTLA-4 ligand)
  "CD86",      # B7-2 — co-stimulatory (CD28/CTLA-4 ligand)
  "ENTPD1",    # CD39 — ectonucleotidase (ATP → AMP)
  "NT5E",      # CD73 — ecto-5'-nucleotidase (AMP → adenosine)
  "FOXP3"      # Forkhead box P3 — Treg master TF
)

# Add to module list
all_modules <- c(base_modules, list(
  IDO1_Tolerogenic_Shield = IDO1_TOLEROGENIC_SHIELD
))

log_msg("  Total modules: ", length(all_modules),
        " (11 canonical + IDO1_Tolerogenic_Shield)")
log_msg("  IDO1 shield genes: ", paste(IDO1_TOLEROGENIC_SHIELD, collapse = ", "))


# =============================================================================
# 2. LOAD SEURAT OBJECT
# =============================================================================

log_msg("Loading SCP2601 Seurat object...")

obj_candidates <- c(
  file.path(DIR_OBJS, "SCP2601_combined_spatial.rds"),
  "output/objects/SCP2601_combined_spatial.rds",
  "SCP2601_combined_spatial.rds"
)

# Allow command-line override
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) obj_candidates <- c(args[1], obj_candidates)

seu <- NULL
for (cand in obj_candidates) {
  if (file.exists(cand) && file.size(cand) > 1000) {
    seu <- tryCatch(readRDS(cand), error = function(e) NULL)
    if (!is.null(seu) && inherits(seu, "Seurat")) {
      log_msg("  Loaded: ", cand)
      break
    }
  }
}

if (is.null(seu)) {
  stop("No SCP2601 Seurat object found. Run 05B_load_and_format_SCP2601.R first.\n",
       "  Searched: ", paste(obj_candidates, collapse = "\n  "))
}

log_msg("  Cells: ", ncol(seu), " | Genes: ", nrow(seu))
log_msg("  Reductions: ", paste(names(seu@reductions), collapse = ", "))
DefaultAssay(seu) <- "RNA"

# --- Detect spatial coordinates ---
has_spatial <- "spatial" %in% names(seu@reductions)
if (has_spatial) {
  xy <- Embeddings(seu, "spatial")[, 1:2]
  log_msg("  ✓ Spatial coords found (", nrow(xy), " cells)")
  log_msg("    X range: [", round(min(xy[, 1], na.rm = TRUE)), ", ",
          round(max(xy[, 1], na.rm = TRUE)), "]")
  log_msg("    Y range: [", round(min(xy[, 2], na.rm = TRUE)), ", ",
          round(max(xy[, 2], na.rm = TRUE)), "]")
} else {
  # Fallback: check metadata
  coord_pairs <- list(c("x", "y"), c("X", "Y"), c("x_um", "y_um"),
                       c("spatial_x", "spatial_y"))
  for (cp in coord_pairs) {
    if (all(cp %in% colnames(seu@meta.data))) {
      xy <- as.matrix(seu@meta.data[, cp])
      has_spatial <- TRUE
      log_msg("  ✓ Spatial coords found in metadata: ", paste(cp, collapse = ", "))
      break
    }
  }
}
if (!has_spatial) {
  log_msg("  ⚠ No spatial coordinates found! Moran's I sweep will be skipped.")
}

# --- Detect metadata columns ---
ct_col <- NULL
for (cand in c("cell_type", "celltype_final_refined", "celltype_final_conservative",
               "Clusters.2", "Clusters", "celltype", "predicted.id")) {
  if (cand %in% colnames(seu@meta.data)) {
    ct_col <- cand
    break
  }
}
if (is.null(ct_col)) {
  seu$cell_type <- "unknown"
  ct_col <- "cell_type"
  log_msg("  WARNING: No cell type column found. Using 'unknown'.")
} else {
  log_msg("  Cell type column: ", ct_col, " (",
          length(unique(seu@meta.data[[ct_col]])), " types)")
}

week_col <- NULL
for (cand in c("week", "stage", "sample_name", "orig.ident")) {
  if (cand %in% colnames(seu@meta.data)) {
    week_col <- cand
    break
  }
}
log_msg("  Week/stage column: ", ifelse(is.null(week_col), "NONE", week_col))


# =============================================================================
# 3. GENE PRESENCE REPORT
# =============================================================================

log_msg("Checking gene presence across all modules...")

gene_presence <- do.call(rbind, lapply(names(all_modules), function(mod) {
  genes <- all_modules[[mod]]
  data.frame(
    module  = mod,
    gene    = genes,
    present = genes %in% rownames(seu),
    stringsAsFactors = FALSE
  )
}))

write.csv(gene_presence, file.path(DIR_TABLES, "SCP2601_misi_v2_gene_presence.csv"),
          row.names = FALSE)

# Summary
gp_summary <- gene_presence %>%
  group_by(module) %>%
  summarise(n = n(), present = sum(present), pct = round(100 * mean(present), 1),
            .groups = "drop")

log_msg("  Gene presence summary:")
for (i in seq_len(nrow(gp_summary))) {
  log_msg("    ", gp_summary$module[i], ": ", gp_summary$present[i], "/",
          gp_summary$n[i], " (", gp_summary$pct[i], "%)")
}


# =============================================================================
# 4. UCELL MODULE SCORING
# =============================================================================

log_msg("Computing UCell module scores for ", length(all_modules), " modules...")

# Install UCell if needed
if (!requireNamespace("UCell", quietly = TRUE)) {
  log_msg("  Installing UCell...")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install("UCell", update = FALSE, ask = FALSE)
}
library(UCell)

# Filter modules to genes present in the dataset
sigs <- lapply(all_modules, function(gs) intersect(gs, rownames(seu)))
sigs <- sigs[vapply(sigs, length, integer(1)) > 0L]
log_msg("  Modules with ≥1 gene: ", length(sigs), " / ", length(all_modules))

# Ensure data layer is available (Seurat v5 compatibility)
tryCatch({
  seu <- JoinLayers(seu)
}, error = function(e) {
  log_msg("  JoinLayers not needed or not available.")
})

# Score with UCell
seu_scored <- UCell::AddModuleScore_UCell(
  seu,
  features = sigs,
  name = "",
  assay = "RNA",
  BPPARAM = BiocParallel::SerialParam()
)

# Extract raw scores into a data.frame
fetch_ucell_col <- function(nm, md) {
  candidates <- c(nm, paste0(nm, "_UCell"), paste0("UCell_", nm))
  hit <- candidates[candidates %in% colnames(md)]
  if (length(hit) == 0) return(rep(NA_real_, nrow(md)))
  md[[hit[1]]]
}

raw_scores <- as.data.frame(
  lapply(names(sigs), function(nm) fetch_ucell_col(nm, seu_scored@meta.data))
)
colnames(raw_scores) <- names(sigs)
rownames(raw_scores) <- colnames(seu_scored)

log_msg("  UCell scoring complete. Modules scored: ",
        paste(colnames(raw_scores), collapse = ", "))


# =============================================================================
# 5. Z-SCORE NORMALIZATION
# =============================================================================

log_msg("Z-score normalizing module scores...")

zscore_cols <- function(df) {
  as.data.frame(lapply(df, function(x) {
    if (all(is.na(x))) return(x)
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(x)))
    as.numeric(scale(x))
  }))
}

z_scores <- zscore_cols(raw_scores)
rownames(z_scores) <- rownames(raw_scores)


# =============================================================================
# 6. COMPOSE MISI SUBINDICES
# =============================================================================

log_msg("Composing MISI v2 subindices...")

# Invert protective modules for vulnerability orientation
z_inv <- z_scores

# Barrier: high = intact = protective → invert
if ("Homing_Barrier" %in% colnames(z_inv))
  z_inv$Homing_Barrier <- -z_inv$Homing_Barrier

# Defense: high = strong antimicrobial = protective → invert
if ("Defense" %in% colnames(z_inv))
  z_inv$Defense <- -z_inv$Defense

# Nutrient Sequestration: high = recaptures EA = protective → invert
if ("Nutrient_Sequestration" %in% colnames(z_inv))
  z_inv$Nutrient_Sequestration <- -z_inv$Nutrient_Sequestration

# IDO1_Tolerogenic_Shield: high = strong tolerance = protective → invert
# (Shield collapse = vulnerability)
if ("IDO1_Tolerogenic_Shield" %in% colnames(z_inv))
  z_inv$IDO1_Tolerogenic_Shield_inv <- -z_inv$IDO1_Tolerogenic_Shield

# Invasion composite: MMP minus TIMP
z_inv$Invasion_MMP_composite <- ifelse(
  !is.na(z_inv$Invasion_MMP) & !is.na(z_inv$Invasion_TIMP),
  z_inv$Invasion_MMP - z_inv$Invasion_TIMP,
  NA_real_
)

# Safe row-mean helper
safe_rowmean <- function(cols, df) {
  present <- cols[cols %in% colnames(df)]
  if (length(present) == 0) return(rep(NA_real_, nrow(df)))
  rowMeans(df[, present, drop = FALSE], na.rm = TRUE)
}

# 4 canonical subindices
subindices <- data.frame(
  MISI_Homing            = safe_rowmean(c("Homing_Glyco", "Homing_Barrier"), z_inv),
  MISI_Nutrient          = safe_rowmean(c("Nutrient_Liberation", "Nutrient_Sequestration"), z_inv),
  MISI_SwitchPressure    = safe_rowmean(c("Switch_Immune", "Switch_Hypoxia"), z_inv),
  MISI_VascularFragility = safe_rowmean(c("Vascular_Angiogenic", "Vascular_EndoDysfx"), z_inv),
  Invasion_MMP_composite = z_inv$Invasion_MMP_composite,
  stringsAsFactors = FALSE
)

# IDO1 Shield score (inverted = collapse vulnerability)
if ("IDO1_Tolerogenic_Shield_inv" %in% colnames(z_inv)) {
  subindices$IDO1_Shield_Collapse <- z_inv$IDO1_Tolerogenic_Shield_inv
  # Also keep the raw (non-inverted) score for interpretability
  subindices$IDO1_Shield_Intact <- z_scores$IDO1_Tolerogenic_Shield
}

# Overall MISI (5 canonical components, equal weights)
components <- cbind(
  subindices$MISI_Homing,
  subindices$MISI_Nutrient,
  subindices$MISI_SwitchPressure,
  subindices$MISI_VascularFragility,
  subindices$Invasion_MMP_composite
)
subindices$MISI <- rowMeans(components, na.rm = TRUE)

# Extended MISI incorporating IDO1 shield collapse
if ("IDO1_Shield_Collapse" %in% colnames(subindices)) {
  components_ext <- cbind(components, subindices$IDO1_Shield_Collapse)
  subindices$MISI_extended <- rowMeans(components_ext, na.rm = TRUE)
}

rownames(subindices) <- colnames(seu)

log_msg("  Subindices computed:")
for (si in colnames(subindices)) {
  vals <- subindices[[si]]
  log_msg("    ", si, ": mean=", round(mean(vals, na.rm = TRUE), 4),
          " sd=", round(sd(vals, na.rm = TRUE), 4),
          " range=[", round(min(vals, na.rm = TRUE), 3), ", ",
          round(max(vals, na.rm = TRUE), 3), "]")
}


# =============================================================================
# 7. ATTACH TO SEURAT OBJECT
# =============================================================================

log_msg("Attaching MISI v2 scores to Seurat metadata...")

# Module z-scores
for (mod in colnames(z_scores)) {
  seu@meta.data[[paste0("MISIv2_", mod, "_z")]] <- z_scores[[mod]]
}

# Subindices
for (si in colnames(subindices)) {
  seu@meta.data[[paste0("MISIv2_", si)]] <- subindices[[si]]
}

# Raw IDO1 shield score for direct visualization
if ("IDO1_Tolerogenic_Shield" %in% colnames(raw_scores)) {
  seu$IDO1_Shield_raw <- raw_scores$IDO1_Tolerogenic_Shield
}

log_msg("  Added ", ncol(z_scores) + ncol(subindices), " MISIv2_ columns")


# =============================================================================
# 8. PER-CELL SCORES TABLE
# =============================================================================

log_msg("Building per-cell scores table...")

cell_table <- tibble(
  cell       = colnames(seu),
  cell_type  = as.character(seu@meta.data[[ct_col]]),
  week       = if (!is.null(week_col)) as.character(seu@meta.data[[week_col]]) else NA_character_,
  has_coords = has_spatial
)

for (mod in colnames(raw_scores))
  cell_table[[paste0(mod, "_raw")]] <- raw_scores[[mod]]
for (mod in colnames(z_scores))
  cell_table[[paste0(mod, "_z")]] <- z_scores[[mod]]
for (si in colnames(subindices))
  cell_table[[si]] <- subindices[[si]]

write.csv(cell_table, file.path(DIR_TABLES, "SCP2601_misi_v2_per_cell_scores.csv"),
          row.names = FALSE)
log_msg("  Saved: SCP2601_misi_v2_per_cell_scores.csv (",
        nrow(cell_table), " cells x ", ncol(cell_table), " columns)")


# =============================================================================
# 9. SUMMARY BY CELLTYPE × WEEK
# =============================================================================

log_msg("Computing summary by celltype × week...")

score_cols_summary <- intersect(
  c("MISI", "MISI_extended", "MISI_Homing", "MISI_Nutrient",
    "MISI_SwitchPressure", "MISI_VascularFragility",
    "Invasion_MMP_composite", "IDO1_Shield_Intact", "IDO1_Shield_Collapse"),
  colnames(subindices)
)

summary_input <- data.frame(
  cell_type = as.character(seu@meta.data[[ct_col]]),
  week      = if (!is.null(week_col)) as.character(seu@meta.data[[week_col]]) else "all",
  stringsAsFactors = FALSE
)
for (sc in score_cols_summary) summary_input[[sc]] <- subindices[[sc]]

summary_tbl <- summary_input %>%
  group_by(cell_type, week) %>%
  summarise(
    n_cells = n(),
    across(all_of(score_cols_summary),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd = ~sd(.x, na.rm = TRUE),
                median = ~median(.x, na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  )

write.csv(summary_tbl, file.path(DIR_TABLES, "SCP2601_misi_v2_summary_by_celltype_week.csv"),
          row.names = FALSE)
log_msg("  Saved: SCP2601_misi_v2_summary_by_celltype_week.csv (",
        nrow(summary_tbl), " strata)")


# =============================================================================
# 10. SPATIAL MORAN'S I RADIUS SWEEP
# =============================================================================

if (RUN_SPATIAL && has_spatial) {
  log_msg("\n=== Spatial Moran's I Radius Sweep ===")

  # Install spdep if needed
  if (!requireNamespace("spdep", quietly = TRUE)) {
    log_msg("  Installing spdep...")
    install.packages("spdep", repos = "https://cloud.r-project.org")
  }
  library(spdep)

  #' Moran's I at a single radius for a single score
  moranI_at_radius <- function(xy_mat, values, radius) {
    valid <- !is.na(values) & complete.cases(xy_mat)
    if (sum(valid) < 30) return(list(I = NA, p = NA, n = sum(valid), mean_nb = NA))

    xy_v  <- xy_mat[valid, , drop = FALSE]
    val_v <- values[valid]

    nb <- tryCatch(spdep::dnearneigh(xy_v, 0, radius, longlat = FALSE),
                   error = function(e) NULL)
    if (is.null(nb)) return(list(I = NA, p = NA, n = length(val_v), mean_nb = NA))

    card_nb <- spdep::card(nb)
    if (all(card_nb == 0)) return(list(I = NA, p = NA, n = length(val_v), mean_nb = 0))

    lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
    mt <- tryCatch(spdep::moran.test(val_v, lw, zero.policy = TRUE),
                   error = function(e) NULL)

    if (is.null(mt)) return(list(I = NA, p = NA, n = length(val_v), mean_nb = mean(card_nb)))

    list(
      I = as.numeric(mt$estimate[["Moran I statistic"]]),
      p = mt$p.value,
      n = length(val_v),
      mean_nb = mean(card_nb)
    )
  }

  # Scores to sweep
  spatial_score_cols <- intersect(
    c("MISI", "MISI_extended", "MISI_Homing", "MISI_Nutrient",
      "MISI_SwitchPressure", "MISI_VascularFragility",
      "Invasion_MMP_composite", "IDO1_Shield_Intact", "IDO1_Shield_Collapse"),
    colnames(subindices)
  )

  # Get XY matrix
  if ("spatial" %in% names(seu@reductions)) {
    xy_mat <- Embeddings(seu, "spatial")[, 1:2]
  } else {
    xy_mat <- xy
  }

  # Run sweep per sample (if multi-sample) or global
  if (!is.null(week_col)) {
    sample_groups <- split(seq_len(ncol(seu)), seu@meta.data[[week_col]])
  } else {
    sample_groups <- list(all = seq_len(ncol(seu)))
  }

  all_spatial_results <- list()

  for (grp_name in names(sample_groups)) {
    idx <- sample_groups[[grp_name]]
    xy_grp <- xy_mat[idx, , drop = FALSE]

    log_msg("  Sample/week: ", grp_name, " (", length(idx), " cells)")

    for (sc in spatial_score_cols) {
      vals <- subindices[[sc]][idx]

      for (r in SPATIAL_RADII) {
        res <- moranI_at_radius(xy_grp, vals, r)
        all_spatial_results[[length(all_spatial_results) + 1]] <- data.frame(
          sample    = grp_name,
          score     = sc,
          radius    = r,
          moran_I   = res$I,
          p_value   = res$p,
          n_cells   = res$n,
          mean_nb   = res$mean_nb,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  spatial_df <- do.call(rbind, all_spatial_results)
  write.csv(spatial_df, file.path(DIR_TABLES, "SCP2601_misi_v2_spatial_radius_sweep.csv"),
            row.names = FALSE)
  log_msg("  Saved: SCP2601_misi_v2_spatial_radius_sweep.csv (",
          nrow(spatial_df), " rows)")

  # Report peaks
  for (sc in spatial_score_cols) {
    sc_data <- spatial_df %>% filter(score == sc, !is.na(moran_I))
    if (nrow(sc_data) > 0) {
      peak <- sc_data[which.max(sc_data$moran_I), ]
      log_msg("    ", sc, ": peak I=", round(peak$moran_I, 4),
              " at radius=", peak$radius,
              " (p=", format(peak$p_value, digits = 3), ")",
              " sample=", peak$sample)
    }
  }

  # --- Separate IDO1 shield sweep for emphasis ---
  ido1_spatial <- spatial_df %>% filter(grepl("IDO1", score))
  if (nrow(ido1_spatial) > 0) {
    write.csv(ido1_spatial,
              file.path(DIR_TABLES, "SCP2601_ido1_shield_spatial_sweep.csv"),
              row.names = FALSE)
    log_msg("  Saved: SCP2601_ido1_shield_spatial_sweep.csv")
  }

} else {
  log_msg("  Skipping Moran's I sweep (RUN_SPATIAL=", RUN_SPATIAL,
          ", has_spatial=", has_spatial, ")")
  spatial_df <- NULL
}


# =============================================================================
# 11. SPATIAL FEATURE MAPS
# =============================================================================

if (has_spatial) {
  log_msg("\n=== Generating Spatial Feature Maps ===")

  spatial_theme <- theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text  = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8, color = "grey40"),
      panel.border = element_rect(linetype = "solid", color = "black",
                                   fill = NA, linewidth = 0.8)
    )

  # Helper: spatial feature map for one score
  plot_spatial_score <- function(seu_obj, score_col, title, subtitle = "") {
    df <- data.frame(
      x = seu_obj$x,
      y = seu_obj$y,
      score = seu_obj@meta.data[[score_col]]
    )
    df <- df[!is.na(df$x) & !is.na(df$y), ]

    ggplot(df, aes(x = x, y = y, color = score)) +
      geom_point(size = 0.5, alpha = 0.8) +
      scale_color_gradient2(low = "navy", mid = "white", high = "firebrick",
                            midpoint = 0, name = "Score") +
      coord_fixed() +
      ggtitle(title, subtitle = subtitle) +
      spatial_theme
  }

  # --- MISI subindex spatial maps ---
  misi_map_features <- c("MISIv2_MISI", "MISIv2_MISI_Homing",
                           "MISIv2_MISI_Nutrient", "MISIv2_MISI_SwitchPressure",
                           "MISIv2_MISI_VascularFragility",
                           "MISIv2_Invasion_MMP_composite")

  if ("MISIv2_IDO1_Shield_Intact" %in% colnames(seu@meta.data)) {
    misi_map_features <- c(misi_map_features,
                            "MISIv2_IDO1_Shield_Intact",
                            "MISIv2_IDO1_Shield_Collapse")
  }

  # Filter to existing columns
  misi_map_features <- intersect(misi_map_features, colnames(seu@meta.data))

  if (length(misi_map_features) > 0 && "x" %in% colnames(seu@meta.data)) {
    plot_list <- lapply(misi_map_features, function(feat) {
      short_name <- gsub("MISIv2_", "", feat)
      plot_spatial_score(seu, feat, short_name,
                         paste0("SCP2601 | n=", ncol(seu)))
    })

    p_combined <- wrap_plots(plot_list, ncol = 3) +
      plot_annotation(
        title = "MISI v2 Subindices — Spatial Feature Maps (SCP2601)",
        subtitle = "Blue = low/protective | Red = high/vulnerable | Physical tissue coordinates",
        theme = theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10, color = "grey40")
        )
      )

    ggsave(file.path(DIR_FIGURES, "SCP2601_misi_spatial_featuremaps.pdf"),
           p_combined, width = 18, height = ceiling(length(plot_list) / 3) * 5,
           limitsize = FALSE)
    ggsave(file.path(DIR_FIGURES, "SCP2601_misi_spatial_featuremaps.png"),
           p_combined, width = 18, height = ceiling(length(plot_list) / 3) * 5,
           dpi = 300, limitsize = FALSE)
    log_msg("  Saved: SCP2601_misi_spatial_featuremaps.pdf/png")
  }

  # --- IDO1 Shield dedicated spatial map ---
  if ("MISIv2_IDO1_Shield_Intact" %in% colnames(seu@meta.data) &&
      "x" %in% colnames(seu@meta.data)) {
    df_ido1 <- data.frame(
      x = seu$x, y = seu$y,
      IDO1_intact  = seu$MISIv2_IDO1_Shield_Intact,
      IDO1_collapse = seu$MISIv2_IDO1_Shield_Collapse
    )
    df_ido1 <- df_ido1[!is.na(df_ido1$x) & !is.na(df_ido1$y), ]

    p_ido1_intact <- ggplot(df_ido1, aes(x = x, y = y, color = IDO1_intact)) +
      geom_point(size = 0.5, alpha = 0.8) +
      scale_color_viridis_c(option = "magma", direction = -1, name = "Shield\nIntact") +
      coord_fixed() + ggtitle("IDO1 Tolerogenic Shield (Intact Score)") +
      spatial_theme

    p_ido1_collapse <- ggplot(df_ido1, aes(x = x, y = y, color = IDO1_collapse)) +
      geom_point(size = 0.5, alpha = 0.8) +
      scale_color_gradient2(low = "navy", mid = "grey90", high = "darkred",
                            midpoint = 0, name = "Shield\nCollapse") +
      coord_fixed() + ggtitle("IDO1 Shield Collapse (Vulnerability)") +
      spatial_theme

    p_ido1_combined <- p_ido1_intact | p_ido1_collapse
    p_ido1_combined <- p_ido1_combined +
      plot_annotation(
        title = "IDO1 Tolerogenic Shield — Spatial Distribution (SCP2601)",
        subtitle = "Left: intact (high = strong tolerance) | Right: collapse (high = vulnerable)",
        theme = theme(
          plot.title = element_text(size = 13, face = "bold"),
          plot.subtitle = element_text(size = 9, color = "grey40")
        )
      )

    ggsave(file.path(DIR_FIGURES, "SCP2601_ido1_shield_spatial_map.pdf"),
           p_ido1_combined, width = 14, height = 6)
    ggsave(file.path(DIR_FIGURES, "SCP2601_ido1_shield_spatial_map.png"),
           p_ido1_combined, width = 14, height = 6, dpi = 300)
    log_msg("  Saved: SCP2601_ido1_shield_spatial_map.pdf/png")
  }
}


# =============================================================================
# 12. MORAN'S I CURVES (VISUALIZATION)
# =============================================================================

if (!is.null(spatial_df) && nrow(spatial_df) > 0) {
  log_msg("Generating Moran's I curve plots...")

  p_moran <- ggplot(spatial_df %>% filter(!is.na(moran_I)),
                    aes(x = radius, y = moran_I, color = score)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    facet_wrap(~sample, scales = "free_y") +
    scale_color_brewer(palette = "Set1") +
    labs(title = "Moran's I Spatial Autocorrelation — MISI v2 Subindices (SCP2601)",
         subtitle = "Peak Moran's I indicates optimal spatial neighborhood scale",
         x = "Radius (µm)", y = "Moran's I", color = "Score") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          strip.text = element_text(face = "bold"))

  ggsave(file.path(DIR_FIGURES, "SCP2601_misi_moranI_curves.pdf"),
         p_moran, width = 14, height = 8)
  ggsave(file.path(DIR_FIGURES, "SCP2601_misi_moranI_curves.png"),
         p_moran, width = 14, height = 8, dpi = 300)
  log_msg("  Saved: SCP2601_misi_moranI_curves.pdf/png")
}


# =============================================================================
# 13. VIOLIN PLOTS BY WEEK
# =============================================================================

if (!is.null(week_col)) {
  log_msg("Generating MISI violin plots by week...")

  violin_features <- intersect(
    c("MISIv2_MISI", "MISIv2_IDO1_Shield_Intact",
      "MISIv2_MISI_Homing", "MISIv2_MISI_VascularFragility"),
    colnames(seu@meta.data)
  )

  if (length(violin_features) > 0) {
    vln_list <- lapply(violin_features, function(feat) {
      VlnPlot(seu, features = feat, group.by = week_col,
              pt.size = 0, cols = brewer.pal(4, "Set2")) +
        geom_hline(yintercept = 0, lty = 2, color = "grey50") +
        ggtitle(gsub("MISIv2_", "", feat)) +
        theme(legend.position = "none",
              plot.title = element_text(face = "bold", size = 10))
    })

    p_vln <- wrap_plots(vln_list, ncol = 2) +
      plot_annotation(
        title = "MISI v2 Scores by Gestational Week (SCP2601)",
        subtitle = "Dashed line = neutral (0) | SCP2601 spatial data"
      )

    ggsave(file.path(DIR_FIGURES, "SCP2601_misi_violins_by_week.pdf"),
           p_vln, width = 12, height = 10)
    log_msg("  Saved: SCP2601_misi_violins_by_week.pdf")
  }
}


# =============================================================================
# 14. MISI vs IDO1 SHIELD SCATTER
# =============================================================================

if (all(c("MISIv2_MISI", "MISIv2_IDO1_Shield_Intact") %in% colnames(seu@meta.data))) {
  log_msg("Generating MISI vs IDO1 Shield scatter...")

  df_scatter <- data.frame(
    MISI = seu$MISIv2_MISI,
    IDO1 = seu$MISIv2_IDO1_Shield_Intact,
    cell_type = as.character(seu@meta.data[[ct_col]]),
    stringsAsFactors = FALSE
  )

  r_val <- cor(df_scatter$MISI, df_scatter$IDO1, use = "complete.obs",
               method = "spearman")

  p_scatter <- ggplot(df_scatter, aes(x = MISI, y = IDO1, color = cell_type)) +
    geom_point(size = 0.5, alpha = 0.4) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
    annotate("text", x = Inf, y = Inf, label = paste0("ρ = ", round(r_val, 3)),
             hjust = 1.1, vjust = 1.3, size = 5, fontface = "bold") +
    labs(title = "MISI Vulnerability vs IDO1 Tolerogenic Shield",
         subtitle = paste0("SCP2601 | Spearman ρ = ", round(r_val, 3),
                           " | Negative = shield collapses with vulnerability"),
         x = "MISI (Vulnerability Index)",
         y = "IDO1 Shield (Intact Score)",
         color = "Cell Type") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "right")

  ggsave(file.path(DIR_FIGURES, "SCP2601_misi_vs_ido1_scatter.pdf"),
         p_scatter, width = 10, height = 8)
  ggsave(file.path(DIR_FIGURES, "SCP2601_misi_vs_ido1_scatter.png"),
         p_scatter, width = 10, height = 8, dpi = 300)
  log_msg("  Saved: SCP2601_misi_vs_ido1_scatter.pdf/png")
  log_msg("  MISI-IDO1 Spearman ρ = ", round(r_val, 3))
}


# =============================================================================
# 15. SAVE SEURAT OBJECT
# =============================================================================

log_msg("\nSaving Seurat object with MISI v2 + IDO1 Shield scores...")
out_path <- file.path(DIR_OBJS, "SCP2601_with_misi_v2.rds")
saveRDS(seu, out_path)
log_msg("  Saved: ", out_path)


# =============================================================================
# 16. FINAL REPORT
# =============================================================================

cat("\n╔══════════════════════════════════════════════════════════════════════╗\n")
cat(  "║  16 — MISI v2 + IDO1 Shield COMPLETE (SCP2601)                    ║\n")
cat(  "╠══════════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Cells: %-57d ║\n", ncol(seu)))
cat(sprintf("║  Genes: %-57d ║\n", nrow(seu)))
cat(sprintf("║  Modules scored: %-48d ║\n", ncol(raw_scores)))
cat(sprintf("║  Spatial sweep: %-49s ║\n",
            ifelse(RUN_SPATIAL && has_spatial, "COMPLETE ✓", "SKIPPED")))
cat(sprintf("║  IDO1 Shield: %-51s ║\n",
            ifelse("IDO1_Tolerogenic_Shield" %in% colnames(raw_scores),
                   "SCORED ✓", "NOT AVAILABLE")))
cat("╠══════════════════════════════════════════════════════════════════════╣\n")
cat("║  Outputs:                                                          ║\n")
cat("║    output/objects/SCP2601_with_misi_v2.rds                         ║\n")
cat("║    output/tables/SCP2601_misi_v2_gene_presence.csv                 ║\n")
cat("║    output/tables/SCP2601_misi_v2_per_cell_scores.csv               ║\n")
cat("║    output/tables/SCP2601_misi_v2_summary_by_celltype_week.csv      ║\n")
cat("║    output/tables/SCP2601_misi_v2_spatial_radius_sweep.csv          ║\n")
cat("║    output/tables/SCP2601_ido1_shield_spatial_sweep.csv             ║\n")
cat("║    output/figures/SCP2601_misi_spatial_featuremaps.pdf              ║\n")
cat("║    output/figures/SCP2601_ido1_shield_spatial_map.pdf               ║\n")
cat("║    output/figures/SCP2601_misi_moranI_curves.pdf                   ║\n")
cat("║    output/figures/SCP2601_misi_violins_by_week.pdf                 ║\n")
cat("║    output/figures/SCP2601_misi_vs_ido1_scatter.pdf                 ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

log_msg("16_compute_misi_SCP2601_spatial done.")