# =============================================================================
# scripts/21_rpac_spatial_routes.R
# rPAC v2 — Physical Spatial Routes on SCP2601 Tissue Coordinates
# =============================================================================
#
# PURPOSE:
#   Adapt the canonical 21_rpac_v2_corrected_routes.R to compute rPAC vectors
#   across PHYSICAL SPATIAL X/Y coordinates (instead of UMAP space).
#
# SEURAT V5 & PERFORMANCE PATCHES:
#   1. JoinLayers() - Required for merged STARmap samples.
#   2. layer="data" - Replaces defunct slot="data".
#   3. sparseMatrixStats - Used for high-speed Z-scoring on large matrices.
#
# AUTHOR:  Shan Kurkcu (PhD Thesis Pipeline)
# DATE:    2026-04-14
# VERSION: 1.2.0 (Full Restoration + v5 Patch)
# =============================================================================

cat("\n╔══════════════════════════════════════════════════════════════════════╗\n")
cat(  "║  21 — rPAC Spatial Routes: Physical Vectors on SCP2601 Tissue      ║\n")
cat(  "╚══════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# 0. SETUP
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(viridis)
  library(sparseMatrixStats) # Crucial for STARmap scale
  library(RANN)              # For fast k-NN
})

# Pipeline infrastructure
tryCatch({
  source("config/config.R")
  source("scripts/R/utils.R")
}, error = function(e) {
  cat("[21] Running in standalone mode.\n")
})

# Directories
DIR_OBJS    <- "output/objects"
DIR_TABLES  <- "output/tables"
DIR_FIGURES <- "output/figures"
DIR_LOGS    <- "output/logs"
for (d in c(DIR_OBJS, DIR_TABLES, DIR_FIGURES, DIR_LOGS)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

log_file <- file.path(DIR_LOGS, "21_rpac_spatial_routes.log")
log_msg <- function(...) {
  msg <- paste0(Sys.time(), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}
log_msg("21_rpac_spatial_routes starting.")

SEED <- 42
set.seed(SEED)

# =============================================================================
# 1. ROUTE DEFINITIONS (v2, HIF-2α corrected)
# =============================================================================

routes_v2_genes <- list(
  NFKB_SR_innate          = list(genes = c("TLR2","TLR4","MYD88","IRAK1","IRAK4","TRAF6","NFKB1","RELA")),
  NFKB_ER_cytokine        = list(genes = c("RELA","IL6","CXCL8","TNF","IL1B","CCL2","PTGS2")),
  HIF1A_ER_metabolic      = list(genes = c("HIF1A","VEGFA","LDHA","PGK1","ENO1","SLC2A1","PDK1")),
  VEGF_SR_angiogenic      = list(genes = c("KDR","FLT1","VEGFA")),
  EA_nutrient_axis        = list(genes = c("ETNK1","PCYT2","SELENOI","EPT1","ETNK2")),
  NK_ER_cytotoxic         = list(genes = c("TBX21","PRF1","GZMB","IFNG","GNLY","NKG7")),
  EPAS1_ER_antiangiogenic = list(genes = c("EPAS1","ARNT","FLT1","PGF")),
  NFKB_EPAS1_bridge       = list(genes = c("RELA","EPAS1","FLT1","PGF")),
  EPAS1_ER_notch          = list(genes = c("EPAS1","DLL4","ANGPT2","NOTCH1"))
)
route_ids <- names(routes_v2_genes)

# =============================================================================
# 2. LOAD & PREPARE SEURAT OBJECT
# =============================================================================

log_msg("Loading Seurat object...")
obj_candidates <- c(
  file.path(DIR_OBJS, "SCP2601_with_misi_v2.rds"),
  file.path(DIR_OBJS, "SCP2601_combined_spatial.rds")
)

seu <- NULL
for (cand in obj_candidates) {
  if (file.exists(cand)) {
    seu <- readRDS(cand)
    log_msg("  Loaded: ", cand)
    break
  }
}

if (is.null(seu)) stop("Seurat object not found. Run scripts 05B and 16 first.")

# SEURAT V5 FIX: Join Layers for merged objects
log_msg("  Ensuring v5 layers are joined...")
seu <- JoinLayers(seu)
DefaultAssay(seu) <- "RNA"

# --- Extract spatial coordinates ---
# --- Extract & Filter Spatial Coordinates ---
log_msg("  Extracting and filtering spatial coordinates...")
xy_full <- NULL
if ("spatial" %in% names(seu@reductions)) {
  xy_full <- Embeddings(seu, "spatial")[, 1:2]
} else {
  for (pair in list(c("x","y"), c("x_um","y_um"), c("spatial_x", "spatial_y"))) {
    if (all(pair %in% colnames(seu@meta.data))) {
      xy_full <- as.matrix(seu@meta.data[, pair])
      break
    }
  }
}

if (is.null(xy_full)) stop("No spatial coordinates found in any expected slot!")

# IDENTIFY VALID CELLS (Those with non-NA coordinates)
valid_spatial_cells <- !is.na(xy_full[,1]) & !is.na(xy_full[,2])
n_missing <- sum(!valid_spatial_cells)

if (n_missing > 0) {
  log_msg("  WARNING: Found ", n_missing, " cells with missing coordinates (likely Sample W9).")
  log_msg("  Subsetting object to include only the ", sum(valid_spatial_cells), " mappable cells.")
  seu <- seu[, valid_spatial_cells]
  xy  <- xy_full[valid_spatial_cells, ]
} else {
  xy <- xy_full
}

# --- Detect MISI/IDO1 columns (Improved Logic) ---
misi_hits <- grep("MISI|misi", colnames(seu@meta.data), value=TRUE)
misi_col  <- if(length(misi_hits) > 0) misi_hits[1] else NULL

ido1_hits <- grep("IDO1|ido1", colnames(seu@meta.data), value=TRUE)
ido1_col  <- if(length(ido1_hits) > 0) ido1_hits[1] else NULL

ct_col    <- "cell_type"

log_msg("  Final Analysis Cells: ", ncol(seu))
log_msg("  MISI col detection: ", ifelse(is.null(misi_col), "MISSING", misi_col))
log_msg("  IDO1 col detection: ", ifelse(is.null(ido1_col), "MISSING", ido1_col))

# =============================================================================
# 3. SPATIAL rPAC SCORING (OPTIMIZED)
# =============================================================================

log_msg("\n=== Computing Spatial rPAC Scores ===")
expr_mat <- GetAssayData(seu, assay = "RNA", layer = "data")

score_route_mean_z <- function(genes, mat) {
  present_genes <- intersect(genes, rownames(mat))
  if (length(present_genes) == 0) return(rep(0, ncol(mat)))
  sub_mat <- mat[present_genes, , drop = FALSE]
  means <- rowMeans(sub_mat)
  sds   <- rowSds(sub_mat)
  sds[sds == 0] <- 1
  z_mat <- (sub_mat - means) / sds
  return(colMeans(z_mat))
}

scores <- data.frame(row.names = colnames(seu))
for (rid in route_ids) {
  scores[[rid]] <- score_route_mean_z(routes_v2_genes[[rid]]$genes, expr_mat)
}

# --- Spatial Gaussian Smoothing ---
log_msg("  Applying spatial Gaussian smoothing (BW = 100 µm)...")
SMOOTH_BANDWIDTH <- 100 
K_SMOOTH <- min(50, ncol(seu) - 1)
nn <- nn2(xy, xy, k = K_SMOOTH + 1)
weights <- exp(-nn$nn.dists[, -1]^2 / (2 * SMOOTH_BANDWIDTH^2))
weights <- weights / rowSums(weights)

scores_smooth <- as.data.frame(lapply(scores, function(s) {
  vapply(1:nrow(nn$nn.idx), function(i) sum(s[nn$nn.idx[i, -1]] * weights[i, ]), numeric(1))
}))

for (rid in route_ids) {
  seu[[paste0("rPAC_", rid)]] <- scores[[rid]]
  seu[[paste0("rPAC_smooth_", rid)]] <- scores_smooth[[rid]]
}

# =============================================================================
# 4. EVT → DECIDUAL VECTORS & INTERFACE ANALYSIS
# =============================================================================

if (!is.na(misi_col) && !is.na(ido1_col)) {
  log_msg("\n=== Computing EVT → Decidual Spatial Vectors ===")
  
  evt_cells <- colnames(seu)[grepl("EVT", seu[[ct_col]], ignore.case=T)]
  dec_cells <- colnames(seu)[grepl("FIB|Stromal|Hofbauer|Macrophage", seu[[ct_col]], ignore.case=T)]
  
  # MISI-High EVT (Q75) vs IDO1-High Decidua (Q75)
  evt_h <- intersect(evt_cells, colnames(seu)[seu[[misi_col]] >= quantile(seu[[misi_col]], 0.75, na.rm=T)])
  dec_h <- intersect(dec_cells, colnames(seu)[seu[[ido1_col]] >= quantile(seu[[ido1_col]], 0.75, na.rm=T)])
  
  if (length(evt_h) > 0 && length(dec_h) > 0) {
    nn_vec <- nn2(xy[dec_h,], xy[evt_h,], k=1)
    vector_table <- data.frame(
      evt_cell = evt_h, dec_cell = dec_h[nn_vec$nn.idx],
      evt_x = xy[evt_h, 1], evt_y = xy[evt_h, 2],
      dec_x = xy[dec_h[nn_vec$nn.idx], 1], dec_y = xy[dec_h[nn_vec$nn.idx], 2],
      distance_um = nn_vec$nn.dists
    )
    write.csv(vector_table, file.path(DIR_TABLES, "SCP2601_rpac_evt_decidual_vectors.csv"))
  }
}

# =============================================================================
# 8. INNOVATIVE: SPATIAL GRADIENT FIELD
# =============================================================================

log_msg("\n=== Computing Spatial Gradient Field ===")
if (!is.na(misi_col)) {
  GRID_N <- 50
  df_field <- data.frame(x = xy[, 1], y = xy[, 2], score = seu[[misi_col]])
  df_field <- df_field[complete.cases(df_field), ]
  
  x_range <- range(df_field$x); y_range <- range(df_field$y)
  x_grid <- seq(x_range[1], x_range[2], length.out = GRID_N)
  y_grid <- seq(y_range[1], y_range[2], length.out = GRID_N)
  
  bw <- max(diff(x_range), diff(y_range)) / 15
  score_grid <- matrix(NA, GRID_N, GRID_N)
  
  for (i in 1:GRID_N) {
    for (j in 1:GRID_N) {
      d <- sqrt((df_field$x - x_grid[i])^2 + (df_field$y - y_grid[j])^2)
      w <- exp(-d^2 / (2 * bw^2))
      if (sum(w) > 1e-10) score_grid[i, j] <- sum(w * df_field$score) / sum(w)
    }
  }
  
  # Quiver plot math
  grid_df <- expand.grid(x = x_grid, y = y_grid)
  grid_df$score <- as.vector(score_grid)
  p_gradient <- ggplot(grid_df, aes(x=x, y=y, fill=score)) + geom_tile() +
    scale_fill_gradient2(low="navy", high="firebrick", midpoint=0) +
    coord_fixed() + theme_minimal() + ggtitle("MISI Vulnerability Gradient Field")
  
  ggsave(file.path(DIR_FIGURES, "SCP2601_rpac_gradient_field.png"), p_gradient, width=10, height=9)
}

# =============================================================================
# 9. INNOVATIVE: IDO1-MISI INTERFACE TENSION MAP
# =============================================================================

if (!is.na(misi_col) && !is.na(ido1_col)) {
  log_msg("Generating Interface Tension Map...")
  seu$interface_zone <- case_when(
    seu[[misi_col]] >= median(seu[[misi_col]]) & seu[[ido1_col]] >= median(seu[[ido1_col]]) ~ "Active Interface",
    seu[[misi_col]] >= median(seu[[misi_col]]) & seu[[ido1_col]] <  median(seu[[ido1_col]]) ~ "Shield Collapsed",
    seu[[misi_col]] <  median(seu[[misi_col]]) & seu[[ido1_col]] >= median(seu[[ido1_col]]) ~ "Protected Zone",
    TRUE ~ "Quiescent"
  )
  
  p_tension <- ggplot(as.data.frame(xy), aes(x=V1, y=V2, color=seu$interface_zone)) +
    geom_point(size=0.4, alpha=0.7) + scale_color_manual(values=c("Active Interface"="#FFB300", "Shield Collapsed"="#D32F2F", "Protected Zone"="#388E3C", "Quiescent"="#1565C0")) +
    coord_fixed() + theme_void() + ggtitle("IDO1-MISI Interface Tension Map")
  
  ggsave(file.path(DIR_FIGURES, "SCP2601_rpac_ido1_interface_tension.png"), p_tension, width=11, height=9)
}

# =============================================================================
# 10. INNOVATIVE: NICHE BOUNDARY ANALYSIS
# =============================================================================

log_msg("Computing Niche Boundary Analysis...")
K_BOUNDARY <- 20
nn_b <- nn2(xy, xy, k = K_BOUNDARY + 1)
ct_vec <- as.character(seu[[ct_col]])
boundary_frac <- vapply(1:nrow(nn_b$nn.idx), function(i) mean(ct_vec[nn_b$nn.idx[i, -1]] != ct_vec[i]), numeric(1))
seu$is_boundary <- boundary_frac >= 0.4

p_boundary <- ggplot(as.data.frame(xy), aes(x=V1, y=V2, color=seu$is_boundary)) +
  geom_point(size=0.3) + scale_color_manual(values=c("FALSE"="grey90", "TRUE"="#E91E63")) +
  coord_fixed() + theme_void() + ggtitle("Tissue Niche Boundaries")

ggsave(file.path(DIR_FIGURES, "SCP2601_rpac_niche_boundary.png"), p_boundary, width=11, height=9)

# =============================================================================
# 11. INNOVATIVE: POLAR ROUTE SUMMARY (RADAR)
# =============================================================================

polar_data <- scores_smooth %>% mutate(cell_type = seu[[ct_col]]) %>%
  pivot_longer(-cell_type, names_to="route", values_to="score") %>%
  group_by(cell_type, route) %>% summarise(mean_score = mean(score, na.rm=T), .groups="drop")

p_polar <- ggplot(polar_data, aes(x=route, y=mean_score, fill=mean_score)) +
  geom_col() + coord_polar() + facet_wrap(~cell_type) + theme_minimal() +
  scale_fill_gradient2(low="navy", high="firebrick")

ggsave(file.path(DIR_FIGURES, "SCP2601_rpac_polar_route_summary.png"), p_polar, width=16, height=12)

# =============================================================================
# 12. FINAL EXPORT
# =============================================================================

saveRDS(seu, file.path(DIR_OBJS, "SCP2601_with_rpac_spatial.rds"))
log_msg("21 complete. Final object saved.")