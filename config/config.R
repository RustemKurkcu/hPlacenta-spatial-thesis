# ======================================================================
# config/config.R
# Global configuration for the Placenta Vicious Cycle analysis pipeline.
#
# Project: Spatiotemporal Dynamics of Placental Nutritional Immunity
# Author: (you)
# Updated: 2026-02-08
# ======================================================================

# -------------------------
# Input objects (edit these paths)
# -------------------------
# If you already built the Seurat objects, point to them here.
PATH_MULTIOME_RDS      <- "data/processed/multiome_rna_seurat.rds"
PATH_SLIDETAGS_RDS     <- "data/processed/slidetags_mapped_to_multiome.rds"  # preferred (already has predicted.id)
PATH_SLIDETAGS_RAW     <- "data/processed/slidetags_spatial_rna_seurat.rds"  # fallback (no mapping)
PATH_STARMAP_RDS       <- "data/processed/starmap_spatial_raw_plus_imputed_seurat.rds"

# Optional: raw inputs (only used by scripts in scripts/01_build_objects/)
DIR_RAW_SLIDETAGS      <- "data/raw/Broad_SCP2601human-placenta-architecture"
DIR_RAW_STARMAP        <- "data/raw/zenodo_spatial"

# -------------------------
# Output directories
# -------------------------
DIR_OBJS    <- "output/objects"
DIR_FIGURES <- "output/figures"
DIR_TABLES  <- "output/tables"
DIR_LOGS    <- "output/logs"
DIR_DOCS    <- "docs"

# -------------------------
# Column names / candidates (pipeline will auto-detect where possible)
# -------------------------
# Week / time
COL_WEEK               <- "week"
COL_WEEK_CANDIDATES    <- c("week", "gestational_week", "gestational.age", "GA", "age", "biosample_id", "sample_id", "orig.ident", "donor_id")

# Cell type annotations
COL_REF_CELLTYPE       <- "cluster"         # multiome reference annotation column (in meta.data)
COL_AUTHOR_CELLTYPE    <- c("cell_type_spatial", "cell_type_cluster", "cell_type", "celltype", "annotation", "cluster")

# Predicted labels from Seurat TransferData/MapQuery
COL_PRED_CELLTYPE      <- "predicted.id"
COL_PRED_SCORE_MAX     <- "prediction.score.max"

# Spatial coordinates
COL_SPATIAL_X_CANDIDATES <- c("spatial_x", "spatial_x_um", "X", "x", "pxl_col_in_fullres", "center_x", "X_centroid")
COL_SPATIAL_Y_CANDIDATES <- c("spatial_y", "spatial_y_um", "Y", "y", "pxl_row_in_fullres", "center_y", "Y_centroid")

# -------------------------
# Analysis parameters
# -------------------------
SEED <- 1
DEFAULT_DIMS <- 1:30

# Cell-type harmonization rules
PRED_SCORE_HIGH   <- 0.70   # use predicted label if author label missing/unknown AND score >= this
PRED_SCORE_REFINE <- 0.85   # optionally refine broad author labels using predicted subtypes if score >= this

# Spatial neighborhood settings
SPATIAL_KNN_K      <- 15
SPATIAL_N_PERM     <- 200   # label permutations for enrichment z-scores (increase for final)

# -------------------------
# Optional features
# -------------------------
# MapQuery projection (projects query into reference UMAP model)
# Default OFF because it requires tight assay/normalization matching.
DO_MAPQUERY <- FALSE
DO_TSNE <- TRUE  # compute tSNE embeddings for query objects (Slide-tags / STARmap)

# If TRUE, will run optional heavy modules if packages available (CellChat, etc.)
RUN_OPTIONAL_HEAVY <- TRUE

# -------------------------
# Gene sets
# -------------------------
# Gene sets are defined in: config/gene_sets.R
source(file.path("config", "gene_sets.R"))

# Ligand-receptor pairs (simple, lightweight interaction scoring)
PATH_LR_PAIRS <- file.path("config", "ligand_receptor_pairs.csv")
