# =============================================================================
# scripts/05B_load_and_format_SCP2601.R
# Load Jian Shu Lab SCP2601 Spatial Seurat Objects & Inject Spatial Coordinates
# =============================================================================
#
# PURPOSE:
#   Load the SCP2601 dataset (Ounadjela et al. 2024, Nature Medicine) STARmap
#   and Slide-tag imputed Seurat objects, extract the physical X/Y spatial
#   coordinates, and inject them as a "spatial" DimReducObject so that the
#   RUN_SPATIAL hook in the downstream MISI/rPAC pipeline triggers automatically.
#
# CONTEXT (Golden Triangle):
#   - ORIGIN DATA: jian-shu-lab/hPlacenta-architecture (GitHub)
#       * Figures.r shows coordinates live in meta.data$x, meta.data$y
#       * add_spatial_metadata() creates CreateDimReducObject(key = "s_")
#       * Cell types in Clusters.2, stage in meta.data$stage
#       * Sample IDs: JS35 (W9), JS36 (W11), JS40 (W8-2)
#   - CORE MATH: RustemKurkcu/placenta-infection-pipeline
#       * scripts/16b_misi_v2_subindices.R (MISI)
#       * scripts/sandbox/21_rpac_v2_corrected_routes.R (rPAC)
#   - SCAFFOLDING: hPlacenta-spatial-thesis zip (target formatting)
#
# DATA SOURCES (SCP2601 on Broad Single Cell Portal):
#   A) STARmap-ISS CSV files:
#      - STARmap-ISS_sample_{W7,W8-2,W9,W11}_imputed_expression.csv
#      - STARmap-ISS_sample_{W7,W8-2,W9,W11}_spots_metadata.csv  (has x,y)
#      - STARmap-ISS_sample_{W7,W8-2,W9,W11}_cell_metadata.csv
#      - STARmap-ISS_sample_{W7,W8-2,W9,W11}_raw_expression.csv
#   B) Slide-tag / Multiome RDS objects (if pre-built by Jian Shu Lab):
#      - multiome_rna_seurat.rds
#      - starmap_spatial_raw_plus_imputed_seurat.rds
#   C) Metadata CSV files:
#      - humanplacenta_cluster.csv   (36,456 barcodes)
#      - humanplacenta_spatial.csv   (spatial coords for Slide-tags)
#      - metadata.csv                (cell annotations)
#
# OUTPUTS:
#   - output/objects/starmap_SCP2601_formatted.rds
#   - output/objects/slidetags_SCP2601_formatted.rds  (if available)
#   - output/objects/SCP2601_combined_spatial.rds
#   - output/tables/SCP2601_spatial_qc_summary.csv
#   - output/tables/SCP2601_gene_inventory.csv
#
# GUARDRAILS:
#   1. MegL (EC 4.4.1.11) cleaves L-methionine/L-cysteine → H₂S + NH₃.
#      MegL does NOT cleave EA.
#   2. EA is the transcriptional trigger (EutV/EutW), not a MegL substrate.
#   3. EPAS1 (HIF-2α), NOT HIF-1α, drives FLT1/sFLT1 in trophoblasts.
#
# EXECUTION:
#   Rscript scripts/05B_load_and_format_SCP2601.R
#
# AUTHOR:  Shan Kurkcu (PhD Thesis Pipeline)
# DATE:    2025-06-19
# VERSION: 1.0.0
# =============================================================================

cat("\n╔══════════════════════════════════════════════════════════════════════╗\n")
cat(  "║  05B — Load & Format SCP2601 Spatial Data (Ounadjela et al. 2024)  ║\n")
cat(  "╚══════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# 0. SETUP
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(Matrix)
  library(tibble)
})

# --- Pipeline infrastructure (if available; else define minimal) ---
tryCatch({
  source("config/config.R")
  source("scripts/R/utils.R")
  cat("[05B] Pipeline infrastructure loaded from config/config.R\n")
}, error = function(e) {
  cat("[05B] Running in standalone mode (no config.R found).\n")
})

# --- Directory setup ---
DIR_OBJS    <- "output/objects"
DIR_TABLES  <- "output/tables"
DIR_FIGURES <- "output/figures"
DIR_LOGS    <- "output/logs"
DATA_DIR    <- "data"  # Where SCP2601 raw CSV files live

for (d in c(DIR_OBJS, DIR_TABLES, DIR_FIGURES, DIR_LOGS, DATA_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# --- Logging ---
log_file <- file.path(DIR_LOGS, "05B_load_and_format_SCP2601.log")
log_msg <- function(...) {
  msg <- paste0(Sys.time(), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}
log_msg("05B starting.")

# --- Constants ---
SEED <- 42
set.seed(SEED)

# SCP2601 sample mapping (from Jian Shu Lab Figures.r):
#   JS35 → W9  (9w_1d)   → STARmap sample W9
#   JS36 → W11 (11w_1d)  → STARmap sample W11
#   JS40 → W8-2 (8w_2d)  → STARmap sample W8-2
#   JS34 → W7            → STARmap sample W7
SAMPLE_MAP <- data.frame(
  js_id       = c("JS34", "JS40", "JS35", "JS36"),
  sample_name = c("W7",   "W8-2", "W9",   "W11"),
  stage       = c("7w",   "8w_2d","9w_1d","11w_1d"),
  week        = c(7,       8,      9,      11),
  stringsAsFactors = FALSE
)

# Cell type palette (from Jian Shu Lab Figures.r, extended)
CELLTYPE_PALETTE <- c(
  Endothelial      = "#D42F2E",
  Erythroblasts    = "#F2AD00",
  EVT              = "darkgreen",
  `EVT-progenitor` = "darkseagreen2",
  FIB1             = "#F28239",
  FIB2             = "#FBE735",
  `Hofbauer cells` = "#8C9ECF",
  STB              = "#92D3E3",
  `STB-progenitor` = "deepskyblue4",
  vCTB             = "#E6A0C4"
)


# =============================================================================
# 1. LOCATE DATA FILES
# =============================================================================

log_msg("Scanning for SCP2601 data files...")

# Define search paths (current directory, data/, and workspace root)
# Define search paths based on the local hPlacenta-architecture directory tree
search_dirs <- unique(c(
  ".", 
  DATA_DIR, 
  "data/raw/zenodo_spatial",
  "data/raw/Broad_SCP2601human-placenta-architecture",
  "data/processed"
))

find_file <- function(pattern, dirs = search_dirs) {
  for (d in dirs) {
    hits <- list.files(d, pattern = pattern, full.names = TRUE, recursive = FALSE)
    # Filter out LFS pointers (< 200 bytes) and PNG files
    for (h in hits) {
      if (file.size(h) > 200) {
        # Check if it's actually a PNG masquerading as CSV
        con <- file(h, "rb")
        magic <- readBin(con, "raw", 4)
        close(con)
        if (!identical(magic, as.raw(c(0x89, 0x50, 0x4e, 0x47)))) {
          return(h)
        }
      }
    }
  }
  return(NULL)
}

# --- STARmap spot metadata (contains X/Y coordinates) ---
starmap_samples <- c("W7", "W8-2", "W9", "W11")
spot_meta_files <- list()
expr_imputed_files <- list()
expr_raw_files <- list()
cell_meta_files <- list()

for (samp in starmap_samples) {
  # Spots metadata (spatial coords)
  spot_meta_files[[samp]] <- find_file(
    paste0("STARmap-ISS_sample_", samp, "_spots_metadata\\.csv$")
  )
  # Imputed expression
  expr_imputed_files[[samp]] <- find_file(
    paste0("STARmap-ISS_sample_", samp, "_imputed_expression\\.csv$")
  )
  # Raw expression
  expr_raw_files[[samp]] <- find_file(
    paste0("STARmap-ISS_sample_", samp, "_raw_expression\\.csv$")
  )
  # Cell metadata
  cell_meta_files[[samp]] <- find_file(
    paste0("STARmap-ISS_sample_", samp, "_cell_metadata\\.csv$")
  )
}

# --- Slide-tag / Multiome pre-built RDS ---
multiome_rds <- find_file("multiome_rna_seurat\\.rds$")
starmap_rds  <- find_file("starmap_spatial_raw_plus_imputed_seurat\\.rds$")

# --- Metadata CSVs ---
spatial_csv  <- find_file("humanplacenta_spatial\\.csv$")
cluster_csv  <- find_file("humanplacenta_cluster\\.csv$")
metadata_csv <- find_file("metadata\\.csv$")

# Report what was found
log_msg("  Data inventory:")
for (samp in starmap_samples) {
  log_msg("  [", samp, "] spots_meta: ",
          ifelse(is.null(spot_meta_files[[samp]]), "MISSING", "FOUND"))
  log_msg("  [", samp, "] imputed_expr: ",
          ifelse(is.null(expr_imputed_files[[samp]]), "MISSING", "FOUND"))
  log_msg("  [", samp, "] raw_expr: ",
          ifelse(is.null(expr_raw_files[[samp]]), "MISSING", "FOUND"))
  log_msg("  [", samp, "] cell_meta: ",
          ifelse(is.null(cell_meta_files[[samp]]), "MISSING", "FOUND"))
}
log_msg("  multiome_rds: ", ifelse(is.null(multiome_rds), "MISSING", "FOUND"))
log_msg("  starmap_rds: ",  ifelse(is.null(starmap_rds), "MISSING", "FOUND"))
log_msg("  spatial_csv: ",  ifelse(is.null(spatial_csv), "MISSING", "FOUND"))
log_msg("  cluster_csv: ",  ifelse(is.null(cluster_csv), "MISSING", "FOUND"))
log_msg("  metadata_csv: ", ifelse(is.null(metadata_csv), "MISSING", "FOUND"))


# =============================================================================
# 2. LOAD FROM CSV (PRIMARY PATH: Build Seurat from CSV files)
# =============================================================================

#' Build a Seurat object from SCP2601 CSV files for one STARmap sample.
#'
#' @param sample_id  Character, e.g. "W7", "W8-2", "W9", "W11"
#' @param spot_meta_path  Path to spots_metadata CSV (must have x, y columns)
#' @param expr_path       Path to imputed expression CSV (genes x spots or spots x genes)
#' @param cell_meta_path  Optional path to cell_metadata CSV
#' @return Seurat object with spatial coordinates injected
build_starmap_seurat <- function(sample_id,
                                 spot_meta_path = NULL,
                                 expr_path = NULL,
                                 cell_meta_path = NULL) {
  log_msg("  Building Seurat for STARmap sample: ", sample_id)
  
  # ---- Read spot metadata (contains spatial coordinates) ----
  if (!is.null(spot_meta_path)) {
    log_msg("    Reading spots metadata: ", basename(spot_meta_path))
    spot_meta <- read.csv(spot_meta_path, stringsAsFactors = FALSE)
    log_msg("    Spots metadata: ", nrow(spot_meta), " spots x ",
            ncol(spot_meta), " columns")
    log_msg("    Columns: ", paste(colnames(spot_meta), collapse = ", "))
    
    # Identify X/Y columns (Jian Shu Lab convention: "x", "y" or "X", "Y")
    xy_pairs <- list(
      c("x", "y"), c("X", "Y"),
      c("x_um", "y_um"), c("X_um", "Y_um"),
      c("x_centroid", "y_centroid"),
      c("spatial_x", "spatial_y"),
      c("x_coord", "y_coord")
    )
    xy_cols <- NULL
    for (pair in xy_pairs) {
      if (all(pair %in% colnames(spot_meta))) {
        xy_cols <- pair
        break
      }
    }
    if (is.null(xy_cols)) {
      log_msg("    WARNING: No X/Y columns found in spots metadata!")
      log_msg("    Available columns: ", paste(colnames(spot_meta), collapse = ", "))
    } else {
      log_msg("    Spatial coordinates found: ", paste(xy_cols, collapse = ", "))
      log_msg("    X range: [", round(min(spot_meta[[xy_cols[1]]], na.rm = TRUE), 1),
              ", ", round(max(spot_meta[[xy_cols[1]]], na.rm = TRUE), 1), "]")
      log_msg("    Y range: [", round(min(spot_meta[[xy_cols[2]]], na.rm = TRUE), 1),
              ", ", round(max(spot_meta[[xy_cols[2]]], na.rm = TRUE), 1), "]")
    }
  } else {
    spot_meta <- NULL
    xy_cols <- NULL
    log_msg("    No spots metadata available.")
  }
  
  # ---- Read expression matrix ----
  if (!is.null(expr_path)) {
    log_msg("    Reading expression: ", basename(expr_path))
    expr_raw <- read.csv(expr_path, row.names = 1, check.names = FALSE,
                         stringsAsFactors = FALSE)
    log_msg("    Expression matrix: ", nrow(expr_raw), " x ", ncol(expr_raw))
    
    # Determine orientation: genes x cells or cells x genes
    # Heuristic: if nrow >> ncol, likely genes x cells
    # If ncol >> nrow, likely cells x genes
    # Also check if row names look like genes
    sample_rownames <- head(rownames(expr_raw), 10)
    sample_colnames <- head(colnames(expr_raw), 10)
    
    # Gene name heuristic: uppercase, short, no # or barcode patterns
    is_gene_like <- function(x) {
      all(nchar(x) < 30) && !any(grepl("#|ATCG", x))
    }
    
    if (is_gene_like(sample_rownames) && !is_gene_like(sample_colnames)) {
      log_msg("    Orientation: genes (rows) x cells (columns)")
      expr_mat <- as.matrix(expr_raw)
    } else if (!is_gene_like(sample_rownames) && is_gene_like(sample_colnames)) {
      log_msg("    Orientation: cells (rows) x genes (columns) — transposing")
      expr_mat <- t(as.matrix(expr_raw))
    } else {
      # Default: assume genes x cells if nrow > ncol
      if (nrow(expr_raw) > ncol(expr_raw)) {
        log_msg("    Assuming genes x cells (nrow > ncol)")
        expr_mat <- as.matrix(expr_raw)
      } else {
        log_msg("    Assuming cells x genes (ncol > nrow) — transposing")
        expr_mat <- t(as.matrix(expr_raw))
      }
    }
    log_msg("    Final matrix: ", nrow(expr_mat), " genes x ", ncol(expr_mat), " cells")
  } else {
    log_msg("    No expression data available for ", sample_id)
    return(NULL)
  }
  
  # ---- Read cell metadata (if available) ----
  cell_meta_df <- NULL
  if (!is.null(cell_meta_path)) {
    log_msg("    Reading cell metadata: ", basename(cell_meta_path))
    cell_meta_df <- read.csv(cell_meta_path, stringsAsFactors = FALSE)
    log_msg("    Cell metadata: ", nrow(cell_meta_df), " cells x ",
            ncol(cell_meta_df), " columns")
  }
  
  # ---- Build Seurat object ----
  log_msg("    Creating Seurat object...")
  expr_sparse <- as(expr_mat, "dgCMatrix")
  seu <- CreateSeuratObject(
    counts = expr_sparse,
    project = paste0("SCP2601_STARmap_", sample_id),
    min.cells = 0,
    min.features = 0
  )
  
  # Add sample-level metadata
  samp_info <- SAMPLE_MAP[SAMPLE_MAP$sample_name == sample_id, ]
  if (nrow(samp_info) > 0) {
    seu$sample_name <- samp_info$sample_name
    seu$stage       <- samp_info$stage
    seu$week        <- samp_info$week
    seu$js_id       <- samp_info$js_id
  } else {
    seu$sample_name <- sample_id
    seu$stage       <- sample_id
    seu$week        <- NA_integer_
    seu$js_id       <- NA_character_
  }
  seu$dataset     <- "STARmap_ISS"
  seu$data_source <- "SCP2601"
  
  # ---- Inject spatial coordinates ----
  if (!is.null(spot_meta) && !is.null(xy_cols)) {
    log_msg("    Injecting spatial coordinates into Seurat...")
    
    # Match spot metadata to cells in the Seurat object
    # Try common ID columns
    id_candidates <- c("cell_id", "spot_id", "barcode", "cell",
                       "CellID", "SpotID", "X1", "V1")
    id_col <- NULL
    if (any(rownames(spot_meta) %in% colnames(seu))) {
      spot_meta$match_id <- rownames(spot_meta)
      id_col <- "match_id"
    } else {
      for (idc in id_candidates) {
        if (idc %in% colnames(spot_meta)) {
          if (any(spot_meta[[idc]] %in% colnames(seu))) {
            id_col <- idc
            break
          }
        }
      }
    }
    
    if (is.null(id_col)) {
      # Try matching by position (row order)
      if (nrow(spot_meta) == ncol(seu)) {
        log_msg("    Matching by position (same count: ", nrow(spot_meta), ")")
        spot_meta$match_id <- colnames(seu)
        id_col <- "match_id"
      } else {
        log_msg("    WARNING: Cannot match spots to cells. Spots: ",
                nrow(spot_meta), " vs Cells: ", ncol(seu))
      }
    }
    
    if (!is.null(id_col)) {
      # Align spot metadata to Seurat cell order
      spot_ordered <- spot_meta[match(colnames(seu), spot_meta[[id_col]]), ]
      
      # Add x, y to metadata (Jian Shu Lab convention)
      seu$x <- spot_ordered[[xy_cols[1]]]
      seu$y <- spot_ordered[[xy_cols[2]]]
      
      # Also store as x_um, y_um for the scaffolding convention
      seu$x_um <- seu$x
      seu$y_um <- seu$y
      
      # Create the "spatial" DimReducObject
      # (This is the KEY step that makes RUN_SPATIAL = TRUE trigger)
      valid_cells <- !is.na(seu$x) & !is.na(seu$y)
      n_valid <- sum(valid_cells)
      log_msg("    Valid spatial coordinates: ", n_valid, " / ", ncol(seu))
      
      if (n_valid > 0) {
        emb <- cbind(seu$x, seu$y)
        colnames(emb) <- c("s_1", "s_2")
        rownames(emb) <- colnames(seu)
        
        # Replace NA with 0 for DimReduc (will be filtered downstream)
        emb[is.na(emb)] <- 0
        
        seu[["spatial"]] <- CreateDimReducObject(
          embeddings = emb,
          key = "s_",
          assay = DefaultAssay(seu)
        )
        log_msg("    ✓ 'spatial' DimReducObject created (key='s_', assay=",
                DefaultAssay(seu), ")")
        
        # Also add coordinate columns expected by the scaffolding
        seu$spatial_x <- seu$x
        seu$spatial_y <- seu$y
        seu$has_spatial <- TRUE
      }
    }
  }
  
  # ---- Merge cell metadata columns ----
  if (!is.null(cell_meta_df)) {
    log_msg("    Merging cell-level metadata...")
    # Try to match by rowname or first column
    if (any(rownames(cell_meta_df) %in% colnames(seu))) {
      for (col in setdiff(colnames(cell_meta_df), colnames(seu@meta.data))) {
        vals <- cell_meta_df[colnames(seu), col]
        if (!all(is.na(vals))) {
          seu[[col]] <- vals
        }
      }
    }
  }
  
  # ---- Detect cell type column (Jian Shu Lab: "Clusters.2") ----
  ct_candidates <- c("Clusters.2", "Clusters", "celltype", "cell_type",
                     "CellType", "cluster", "predicted.id", "predicted.celltype")
  for (ctc in ct_candidates) {
    if (ctc %in% colnames(seu@meta.data)) {
      seu$cell_type <- as.character(seu@meta.data[[ctc]])
      log_msg("    Cell type column: ", ctc, " → cell_type")
      break
    }
  }
  
  # ---- Normalize ----
  log_msg("    Normalizing...")
  seu <- NormalizeData(seu, verbose = FALSE)
  
  log_msg("    Done. Cells: ", ncol(seu), " | Genes: ", nrow(seu),
          " | Has spatial: ", "spatial" %in% names(seu@reductions))
  
  
  return(seu)
}


# =============================================================================
# 3. LOAD FROM PRE-BUILT RDS (FALLBACK PATH)
# =============================================================================

#' Load a pre-built Seurat RDS and inject spatial coordinates if needed.
#'
#' @param rds_path Path to .rds file
#' @param dataset_label Character label
#' @return Seurat object with spatial DimReduc, or NULL
load_and_format_rds <- function(rds_path, dataset_label = "SCP2601") {
  if (is.null(rds_path) || !file.exists(rds_path) || file.size(rds_path) < 500) {
    log_msg("  Skipping ", dataset_label, " RDS (missing or too small)")
    return(NULL)
  }
  
  log_msg("  Loading pre-built RDS: ", basename(rds_path))
  seu <- tryCatch(readRDS(rds_path), error = function(e) {
    log_msg("  ERROR reading RDS: ", conditionMessage(e))
    return(NULL)
  })
  
  if (is.null(seu) || !inherits(seu, "Seurat")) {
    log_msg("  Not a valid Seurat object.")
    return(NULL)
  }
  
  log_msg("  Loaded: ", ncol(seu), " cells x ", nrow(seu), " genes")
  log_msg("  Reductions: ", paste(names(seu@reductions), collapse = ", "))
  log_msg("  Metadata cols: ", paste(head(colnames(seu@meta.data), 20), collapse = ", "))
  
  # Check if spatial coordinates already exist
  if ("spatial" %in% names(seu@reductions)) {
    log_msg("  ✓ 'spatial' reduction already present.")
    return(seu)
  }
  
  # Try to find x, y in metadata (Jian Shu Lab convention)
  xy_pairs <- list(
    c("x", "y"), c("X", "Y"),
    c("x_um", "y_um"), c("x_centroid", "y_centroid"),
    c("spatial_x", "spatial_y"), c("x_coord", "y_coord")
  )
  
  for (pair in xy_pairs) {
    if (all(pair %in% colnames(seu@meta.data))) {
      log_msg("  Found coords: ", paste(pair, collapse = ", "),
              " → creating 'spatial' reduction")
      emb <- as.matrix(seu@meta.data[, pair])
      colnames(emb) <- c("s_1", "s_2")
      rownames(emb) <- colnames(seu)
      seu[["spatial"]] <- CreateDimReducObject(
        embeddings = emb, key = "s_", assay = DefaultAssay(seu)
      )
      seu$has_spatial <- TRUE
      log_msg("  ✓ 'spatial' DimReducObject created")
      break
    }
  }
  
  if (!"spatial" %in% names(seu@reductions)) {
    log_msg("  WARNING: No spatial coordinates found in pre-built RDS.")
    seu$has_spatial <- FALSE
  }
  
  return(seu)
}


# =============================================================================
# 4. EXECUTE: Build/Load all available objects
# =============================================================================

log_msg("\n=== Loading SCP2601 objects ===")

starmap_objects <- list()

# --- Strategy A: Build from CSV files ---
for (samp in starmap_samples) {
  if (!is.null(expr_imputed_files[[samp]]) || !is.null(expr_raw_files[[samp]])) {
    # Prefer imputed expression
    expr_path <- if (!is.null(expr_imputed_files[[samp]])) {
      expr_imputed_files[[samp]]
    } else {
      expr_raw_files[[samp]]
    }
    
    obj <- tryCatch(
      build_starmap_seurat(
        sample_id       = samp,
        spot_meta_path  = spot_meta_files[[samp]],
        expr_path       = expr_path,
        cell_meta_path  = cell_meta_files[[samp]]
      ),
      error = function(e) {
        log_msg("  ERROR building ", samp, ": ", conditionMessage(e))
        NULL
      }
    )
    
    if (!is.null(obj)) {
      starmap_objects[[samp]] <- obj
      log_msg("  ✓ ", samp, ": ", ncol(obj), " cells loaded from CSV")
    }
  }
}

# --- Strategy B: Load pre-built RDS (if CSV path yielded nothing) ---
if (length(starmap_objects) == 0) {
  log_msg("  No CSV-built objects. Trying pre-built RDS...")
  
  starmap_seu <- load_and_format_rds(starmap_rds, "STARmap RDS")
  if (!is.null(starmap_seu)) {
    # Split by sample if multi-sample
    if ("stage" %in% colnames(starmap_seu@meta.data)) {
      stages <- unique(starmap_seu$stage)
      for (stg in stages) {
        cells <- colnames(starmap_seu)[starmap_seu$stage == stg]
        if (length(cells) > 10) {
          sub_obj <- subset(starmap_seu, cells = cells)
          samp_match <- SAMPLE_MAP$sample_name[SAMPLE_MAP$stage == stg]
          if (length(samp_match) > 0) {
            starmap_objects[[samp_match[1]]] <- sub_obj
            log_msg("  ✓ ", samp_match[1], " (", stg, "): ",
                    ncol(sub_obj), " cells from RDS")
          }
        }
      }
    } else {
      starmap_objects[["combined"]] <- starmap_seu
    }
  }
}

# --- Multiome reference ---
multiome_seu <- load_and_format_rds(multiome_rds, "Multiome RDS")


# =============================================================================
# 5. MERGE INTO COMBINED SPATIAL OBJECT
# =============================================================================

log_msg("\n=== Merging SCP2601 spatial objects ===")

if (length(starmap_objects) > 1) {
  # Merge multiple STARmap samples
  base_obj <- starmap_objects[[1]]
  merge_list <- starmap_objects[-1]
  
  combined <- tryCatch({
    merged <- merge(base_obj, y = merge_list,
                    add.cell.ids = names(starmap_objects),
                    project = "SCP2601_combined")
    log_msg("  Merged ", length(starmap_objects), " samples: ",
            ncol(merged), " total cells")
    merged
  }, error = function(e) {
    log_msg("  Merge failed: ", conditionMessage(e))
    log_msg("  Using first sample only.")
    base_obj
  })
} else if (length(starmap_objects) == 1) {
  combined <- starmap_objects[[1]]
  log_msg("  Single sample: ", ncol(combined), " cells")
} else {
  log_msg("  WARNING: No SCP2601 spatial objects could be loaded.")
  log_msg("  Please ensure data files are downloaded from:")
  log_msg("    https://singlecell.broadinstitute.org/single_cell/study/SCP2601")
  log_msg("  Place CSV files in: ", DATA_DIR)
  combined <- NULL
}


# =============================================================================
# 6. POST-PROCESSING & QC
# =============================================================================

if (!is.null(combined)) {
  log_msg("\n=== Post-processing ===")
  
  # Normalize if not already done
  combined <- NormalizeData(combined, verbose = FALSE)
  
  # Find variable features (Using dispersion because data is imputed/continuous)
  combined <- FindVariableFeatures(combined, selection.method = "dispersion", nfeatures = 3000, verbose = FALSE)  
  # Scale data
  combined <- ScaleData(combined, verbose = FALSE)
  
  # PCA
  combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
  log_msg("  PCA complete (30 PCs)")
  
  # UMAP
  combined <- RunUMAP(combined, dims = 1:30, verbose = FALSE)
  log_msg("  UMAP complete")
  
  # Verify spatial reduction exists
  has_spatial_reduction <- "spatial" %in% names(combined@reductions)
  log_msg("  Has 'spatial' reduction: ", has_spatial_reduction)
  
  if (has_spatial_reduction) {
    # Verify coordinate range makes sense
    xy <- Embeddings(combined, "spatial")
    log_msg("  Spatial coord range: x=[", round(min(xy[, 1], na.rm = TRUE), 1),
            ", ", round(max(xy[, 1], na.rm = TRUE), 1),
            "] y=[", round(min(xy[, 2], na.rm = TRUE), 1),
            ", ", round(max(xy[, 2], na.rm = TRUE), 1), "]")
  }
  
  # ---- Gene inventory for thesis-relevant genes ----
  log_msg("  Checking thesis-relevant gene availability...")
  
  thesis_genes <- list(
    # IDO1 Tolerogenic Shield (12 genes — from Part 3)
    IDO1_Tolerogenic_Shield = c("IDO1", "TGFB1", "IL10", "CD274", "PDCD1LG2",
                                "HAVCR2", "LGALS9", "CD80", "CD86",
                                "ENTPD1", "NT5E", "FOXP3"),
    # MISI modules (key markers)
    MISI_Homing     = c("GALNT1", "GALNT2", "CDH1", "CDH5", "OCLN"),
    MISI_Nutrient   = c("PLD1", "PLD2", "PCYT2", "ETNK1"),
    MISI_Switch     = c("TLR2", "TLR4", "TNF", "IL6", "HIF1A"),
    MISI_Vascular   = c("FLT1", "PGF", "KDR", "VEGFA", "ENG"),
    MISI_Invasion   = c("MMP1", "MMP2", "MMP9", "TIMP1", "TIMP2"),
    # rPAC routes (key TFs and targets)
    rPAC_EPAS1_axis = c("EPAS1", "ARNT", "FLT1", "PGF", "DLL4", "ANGPT2"),
    rPAC_NFKB       = c("RELA", "NFKB1", "MYD88", "TRAF6"),
    # Fn-specific response markers
    Fn_response     = c("NFKB1", "RELA", "EPAS1", "IL6", "CXCL8", "TNF",
                        "PTGS2", "CCL2", "IDO1")
  )
  
  gene_inventory <- do.call(rbind, lapply(names(thesis_genes), function(mod) {
    genes <- thesis_genes[[mod]]
    present <- genes %in% rownames(combined)
    data.frame(
      module = mod,
      gene = genes,
      present = present,
      stringsAsFactors = FALSE
    )
  }))
  
  write.csv(gene_inventory, file.path(DIR_TABLES, "SCP2601_gene_inventory.csv"),
            row.names = FALSE)
  
  n_present <- sum(gene_inventory$present)
  n_total   <- nrow(gene_inventory)
  log_msg("  Gene inventory: ", n_present, "/", n_total, " thesis genes detected")
  
  # Per-module summary
  gene_summary <- gene_inventory %>%
    group_by(module) %>%
    summarise(
      n_genes = n(),
      n_present = sum(present),
      pct_present = round(100 * mean(present), 1),
      missing = paste(gene[!present], collapse = ", "),
      .groups = "drop"
    )
  log_msg("  Module detection rates:")
  for (i in seq_len(nrow(gene_summary))) {
    log_msg("    ", gene_summary$module[i], ": ",
            gene_summary$n_present[i], "/", gene_summary$n_genes[i],
            " (", gene_summary$pct_present[i], "%)",
            if (nzchar(gene_summary$missing[i])) paste0(" [missing: ", gene_summary$missing[i], "]") else "")
  }
  
  
  # ---- QC Summary table ----
  qc_df <- data.frame(
    metric = c("n_cells", "n_genes", "n_samples",
               "has_spatial_reduction", "n_thesis_genes_detected",
               "median_nCount_RNA", "median_nFeature_RNA"),
    value = c(
      ncol(combined),
      nrow(combined),
      length(unique(combined$sample_name)),
      has_spatial_reduction,
      n_present,
      median(combined$nCount_RNA, na.rm = TRUE),
      median(combined$nFeature_RNA, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  write.csv(qc_df, file.path(DIR_TABLES, "SCP2601_spatial_qc_summary.csv"),
            row.names = FALSE)
  log_msg("  QC summary saved.")
  
  
  # =============================================================================
  # 7. SAVE OUTPUTS
  # =============================================================================
  
  log_msg("\n=== Saving outputs ===")
  
  # Save individual sample objects
  for (samp in names(starmap_objects)) {
    out_path <- file.path(DIR_OBJS, paste0("starmap_", samp, "_SCP2601.rds"))
    saveRDS(starmap_objects[[samp]], out_path)
    log_msg("  Saved: ", basename(out_path))
  }
  
  # Save combined object
  combined_path <- file.path(DIR_OBJS, "SCP2601_combined_spatial.rds")
  saveRDS(combined, combined_path)
  log_msg("  Saved: ", basename(combined_path), " (",
          ncol(combined), " cells)")
  
  # Save multiome reference if available
  if (!is.null(multiome_seu)) {
    out_path <- file.path(DIR_OBJS, "multiome_SCP2601_formatted.rds")
    saveRDS(multiome_seu, out_path)
    log_msg("  Saved: ", basename(out_path))
  }
  
} else {
  log_msg("\n  ⚠ No data loaded. Cannot proceed.")
  log_msg("  Download data from SCP2601 and re-run this script.")
}


# =============================================================================
# 8. FINAL REPORT
# =============================================================================

cat("\n╔══════════════════════════════════════════════════════════════════════╗\n")
cat(  "║  05B — SCP2601 Loading Complete                                    ║\n")
cat(  "╠══════════════════════════════════════════════════════════════════════╣\n")

if (!is.null(combined)) {
  cat(sprintf("║  Total cells: %-50d ║\n", ncol(combined)))
  cat(sprintf("║  Total genes: %-50d ║\n", nrow(combined)))
  cat(sprintf("║  Samples: %-54s ║\n",
              paste(unique(combined$sample_name), collapse = ", ")))
  cat(sprintf("║  Spatial reduction: %-44s ║\n",
              ifelse("spatial" %in% names(combined@reductions), "YES ✓", "NO ✗")))
  cat(sprintf("║  RUN_SPATIAL trigger: %-42s ║\n",
              ifelse("spatial" %in% names(combined@reductions),
                     "READY (will auto-fire)", "NOT AVAILABLE")))
  cat(sprintf("║  Thesis genes: %d / %d detected%-24s ║\n",
              n_present, n_total, ""))
} else {
  cat("║  WARNING: No data loaded.                                         ║\n")
  cat("║  Download SCP2601 CSV files and re-run.                           ║\n")
}

cat("╠══════════════════════════════════════════════════════════════════════╣\n")
cat("║  Outputs:                                                          ║\n")
cat("║    output/objects/SCP2601_combined_spatial.rds                      ║\n")
cat("║    output/objects/starmap_{sample}_SCP2601.rds                     ║\n")
cat("║    output/tables/SCP2601_gene_inventory.csv                        ║\n")
cat("║    output/tables/SCP2601_spatial_qc_summary.csv                    ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

log_msg("05B complete.")

