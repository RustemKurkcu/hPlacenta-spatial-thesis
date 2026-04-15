# =============================================================================
# Script: 01_preprocess_harmony_embeddings.R
# Purpose: Preprocess and integrate 4-week STARmap data (W7, W8-2, W9, W11)
#          using Seurat v5 + SCTransform + Harmony, then produce UMAP/t-SNE.
#
# Thesis Context:
#   - This script is Script 01 of the active spatial thesis pipeline.
#   - It is designed for defense-ready reproducibility and full auditability.
#
# Sections (A-J):
#   A) Header + Run Manifest
#   B) Logging + Config
#   C) Data Ingestion + Metadata Harmonization
#   D) Per-week Seurat Object Construction + QC
#   E) Merge + Seurat v5 Hygiene
#   F) SCTransform Workflow
#   G) PCA + Harmony Integration
#   H) UMAP + t-SNE
#   I) Clustering + Diagnostics
#   J) Export Artifacts
#
# Author: Thesis Pipeline Team
# Version: 1.0.1
# Date: 2026-04-15
# =============================================================================

# -----------------------------------------------------------------------------
# Dependency checks (self-contained preflight)
# -----------------------------------------------------------------------------
required_pkgs <- c(
  "Seurat",
  "SeuratObject",
  "harmony",
  "dplyr",
  "ggplot2",
  "patchwork",
  "Matrix",
  "jsonlite",
  "readr",
  "stringr",
  "purrr",
  "tibble"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall these packages before running this script."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(jsonlite)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
})

# =============================================================================
# A) Header + Run Manifest
# =============================================================================

PIPELINE_NAME <- "01_preprocess_harmony_embeddings"
PIPELINE_VERSION <- "1.0.1"
RUN_TIMESTAMP <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
RUN_SEED <- 42L
set.seed(RUN_SEED)

# =============================================================================
# B) Logging + Config
# =============================================================================

# ---- Output layout ----
OUT_ROOT <- "output"
DIR_LOGS <- file.path(OUT_ROOT, "logs")
DIR_OBJECTS <- file.path(OUT_ROOT, "objects")
DIR_TABLES <- file.path(OUT_ROOT, "tables")
DIR_FIGURES <- file.path(OUT_ROOT, "figures")
DIR_REPORTS <- file.path(OUT_ROOT, "reports")

for (d in c(DIR_LOGS, DIR_OBJECTS, DIR_TABLES, DIR_FIGURES, DIR_REPORTS)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG_FILE <- file.path(DIR_LOGS, "01_preprocess_harmony_embeddings.log")

# ---- Robust logging helper ----
log_msg <- function(..., .level = "INFO") {
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  txt <- paste0(..., collapse = "")
  line <- paste0(stamp, " [", .level, "] ", txt)
  cat(line, "\n")
  cat(line, "\n", file = LOG_FILE, append = TRUE)
}

log_msg("============================================================")
log_msg("Starting ", PIPELINE_NAME, " v", PIPELINE_VERSION)
log_msg("Run timestamp: ", RUN_TIMESTAMP)
log_msg("Seed: ", RUN_SEED)
log_msg("============================================================")

# ---- Config block (tunable, fully logged) ----
cfg <- list(
  # Data roots (course-corrected to match local project layout)
  data_roots = list(
    broad = "data/raw/Broad_SCP2601human-placenta-architecture",
    zenodo = "data/raw/zenodo_spatial"
  ),
  week_ids = c("W7", "W8-2", "W9", "W11"),

  # STARmap file conventions used by Jian Shu lab exports
  expression_patterns = c(
    "STARmap-ISS_sample_{week}_imputed_expression\\.csv$",
    "STARmap-ISS_sample_{week}_raw_expression\\.csv$"
  ),
  spots_metadata_patterns = c(
    "STARmap-ISS_sample_{week}_spots_metadata\\.csv$"
  ),
  cell_metadata_patterns = c(
    "STARmap-ISS_sample_{week}_cell_metadata\\.csv$"
  ),

  # QC thresholds (selected to balance STARmap sparsity vs biological retention)
  qc = list(
    min_features = 150L,
    max_features = 9000L,
    min_counts = 300L,
    max_counts = 80000L,
    max_percent_mt = 20
  ),

  # SCTransform settings
  sct = list(
    vst_flavor = "v2",
    variable_features_n = 3000L,
    regress_vars = c("percent.mt"),
    conserve_memory = TRUE,
    return_only_var_genes = FALSE
  ),

  # Dimensionality / integration settings
  dims_use = 1:30,
  harmony = list(
    group_by = "week",
    reduction = "pca",
    assay_use = "SCT",
    max_iter_harmony = 20L,
    theta = 2,
    lambda = 1,
    sigma = 0.1
  ),

  # Embeddings
  umap = list(
    n_neighbors = 30L,
    min_dist = 0.3,
    spread = 1.0,
    metric = "cosine"
  ),
  tsne = list(
    perplexity = 30,
    check_duplicates = FALSE
  ),

  # Clustering
  clustering = list(
    resolution = 0.6,
    algorithm = 1L
  )
)

log_msg("Config loaded. Weeks: ", paste(cfg$week_ids, collapse = ", "))
log_msg("Data roots: broad=", cfg$data_roots$broad, " ; zenodo=", cfg$data_roots$zenodo)
log_msg("QC thresholds: min_features=", cfg$qc$min_features,
        ", max_features=", cfg$qc$max_features,
        ", min_counts=", cfg$qc$min_counts,
        ", max_counts=", cfg$qc$max_counts,
        ", max_percent_mt=", cfg$qc$max_percent_mt)

append_methods_provenance <- function(report_path, week, expression_file, spots_file, cellmeta_file, notes = NULL) {
  if (!file.exists(report_path)) {
    header <- c(
      "# 01 Methods and Data Provenance Log",
      "",
      "This file is auto-appended by `01_preprocess_harmony_embeddings.R`.",
      "Each run records the exact modalities/files used for each week.",
      ""
    )
    writeLines(header, con = report_path)
  }

  lines <- c(
    paste0("## Run: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("- Script: ", PIPELINE_NAME, " v", PIPELINE_VERSION),
    paste0("- Week: ", week),
    "- Modalities used:",
    paste0("  - Expression matrix: ", expression_file),
    paste0("  - Spatial coordinates (spots metadata): ", spots_file),
    paste0("  - Cell metadata (optional): ", ifelse(is.na(cellmeta_file), "MISSING", cellmeta_file)),
    "- Method notes:",
    "  - Expression + spots metadata were sourced from split raw directories (Broad + Zenodo).",
    "  - Spatial coordinates are harmonized to x_um/y_um before Seurat construction.",
    "  - SCTransform + Harmony embeddings are generated downstream in this script."
  )
  if (!is.null(notes) && length(notes) > 0) lines <- c(lines, paste0("- Additional notes: ", notes))
  lines <- c(lines, "")
  write(lines, file = report_path, append = TRUE)
}

# =============================================================================
# C) Data Ingestion + Metadata Harmonization
# =============================================================================

find_first_regex <- function(roots, patterns, week, prefer_root = NULL) {
  roots_order <- roots
  if (!is.null(prefer_root) && prefer_root %in% roots) {
    roots_order <- c(prefer_root, setdiff(roots, prefer_root))
  }
  for (root in roots_order) {
    if (!dir.exists(root)) next
    for (pat in patterns) {
      rx <- stringr::str_replace_all(pat, "\\{week\\}", week)
      hits <- list.files(root, pattern = rx, full.names = TRUE, recursive = TRUE)
      if (length(hits) > 0) return(hits[[1]])
    }
  }
  NA_character_
}

read_csv_flexible <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

standardize_spots_metadata <- function(md_raw, week_id, barcodes) {
  md <- md_raw

  barcode_candidates <- c("barcode", "barcodes", "cell", "cell_id", "spot_id", "X")
  barcode_col <- intersect(barcode_candidates, names(md))
  if (length(barcode_col) == 0) {
    stop("No barcode-like column found in spots metadata for week ", week_id)
  }
  barcode_col <- barcode_col[[1]]

  x_candidates <- c("x", "x_um", "spatial_x", "X_coord", "xcoord")
  y_candidates <- c("y", "y_um", "spatial_y", "Y_coord", "ycoord")
  x_col <- intersect(x_candidates, names(md))
  y_col <- intersect(y_candidates, names(md))
  if (length(x_col) == 0 || length(y_col) == 0) {
    stop("Missing coordinate columns in metadata for week ", week_id)
  }
  x_col <- x_col[[1]]
  y_col <- y_col[[1]]

  out <- md %>%
    mutate(
      barcode = as.character(.data[[barcode_col]]),
      x_um = as.numeric(.data[[x_col]]),
      y_um = as.numeric(.data[[y_col]]),
      week = week_id,
      sample_id = week_id
    )

  out <- out %>% distinct(barcode, .keep_all = TRUE)

  if (length(barcodes) > 0) {
    # Keep only cells present in counts matrix
    out <- out %>% filter(barcode %in% barcodes)
    # Make rownames-compatible ordering
    out <- out[match(barcodes, out$barcode), , drop = FALSE]

    missing_md <- sum(is.na(out$barcode))
    if (missing_md > 0) {
      stop("Metadata missing for ", missing_md, " barcodes in week ", week_id)
    }
  }

  out
}

coerce_expression_with_spots <- function(expr_path, spots_md, week_id) {
  expr <- read_csv_flexible(expr_path)
  if (ncol(expr) < 2) {
    stop("Expression CSV has too few columns for week ", week_id, ": ", expr_path)
  }

  first_name <- names(expr)[1]
  if (tolower(first_name) %in% c("x", "gene", "gene_name", "genes", "feature", "feature_id", "")) {
    expr <- expr %>% rename(gene = 1)
  } else {
    expr <- expr %>% rename(gene = 1)
  }
  expr$gene <- as.character(expr$gene)
  expr <- expr %>% filter(!is.na(gene), gene != "") %>% distinct(gene, .keep_all = TRUE)

  candidate_cols <- setdiff(colnames(expr), "gene")
  spot_barcodes <- as.character(spots_md$barcode)
  col_overlap <- sum(candidate_cols %in% spot_barcodes)
  row_overlap <- sum(expr$gene %in% spot_barcodes)

  if (row_overlap > col_overlap) {
    mat_raw <- as.matrix(expr[, setdiff(colnames(expr), "gene"), drop = FALSE])
    suppressWarnings(storage.mode(mat_raw) <- "numeric")
    rownames(mat_raw) <- expr$gene
    mat_raw <- t(mat_raw)
    colnames(mat_raw) <- make.unique(colnames(mat_raw))
    rownames(mat_raw) <- make.unique(rownames(mat_raw))
    mat <- mat_raw
  } else {
    mat_raw <- as.matrix(expr[, candidate_cols, drop = FALSE])
    suppressWarnings(storage.mode(mat_raw) <- "numeric")
    rownames(mat_raw) <- expr$gene
    colnames(mat_raw) <- candidate_cols
    mat <- mat_raw
  }

  keep_barcodes <- intersect(colnames(mat), spot_barcodes)
  if (length(keep_barcodes) == 0) {
    stop("No overlapping barcodes between expression and spots metadata for week ", week_id)
  }
  mat <- mat[, keep_barcodes, drop = FALSE]
  mat <- as(Matrix(mat, sparse = TRUE), "dgCMatrix")
  mat
}

week_inputs <- lapply(cfg$week_ids, function(wk) {
  list(
    week = wk,
    expression = find_first_regex(
      roots = c(cfg$data_roots$broad, cfg$data_roots$zenodo),
      patterns = cfg$expression_patterns,
      week = wk,
      prefer_root = cfg$data_roots$broad
    ),
    spots_metadata = find_first_regex(
      roots = c(cfg$data_roots$zenodo, cfg$data_roots$broad),
      patterns = cfg$spots_metadata_patterns,
      week = wk,
      prefer_root = cfg$data_roots$zenodo
    ),
    cell_metadata = find_first_regex(
      roots = c(cfg$data_roots$zenodo, cfg$data_roots$broad),
      patterns = cfg$cell_metadata_patterns,
      week = wk,
      prefer_root = cfg$data_roots$zenodo
    )
  )
})
names(week_inputs) <- cfg$week_ids

for (wk in cfg$week_ids) {
  wi <- week_inputs[[wk]]
  log_msg("Week ", wk, " file discovery:")
  log_msg("  expression      : ", wi$expression)
  log_msg("  spots_metadata  : ", wi$spots_metadata)
  log_msg("  cell_metadata   : ", wi$cell_metadata)
  if (is.na(wi$expression) || is.na(wi$spots_metadata)) {
    stop("Required files missing for week ", wk, ". Check Broad/Zenodo data roots and STARmap naming patterns.")
  }
}

# =============================================================================
# D) Per-week Seurat Object Construction + QC
# =============================================================================

week_objects <- list()
qc_summary <- list()

for (wk in cfg$week_ids) {
  log_msg("--- Processing week ", wk, " ---")
  wi <- week_inputs[[wk]]

  spots_md_raw <- read_csv_flexible(wi$spots_metadata)
  names(spots_md_raw) <- make.names(names(spots_md_raw), unique = TRUE)
  spots_md <- standardize_spots_metadata(spots_md_raw, wk, barcodes = character(0))
  counts <- coerce_expression_with_spots(wi$expression, spots_md, wk)
  log_msg("  Counts loaded from STARmap CSV: genes=", nrow(counts), ", cells=", ncol(counts))

  md <- standardize_spots_metadata(spots_md_raw, wk, colnames(counts))

  seu_w <- CreateSeuratObject(counts = counts, project = paste0("STARmap_", wk))

  # Add harmonized metadata in cell order
  rownames(md) <- md$barcode
  md <- md[colnames(seu_w), , drop = FALSE]
  seu_w <- AddMetaData(seu_w, metadata = md)

  # Mito percent; if MT genes absent, set 0 and log warning
  mt_pattern <- "^MT-"
  mt_features <- grep(mt_pattern, rownames(seu_w), value = TRUE)
  if (length(mt_features) > 0) {
    seu_w[["percent.mt"]] <- PercentageFeatureSet(seu_w, pattern = mt_pattern)
  } else {
    seu_w[["percent.mt"]] <- 0
    log_msg("  No MT- genes detected for week ", wk, ". percent.mt set to 0.", .level = "WARN")
  }

  # QC filters
  n_before <- ncol(seu_w)
  keep <- (
    seu_w$nFeature_RNA >= cfg$qc$min_features &
      seu_w$nFeature_RNA <= cfg$qc$max_features &
      seu_w$nCount_RNA >= cfg$qc$min_counts &
      seu_w$nCount_RNA <= cfg$qc$max_counts &
      seu_w$percent.mt <= cfg$qc$max_percent_mt &
      !is.na(seu_w$x_um) & !is.na(seu_w$y_um)
  )

  seu_w <- subset(seu_w, cells = colnames(seu_w)[keep])
  n_after <- ncol(seu_w)

  log_msg("  QC retained ", n_after, "/", n_before, " cells (", round(100 * n_after / max(n_before, 1), 2), "%).")

  qc_summary[[wk]] <- tibble(
    week = wk,
    cells_before_qc = n_before,
    cells_after_qc = n_after,
    pct_retained = round(100 * n_after / max(n_before, 1), 2),
    median_nFeature = median(seu_w$nFeature_RNA),
    median_nCount = median(seu_w$nCount_RNA),
    median_percent_mt = median(seu_w$percent.mt)
  )

  append_methods_provenance(
    report_path = file.path(DIR_REPORTS, "01_methods_and_provenance.md"),
    week = wk,
    expression_file = wi$expression,
    spots_file = wi$spots_metadata,
    cellmeta_file = wi$cell_metadata
  )

  week_objects[[wk]] <- seu_w
}

qc_summary_tbl <- bind_rows(qc_summary)
write_csv(qc_summary_tbl, file.path(DIR_TABLES, "01_qc_summary_by_week.csv"))
log_msg("Per-week QC summary saved.")

# =============================================================================
# E) Merge + Seurat v5 Hygiene
# =============================================================================

log_msg("Merging week-level Seurat objects...")
seu <- Reduce(function(x, y) merge(x, y), week_objects)

# Seurat v5 required for merged-layer workflows
log_msg("Applying JoinLayers() for Seurat v5 compatibility.")
seu <- JoinLayers(seu)
DefaultAssay(seu) <- "RNA"

# Ensure harmonized week column exists
if (!"week" %in% colnames(seu@meta.data)) {
  stop("Merged object missing 'week' metadata column.")
}

# Explicit check for W9 inclusion
wk_counts <- table(seu$week)
log_msg("Cells by week after merge+QC: ", paste(names(wk_counts), as.integer(wk_counts), sep = "=", collapse = "; "))
if (!"W9" %in% names(wk_counts)) {
  stop("W9 is missing after merge/QC. Abort to protect 4-week integration objective.")
}

# =============================================================================
# F) SCTransform Workflow
# =============================================================================

log_msg("Running SCTransform on merged object...")
seu <- SCTransform(
  object = seu,
  assay = "RNA",
  new.assay.name = "SCT",
  vst.flavor = cfg$sct$vst_flavor,
  variable.features.n = cfg$sct$variable_features_n,
  vars.to.regress = cfg$sct$regress_vars,
  conserve.memory = cfg$sct$conserve_memory,
  return.only.var.genes = cfg$sct$return_only_var_genes,
  verbose = TRUE
)
DefaultAssay(seu) <- "SCT"
log_msg("SCTransform complete. Default assay set to SCT.")

# =============================================================================
# G) PCA + Harmony Integration
# =============================================================================

log_msg("Running PCA...")
seu <- RunPCA(seu, assay = "SCT", npcs = max(cfg$dims_use), verbose = FALSE)

log_msg("Running Harmony integration by: ", cfg$harmony$group_by)
seu <- RunHarmony(
  object = seu,
  group.by.vars = cfg$harmony$group_by,
  reduction = cfg$harmony$reduction,
  assay.use = cfg$harmony$assay_use,
  max.iter.harmony = cfg$harmony$max_iter_harmony,
  theta = cfg$harmony$theta,
  lambda = cfg$harmony$lambda,
  sigma = cfg$harmony$sigma,
  verbose = TRUE
)

# =============================================================================
# H) UMAP + t-SNE
# =============================================================================

log_msg("Computing UMAP (global structure) on Harmony dims...")
seu <- RunUMAP(
  object = seu,
  reduction = "harmony",
  dims = cfg$dims_use,
  reduction.name = "umap_harmony",
  n.neighbors = cfg$umap$n_neighbors,
  min.dist = cfg$umap$min_dist,
  spread = cfg$umap$spread,
  metric = cfg$umap$metric,
  verbose = TRUE
)

log_msg("Computing t-SNE (local neighborhoods) on Harmony dims...")
seu <- RunTSNE(
  object = seu,
  reduction = "harmony",
  dims = cfg$dims_use,
  reduction.name = "tsne_harmony",
  perplexity = cfg$tsne$perplexity,
  check_duplicates = cfg$tsne$check_duplicates,
  verbose = TRUE
)

# Core embedding figures
p_umap_week <- DimPlot(seu, reduction = "umap_harmony", group.by = "week", pt.size = 0.2) +
  ggtitle("UMAP (Harmony) by Week")
p_tsne_week <- DimPlot(seu, reduction = "tsne_harmony", group.by = "week", pt.size = 0.2) +
  ggtitle("t-SNE (Harmony) by Week")

pdf(file.path(DIR_FIGURES, "01_embeddings_week_umap_tsne.pdf"), width = 14, height = 6)
print(p_umap_week + p_tsne_week)
dev.off()

ggsave(file.path(DIR_FIGURES, "01_umap_harmony_by_week.png"), p_umap_week, width = 8, height = 6, dpi = 300)
ggsave(file.path(DIR_FIGURES, "01_tsne_harmony_by_week.png"), p_tsne_week, width = 8, height = 6, dpi = 300)
log_msg("Embedding figures saved.")

# =============================================================================
# I) Clustering + Diagnostics
# =============================================================================

log_msg("Running neighbors + clusters on Harmony space...")
seu <- FindNeighbors(seu, reduction = "harmony", dims = cfg$dims_use, verbose = FALSE)
seu <- FindClusters(
  seu,
  resolution = cfg$clustering$resolution,
  algorithm = cfg$clustering$algorithm,
  verbose = FALSE
)

# Basic diagnostics table
diag_tbl <- seu@meta.data %>%
  as_tibble(rownames = "cell_id") %>%
  count(week, seurat_clusters, name = "n_cells") %>%
  arrange(week, seurat_clusters)

write_csv(diag_tbl, file.path(DIR_TABLES, "01_cluster_composition_by_week.csv"))
log_msg("Cluster composition diagnostics saved.")

# Cluster visualization
p_umap_cluster <- DimPlot(seu, reduction = "umap_harmony", group.by = "seurat_clusters", label = TRUE, pt.size = 0.2) +
  ggtitle("UMAP (Harmony) by Cluster")
p_tsne_cluster <- DimPlot(seu, reduction = "tsne_harmony", group.by = "seurat_clusters", label = TRUE, pt.size = 0.2) +
  ggtitle("t-SNE (Harmony) by Cluster")

pdf(file.path(DIR_FIGURES, "01_embeddings_clusters_umap_tsne.pdf"), width = 14, height = 6)
print(p_umap_cluster + p_tsne_cluster)
dev.off()

# =============================================================================
# J) Export Artifacts
# =============================================================================

obj_path <- file.path(DIR_OBJECTS, "01_integrated_harmony_sct_4weeks.rds")
saveRDS(seu, obj_path)
log_msg("Saved integrated Seurat object: ", obj_path)

# Run manifest for reproducibility
manifest <- list(
  pipeline = PIPELINE_NAME,
  version = PIPELINE_VERSION,
  run_timestamp = RUN_TIMESTAMP,
  seed = RUN_SEED,
  weeks_requested = cfg$week_ids,
  weeks_observed = sort(unique(as.character(seu$week))),
  n_cells_final = ncol(seu),
  n_genes_final = nrow(seu),
  qc = cfg$qc,
  sct = cfg$sct,
  harmony = cfg$harmony,
  umap = cfg$umap,
  tsne = cfg$tsne,
  clustering = cfg$clustering,
  object_output = obj_path,
  tables_output = c(
    file.path(DIR_TABLES, "01_qc_summary_by_week.csv"),
    file.path(DIR_TABLES, "01_cluster_composition_by_week.csv")
  ),
  reports_output = c(
    file.path(DIR_REPORTS, "01_methods_and_provenance.md"),
    file.path(DIR_REPORTS, "01_preprocess_harmony_manifest.json")
  ),
  figures_output = c(
    file.path(DIR_FIGURES, "01_embeddings_week_umap_tsne.pdf"),
    file.path(DIR_FIGURES, "01_umap_harmony_by_week.png"),
    file.path(DIR_FIGURES, "01_tsne_harmony_by_week.png"),
    file.path(DIR_FIGURES, "01_embeddings_clusters_umap_tsne.pdf")
  )
)
manifest_path <- file.path(DIR_REPORTS, "01_preprocess_harmony_manifest.json")
write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
log_msg("Saved run manifest: ", manifest_path)

log_msg("============================================================")
log_msg("Script complete: ", PIPELINE_NAME)
log_msg("Final dimensions: ", nrow(seu), " genes x ", ncol(seu), " cells")
log_msg("============================================================")
