#!/usr/bin/env Rscript

# =============================================================================
# Script: 03_spatial_cellchat.R
# Purpose: Compute niche-aware spatial CellChat communication from scored object.
# Notes:   Compute-only script (no figures); stores CellChat object + metadata.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(jsonlite)
})

source("R/celltype_dictionary.R")

PIPELINE_NAME <- "03_spatial_cellchat"
PIPELINE_VERSION <- "1.0.0"
RUN_TIMESTAMP <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
RUN_SEED <- 42L
set.seed(RUN_SEED)

OUT_ROOT <- "output"
DIR_LOGS <- file.path(OUT_ROOT, "logs")
DIR_OBJECTS <- file.path(OUT_ROOT, "objects")
DIR_REPORTS <- file.path(OUT_ROOT, "reports")
for (d in c(DIR_LOGS, DIR_OBJECTS, DIR_REPORTS)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
LOG_FILE <- file.path(DIR_LOGS, "03_spatial_cellchat.log")

log_msg <- function(..., .level = "INFO") {
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  txt <- paste0(..., collapse = "")
  line <- paste0(stamp, " [", .level, "] ", txt)
  cat(line, "\n")
  cat(line, "\n", file = LOG_FILE, append = TRUE)
}

configure_future_runtime <- function(max_size_gb = 8, workers = 10) {
  old_plan <- NULL
  old_max <- getOption("future.globals.maxSize")
  safe_workers <- as.integer(workers)
  if (as.numeric(max_size_gb) > 10 && safe_workers > 6) {
    safe_workers <- 6L
    log_msg(
      "Safe worker cap applied: max_size_gb=", max_size_gb,
      " => workers capped to ", safe_workers, " to avoid swap pressure."
    )
  }

  if (requireNamespace("future", quietly = TRUE)) {
    old_plan <- future::plan()
    options(future.globals.maxSize = as.numeric(max_size_gb) * 1024^3)
    future::plan(future::multisession, workers = safe_workers)
    log_msg(
      "Configured future runtime: plan=multisession (workers=", safe_workers,
      "), future.globals.maxSize=", round(getOption("future.globals.maxSize") / 1024^3, 2), " GiB."
    )
  } else {
    log_msg("Package 'future' not installed; skipping future runtime configuration.", .level = "WARN")
  }

  list(old_plan = old_plan, old_max = old_max)
}

restore_future_runtime <- function(state) {
  if (!requireNamespace("future", quietly = TRUE)) return(invisible(NULL))
  if (!is.null(state$old_max)) options(future.globals.maxSize = state$old_max)
  if (!is.null(state$old_plan)) future::plan(state$old_plan)
  invisible(NULL)
}

compute_min_nonzero_distance <- function(coords_mat) {
  if (nrow(coords_mat) < 2) stop("Need at least two cells to compute spatial distances.")
  xy <- as.matrix(coords_mat[, c("x_cent", "y_cent"), drop = FALSE])

  if (requireNamespace("RANN", quietly = TRUE)) {
    nn <- RANN::nn2(data = xy, query = xy, k = 2)
    d <- nn$nn.dists[, 2]
    d <- d[is.finite(d) & d > 0]
    if (length(d) > 0) return(min(d))
  }
  if (requireNamespace("FNN", quietly = TRUE)) {
    nn <- FNN::get.knn(xy, k = 2)
    d <- nn$nn.dist[, 2]
    d <- d[is.finite(d) & d > 0]
    if (length(d) > 0) return(min(d))
  }

  # Fallback: brute force on a capped subset to avoid O(N^2) memory blowups.
  set.seed(42L)
  idx <- sample(seq_len(nrow(xy)), size = min(5000, nrow(xy)))
  dmat <- as.matrix(stats::dist(xy[idx, , drop = FALSE]))
  dvec <- dmat[upper.tri(dmat)]
  dvec <- dvec[is.finite(dvec) & dvec > 0]
  if (length(dvec) == 0) stop("Could not compute a non-zero inter-cell distance.")
  min(dvec)
}

resolve_object_path <- function(path_rds) {
  path_qs <- sub("\\.rds$", ".qs", path_rds)
  if (file.exists(path_qs)) return(path_qs)
  if (file.exists(path_rds)) return(path_rds)
  stop("Missing object: ", path_rds)
}

read_object <- function(path) {
  if (grepl("\\.qs$", path)) {
    if (!requireNamespace("qs", quietly = TRUE)) stop("Need package 'qs' to read: ", path)
    return(qs::qread(path))
  }
  readRDS(path)
}

get_normalized_assay_data <- function(seu, assay = "RNA", preferred_layers = c("data", "lognorm", "counts")) {
  is_valid_matrix <- function(x) !is.null(x) && (is.matrix(x) || inherits(x, "Matrix"))
  has_shape <- function(x) is_valid_matrix(x) && nrow(x) > 0 && ncol(x) > 0

  # SeuratObject v5 uses `layer`; older objects/scripts may still rely on slot.
  for (ly in preferred_layers) {
    mat <- tryCatch(
      SeuratObject::GetAssayData(seu, assay = assay, layer = ly),
      error = function(e) NULL
    )
    if (has_shape(mat)) return(list(mat = mat, source = paste0("layer:", ly)))
  }
  # Legacy Seurat (v4 and older) fallback.
  legacy_slots <- unique(c(preferred_layers, "data", "counts"))
  for (sl in legacy_slots) {
    mat <- tryCatch(
      SeuratObject::GetAssayData(seu, assay = assay, slot = sl),
      error = function(e) NULL
    )
    if (has_shape(mat)) return(list(mat = mat, source = paste0("slot:", sl)))
  }
  stop(
    "Unable to extract assay data from assay '", assay, "'. ",
    "Tried layers/slots: ", paste(unique(c(preferred_layers, legacy_slots)), collapse = ", "),
    ". Please confirm assay/layer names in this Seurat object."
  )
}

record_artifact_manifest <- function(
  manifest_path,
  pipeline,
  version,
  run_timestamp,
  seed,
  source_data,
  output_object,
  script_path,
  notes,
  hypothesis,
  methods_blurb,
  thesis_aim
) {
  manifest <- list(
    pipeline = pipeline,
    version = version,
    run_timestamp = run_timestamp,
    seed = seed,
    source_data = source_data,
    output_object = output_object,
    script_path = script_path,
    notes = notes,
    hypothesis = hypothesis,
    methods_blurb = methods_blurb,
    thesis_aim = thesis_aim
  )
  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
}

input_obj <- resolve_object_path(file.path(DIR_OBJECTS, "02_scored_misi_ido1.rds"))
output_obj <- file.path(DIR_OBJECTS, "03_spatial_cellchat.rds")
manifest_path <- file.path(DIR_REPORTS, "03_spatial_cellchat_manifest.json")

log_msg("Loading scored object: ", input_obj)
seu <- read_object(input_obj)
DefaultAssay(seu) <- "RNA"

seu$week <- factor(as.character(seu$week), levels = c("W7", "W8-2", "W9", "W11"))
celltype_col <- if ("predicted.celltype" %in% colnames(seu@meta.data)) "predicted.celltype" else pick_celltype_source_column(seu@meta.data)
if (is.na(celltype_col)) celltype_col <- "seurat_clusters"
seu$celltype_plot <- rename_with_true_names(as.character(seu@meta.data[[celltype_col]]))
Idents(seu) <- "celltype_plot"

if (!all(c("x_um", "y_um") %in% colnames(seu@meta.data))) {
  stop("Spatial columns x_um/y_um are required for spatial CellChat.")
}

coords <- as.matrix(seu@meta.data[, c("x_um", "y_um"), drop = FALSE])
colnames(coords) <- c("x_cent", "y_cent")
rownames(coords) <- colnames(seu)

meta <- seu@meta.data %>%
  dplyr::mutate(group = as.character(celltype_plot))
rownames(meta) <- colnames(seu)

# --- Bypass Seurat v5 NormalizeData bug and prefer SCT ---
if ("SCT" %in% Assays(seu)) {
  log_msg("Detected SCT assay. Using SCT normalized data for CellChat.")
  DefaultAssay(seu) <- "SCT"
  assay_payload <- get_normalized_assay_data(seu, assay = "SCT", preferred_layers = c("data"))
  data.input <- assay_payload$mat
} else {
  DefaultAssay(seu) <- "RNA"
  assay_payload <- get_normalized_assay_data(seu, assay = "RNA")
  data.input <- assay_payload$mat

  if (!grepl("data|lognorm", assay_payload$source)) {
    log_msg("RNA layer:data is empty. Applying manual CPM + log1p normalization directly.", .level = "WARN")
    # Bypass the buggy Seurat::NormalizeData call entirely
    if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package required.")
    lib_size <- Matrix::colSums(data.input)
    lib_size[lib_size == 0] <- 1
    data.input <- Matrix::t(Matrix::t(data.input) / lib_size) * 1e4
    data.input <- log1p(data.input)
  }
}
log_msg("Final CellChat input matrix: ", nrow(data.input), " genes x ", ncol(data.input), " cells.")
scale.factors <- list(spot = 1, spot.diameter = 1)
spatial.factors <- list(ratio = 1, tol = 0)

using_spatial_cellchat <- requireNamespace("SpatialCellChat", quietly = TRUE)
db_human <- tryCatch(
  SpatialCellChat::CellChatDB.human,
  error = function(e) CellChat::CellChatDB.human
)

if (using_spatial_cellchat) {
  log_msg("Detected SpatialCellChat package; using SpatialCellChat API (CellChat v3).")
  create_fn <- if ("createSpatialCellChat" %in% getNamespaceExports("SpatialCellChat")) {
    SpatialCellChat::createSpatialCellChat
  } else {
    stop("SpatialCellChat is installed but `createSpatialCellChat` is not exported.")
  }
  subset_fn <- SpatialCellChat::subsetData
  over_gene_fn <- SpatialCellChat::identifyOverExpressedGenes
  over_inter_fn <- SpatialCellChat::identifyOverExpressedInteractions
  commprob_fn <- SpatialCellChat::computeCommunProb
  filter_fn <- SpatialCellChat::filterCommunication
  pathway_fn <- SpatialCellChat::computeCommunProbPathway
  aggregate_fn <- SpatialCellChat::aggregateNet
  subset_comm_fn <- SpatialCellChat::subsetCommunication
} else {
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop("Neither 'SpatialCellChat' nor 'CellChat' is installed.")
  }
  log_msg("Using CellChat package (v1/v2 API). Install SpatialCellChat for CellChat v3 features.", .level = "WARN")
  create_fn <- CellChat::createCellChat
  subset_fn <- CellChat::subsetData
  over_gene_fn <- CellChat::identifyOverExpressedGenes
  over_inter_fn <- CellChat::identifyOverExpressedInteractions
  commprob_fn <- CellChat::computeCommunProb
  filter_fn <- CellChat::filterCommunication
  pathway_fn <- CellChat::computeCommunProbPathway
  aggregate_fn <- CellChat::aggregateNet
  subset_comm_fn <- CellChat::subsetCommunication
}

create_args <- list(
  object = data.input,
  meta = meta,
  group.by = "group",
  datatype = "spatial",
  coordinates = coords
)
create_formals <- names(formals(create_fn))
if ("spatial.factors" %in% create_formals) {
  create_args$spatial.factors <- spatial.factors
} else if ("scale.factors" %in% create_formals) {
  create_args$scale.factors <- scale.factors
}

log_msg("Creating spatial CellChat object using physical coordinates.")
cellchat_base <- do.call(create_fn, create_args)
cellchat_base@DB <- db_human
if (is.null(cellchat_base@DB) || is.null(cellchat_base@DB$interaction)) {
  stop(
    "CellChat DB is missing. Please ensure CellChat ligand-receptor DB is available ",
    "(e.g., install/load CellChat package data)."
  )
}

future_state <- configure_future_runtime(max_size_gb = 80, workers = 10)
on.exit(restore_future_runtime(future_state), add = TRUE)

subset_db_fn <- if (using_spatial_cellchat) SpatialCellChat::subsetDB else CellChat::subsetDB
comm_formals <- names(formals(commprob_fn))
min_dist <- compute_min_nonzero_distance(coords)
dynamic_scale_distance <- 1 / min_dist
log_msg(
  "Dynamic scale.distance computed from min non-zero physical distance: ",
  signif(min_dist, 6), " um -> scale.distance=", signif(dynamic_scale_distance, 6), "."
)

cluster_tbl <- seu@meta.data %>%
  dplyr::mutate(seurat_clusters = as.character(seurat_clusters)) %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::summarise(mean_misi = mean(MISI_Vulnerability, na.rm = TRUE), n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n >= 20) %>%
  dplyr::arrange(dplyr::desc(mean_misi), dplyr::desc(n))
top3_vulnerable_clusters <- head(cluster_tbl$seurat_clusters, 3)
top3_celltypes <- unique(as.character(seu$celltype_plot[seu$seurat_clusters %in% top3_vulnerable_clusters]))

run_cellchat_architecture <- function(architecture_name, db_search, interaction_range, output_object, manifest_out, is_juxtacrine = FALSE, pathway_focus = NULL) {
  log_msg("Running architecture: ", architecture_name, " (interaction.range=", interaction_range, ")")
  cellchat <- cellchat_base
  cellchat@DB <- subset_db_fn(db_human, search = db_search)

  if (!is.null(pathway_focus)) {
    log_msg("  FAST-TRACK ENABLED: Filtering to target pathways.")
    cellchat@DB$interaction <- cellchat@DB$interaction[cellchat@DB$interaction$pathway_name %in% pathway_focus, ]
  }
  if (nrow(cellchat@DB$interaction) == 0) stop("Filtered DB is empty.")

  # THE BIOPHYSICAL WORKAROUND: Force continuous coordinate math by masking the annotation
  cellchat@DB$interaction$annotation <- "Secreted Signaling"

  cellchat <- subset_fn(cellchat)
  cellchat <- over_gene_fn(cellchat)
  cellchat <- over_inter_fn(cellchat)

  # THE MATRIX FIX: contact.dependent=FALSE prevents dimension collisions.
  # (raw.use is intentionally omitted to allow matching subsetted matrices)
  comm_args <- list(
    object = cellchat,
    distance.use = TRUE,
    interaction.range = interaction_range,
    contact.dependent = FALSE
  )

  # API Patch: Only add legacy arguments if the specific API version accepts them
  if ("scale.distance" %in% comm_formals) comm_args$scale.distance <- dynamic_scale_distance
  if ("type" %in% comm_formals) comm_args$type <- "truncatedMean"
  if ("trim" %in% comm_formals) comm_args$trim <- 0

  log_msg("Computing probabilities...")
  cellchat <- do.call(commprob_fn, comm_args)

  tryCatch({
    cellchat <- filter_fn(cellchat, min.cells = 4)
    cellchat <- pathway_fn(cellchat)
    cellchat <- aggregate_fn(cellchat)
  }, error = function(e) log_msg("  WARN during aggregation: ", conditionMessage(e)))

  bundle <- list(
    architecture = architecture_name,
    cellchat = cellchat,
    source_scored_object = input_obj,
    db_search = db_search,
    interaction_range = interaction_range,
    vulnerable_clusters = top3_vulnerable_clusters,
    vulnerable_celltypes = top3_celltypes,
    pathway_focus = pathway_focus
  )
  saveRDS(bundle, output_object)
  log_msg("Saved bundle: ", output_object)

  record_artifact_manifest(
    manifest_path = manifest_out,
    pipeline = PIPELINE_NAME,
    version = PIPELINE_VERSION,
    run_timestamp = RUN_TIMESTAMP,
    seed = RUN_SEED,
    source_data = input_obj,
    output_object = output_object,
    script_path = "scripts/01_active_pipeline/03_spatial_cellchat.R",
    notes = c(
      paste0("Architecture = ", architecture_name),
      paste0("DB search = ", paste(db_search, collapse = ", ")),
      paste0("interaction.range = ", interaction_range, " um"),
      paste0("pathway_focus = ", ifelse(is.null(pathway_focus), "FULL", paste(pathway_focus, collapse = ",")))
    ),
    hypothesis = "Vulnerable spatial niches exhibit enhanced tolerogenic/exhaustion signaling to break down the IDO1 shield.",
    methods_blurb = "Dual-architecture SpatialCellChat run with fast-track targeted pathways and full-track comprehensive signaling.",
    thesis_aim = "Disentangle local synaptic signaling from longer-range cytokine diffusion across vulnerability niches."
  )
}

thesis_pathways <- c("MMP", "WNT", "TGFb", "CXCL", "CCL", "SPP1", "NOTCH", "FN1", "MHC-I", "MHC-II", "CD45", "TIGIT", "PD-L1")

# 1. FAST-TRACK
run_cellchat_architecture("juxtacrine_fast", c("Cell-Cell Contact", "ECM-Receptor"), 10, file.path(DIR_OBJECTS, "03_spatial_juxtacrine_fast.rds"), file.path(DIR_REPORTS, "03_manifest_juxt_fast.json"), TRUE, thesis_pathways)
run_cellchat_architecture("paracrine_fast", "Secreted Signaling", 100, file.path(DIR_OBJECTS, "03_spatial_paracrine_fast.rds"), file.path(DIR_REPORTS, "03_manifest_para_fast.json"), FALSE, thesis_pathways)

# 2. FULL-TRACK
run_cellchat_architecture("juxtacrine_full", c("Cell-Cell Contact", "ECM-Receptor"), 10, file.path(DIR_OBJECTS, "03_spatial_juxtacrine_full.rds"), file.path(DIR_REPORTS, "03_manifest_juxt_full.json"), TRUE, NULL)
run_cellchat_architecture("paracrine_full", "Secreted Signaling", 100, file.path(DIR_OBJECTS, "03_spatial_paracrine_full.rds"), file.path(DIR_REPORTS, "03_manifest_para_full.json"), FALSE, NULL)
