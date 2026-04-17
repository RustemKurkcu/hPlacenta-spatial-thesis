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

  if (requireNamespace("future", quietly = TRUE)) {
    old_plan <- future::plan()
    options(future.globals.maxSize = as.numeric(max_size_gb) * 1024^3)
    future::plan(future::multisession, workers = workers)
    log_msg(
      "Configured future runtime: plan=multisession (workers=", workers,
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

assay_payload <- get_normalized_assay_data(seu, assay = "RNA")
data.input <- assay_payload$mat
log_msg(
  "Assay matrix extracted for CellChat from ", assay_payload$source, ": ",
  nrow(data.input), " genes x ", ncol(data.input), " cells."
)
if (!grepl("data|lognorm", assay_payload$source)) {
  log_msg(
    "Expression matrix came from counts-like layer/slot (", assay_payload$source, "). ",
    "Trying NormalizeData() on RNA assay to populate layer:data.",
    .level = "WARN"
  )
  normalized_from_seurat <- FALSE
  seu <- tryCatch(
    Seurat::NormalizeData(seu, assay = "RNA", verbose = FALSE),
    error = function(e) {
      log_msg("NormalizeData failed: ", conditionMessage(e), .level = "WARN")
      seu
    }
  )
  norm_layer <- tryCatch(
    SeuratObject::GetAssayData(seu, assay = "RNA", layer = "data"),
    error = function(e) NULL
  )
  if (!is.null(norm_layer) && nrow(norm_layer) > 0 && ncol(norm_layer) > 0) {
    data.input <- norm_layer
    normalized_from_seurat <- TRUE
    log_msg("Using RNA layer:data generated by NormalizeData().")
  }
  if (!normalized_from_seurat) {
    log_msg("Falling back to manual CPM + log1p normalization for CellChat input.", .level = "WARN")
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Package 'Matrix' is required to normalize count matrices for CellChat.")
    }
    lib_size <- Matrix::colSums(data.input)
    lib_size[lib_size == 0] <- 1
    data.input <- Matrix::t(Matrix::t(data.input) / lib_size) * 1e4
    data.input <- log1p(data.input)
  }
}
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
if (!"annotation" %in% colnames(cellchat_base@DB$interaction)) {
  cellchat_base@DB$interaction$annotation <- "Secreted Signaling"
}

future_state <- configure_future_runtime(max_size_gb = 32, workers = 10)
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
top2_vulnerable_clusters <- head(cluster_tbl$seurat_clusters, 2)
top2_celltypes <- unique(as.character(seu$celltype_plot[seu$seurat_clusters %in% top2_vulnerable_clusters]))

run_cellchat_architecture <- function(
  architecture_name,
  db_search,
  interaction_range,
  output_object,
  manifest_out
) {
  log_msg("Running architecture: ", architecture_name, " (interaction.range=", interaction_range, ")")
  cellchat <- cellchat_base
  cellchat@DB <- subset_db_fn(db_human, search = db_search)
  if (is.null(cellchat@DB) || is.null(cellchat@DB$interaction)) {
    stop("Filtered DB is missing interaction table for architecture: ", architecture_name)
  }
  if (!"annotation" %in% colnames(cellchat@DB$interaction)) {
    cellchat@DB$interaction$annotation <- "Secreted Signaling"
  }

  cellchat <- subset_fn(cellchat)
  cellchat <- over_gene_fn(cellchat)
  cellchat <- over_inter_fn(cellchat)

  comm_args <- list(
    object = cellchat,
    distance.use = TRUE,
    interaction.range = interaction_range,
    contact.dependent = FALSE
  )
  if ("scale.distance" %in% comm_formals) comm_args$scale.distance <- dynamic_scale_distance
  if ("type" %in% comm_formals) comm_args$type <- "truncatedMean"
  if ("trim" %in% comm_formals) comm_args$trim <- 0.1
  if ("raw.use" %in% comm_formals) comm_args$raw.use <- TRUE
  if ("contact.dependent.forced" %in% comm_formals) comm_args$contact.dependent.forced <- FALSE

  log_msg("Computing communication probabilities for ", architecture_name, ".")
  cellchat <- do.call(commprob_fn, comm_args)
  cellchat <- filter_fn(cellchat, min.cells = 10)
  cellchat <- pathway_fn(cellchat)
  cellchat <- aggregate_fn(cellchat)

  comm_df <- subset_comm_fn(cellchat)
  comm_vulnerable <- comm_df %>%
    dplyr::filter(source %in% top2_celltypes | target %in% top2_celltypes)
  pathway_rank <- comm_vulnerable %>%
    dplyr::group_by(pathway_name) %>%
    dplyr::summarise(vulnerable_prob = sum(prob, na.rm = TRUE), n_edges = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(vulnerable_prob), dplyr::desc(n_edges))

  bundle <- list(
    architecture = architecture_name,
    cellchat = cellchat,
    source_scored_object = input_obj,
    db_search = db_search,
    interaction_range = interaction_range,
    vulnerable_clusters = top2_vulnerable_clusters,
    vulnerable_celltypes = top2_celltypes,
    vulnerable_pathway_rank = pathway_rank,
    spatial_meta = seu@meta.data %>%
      dplyr::select(x_um, y_um, week, seurat_clusters, celltype_plot, MISI_Vulnerability, IDO1_Tolerogenic_Shield),
    run_metadata = list(
      pipeline = PIPELINE_NAME,
      version = PIPELINE_VERSION,
      run_timestamp = RUN_TIMESTAMP,
      seed = RUN_SEED
    )
  )
  saveRDS(bundle, output_object)
  log_msg("Saved ", architecture_name, " bundle: ", output_object)

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
      paste0("Dual-architecture run: ", architecture_name),
      paste0("DB search = ", paste(db_search, collapse = ", ")),
      paste0("interaction.range = ", interaction_range, " um")
    ),
    hypothesis = "Vulnerable spatial niches exhibit enhanced tolerogenic/exhaustion signaling to break down the IDO1 shield.",
    methods_blurb = "Spatial CellChat with dual-architecture modeling to separate contact-mediated vs secreted diffusion signaling.",
    thesis_aim = "Quantify and compare juxtacrine and paracrine communication programs in high-vulnerability tissue microdomains."
  )
  log_msg("Saved manifest: ", manifest_out)
}

juxtacrine_obj <- file.path(DIR_OBJECTS, "03_spatial_cellchat_juxtacrine.rds")
juxtacrine_manifest <- file.path(DIR_REPORTS, "03_spatial_cellchat_juxtacrine_manifest.json")
run_cellchat_architecture(
  architecture_name = "juxtacrine",
  db_search = c("Cell-Cell Contact", "ECM-Receptor"),
  interaction_range = 10,
  output_object = juxtacrine_obj,
  manifest_out = juxtacrine_manifest
)

paracrine_obj <- file.path(DIR_OBJECTS, "03_spatial_cellchat_paracrine.rds")
paracrine_manifest <- file.path(DIR_REPORTS, "03_spatial_cellchat_paracrine_manifest.json")
run_cellchat_architecture(
  architecture_name = "paracrine",
  db_search = "Secreted Signaling",
  interaction_range = 100,
  output_object = paracrine_obj,
  manifest_out = paracrine_manifest
)
