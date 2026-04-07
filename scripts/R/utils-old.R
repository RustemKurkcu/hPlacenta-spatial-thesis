# ======================================================================
# scripts/R/utils.R
# Unified utility functions -- Seurat v5.4.0 tested
# ======================================================================

suppressPackageStartupMessages({
  library(stringr)
  library(ggplot2)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

log_msg <- function(msg, logfile = NULL) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", ts, "] ", msg)
  message(line)
  if (!is.null(logfile)) {
    ensure_dir(dirname(logfile))
    cat(line, "\n", file = logfile, append = TRUE)
  }
  invisible(line)
}

pick_first_present <- function(df, candidates) {
  for (nm in candidates) {
    if (nm %in% colnames(df)) return(nm)
  }
  NULL
}

# ---- Week parsing ----
parse_week <- function(x) {
  if (is.null(x)) return(NA_real_)
  x <- as.character(x)
  w <- stringr::str_match(x, "(?i)\\bW\\s*([0-9]{1,2})\\b")[,2]
  if (all(is.na(w))) {
    w <- stringr::str_match(x, "(?i)Week\\s*([0-9]{1,2})")[,2]
  }
  as.numeric(w)
}

ensure_week_column <- function(obj, candidates) {
  md <- obj@meta.data
  if ("week" %in% colnames(md) && !all(is.na(md$week))) return(obj)
  for (nm in candidates) {
    if (nm %in% colnames(md)) {
      w <- parse_week(md[[nm]])
      if (!all(is.na(w))) { obj$week <- w; return(obj) }
    }
  }
  w <- parse_week(rownames(md))
  if (!all(is.na(w))) { obj$week <- w; return(obj) }
  obj$week <- NA_real_
  obj
}

# ---- Spatial coordinates ----
ensure_spatial_coords <- function(obj, x_candidates, y_candidates) {
  md <- obj@meta.data
  xcol <- pick_first_present(md, x_candidates)
  ycol <- pick_first_present(md, y_candidates)
  if (is.null(xcol) || is.null(ycol)) {
    obj$spatial_x_use <- NA_real_
    obj$spatial_y_use <- NA_real_
    return(obj)
  }
  obj$spatial_x_use <- as.numeric(md[[xcol]])
  obj$spatial_y_use <- as.numeric(md[[ycol]])
  obj
}

ensure_celltype <- function(obj, col_celltype) {
  if (!(col_celltype %in% colnames(obj@meta.data))) {
    stop("Expected cell-type column not found: ", col_celltype)
  }
  obj
}

# ======================================================================
# FIX #1: Rename metadata columns that clash with reductions
# ======================================================================
fix_metadata_reduction_clash <- function(obj) {
  md_cols <- colnames(obj@meta.data)
  clash_patterns <- c("^umap_[0-9]+$", "^tsne_[0-9]+$",
                      "^UMAP_[0-9]+$", "^tSNE_[0-9]+$")
  for (pat in clash_patterns) {
    hits <- grep(pat, md_cols, value = TRUE)
    for (h in hits) {
      new_name <- paste0("meta_", h)
      obj@meta.data[[new_name]] <- obj@meta.data[[h]]
      obj@meta.data[[h]] <- NULL
    }
  }
  obj
}

# ======================================================================
# FIX #2: Join split layers (Seurat v5 compatible)
# ======================================================================
safe_join_layers <- function(obj, assay = NULL) {
  assay <- assay %||% Seurat::DefaultAssay(obj)
  tryCatch({
    current_layers <- Seurat::Layers(obj[[assay]])
    split_counts <- grep("^counts\\.", current_layers, value = TRUE)
    split_data   <- grep("^data\\.",   current_layers, value = TRUE)
    if (length(split_counts) > 1 || length(split_data) > 1) {
      message("  Joining split layers in assay '", assay, "'...")
      obj[[assay]] <- Seurat::JoinLayers(obj[[assay]])
    }
  }, error = function(e) {
    message("  Note: JoinLayers skipped (", e$message, ")")
  })
  obj
}

# ---- Safe assay getter ----
get_assay_matrix <- function(obj, assay = NULL, layer = "data") {
  assay <- assay %||% Seurat::DefaultAssay(obj)
  m <- tryCatch(Seurat::GetAssayData(obj, assay = assay, layer = layer),
                error = function(e) NULL)
  if (is.null(m)) {
    m <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = layer),
                  error = function(e) NULL)
  }
  m
}

set_scalar_meta <- function(obj, key, value) {
  obj@meta.data[[key]] <- rep(value, ncol(obj))
  obj
}

# ---- Plot helpers ----
save_plot <- function(p, path, w = 8, h = 6, dpi = 300) {
  ensure_dir(dirname(path))
  ggplot2::ggsave(path, plot = p, width = w, height = h, dpi = dpi)
  invisible(path)
}

# ---- Count-type inference ----
infer_counts_type <- function(obj, assay = "RNA", layer = "counts", sample_n = 5000) {
  mat <- tryCatch(get_assay_matrix(obj, assay = assay, layer = layer),
                  error = function(e) NULL)
  if (is.null(mat)) return("unknown")
  vals <- if (inherits(mat, "dgCMatrix") || inherits(mat, "dgTMatrix")) mat@x else as.numeric(mat)
  if (length(vals) == 0) return("unknown")
  if (length(vals) > sample_n) vals <- sample(vals, sample_n)
  frac_int <- mean(abs(vals - round(vals)) < 1e-6)
  any_neg  <- any(vals < -1e-9)
  if (any_neg) return("normalized_like")
  if (frac_int >= 0.99) return("UMI_like")
  "normalized_like"
}

# ======================================================================
# FIX #3: Adaptive AddModuleScore
# ======================================================================
add_module_score_safe <- function(obj, genes, score_name,
                                  assay = NULL, layer = "data",
                                  min_genes = 3, seed = 1) {
  assay <- assay %||% Seurat::DefaultAssay(obj)
  feats <- rownames(obj[[assay]])
  found <- intersect(genes, feats)
  
  if (length(found) < min_genes) {
    warning(sprintf("Only %d/%d genes found for '%s'; skipping.",
                    length(found), length(genes), score_name))
    obj@meta.data[[score_name]] <- NA_real_
    return(obj)
  }
  
  ncells <- ncol(obj)
  ngenes <- length(found)
  
  # Adaptive nbin
  nbin <- if (ncells < 50) 3
  else if (ncells < 100) 5
  else if (ncells < 500) 10
  else if (ncells < 2000) 15
  else 24
  nbin <- min(nbin, max(3, floor(ncells / 2)))
  
  # Adaptive ctrl
  ctrl_size <- min(100, max(5, ngenes, floor(ncells / 10)))
  
  tryCatch({
    obj <- Seurat::AddModuleScore(
      obj, features = list(found), name = score_name,
      assay = assay, nbin = nbin, ctrl = ctrl_size, seed = seed
    )
    appended <- paste0(score_name, "1")
    if (appended %in% colnames(obj@meta.data)) {
      obj@meta.data[[score_name]] <- obj@meta.data[[appended]]
      obj@meta.data[[appended]] <- NULL
    }
  }, error = function(e) {
    warning(sprintf("AddModuleScore failed for '%s': %s. Setting to NA.",
                    score_name, e$message))
    obj@meta.data[[score_name]] <<- NA_real_
  })
  obj
}

# ---- Enhanced module scoring with logging ----
add_module_score_enhanced <- function(obj, genes, score_name,
                                      assay = NULL, layer = "data",
                                      min_genes = 3, seed = 1,
                                      verbose = TRUE) {
  assay <- assay %||% Seurat::DefaultAssay(obj)
  feats <- rownames(obj[[assay]])
  found <- intersect(genes, feats)
  missing <- setdiff(genes, feats)
  if (verbose && length(missing) > 0) {
    log_msg(sprintf("  %s: %d/%d genes found (%d missing)",
                    score_name, length(found), length(genes), length(missing)))
  }
  add_module_score_safe(obj, found, score_name, assay, layer, min_genes, seed)
}

# ---- Check if data layer exists (v5 compatible) ----
has_data_layer <- function(obj, assay = NULL) {
  assay <- assay %||% Seurat::DefaultAssay(obj)
  tryCatch({
    lyrs <- Seurat::Layers(obj[[assay]])
    "data" %in% lyrs
  }, error = function(e) {
    tryCatch({
      d <- Seurat::GetAssayData(obj, assay = assay, slot = "data")
      !is.null(d) && (if (inherits(d, "dgCMatrix")) length(d@x) > 0 else length(d) > 0)
    }, error = function(e2) FALSE)
  })
}

# ---- ensure_pca_umap (metadata fix + layer join) ----
ensure_pca_umap <- function(obj, assay = NULL, dims = 1:30,
                            umap_return_model = FALSE, force = FALSE) {
  assay <- assay %||% Seurat::DefaultAssay(obj)
  
  if (!force && "umap" %in% names(obj@reductions)) {
    emb <- tryCatch(Seurat::Embeddings(obj, "umap"), error = function(e) NULL)
    if (!is.null(emb) && is.numeric(emb) && nrow(emb) > 0) return(obj)
    force <- TRUE
  }
  
  obj <- fix_metadata_reduction_clash(obj)
  Seurat::DefaultAssay(obj) <- assay
  obj <- safe_join_layers(obj, assay)
  
  if (assay != "SCT") {
    if (!has_data_layer(obj, assay)) {
      obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    }
    obj <- Seurat::FindVariableFeatures(obj,nfeatures = 5000, verbose = FALSE)
    obj <- Seurat::ScaleData(obj, features = Seurat::VariableFeatures(obj), verbose = FALSE)
  }
  
  obj <- Seurat::RunPCA(obj, verbose = FALSE)
  obj <- Seurat::RunUMAP(obj, dims = dims, return.model = umap_return_model, verbose = FALSE)
  obj
}
