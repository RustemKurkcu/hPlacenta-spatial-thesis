# ======================================================================
# scripts/R/utils.R
# Unified utility functions -- Seurat v5.4.0 tested
# VariableFeatures set to 5000 by default
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
    # Seurat v5: Layers/JoinLayers are in SeuratObject (not exported by Seurat)
    current_layers <- SeuratObject::Layers(obj[[assay]])
    split_counts <- grep("^counts\\.", current_layers, value = TRUE)
    split_data   <- grep("^data\\.",   current_layers, value = TRUE)
    if (length(split_counts) > 1 || length(split_data) > 1) {
      message("  Joining split layers in assay '", assay, "'...")
      obj[[assay]] <- SeuratObject::JoinLayers(obj[[assay]])
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
  min_genes = 3, seed = 1,
  verbose = FALSE) {

  assay <- assay %||% Seurat::DefaultAssay(obj)
  feats <- rownames(obj[[assay]])
  found <- intersect(genes, feats)

  if (length(found) < min_genes) {
    warning(sprintf("Only %d/%d genes found for '%s'; skipping.",
      length(found), length(genes), score_name))
    obj@meta.data[[score_name]] <- NA_real_
    return(obj)
  }

  # ---- Choose stable AddModuleScore parameters for Seurat v5 ----
  # AddModuleScore bins features by average expression using `cut_number`.
  # If there are too few unique average-expression values (common in sparse / panel data
  # or when the 'data' layer is empty), large nbin will fail.
  total_features <- length(feats)

  # Compute average expression to estimate number of unique bins available
  avg_expr <- tryCatch({
    m <- get_assay_matrix(obj, assay = assay, layer = layer)
    Matrix::rowMeans(m)
  }, error = function(e) NULL)

  n_unique <- if (!is.null(avg_expr)) {
    length(unique(round(as.numeric(avg_expr), 6)))
  } else {
    # Conservative fallback
    max(10L, min(24L, total_features))
  }

  # Start with a gene-count based heuristic and clamp to what the data can support
  nbin <- min(24L, max(5L, floor(total_features / 100)))
  nbin <- min(nbin, max(3L, n_unique - 1L))

  # ctrl must be <= (approx bin size - 1) to avoid "sample larger than population"
  approx_bin_size <- max(1L, floor(total_features / max(1L, nbin)))
  ctrl_size <- min(100L, max(5L, floor(total_features / 50)))
  ctrl_size <- min(ctrl_size, max(1L, approx_bin_size - 1L))

  if (verbose) {
    message(sprintf("AddModuleScore('%s'): found=%d/%d, total_features=%d, nbin=%d, ctrl=%d, n_unique(avg)=%d",
      score_name, length(found), length(genes), total_features, nbin, ctrl_size, n_unique))
  }

  tryCatch({
    obj <- Seurat::AddModuleScore(
      obj,
      features = list(found),
      name = paste0(score_name, "_"),
      assay = assay,
      nbin = nbin,
      ctrl = ctrl_size,
      seed = seed
    )
    tmp <- paste0(score_name, "_1")
    if (tmp %in% colnames(obj@meta.data)) {
      obj@meta.data[[score_name]] <- obj@meta.data[[tmp]]
      obj@meta.data[[tmp]] <- NULL
    }
  }, error = function(e) {
    warning(sprintf("AddModuleScore failed for '%s' (nbin=%d, ctrl=%d): %s. Setting to NA.",
      score_name, nbin, ctrl_size, e$message))
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
    lyrs <- SeuratObject::Layers(obj[[assay]])
    "data" %in% lyrs
  }, error = function(e) {  
    tryCatch({
      # SeuratObject v5 uses `layer=` (slot is deprecated)
      d <- Seurat::GetAssayData(obj, assay = assay, layer = "data")
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
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = 5000, verbose = FALSE)
    v_feats <- Seurat::VariableFeatures(obj)
    if (length(v_feats) == 0) v_feats <- rownames(obj[[assay]])
    obj <- Seurat::ScaleData(obj, features = v_feats, verbose = FALSE)
  }
  
  obj <- Seurat::RunPCA(obj, verbose = FALSE)
  obj <- Seurat::RunUMAP(obj, dims = dims, return.model = umap_return_model, verbose = FALSE)

  # Standardize UMAP dim names to avoid clashes with metadata columns like `umap_1/umap_2`
  # (these clashes can break FeaturePlot/FetchData in Seurat v5)
  if ("umap" %in% names(obj@reductions)) {
    try({
      red <- obj@reductions$umap
      emb <- red@cell.embeddings
      emb <- as.matrix(emb)
      storage.mode(emb) <- "double"
      if (ncol(emb) >= 2) {
        colnames(emb)[1:2] <- c("UMAP_1", "UMAP_2")
      }
      red@cell.embeddings <- emb
      red@key <- "UMAP_"
      obj@reductions$umap <- red
    }, silent = TRUE)
  }
  obj
}
