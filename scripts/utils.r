# Core utility fallbacks for pipeline scripts.
# This file is intentionally lightweight so updated script bundles can run
# even if a larger historical utils.R is missing pieces.

`%||%` <- function(a, b) if (!is.null(a)) a else b

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

log_msg <- function(msg, logfile = NULL) {
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  cat(line, "\n")
  if (!is.null(logfile)) cat(line, "\n", file = logfile, append = TRUE)
  invisible(line)
}

pick_first_present <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

ensure_week_column <- function(obj, candidates = c("week", "gestational_week", "donor_id", "biosample_id", "sample_id")) {
  if (!inherits(obj, "Seurat")) return(obj)
  md <- obj@meta.data
  if ("week" %in% colnames(md) && any(!is.na(md$week))) return(obj)
  
  src <- pick_first_present(md, candidates)
  if (is.null(src)) {
    obj$week <- NA_real_
    return(obj)
  }
  
  vals <- as.character(md[[src]])
  wk <- suppressWarnings(as.numeric(vals))
  if (all(is.na(wk))) {
    wk <- suppressWarnings(as.numeric(sub(".*?(\\d+).*", "\\1", vals)))
  }
  obj$week <- wk
  obj
}

ensure_spatial_coords <- function(obj, x_candidates, y_candidates) {
  if (!inherits(obj, "Seurat")) return(obj)
  if (!inherits(obj, "Seurat")) return(obj)
  md <- obj@meta.data
  xcol <- pick_first_present(md, x_candidates)
  ycol <- pick_first_present(md, y_candidates)
  obj$spatial_x_use <- if (!is.null(xcol)) suppressWarnings(as.numeric(md[[xcol]])) else NA_real_
  obj$spatial_y_use <- if (!is.null(ycol)) suppressWarnings(as.numeric(md[[ycol]])) else NA_real_
  obj
}

print_seurat_diagnostic <- function(obj, label = "Seurat object") {
  cat(strrep("=", 70), "\n")
  cat(label, "Diagnostic\n")
  cat(strrep("=", 70), "\n")
  cat("Cells:", ncol(obj), "\n")
  cat("Assays:", paste(names(obj@assays), collapse = ", "), "\n")
  cat("Reductions:", paste(names(obj@reductions), collapse = ", "), "\n")
  cat("Metadata columns:", ncol(obj@meta.data), "\n")
  if ("week" %in% colnames(obj@meta.data)) {
    w <- sort(unique(stats::na.omit(obj@meta.data$week)))
    if (length(w) > 0) cat("Weeks:", paste(w, collapse = ", "), "\n")
  }
  cat(strrep("=", 70), "\n\n")
  invisible(NULL)
}

save_plot <- function(plot_obj, path, w = 8, h = 6, dpi = 300) {
  ggplot2::ggsave(filename = path, plot = plot_obj, width = w, height = h, dpi = dpi)
}
save_plot <- function(plot_obj, path, w = 8, h = 6, dpi = 300) {
  ggplot2::ggsave(filename = path, plot = plot_obj, width = w, height = h, dpi = dpi)
}

# Seurat v5-compatible helpers -------------------------------------------------

safe_join_layers <- function(obj, assay = NULL) {
  assay <- assay %||% Seurat::DefaultAssay(obj)
  tryCatch({
    current_layers <- SeuratObject::Layers(obj[[assay]])
    split_counts <- grep("^counts\\.", current_layers, value = TRUE)
    split_data <- grep("^data\\.", current_layers, value = TRUE)
    if (length(split_counts) > 1 || length(split_data) > 1) {
      message("  Joining split layers in assay '", assay, "'...")
      obj[[assay]] <- SeuratObject::JoinLayers(obj[[assay]])
    }
  }, error = function(e) {
    message("  Note: JoinLayers skipped (", e$message, ")")
  })
  obj
}

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

has_data_layer <- function(obj, assay = NULL) {
  assay <- assay %||% Seurat::DefaultAssay(obj)
  tryCatch({
    lyrs <- SeuratObject::Layers(obj[[assay]])
    "data" %in% lyrs
  }, error = function(e) {
    tryCatch({
      d <- Seurat::GetAssayData(obj, assay = assay, layer = "data")
      !is.null(d) && (if (inherits(d, "dgCMatrix")) length(d@x) > 0 else length(d) > 0)
    }, error = function(e2) FALSE)
  })
}

# STARmap helper: prefer imputed assay (more genes), otherwise best available
select_starmap_assay <- function(obj, prefer_imputed = TRUE) {
  stopifnot(inherits(obj, "Seurat"))
  assay_names <- names(obj@assays)
  if (length(assay_names) == 0) stop("No assays found on STARmap object")
  
  if (prefer_imputed) {
    imputed_hits <- grep("imput|magic|alra|denois", assay_names,
                         ignore.case = TRUE, value = TRUE)
    if (length(imputed_hits) > 0) return(imputed_hits[[1]])
  }
  
  if ("SCT" %in% assay_names) return("SCT")
  if ("RNA" %in% assay_names) return("RNA")
  
  # Fallback: select assay with the largest feature space
  n_feats <- vapply(assay_names, function(a) nrow(obj[[a]]), numeric(1))
  assay_names[[which.max(n_feats)]]
}

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
  
  total_features <- length(feats)
  avg_expr <- tryCatch({
    m <- get_assay_matrix(obj, assay = assay, layer = layer)
    Matrix::rowMeans(m)
  }, error = function(e) NULL)
  
  n_unique <- if (!is.null(avg_expr)) {
    length(unique(round(as.numeric(avg_expr), 6)))
  } else {
    max(10L, min(24L, total_features))
  }
  
  nbin_init <- min(24L, max(5L, floor(total_features / 100)))
  nbin_init <- min(nbin_init, max(3L, n_unique - 1L))
  nbin_try <- unique(c(nbin_init, 12L, 8L, 6L, 4L, 3L))
  nbin_try <- nbin_try[nbin_try >= 3L]
  
  score_ok <- FALSE
  last_err <- NULL
  for (nbin in nbin_try) {
    approx_bin_size <- max(1L, floor(total_features / max(1L, nbin)))
    ctrl_size <- min(100L, max(5L, floor(total_features / 50)))
    ctrl_size <- min(ctrl_size, max(1L, approx_bin_size - 1L))
    
    if (verbose) {
      message(sprintf("AddModuleScore('%s'): found=%d/%d, nbin=%d, ctrl=%d",
                      score_name, length(found), length(genes), nbin, ctrl_size))
    }
    
    scored <- tryCatch({
      obj <- Seurat::AddModuleScore(
        object = obj,
        features = list(found),
        name = paste0(score_name, "_"),
        assay = assay,
        nbin = nbin,
        ctrl = ctrl_size,
        seed = seed
      )
      TRUE
    }, error = function(e) {
      last_err <<- e$message
      FALSE
    })
    
    if (scored) {
      tmp <- paste0(score_name, "_1")
      if (tmp %in% colnames(obj@meta.data)) {
        obj@meta.data[[score_name]] <- obj@meta.data[[tmp]]
        obj@meta.data[[tmp]] <- NULL
      }
      score_ok <- TRUE
      break
    }
  }
  
  if (!score_ok) {
    warning(sprintf("AddModuleScore failed for '%s' after adaptive retries: %s. Setting to NA.",
                    score_name, last_err %||% "unknown error"))
    obj@meta.data[[score_name]] <- NA_real_
  }
  
  obj
}
