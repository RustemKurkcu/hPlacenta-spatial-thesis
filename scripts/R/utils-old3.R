<<<<<<< HEAD
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

# Enhanced-compatibility wrappers ----------------------------------------------

log_msg_enhanced <- function(msg, logfile = NULL, verbose = TRUE) {
  if (isTRUE(verbose)) {
    log_msg(msg, logfile = logfile)
  } else if (!is.null(logfile)) {
    line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
    cat(line, "\n", file = logfile, append = TRUE)
  }
  invisible(msg)
}

check_gene_availability <- function(obj, genes, assay = NULL) {
  assay <- assay %||% Seurat::DefaultAssay(obj)
  feats <- rownames(obj[[assay]])
  present <- intersect(unique(genes), feats)
  total <- length(unique(genes))
  list(
    genes_present = present,
    genes_missing = setdiff(unique(genes), feats),
    n_present = length(present),
    total_requested = total,
    pct_present = if (total > 0) 100 * length(present) / total else 0
  )
}

add_modules_from_list <- function(obj, modules, assay = NULL, prefix = "score_", seed = 1,
                                  min_genes = 3, verbose = FALSE) {
  assay <- assay %||% Seurat::DefaultAssay(obj)
  Seurat::DefaultAssay(obj) <- assay
  obj <- safe_join_layers(obj, assay = assay)
  if (!has_data_layer(obj, assay = assay)) {
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  }
  
  for (nm in names(modules)) {
    genes <- modules[[nm]]
    score_name <- paste0(prefix, nm)
    if (isTRUE(verbose)) {
      log_msg(sprintf("Scoring module '%s' in assay '%s'", score_name, assay))
    }
    obj <- add_module_score_safe(
      obj = obj,
      genes = genes,
      score_name = score_name,
      assay = assay,
      seed = seed,
      min_genes = min_genes,
      verbose = verbose
    )
  }
  obj
}

add_modules_from_list_enhanced <- function(obj, modules, assay = NULL, prefix = "score_",
                                           seed = 1, min_genes = 3, verbose = TRUE) {
  assay <- assay %||% Seurat::DefaultAssay(obj)
  for (nm in names(modules)) {
    av <- check_gene_availability(obj, modules[[nm]], assay = assay)
    log_msg_enhanced(sprintf("  %s: %d/%d genes available (%.1f%%)",
                             nm, av$n_present, av$total_requested, av$pct_present),
                     verbose = verbose)
  }
  add_modules_from_list(obj, modules, assay = assay, prefix = prefix,
                        seed = seed, min_genes = min_genes, verbose = verbose)
}

# Spatial helpers ----------------------------------------------------------------

check_required_packages <- function(pkgs, context = NULL) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    where <- if (!is.null(context)) paste0(" for ", context) else ""
    stop("Missing required package(s)", where, ": ", paste(missing, collapse = ", "),
         ". Please install them in R (e.g., install.packages(...)).")
  }
  invisible(TRUE)
}

scatter_spatial <- function(df, x = "spatial_x_use", y = "spatial_y_use", color = "celltype",
                            title = NULL, point_size = 0.6, alpha = 0.9) {
  ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]], color = .data[[color]])) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic() +
    ggplot2::labs(title = title, x = "spatial x", y = "spatial y", color = color)
}

compute_spatial_knn_graph <- function(x, y, k = 15) {
  check_required_packages(c("RANN"), context = "compute_spatial_knn_graph")
  coords <- cbind(as.numeric(x), as.numeric(y))
  keep <- is.finite(coords[, 1]) & is.finite(coords[, 2])
  coords <- coords[keep, , drop = FALSE]
  if (nrow(coords) < 3) stop("Too few finite spatial coordinates for KNN")
  k_use <- max(1L, min(as.integer(k), nrow(coords) - 1L))
  nn <- RANN::nn2(coords, k = k_use + 1L)
  list(idx = nn$nn.idx[, -1, drop = FALSE],
       dist = nn$nn.dists[, -1, drop = FALSE],
       k = k_use)
}

compute_spatial_knn <- function(x, y, k = 15) {
  compute_spatial_knn_graph(x, y, k = k)$idx
}

neighbor_enrichment <- function(labels, knn_idx, n_perm = 200, seed = 1) {
  labels <- as.character(labels)
  lev <- sort(unique(labels))
  n <- length(labels)
  if (nrow(knn_idx) != n) stop("labels and knn_idx have incompatible sizes")
  
  count_edges <- function(lab) {
    from <- rep(lab, each = ncol(knn_idx))
    to <- lab[as.vector(knn_idx)]
    tab <- table(factor(from, levels = lev), factor(to, levels = lev))
    as.matrix(tab)
  }
  
  observed <- count_edges(labels)
  set.seed(seed)
  perms <- array(0, dim = c(length(lev), length(lev), n_perm),
                 dimnames = list(lev, lev, NULL))
  for (i in seq_len(n_perm)) {
    perms[, , i] <- count_edges(sample(labels, replace = FALSE))
  }
  expected <- apply(perms, c(1, 2), mean)
  sd_null <- apply(perms, c(1, 2), stats::sd)
  z <- (observed - expected) / (sd_null + 1e-9)
  
  list(observed = observed, expected = expected, z = z, sd_null = sd_null)
}
=======
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

# Enhanced-compatibility wrappers ----------------------------------------------

log_msg_enhanced <- function(msg, logfile = NULL, verbose = TRUE) {
  if (isTRUE(verbose)) {
    log_msg(msg, logfile = logfile)
  } else if (!is.null(logfile)) {
    line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
    cat(line, "\n", file = logfile, append = TRUE)
  }
  invisible(msg)
}

check_gene_availability <- function(obj, genes, assay = NULL) {
  assay <- assay %||% Seurat::DefaultAssay(obj)
  feats <- rownames(obj[[assay]])
  present <- intersect(unique(genes), feats)
  total <- length(unique(genes))
  list(
    genes_present = present,
    genes_missing = setdiff(unique(genes), feats),
    n_present = length(present),
    total_requested = total,
    pct_present = if (total > 0) 100 * length(present) / total else 0
  )
}

add_modules_from_list <- function(obj, modules, assay = NULL, prefix = "score_", seed = 1,
                                  min_genes = 3, verbose = FALSE) {
  assay <- assay %||% Seurat::DefaultAssay(obj)
  Seurat::DefaultAssay(obj) <- assay
  obj <- safe_join_layers(obj, assay = assay)
  if (!has_data_layer(obj, assay = assay)) {
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  }
  
  for (nm in names(modules)) {
    genes <- modules[[nm]]
    score_name <- paste0(prefix, nm)
    if (isTRUE(verbose)) {
      log_msg(sprintf("Scoring module '%s' in assay '%s'", score_name, assay))
    }
    obj <- add_module_score_safe(
      obj = obj,
      genes = genes,
      score_name = score_name,
      assay = assay,
      seed = seed,
      min_genes = min_genes,
      verbose = verbose
    )
  }
  obj
}

add_modules_from_list_enhanced <- function(obj, modules, assay = NULL, prefix = "score_",
                                           seed = 1, min_genes = 3, verbose = TRUE) {
  assay <- assay %||% Seurat::DefaultAssay(obj)
  for (nm in names(modules)) {
    av <- check_gene_availability(obj, modules[[nm]], assay = assay)
    log_msg_enhanced(sprintf("  %s: %d/%d genes available (%.1f%%)",
                             nm, av$n_present, av$total_requested, av$pct_present),
                     verbose = verbose)
  }
  add_modules_from_list(obj, modules, assay = assay, prefix = prefix,
                        seed = seed, min_genes = min_genes, verbose = verbose)
}

# Spatial helpers ----------------------------------------------------------------

check_required_packages <- function(pkgs, context = NULL) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    where <- if (!is.null(context)) paste0(" for ", context) else ""
    stop("Missing required package(s)", where, ": ", paste(missing, collapse = ", "),
         ". Please install them in R (e.g., install.packages(...)).")
  }
  invisible(TRUE)
}

scatter_spatial <- function(df, x = "spatial_x_use", y = "spatial_y_use", color = "celltype",
                            title = NULL, point_size = 0.6, alpha = 0.9) {
  ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]], color = .data[[color]])) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic() +
    ggplot2::labs(title = title, x = "spatial x", y = "spatial y", color = color)
}

compute_spatial_knn_graph <- function(x, y, k = 15) {
  check_required_packages(c("RANN"), context = "compute_spatial_knn_graph")
  coords <- cbind(as.numeric(x), as.numeric(y))
  keep <- is.finite(coords[, 1]) & is.finite(coords[, 2])
  coords <- coords[keep, , drop = FALSE]
  if (nrow(coords) < 3) stop("Too few finite spatial coordinates for KNN")
  k_use <- max(1L, min(as.integer(k), nrow(coords) - 1L))
  nn <- RANN::nn2(coords, k = k_use + 1L)
  list(idx = nn$nn.idx[, -1, drop = FALSE],
       dist = nn$nn.dists[, -1, drop = FALSE],
       k = k_use)
}

compute_spatial_knn <- function(x, y, k = 15) {
  compute_spatial_knn_graph(x, y, k = k)$idx
}

neighbor_enrichment <- function(labels, knn_idx, n_perm = 200, seed = 1) {
  labels <- as.character(labels)
  lev <- sort(unique(labels))
  n <- length(labels)
  if (nrow(knn_idx) != n) stop("labels and knn_idx have incompatible sizes")
  
  count_edges <- function(lab) {
    from <- rep(lab, each = ncol(knn_idx))
    to <- lab[as.vector(knn_idx)]
    tab <- table(factor(from, levels = lev), factor(to, levels = lev))
    as.matrix(tab)
  }
  
  observed <- count_edges(labels)
  set.seed(seed)
  perms <- array(0, dim = c(length(lev), length(lev), n_perm),
                 dimnames = list(lev, lev, NULL))
  for (i in seq_len(n_perm)) {
    perms[, , i] <- count_edges(sample(labels, replace = FALSE))
  }
  expected <- apply(perms, c(1, 2), mean)
  sd_null <- apply(perms, c(1, 2), stats::sd)
  z <- (observed - expected) / (sd_null + 1e-9)
  
  list(observed = observed, expected = expected, z = z, sd_null = sd_null)
}
>>>>>>> 54c1a5c675ca223d183141f8a75de090237ce902
