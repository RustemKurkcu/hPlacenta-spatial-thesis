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
# Version: 1.1.1
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
  "tibble",
  "future"
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
  library(future)
})


record_artifact_manifest <- function(
  manifest_path,
  pipeline,
  version,
  run_timestamp,
  seed,
  weeks_requested,
  weeks_observed,
  n_cells_final,
  n_genes_final,
  qc,
  sct,
  harmony,
  umap,
  tsne,
  clustering,
  object_output,
  tables_output,
  reports_output,
  figures_output,
  source_data = NULL,
  compute_script = NULL,
  plotting_script = NULL
) {
  manifest <- list(
    pipeline = pipeline,
    version = version,
    run_timestamp = run_timestamp,
    seed = seed,
    weeks_requested = weeks_requested,
    weeks_observed = weeks_observed,
    n_cells_final = n_cells_final,
    n_genes_final = n_genes_final,
    qc = qc,
    sct = sct,
    harmony = harmony,
    umap = umap,
    tsne = tsne,
    clustering = clustering,
    object_output = object_output,
    tables_output = tables_output,
    reports_output = reports_output,
    figures_output = figures_output,
    source_data = source_data,
    compute_script = compute_script,
    plotting_script = plotting_script
  )
  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
  invisible(manifest)
}

# =============================================================================
# A) Header + Run Manifest
# =============================================================================

PIPELINE_NAME <- "01_preprocess_harmony_embeddings"
PIPELINE_VERSION <- "1.1.1"
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
DIR_TROUBLE <- file.path(OUT_ROOT, "troubleshooting", "01_preprocess_harmony_embeddings")

for (d in c(DIR_LOGS, DIR_OBJECTS, DIR_TABLES, DIR_FIGURES, DIR_REPORTS, DIR_TROUBLE)) {
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
  io = list(
    vroom_connection_size = 10485760L, # 10 MB; helps very long header/line STARmap CSVs
    prefer_qs_checkpoints = TRUE
  ),

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
    max_features = Inf,
    min_counts = 300L,
    max_counts = Inf,
    max_percent_mt = 20
  ),

  # SCTransform settings
  sct = list(
    vst_flavor = "v2",
    variable_features_n = 3000L,
    regress_vars = c("percent.mt"),
    conserve_memory = TRUE,
    return_only_var_genes = FALSE,
    future_maxsize_gb = 8,
    future_plan = "sequential"
  ),

  # Resume/checkpoint controls
  resume = list(
    use_checkpoints = TRUE,
    force_recompute = FALSE,
    resume_from = "week_qc" # one of: week_qc, merged, sct, pca
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
Sys.setenv(VROOM_CONNECTION_SIZE = as.character(cfg$io$vroom_connection_size))
log_msg("VROOM_CONNECTION_SIZE set to ", cfg$io$vroom_connection_size, " bytes")
options(future.globals.maxSize = as.numeric(cfg$sct$future_maxsize_gb) * 1024^3)
if (identical(cfg$sct$future_plan, "sequential")) {
  future::plan(future::sequential)
} else if (identical(cfg$sct$future_plan, "multisession")) {
  future::plan(future::multisession)
} else {
  future::plan(future::sequential)
  log_msg("Unknown future plan '", cfg$sct$future_plan, "'. Falling back to sequential.", .level = "WARN")
}
log_msg("future.globals.maxSize set to ", cfg$sct$future_maxsize_gb, " GiB")
log_msg("future::plan set to ", cfg$sct$future_plan)
log_msg("Resume mode: ", cfg$resume$resume_from, " (use_checkpoints=", cfg$resume$use_checkpoints, ")")

# Prefer qs for large object checkpoints when available, with .rds fallback for portability.
qs_available <- isTRUE(cfg$io$prefer_qs_checkpoints) && requireNamespace("qs", quietly = TRUE)
if (isTRUE(cfg$io$prefer_qs_checkpoints) && !qs_available) {
  log_msg("Package 'qs' not installed; checkpoint IO will use .rds fallback.", .level = "WARN")
}

ckpt_qs_path <- function(path_rds) sub("\\.rds$", ".qs", path_rds)

ckpt_exists <- function(path_rds) {
  path_qs <- ckpt_qs_path(path_rds)
  file.exists(path_qs) || file.exists(path_rds)
}

ckpt_read <- function(path_rds) {
  path_qs <- ckpt_qs_path(path_rds)
  if (qs_available && file.exists(path_qs)) {
    log_msg("Loading qs checkpoint: ", path_qs)
    return(qs::qread(path_qs))
  }
  if (file.exists(path_rds)) {
    log_msg("Loading rds checkpoint: ", path_rds)
    return(readRDS(path_rds))
  }
  if (file.exists(path_qs)) {
    log_msg("Loading qs checkpoint: ", path_qs)
    return(qs::qread(path_qs))
  }
  stop("Checkpoint not found (.rds or .qs): ", path_rds)
}

ckpt_write <- function(object, path_rds) {
  if (qs_available) {
    path_qs <- ckpt_qs_path(path_rds)
    qs::qsave(object, path_qs, preset = "high")
    return(path_qs)
  }
  saveRDS(object, path_rds)
  path_rds
}

assert_no_spot_level_metadata <- function(seu_obj, context_label = "object") {
  forbidden <- intersect(c("gene", "spot_id"), colnames(seu_obj@meta.data))
  if (length(forbidden) > 0) {
    stop(
      "Detected spot-level metadata columns in ", context_label, ": ",
      paste(forbidden, collapse = ", "),
      ". Remove stale checkpoints and rerun from resume_from='week_qc' with force_recompute=TRUE."
    )
  }
}

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

is_git_lfs_pointer_file <- function(path) {
  if (is.na(path) || !file.exists(path)) return(FALSE)
  lines <- tryCatch(readLines(path, n = 3, warn = FALSE), error = function(e) character(0))
  if (length(lines) == 0) return(FALSE)
  any(grepl("^version https://git-lfs.github.com/spec/v1$", lines))
}

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
      if (length(hits) > 0) {
        non_pointer_hits <- hits[!vapply(hits, is_git_lfs_pointer_file, logical(1))]
        if (length(non_pointer_hits) > 0) return(non_pointer_hits[[1]])
        return(hits[[1]])
      }
    }
  }
  NA_character_
}

detect_delimiter <- function(path) {
  lines <- readLines(path, n = 100, warn = FALSE, encoding = "UTF-8")
  if (length(lines) == 0) return(list(delim = ",", skip = 0L))

  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) return(list(delim = ",", skip = 0L))

  header_line <- lines[[1]]
  cands <- c("," = ",", "tab" = "	", ";" = ";", "|" = "|")
  scores <- vapply(cands, function(d) stringr::str_count(header_line, fixed(d)), numeric(1))

  if (all(scores <= 0)) {
    return(list(delim = ",", skip = 0L))
  }

  best_delim <- unname(cands[[which.max(scores)]])
  list(delim = best_delim, skip = 0L)
}

read_csv_flexible <- function(path) {
  if (is_git_lfs_pointer_file(path)) {
    stop(
      "Input file is a Git LFS pointer, not real CSV content: ", path,
      "\nRun `git lfs pull` in the repository (or replace pointer with actual file) and retry."
    )
  }

  probe <- detect_delimiter(path)
  delim_guess <- probe$delim
  skip_guess <- probe$skip

  safe_read_delim <- function(file, delim, skip) {
    tryCatch(
      readr::read_delim(
        file = file,
        delim = delim,
        skip = skip,
        show_col_types = FALSE,
        progress = FALSE,
        guess_max = 10000,
        name_repair = "minimal"
      ),
      error = function(e) {
        msg <- conditionMessage(e)
        if (grepl("connection buffer|VROOM_CONNECTION_SIZE", msg, ignore.case = TRUE)) {
          Sys.setenv(VROOM_CONNECTION_SIZE = as.character(52428800L)) # 50 MB retry
          return(readr::read_delim(
            file = file,
            delim = delim,
            skip = skip,
            show_col_types = FALSE,
            progress = FALSE,
            guess_max = 10000
          ))
        }
        stop(e)
      }
    )
  }

  df <- safe_read_delim(
    file = path,
    delim = delim_guess,
    skip = skip_guess
  )

  # Fallback attempts if file parsed into one column
  if (ncol(df) < 2) {
    for (sk in c(0L, 1L, 2L)) {
      for (d in c(",", "	", ";", "|")) {
        trial <- safe_read_delim(file = path, delim = d, skip = sk)
        if (ncol(trial) > ncol(df)) {
          df <- trial
          delim_guess <- d
          skip_guess <- sk
        }
        if (ncol(df) >= 2) break
      }
      if (ncol(df) >= 2) break
    }
  }

  # Last fallback: base reader (handles some odd quoting/encoding edge cases)
  if (ncol(df) < 2) {
    trial <- tryCatch(
      utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, skip = skip_guess),
      error = function(e) NULL
    )
    if (!is.null(trial) && ncol(trial) > ncol(df)) {
      df <- tibble::as_tibble(trial, .name_repair = "minimal")
    }
  }

  attr(df, "delim_used") <- delim_guess
  attr(df, "skip_used") <- skip_guess
  df
}

write_troubleshooting_bundle <- function(week_id, expr_path, spots_path, expr_df = NULL, spots_df = NULL, err_msg = NULL) {
  wk_dir <- file.path(DIR_TROUBLE, paste0("week_", week_id))
  dir.create(wk_dir, recursive = TRUE, showWarnings = FALSE)

  diag_lines <- c(
    paste0("timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("week: ", week_id),
    paste0("expression_path: ", expr_path),
    paste0("spots_path: ", spots_path),
    paste0("expression_file_size_bytes: ", ifelse(file.exists(expr_path), file.info(expr_path)$size, NA)),
    paste0("error: ", ifelse(is.null(err_msg), "none", err_msg))
  )

  # Keep diagnostics compact (<100MB by design: only heads/previews)
  if (!is.null(expr_df)) {
    diag_lines <- c(diag_lines, paste0("expr_dims: ", nrow(expr_df), " x ", ncol(expr_df)))
    diag_lines <- c(diag_lines,
                    paste0("expr_delim_used: ", attr(expr_df, "delim_used")),
                    paste0("expr_skip_used: ", attr(expr_df, "skip_used")))
    utils::write.csv(utils::head(expr_df[, seq_len(min(50, ncol(expr_df))), drop = FALSE], 200),
                     file = file.path(wk_dir, "expression_preview_head200x50.csv"),
                     row.names = FALSE)
    first_lines <- readLines(expr_path, n = 50, warn = FALSE)
    writeLines(utils::head(first_lines, 20),
               con = file.path(wk_dir, "expression_raw_first20lines.txt"))
    # Binary preview for cases where text read appears empty
    con <- file(expr_path, "rb")
    on.exit(close(con), add = TRUE)
    raw_head <- readBin(con, what = "raw", n = 256)
    writeLines(paste(raw_head, collapse = " "), con = file.path(wk_dir, "expression_raw_first256bytes_hex.txt"))
  }
  if (!is.null(spots_df)) {
    diag_lines <- c(diag_lines, paste0("spots_dims: ", nrow(spots_df), " x ", ncol(spots_df)))
    utils::write.csv(utils::head(spots_df, 500),
                     file = file.path(wk_dir, "spots_preview_head500.csv"),
                     row.names = FALSE)
  }
  writeLines(diag_lines, con = file.path(wk_dir, "diagnostics.txt"))
}

standardize_spots_metadata <- function(md_raw, week_id, barcodes, cell_md_raw = NULL) {
  resolve_col <- function(df, candidates) {
    if (is.null(df) || ncol(df) == 0) return(NULL)
    nm <- names(df)
    nm_norm <- tolower(gsub("[^a-z0-9]", "", nm))
    cand_norm <- tolower(gsub("[^a-z0-9]", "", candidates))
    hit <- which(nm_norm %in% cand_norm)
    if (length(hit) == 0) return(NULL)
    nm[[hit[[1]]]]
  }

  make_harmonized_md <- function(df, source_label) {
    barcode_col <- resolve_col(df, c("barcode", "barcodes", "cell", "cell_id", "spot_id", "id", "X"))
    x_col <- resolve_col(df, c("x", "x_um", "spatial_x", "X_coord", "xcoord", "x_pixel", "xpixel"))
    y_col <- resolve_col(df, c("y", "y_um", "spatial_y", "Y_coord", "ycoord", "y_pixel", "ypixel"))

    if (is.null(barcode_col) || is.null(x_col) || is.null(y_col)) return(NULL)

    out <- df %>%
      dplyr::mutate(
        barcode = as.character(.data[[barcode_col]]),
        x_um = suppressWarnings(as.numeric(.data[[x_col]])),
        y_um = suppressWarnings(as.numeric(.data[[y_col]])),
        week = week_id,
        sample_id = week_id,
        metadata_source = source_label
      )
    out$barcode <- trimws(out$barcode)
    out <- out %>% dplyr::filter(!is.na(barcode), barcode != "") %>% dplyr::distinct(barcode, .keep_all = TRUE)
    out
  }

  md_primary <- make_harmonized_md(md_raw, "spots_metadata")
  md_cell <- make_harmonized_md(cell_md_raw, "cell_metadata")

  # Prefer cell-level metadata when available; spots metadata is transcript-level and can carry per-spot fields.
  out <- md_cell
  if (is.null(out)) {
    out <- md_primary
  } else {
    log_msg("  Using cell metadata for harmonization in week ", week_id, ".")
  }

  if (is.null(out)) {
    spot_cols <- if (!is.null(md_raw)) paste(names(md_raw), collapse = ", ") else "<none>"
    cell_cols <- if (!is.null(cell_md_raw)) paste(names(cell_md_raw), collapse = ", ") else "<none>"
    stop(
      "Could not resolve barcode/x/y columns for week ", week_id,
      ". Checked spots metadata and cell metadata.",
      "\nspots columns: ", spot_cols,
      "\ncell columns: ", cell_cols
    )
  }

  if (length(barcodes) > 0) {
    out <- out %>% dplyr::filter(barcode %in% barcodes)
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
    write_troubleshooting_bundle(
      week_id = week_id,
      expr_path = expr_path,
      spots_path = "N/A",
      expr_df = expr,
      spots_df = spots_md,
      err_msg = "Expression parsed with <2 columns; likely delimiter/format issue."
    )
    stop("Expression CSV has too few columns for week ", week_id, ": ", expr_path)
  }

  first_name <- names(expr)[1]
  if (tolower(first_name) %in% c("x", "gene", "gene_name", "genes", "feature", "feature_id", "")) {
    expr <- expr %>% rename(gene = 1)
  } else {
    expr <- expr %>% rename(gene = 1)
  }
  expr$gene <- as.character(expr$gene)
  expr$gene <- trimws(expr$gene)
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

  keep_barcodes <- intersect(spot_barcodes, colnames(mat))
  if (length(keep_barcodes) == 0) {
    stop("No overlapping barcodes between expression and spots metadata for week ", week_id)
  }
  mat <- mat[, keep_barcodes, drop = FALSE]
  mat <- as(Matrix(mat, sparse = TRUE), "dgCMatrix")
  mat
}

coerce_expression_without_metadata <- function(expr_path, week_id) {
  expr <- read_csv_flexible(expr_path)
  if (ncol(expr) < 2) {
    stop("Expression CSV has too few columns for week ", week_id, ": ", expr_path)
  }

  expr <- expr %>% dplyr::rename(id_col = 1)
  expr$id_col <- trimws(as.character(expr$id_col))
  expr <- expr %>% dplyr::filter(!is.na(id_col), id_col != "") %>% dplyr::distinct(id_col, .keep_all = TRUE)

  candidate_cols <- setdiff(colnames(expr), "id_col")
  id_vals <- expr$id_col

  barcode_like <- mean(grepl("^(W[0-9]|[A-Za-z]+_tile_|cell|spot)", id_vals), na.rm = TRUE)

  if (is.finite(barcode_like) && barcode_like >= 0.05) {
    mat_raw <- as.matrix(expr[, candidate_cols, drop = FALSE])
    suppressWarnings(storage.mode(mat_raw) <- "numeric")
    rownames(mat_raw) <- expr$id_col
    mat_raw <- t(mat_raw)
    colnames(mat_raw) <- make.unique(colnames(mat_raw))
    rownames(mat_raw) <- make.unique(rownames(mat_raw))
    mat <- mat_raw
  } else {
    mat_raw <- as.matrix(expr[, candidate_cols, drop = FALSE])
    suppressWarnings(storage.mode(mat_raw) <- "numeric")
    rownames(mat_raw) <- expr$id_col
    colnames(mat_raw) <- candidate_cols
    mat <- mat_raw
  }

  mat <- as(Matrix(mat, sparse = TRUE), "dgCMatrix")
  md <- tibble::tibble(
    barcode = colnames(mat),
    x_um = NA_real_,
    y_um = NA_real_,
    week = week_id,
    sample_id = week_id,
    metadata_source = "placeholder_from_expression"
  )

  list(counts = mat, metadata = md)
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
  if (!is.na(wi$expression) && is_git_lfs_pointer_file(wi$expression)) log_msg("  WARN expression is Git LFS pointer: ", wi$expression, .level = "WARN")
  if (!is.na(wi$spots_metadata) && is_git_lfs_pointer_file(wi$spots_metadata)) log_msg("  WARN spots_metadata is Git LFS pointer: ", wi$spots_metadata, .level = "WARN")
  if (!is.na(wi$cell_metadata) && is_git_lfs_pointer_file(wi$cell_metadata)) log_msg("  WARN cell_metadata is Git LFS pointer: ", wi$cell_metadata, .level = "WARN")
  if (is.na(wi$expression)) {
    stop("Expression file missing for week ", wk, ". Check data roots and STARmap naming patterns.")
  }
  if (is.na(wi$spots_metadata) && is.na(wi$cell_metadata)) {
    stop("Both spots and cell metadata are missing for week ", wk, ". Need at least one metadata source with coordinates.")
  }
}

# =============================================================================
# D) Per-week Seurat Object Construction + QC
# =============================================================================

week_objects <- list()
qc_summary <- list()
resume_stage <- cfg$resume$resume_from

if (resume_stage == "week_qc") for (wk in cfg$week_ids) {
  log_msg("--- Processing week ", wk, " ---")
  wi <- week_inputs[[wk]]
  wk_ckpt <- file.path(DIR_OBJECTS, paste0("01_week_", wk, "_post_qc.rds"))

  if (isTRUE(cfg$resume$use_checkpoints) && isFALSE(cfg$resume$force_recompute) && ckpt_exists(wk_ckpt)) {
    seu_w <- ckpt_read(wk_ckpt)
    assert_no_spot_level_metadata(seu_w, context_label = paste0("week checkpoint ", wk))
    week_objects[[wk]] <- seu_w
    qc_summary[[wk]] <- tibble(
      week = wk,
      cells_before_qc = ncol(seu_w),
      cells_after_qc = ncol(seu_w),
      pct_retained = 100,
      median_nFeature = median(seu_w$nFeature_RNA),
      median_nCount = median(seu_w$nCount_RNA),
      median_percent_mt = median(seu_w$percent.mt)
    )
    next
  }

  spots_ptr <- !is.na(wi$spots_metadata) && is_git_lfs_pointer_file(wi$spots_metadata)
  cell_ptr <- !is.na(wi$cell_metadata) && is_git_lfs_pointer_file(wi$cell_metadata)

  if (spots_ptr && (is.na(wi$cell_metadata) || cell_ptr)) {
    log_msg("  Spots/cell metadata unavailable as real CSVs for week ", wk,
            ". Building placeholder metadata from expression matrix.", .level = "WARN")
    log_msg("  Reading expression (fallback mode): ", wi$expression)
    inferred <- coerce_expression_without_metadata(wi$expression, wk)
    counts <- inferred$counts
    md <- inferred$metadata
  } else {
    # Prefer cell_metadata as primary coordinate source when available;
    # spots_metadata can be very large (gene-level spots) and slow to parse.
    use_cell_primary <- !is.na(wi$cell_metadata) && !cell_ptr

    if (use_cell_primary) {
      log_msg("  Using cell metadata as primary coordinate source for week ", wk, ".")
      cell_md_raw <- read_csv_flexible(wi$cell_metadata)
      spots_md_raw <- NULL
    } else {
      spots_md_raw <- if (!is.na(wi$spots_metadata) && !spots_ptr) read_csv_flexible(wi$spots_metadata) else NULL
      cell_md_raw <- if (!is.na(wi$cell_metadata) && !cell_ptr) read_csv_flexible(wi$cell_metadata) else NULL
    }

    if (!is.null(spots_md_raw)) names(spots_md_raw) <- make.names(names(spots_md_raw), unique = TRUE)
    if (!is.null(cell_md_raw)) names(cell_md_raw) <- make.names(names(cell_md_raw), unique = TRUE)

    spots_md <- standardize_spots_metadata(
      md_raw = spots_md_raw,
      week_id = wk,
      barcodes = character(0),
      cell_md_raw = cell_md_raw
    )

    counts <- tryCatch(
      {
        log_msg("  Reading expression with metadata alignment: ", wi$expression)
        coerce_expression_with_spots(wi$expression, spots_md, wk)
      },
      error = function(e) {
        write_troubleshooting_bundle(
          week_id = wk,
          expr_path = wi$expression,
          spots_path = wi$spots_metadata,
          expr_df = tryCatch(read_csv_flexible(wi$expression), error = function(...) NULL),
          spots_df = if (!is.null(spots_md_raw)) spots_md_raw else cell_md_raw,
          err_msg = conditionMessage(e)
        )
        stop(e)
      }
    )

    md <- standardize_spots_metadata(
      md_raw = spots_md_raw,
      week_id = wk,
      barcodes = colnames(counts),
      cell_md_raw = cell_md_raw
    )
  }

  log_msg("  Counts loaded from STARmap CSV: genes=", nrow(counts), ", cells=", ncol(counts))

  seu_w <- CreateSeuratObject(counts = counts, project = paste0("STARmap_", wk))

  # Add harmonized metadata in cell order
  md <- as.data.frame(md, stringsAsFactors = FALSE)
  keep_md_cols <- c("barcode", "cell_id", "tile", "x_pixel", "y_pixel", "z_pixel", "x_um", "y_um", "z_um", "week", "sample_id", "metadata_source")
  keep_md_cols <- intersect(keep_md_cols, colnames(md))
  md <- md[, keep_md_cols, drop = FALSE]
  rownames(md) <- md$barcode
  md <- md[colnames(seu_w), , drop = FALSE]
  seu_w <- AddMetaData(seu_w, metadata = md)
  assert_no_spot_level_metadata(seu_w, context_label = paste0("week object ", wk, " after AddMetaData"))

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
  keep_primary <- (
    seu_w$nFeature_RNA >= cfg$qc$min_features &
      seu_w$nFeature_RNA <= cfg$qc$max_features &
      seu_w$nCount_RNA >= cfg$qc$min_counts &
      seu_w$nCount_RNA <= cfg$qc$max_counts &
      seu_w$percent.mt <= cfg$qc$max_percent_mt &
      !is.na(seu_w$x_um) & !is.na(seu_w$y_um)
  )

  n_keep_primary <- sum(keep_primary)
  if (n_keep_primary == 0) {
    log_msg(
      "  Primary QC retained 0 cells for ", wk,
      ". Applying adaptive fallback for imputed STARmap (coordinate-valid + counts>0).",
      .level = "WARN"
    )
    if (all(is.na(seu_w$x_um)) || all(is.na(seu_w$y_um))) {
      keep <- (seu_w$nCount_RNA > 0)
    } else {
      keep <- (!is.na(seu_w$x_um) & !is.na(seu_w$y_um) & seu_w$nCount_RNA > 0)
    }
  } else {
    keep <- keep_primary
  }

  if (sum(keep) == 0) {
    stop(
      "No cells retained after fallback QC for week ", wk,
      ". Check barcode alignment and coordinate columns in spots metadata."
    )
  }

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

  wk_ckpt_saved <- ckpt_write(seu_w, wk_ckpt)
  log_msg("  Saved week checkpoint: ", wk_ckpt_saved)
  week_objects[[wk]] <- seu_w
}

if (length(qc_summary) > 0) {
  qc_summary_tbl <- bind_rows(qc_summary)
  write_csv(qc_summary_tbl, file.path(DIR_TABLES, "01_qc_summary_by_week.csv"))
  log_msg("Per-week QC summary saved.")
} else {
  log_msg("QC summary not regenerated in this run (resume_from=", resume_stage, ").", .level = "WARN")
}

# =============================================================================
# E) Merge + Seurat v5 Hygiene
# =============================================================================

log_msg("Merging week-level Seurat objects...")
merged_ckpt <- file.path(DIR_OBJECTS, "01_merged_post_qc.rds")
if (resume_stage %in% c("merged", "sct", "pca")) {
  if (!ckpt_exists(merged_ckpt)) {
    stop("resume_from='", resume_stage, "' requires merged checkpoint: ", merged_ckpt)
  }
  seu <- ckpt_read(merged_ckpt)
  assert_no_spot_level_metadata(seu, context_label = "merged checkpoint")
} else if (isTRUE(cfg$resume$use_checkpoints) && isFALSE(cfg$resume$force_recompute) && ckpt_exists(merged_ckpt)) {
  seu <- ckpt_read(merged_ckpt)
  assert_no_spot_level_metadata(seu, context_label = "merged checkpoint")
} else {
  seu <- Reduce(function(x, y) merge(x, y), week_objects)
  assert_no_spot_level_metadata(seu, context_label = "merged object")
  merged_ckpt_saved <- ckpt_write(seu, merged_ckpt)
  log_msg("Saved merged checkpoint: ", merged_ckpt_saved)
}

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
post_sct_ckpt <- file.path(DIR_OBJECTS, "01_post_sct.rds")
if (resume_stage %in% c("sct", "pca")) {
  if (!ckpt_exists(post_sct_ckpt)) {
    stop("resume_from='", resume_stage, "' requires SCT checkpoint: ", post_sct_ckpt)
  }
  seu <- ckpt_read(post_sct_ckpt)
  assert_no_spot_level_metadata(seu, context_label = "SCTransform checkpoint")
} else if (isTRUE(cfg$resume$use_checkpoints) && isFALSE(cfg$resume$force_recompute) && ckpt_exists(post_sct_ckpt)) {
  seu <- ckpt_read(post_sct_ckpt)
  assert_no_spot_level_metadata(seu, context_label = "SCTransform checkpoint")
} else {
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
  post_sct_ckpt_saved <- ckpt_write(seu, post_sct_ckpt)
  log_msg("Saved SCTransform checkpoint: ", post_sct_ckpt_saved)
}
assert_no_spot_level_metadata(seu, context_label = "post-SCT object")
DefaultAssay(seu) <- "SCT"
log_msg("SCTransform complete. Default assay set to SCT.")

# =============================================================================
# G) PCA + Harmony Integration
# =============================================================================

log_msg("Running PCA...")
seu <- RunPCA(seu, assay = "SCT", npcs = max(cfg$dims_use), verbose = FALSE)

log_msg("Running Harmony integration by: ", cfg$harmony$group_by)
seu <- harmony::RunHarmony(
  object = seu,
  group.by.vars = cfg$harmony$group_by,
  reduction.use = cfg$harmony$reduction,
  assay.use = cfg$harmony$assay_use,
  max.iter.harmony = cfg$harmony$max_iter_harmony,
  theta = cfg$harmony$theta,
  lambda = cfg$harmony$lambda,
  sigma = cfg$harmony$sigma,
  verbose = TRUE
)
assert_no_spot_level_metadata(seu, context_label = "post-Harmony object")
post_harmony_saved <- ckpt_write(seu, file.path(DIR_OBJECTS, "01_post_harmony.rds"))
log_msg("Saved Harmony checkpoint: ", post_harmony_saved)

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
diag_tbl <- as.data.frame(seu@meta.data[, c("week", "seurat_clusters")]) %>%
  dplyr::count(week, seurat_clusters, name = "n_cells") %>%
  dplyr::arrange(week, seurat_clusters)

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
obj_path_saved <- ckpt_write(seu, obj_path)
log_msg("Saved integrated Seurat object: ", obj_path_saved)

# Run manifest for reproducibility
manifest_path <- file.path(DIR_REPORTS, "01_preprocess_harmony_manifest.json")
record_artifact_manifest(
  manifest_path = manifest_path,
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
  object_output = obj_path_saved,
  tables_output = c(
    file.path(DIR_TABLES, "01_qc_summary_by_week.csv"),
    file.path(DIR_TABLES, "01_cluster_composition_by_week.csv")
  ),
  reports_output = c(
    file.path(DIR_REPORTS, "01_methods_and_provenance.md"),
    file.path(DIR_REPORTS, "01_preprocess_harmony_manifest.json"),
    DIR_TROUBLE
  ),
  figures_output = c(
    file.path(DIR_FIGURES, "01_embeddings_week_umap_tsne.pdf"),
    file.path(DIR_FIGURES, "01_umap_harmony_by_week.png"),
    file.path(DIR_FIGURES, "01_tsne_harmony_by_week.png"),
    file.path(DIR_FIGURES, "01_embeddings_clusters_umap_tsne.pdf")
  ),
  source_data = unname(vapply(week_inputs, function(x) x$expression, character(1))),
  compute_script = "scripts/01_active_pipeline/01_preprocess_harmony_embeddings.R",
  plotting_script = "scripts/01_active_pipeline/01c_plot_spatial_embeddings.R"
)
log_msg("Saved run manifest: ", manifest_path)

log_msg("============================================================")
log_msg("Script complete: ", PIPELINE_NAME)
log_msg("Final dimensions: ", nrow(seu), " genes x ", ncol(seu), " cells")
log_msg("============================================================")
