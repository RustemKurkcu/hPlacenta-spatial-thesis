source("R/config.R")
source("R/helpers_io.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(qs)
  library(ggplot2)
})

log_file <- file.path(CFG$dirs$logs, "pipeline.log")
log_msg("17_cellchat_by_condition starting.", log_file = log_file)

ensure_dir(CFG$dirs$tables)
ensure_dir(CFG$dirs$figures)
ensure_dir(CFG$dirs$objects)

get_expr_data <- function(seu_obj, assay = "RNA", layer_name = "data") {
  m <- tryCatch(SeuratObject::LayerData(seu_obj, assay = assay, layer = layer_name), error = function(e) NULL)
  if (!is.null(m) && ncol(m) > 0) return(m)
  m <- tryCatch(SeuratObject::GetAssayData(seu_obj, assay = assay, layer = layer_name), error = function(e) NULL)
  if (!is.null(m) && ncol(m) > 0) return(m)
  m <- tryCatch(SeuratObject::GetAssayData(seu_obj, assay = assay, slot = layer_name), error = function(e) NULL)
  if (is.null(m) || ncol(m) == 0) stop("Could not retrieve expression layer '", layer_name, "' from assay '", assay, "'.")
  m
}

obj_candidates <- c(
  file.path(CFG$dirs$objects, "seu_with_misi.qs")
)
obj_path <- obj_candidates[1]
if (!file.exists(obj_path)) stop("Missing required Seurat object for script 17_cellchat_by_condition: ", obj_path)

seu <- tryCatch(qs::qread(obj_path, nthreads = 1), error = function(e) NULL)
if (is.null(seu) || !inherits(seu, "Seurat")) stop("Could not load required Seurat object: ", obj_path)
log_msg("Using object: ", obj_path, log_file = log_file)

if (!requireNamespace("CellChat", quietly = TRUE)) {
  note <- data.frame(
    status = "skipped",
    reason = "CellChat not installed",
    object_used = obj_path,
    stringsAsFactors = FALSE
  )
  write.csv(note, file.path(CFG$dirs$tables, "cellchat_by_condition_status.csv"), row.names = FALSE)
  log_msg("17_cellchat_by_condition skipped: CellChat not installed.", log_file = log_file)
  quit(save = "no")
}

suppressPackageStartupMessages(library(CellChat))

if (!"RNA" %in% Assays(seu)) stop("RNA assay is required for CellChat script")
DefaultAssay(seu) <- "RNA"

if (!CFG$cols$cell_type %in% colnames(seu@meta.data)) {
  if ("celltype_refined" %in% colnames(seu@meta.data)) {
    seu[[CFG$cols$cell_type]] <- seu$celltype_refined
  } else if ("celltype_author" %in% colnames(seu@meta.data)) {
    seu[[CFG$cols$cell_type]] <- seu$celltype_author
  } else {
    stop("No cell type grouping column found.")
  }
}
if (!"condition" %in% colnames(seu@meta.data)) {
  if (all(c(CFG$cols$infection, CFG$cols$hpi) %in% colnames(seu@meta.data))) {
    seu$condition <- paste0(seu[[CFG$cols$infection]][, 1], "_", seu[[CFG$cols$hpi]][, 1])
  } else {
    stop("condition (or infection+hpi) metadata is required")
  }
}

# Expression matrix for CellChat (Seurat v5/compat fallback)
expr <- get_expr_data(seu, assay = "RNA", layer_name = "data")
meta <- seu@meta.data[colnames(expr), , drop = FALSE]
meta$group <- as.character(meta[[CFG$cols$cell_type]])

# optional downsample to cap memory per condition
max_cells_per_condition <- 30000
set.seed(42)

RUN_PER_CONDITION <- TRUE
RUN_JOINT_CELLCHAT <- TRUE
JOINT_CONDITION_PAIRS <- list(
  c("UI_24h", "Lm_24h"),
  c("UI_24h", "Pf_24h"),
  c("UI_24h", "Tg_24h")
)

compute_prob_backoff <- function(cc_obj, cond_label) {
  attempts <- list(
    list(label = "raw_counts", raw.use = TRUE, type = "triMean", trim = NULL),
    list(label = "norm_triMean", raw.use = FALSE, type = "triMean", trim = NULL),
    list(label = "norm_trunc_0.1", raw.use = FALSE, type = "truncatedMean", trim = 0.1),
    list(label = "norm_trunc_0.2", raw.use = FALSE, type = "truncatedMean", trim = 0.2)
  )
  for (i in seq_along(attempts)) {
    a <- attempts[[i]]
    cc_try <- tryCatch({
      if (is.null(a$trim)) {
        CellChat::computeCommunProb(cc_obj, raw.use = a$raw.use, population.size = TRUE, type = a$type)
      } else {
        CellChat::computeCommunProb(cc_obj, raw.use = a$raw.use, population.size = TRUE, type = a$type, trim = a$trim)
      }
    }, error = function(e) NULL)

    if (is.null(cc_try)) next
    comm_try <- tryCatch(CellChat::subsetCommunication(cc_try), error = function(e) data.frame())
    if (nrow(comm_try) > 0) {
      diag <- data.frame(step = i, method = a$label, n_comm = nrow(comm_try), stringsAsFactors = FALSE)
      write.csv(diag, file.path(CFG$dirs$tables, paste0("cellchat_prob_method_", cond_label, ".csv")), row.names = FALSE)
      return(cc_try)
    }
  }
  cc_obj
}

run_cellchat_for <- function(cond_label, cells, min_cells = 10) {
  if (length(cells) < 50) {
    write.csv(data.frame(note = "too_few_cells", n_cells = length(cells)),
              file.path(CFG$dirs$tables, paste0("cellchat_comm_", cond_label, ".csv")), row.names = FALSE)
    return(NULL)
  }

  if (length(cells) > max_cells_per_condition) {
    cells <- sample(cells, max_cells_per_condition)
  }
  data.sub <- expr[, cells, drop = FALSE]
  meta.sub <- meta[cells, , drop = FALSE]
  rownames(meta.sub) <- colnames(data.sub)

  group_counts <- meta.sub %>% count(group, name = "n_cells") %>% arrange(desc(n_cells))
  write.csv(group_counts, file.path(CFG$dirs$tables, paste0("cellchat_group_counts_", cond_label, ".csv")), row.names = FALSE)

  cc <- CellChat::createCellChat(object = data.sub, meta = meta.sub, group.by = "group")
  cc@DB <- CellChat::CellChatDB.human
  cc <- CellChat::subsetData(cc)
  cc <- CellChat::identifyOverExpressedGenes(cc)
  cc <- CellChat::identifyOverExpressedInteractions(cc)

  cc <- compute_prob_backoff(cc, cond_label)

  # adaptive min.cells search so users can see where interactions appear/disappear
  min_cells_grid <- sort(unique(c(3, 5, 10, 20, min_cells)))
  grid_results <- lapply(min_cells_grid, function(m) {
    cc_tmp <- tryCatch(CellChat::filterCommunication(cc, min.cells = m), error = function(e) NULL)
    if (is.null(cc_tmp)) return(data.frame(min_cells = m, n_comm = NA_integer_))
    comm_tmp <- tryCatch(CellChat::subsetCommunication(cc_tmp), error = function(e) data.frame())
    data.frame(min_cells = m, n_comm = nrow(comm_tmp))
  })
  grid_df <- do.call(rbind, grid_results)
  write.csv(grid_df, file.path(CFG$dirs$tables, paste0("cellchat_min_cells_grid_", cond_label, ".csv")), row.names = FALSE)

  viable <- grid_df$min_cells[!is.na(grid_df$n_comm) & grid_df$n_comm > 0]
  if (length(viable) == 0) {
    write.csv(data.frame(note = "no_interactions_after_compute_prob", condition = cond_label),
              file.path(CFG$dirs$tables, paste0("cellchat_comm_", cond_label, ".csv")), row.names = FALSE)
    qs::qsave(cc, file.path(CFG$dirs$objects, paste0("cellchat_", cond_label, ".qs")), preset = "high")
    return(TRUE)
  }

  min_cells_used <- max(viable)
  cc <- CellChat::filterCommunication(cc, min.cells = min_cells_used)
  comm <- tryCatch(CellChat::subsetCommunication(cc), error = function(e) data.frame())
  write.csv(data.frame(min_cells_used = min_cells_used, n_comm = nrow(comm)),
            file.path(CFG$dirs$tables, paste0("cellchat_selected_threshold_", cond_label, ".csv")), row.names = FALSE)

  comm <- comm %>%
    left_join(group_counts %>% rename(source = group, n_source = n_cells), by = "source") %>%
    left_join(group_counts %>% rename(target = group, n_target = n_cells), by = "target") %>%
    mutate(condition = cond_label)
  write.csv(comm, file.path(CFG$dirs$tables, paste0("cellchat_comm_", cond_label, ".csv")), row.names = FALSE)

  focus <- comm %>%
    filter(
      ligand %in% c("AREG", "HGF", "VEGFA", "PGF", "IL1B", "TNF", "CCL4", "CCL20", "CXCL8", "LGALS3") |
        receptor %in% c("EGFR", "MET", "KDR", "FLT1", "IL1R1", "ITGB1")
    )
  write.csv(focus, file.path(CFG$dirs$tables, paste0("cellchat_targeted_axes_", cond_label, ".csv")), row.names = FALSE)

  # Pathway aggregation can fail when communication matrix is sparse; handle gracefully.
  cc <- tryCatch(CellChat::computeCommunProbPathway(cc), error = function(e) {
    log_msg("computeCommunProbPathway failed for ", cond_label, ": ", conditionMessage(e), log_file = log_file)
    cc
  })
  cc <- tryCatch(CellChat::aggregateNet(cc), error = function(e) {
    log_msg("aggregateNet failed for ", cond_label, ": ", conditionMessage(e), log_file = log_file)
    cc
  })

  mat <- tryCatch(cc@net$count, error = function(e) NULL)
  if (!is.null(mat) && length(mat) > 0 && sum(mat, na.rm = TRUE) > 0) {
    mat_df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
    colnames(mat_df) <- c("source", "target", "n_interactions")
    p <- ggplot(mat_df, aes(x = source, y = target, fill = n_interactions)) +
      geom_tile(color = "grey85") +
      scale_fill_viridis_c(trans = "log1p") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste0("CellChat interaction counts: ", cond_label), x = "Source", y = "Target", fill = "count")
    save_plot(file.path(CFG$dirs$figures, paste0("FigCellChat_heatmap_", cond_label, ".pdf")), p, w = 10, h = 8)
  }

  qs::qsave(cc, file.path(CFG$dirs$objects, paste0("cellchat_", cond_label, ".qs")), preset = "high")
  invisible(TRUE)
}

run_joint_cellchat <- function(cond_a, cond_b, min_cells = 10) {
  label <- paste0(gsub("[^A-Za-z0-9_.-]", "_", cond_a), "_vs_", gsub("[^A-Za-z0-9_.-]", "_", cond_b))
  cells <- rownames(meta)[meta$condition %in% c(cond_a, cond_b)]
  if (length(cells) < 100) {
    write.csv(data.frame(note = "too_few_cells_for_joint", cond_a = cond_a, cond_b = cond_b, n_cells = length(cells)),
              file.path(CFG$dirs$tables, paste0("cellchat_joint_", label, ".csv")), row.names = FALSE)
    return(NULL)
  }
  if (length(cells) > max_cells_per_condition * 2) cells <- sample(cells, max_cells_per_condition * 2)

  data.sub <- expr[, cells, drop = FALSE]
  meta.sub <- meta[cells, , drop = FALSE] %>%
    mutate(
      condition = factor(condition, levels = c(cond_a, cond_b)),
      group = as.character(group)
    ) %>%
    arrange(condition, group)
  data.sub <- data.sub[, rownames(meta.sub), drop = FALSE]
  ordered_joint_levels <- unique(paste0(as.character(meta.sub$condition), "__", meta.sub$group))
  meta.sub$joint_group <- factor(
    paste0(as.character(meta.sub$condition), "__", meta.sub$group),
    levels = ordered_joint_levels
  )
  rownames(meta.sub) <- colnames(data.sub)

  group_counts <- meta.sub %>% count(joint_group, name = "n_cells") %>% arrange(desc(n_cells))
  write.csv(group_counts, file.path(CFG$dirs$tables, paste0("cellchat_joint_group_counts_", label, ".csv")), row.names = FALSE)

  cc <- CellChat::createCellChat(object = data.sub, meta = meta.sub, group.by = "joint_group")
  cc@DB <- CellChat::CellChatDB.human
  cc <- CellChat::subsetData(cc)
  cc <- CellChat::identifyOverExpressedGenes(cc)
  cc <- CellChat::identifyOverExpressedInteractions(cc)
  cc <- compute_prob_backoff(cc, paste0("joint_", label))

  min_cells_grid <- sort(unique(c(3, 5, 10, 20, min_cells)))
  grid_results <- lapply(min_cells_grid, function(m) {
    cc_tmp <- tryCatch(CellChat::filterCommunication(cc, min.cells = m), error = function(e) NULL)
    if (is.null(cc_tmp)) return(data.frame(min_cells = m, n_comm = NA_integer_))
    comm_tmp <- tryCatch(CellChat::subsetCommunication(cc_tmp), error = function(e) data.frame())
    data.frame(min_cells = m, n_comm = nrow(comm_tmp))
  })
  grid_df <- do.call(rbind, grid_results)
  write.csv(grid_df, file.path(CFG$dirs$tables, paste0("cellchat_joint_min_cells_grid_", label, ".csv")), row.names = FALSE)

  viable <- grid_df$min_cells[!is.na(grid_df$n_comm) & grid_df$n_comm > 0]
  if (length(viable) == 0) {
    write.csv(data.frame(note = "no_joint_interactions_after_compute_prob", cond_a = cond_a, cond_b = cond_b),
              file.path(CFG$dirs$tables, paste0("cellchat_joint_", label, ".csv")), row.names = FALSE)
    qs::qsave(cc, file.path(CFG$dirs$objects, paste0("cellchat_joint_", label, ".qs")), preset = "high")
    return(TRUE)
  }

  min_cells_used <- max(viable)
  cc <- CellChat::filterCommunication(cc, min.cells = min_cells_used)
  cc <- tryCatch(CellChat::computeCommunProbPathway(cc), error = function(e) cc)
  cc <- tryCatch(CellChat::aggregateNet(cc), error = function(e) cc)
  comm <- tryCatch(CellChat::subsetCommunication(cc), error = function(e) data.frame())
  write.csv(data.frame(min_cells_used = min_cells_used, n_comm = nrow(comm)),
            file.path(CFG$dirs$tables, paste0("cellchat_joint_selected_threshold_", label, ".csv")), row.names = FALSE)
  write.csv(comm, file.path(CFG$dirs$tables, paste0("cellchat_joint_", label, ".csv")), row.names = FALSE)

  if (nrow(comm) > 0 && all(c("source", "target", "prob") %in% colnames(comm))) {
    comm_collapsed <- comm %>%
      mutate(
        source_condition = sub("__.*$", "", source),
        source_group = sub("^.*__", "", source),
        target_condition = sub("__.*$", "", target),
        target_group = sub("^.*__", "", target)
      ) %>%
      group_by(source_condition, source_group, target_condition, target_group) %>%
      summarise(prob_sum = sum(prob, na.rm = TRUE), n_edges = dplyr::n(), .groups = "drop") %>%
      arrange(desc(prob_sum))
    write.csv(comm_collapsed, file.path(CFG$dirs$tables, paste0("cellchat_joint_collapsed_", label, ".csv")), row.names = FALSE)
  }

  mat <- tryCatch(cc@net$count, error = function(e) NULL)
  if (!is.null(mat) && length(mat) > 0 && sum(mat, na.rm = TRUE) > 0) {
    mat_df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
    colnames(mat_df) <- c("source", "target", "n_interactions")
    p <- ggplot(mat_df, aes(x = source, y = target, fill = n_interactions)) +
      geom_tile(color = "grey85") +
      scale_fill_viridis_c(trans = "log1p") +
      theme_bw(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste0("Joint CellChat: ", cond_a, " vs ", cond_b), x = "Source", y = "Target", fill = "count")
    save_plot(file.path(CFG$dirs$figures, paste0("FigCellChat_joint_heatmap_", label, ".pdf")), p, w = 11, h = 9)
  }

  qs::qsave(cc, file.path(CFG$dirs$objects, paste0("cellchat_joint_", label, ".qs")), preset = "high")
  TRUE
}

conditions <- sort(unique(as.character(meta$condition)))
ran <- c()
if (isTRUE(RUN_PER_CONDITION)) {
  for (cond in conditions) {
    cond_label <- gsub("[^A-Za-z0-9_.-]", "_", cond)
    cells <- rownames(meta)[meta$condition == cond]
    ok <- tryCatch(run_cellchat_for(cond_label, cells, min_cells = 10), error = function(e) {
      log_msg("Condition failed ", cond_label, ": ", conditionMessage(e), log_file = log_file)
      NULL
    })
    if (!is.null(ok)) ran <- c(ran, cond_label)
  }
}

if (isTRUE(RUN_JOINT_CELLCHAT)) {
  for (pair in JOINT_CONDITION_PAIRS) {
    if (length(pair) != 2) next
    if (!all(pair %in% conditions)) next
    tryCatch(
      run_joint_cellchat(pair[1], pair[2], min_cells = 10),
      error = function(e) log_msg("Joint CellChat failed for ", pair[1], " vs ", pair[2], ": ", conditionMessage(e), log_file = log_file)
    )
  }
}

# Pairwise comparison helper
compare_conditions <- function(a, b, min_prob = 1e-4, min_cells = 10) {
  fa <- file.path(CFG$dirs$tables, paste0("cellchat_comm_", a, ".csv"))
  fb <- file.path(CFG$dirs$tables, paste0("cellchat_comm_", b, ".csv"))
  if (!file.exists(fa) || !file.exists(fb)) return(NULL)
  ca <- read.csv(fa, stringsAsFactors = FALSE)
  cb <- read.csv(fb, stringsAsFactors = FALSE)
  req_cols <- c("source", "target", "ligand", "receptor", "prob", "n_source", "n_target")
  if (!all(req_cols %in% colnames(ca)) || !all(req_cols %in% colnames(cb))) return(NULL)

  joined <- dplyr::full_join(ca, cb, by = c("source", "target", "ligand", "receptor"), suffix = c("_a", "_b")) %>%
    mutate(
      prob_a = ifelse(is.na(prob_a), 0, prob_a),
      prob_b = ifelse(is.na(prob_b), 0, prob_b),
      n_source_a = ifelse(is.na(n_source_a), 0, n_source_a),
      n_target_a = ifelse(is.na(n_target_a), 0, n_target_a),
      n_source_b = ifelse(is.na(n_source_b), 0, n_source_b),
      n_target_b = ifelse(is.na(n_target_b), 0, n_target_b),
      support_ok = pmax(n_source_a, n_source_b) >= min_cells & pmax(n_target_a, n_target_b) >= min_cells,
      signal_ok = pmax(prob_a, prob_b) >= min_prob,
      log2FC = log2((prob_a + 1e-4) / (prob_b + 1e-4))
    ) %>%
    filter(support_ok, signal_ok) %>%
    arrange(desc(abs(log2FC)))
  write.csv(joined, file.path(CFG$dirs$tables, paste0("cellchat_compare_", a, "_vs_", b, ".csv")), row.names = FALSE)
}

make_native_bubble <- function(a, b) {
  qa <- file.path(CFG$dirs$objects, paste0("cellchat_", a, ".qs"))
  qb <- file.path(CFG$dirs$objects, paste0("cellchat_", b, ".qs"))
  if (!file.exists(qa) || !file.exists(qb)) return(invisible(NULL))
  cca <- tryCatch(qs::qread(qa), error = function(e) NULL)
  ccb <- tryCatch(qs::qread(qb), error = function(e) NULL)
  if (is.null(cca) || is.null(ccb)) return(invisible(NULL))

  merged <- tryCatch(CellChat::mergeCellChat(list(condA = cca, condB = ccb), add.names = c(a, b)), error = function(e) NULL)
  if (is.null(merged)) return(invisible(NULL))

  p <- tryCatch(
    CellChat::netVisual_bubble(merged, comparison = c(1, 2), angle.x = 45),
    error = function(e) NULL
  )
  if (!is.null(p)) {
    save_plot(file.path(CFG$dirs$figures, paste0("FigCellChat_bubble_", a, "_vs_", b, ".pdf")), p, w = 10, h = 10)
  }
}

# default comparisons with flexible naming (handles e.g., Lm_24h and L.monocytogenes_24)
pick_cond <- function(x, infection_regex, hpi_regex) {
  x[grepl(infection_regex, x, ignore.case = TRUE) & grepl(hpi_regex, x, ignore.case = TRUE)][1]
}
lm24 <- pick_cond(ran, "lm|l\\.monocytogenes|listeria", "24")
ui24 <- pick_cond(ran, "^ui", "24")
lm48 <- pick_cond(ran, "lm|l\\.monocytogenes|listeria", "48")
ui48 <- pick_cond(ran, "^ui", "48")

if (!is.na(lm24) && !is.na(ui24)) {
  compare_conditions(lm24, ui24)
  make_native_bubble(lm24, ui24)
}
if (!is.na(lm48) && !is.na(ui48)) {
  compare_conditions(lm48, ui48)
  make_native_bubble(lm48, ui48)
}

run_status <- data.frame(
  object_used = obj_path,
  n_conditions = length(conditions),
  n_conditions_ran = length(ran),
  run_joint_cellchat = RUN_JOINT_CELLCHAT,
  max_cells_per_condition = max_cells_per_condition,
  stringsAsFactors = FALSE
)
write.csv(run_status, file.path(CFG$dirs$tables, "cellchat_by_condition_run_status.csv"), row.names = FALSE)

write_legend(
  "FigCellChatByCondition", "Condition-specific CellChat (non-spatial)",
  hypothesis = "Condition-specific communication networks reveal pathogen/time-specific sender-receiver changes not visible in pooled analysis.",
  methods = "Ran CellChat per condition and (optionally) in joint UI-vs-infection mode using condition+celltype groups; exported threshold grids, full/targeted interaction tables, and heatmaps.",
  readout = "Use per-condition and joint min-cells grids to identify viable thresholds; then inspect joint/collapsed CSVs and compare tables for robust shifts.",
  interpretation_template = "- Interpret low-cell-count edges cautiously using n_source/n_target.\n- Prioritize targeted axes and validate with independent imaging/orthogonal assays.",
  outfile = file.path(CFG$dirs$legends, "FigCellChatByCondition_legend.md")
)

log_msg("17_cellchat_by_condition done.", log_file = log_file)
