source("R/config.R")
source("R/helpers_io.R")
source("R/helpers_methods.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(qs)
  library(ggplot2)
  library(tidyr)
  library(readr)
})

log_file <- file.path(CFG$dirs$logs, "pipeline.log")
log_msg("16_compute_misi (v2) starting.", log_file = log_file)

ensure_dir(CFG$dirs$tables)
ensure_dir(CFG$dirs$figures)
ensure_dir(CFG$dirs$objects)

safe_write_csv <- function(df, path) {
  tryCatch(readr::write_csv(df, path), error = function(e) {
    write.csv(df, path, row.names = FALSE)
  })
}

get_expr_data <- function(seu_obj, assay = "RNA", layer_name = "data") {
  m <- tryCatch(SeuratObject::LayerData(seu_obj, assay = assay, layer = layer_name), error = function(e) NULL)
  if (!is.null(m) && ncol(m) > 0) return(m)
  m <- tryCatch(SeuratObject::GetAssayData(seu_obj, assay = assay, layer = layer_name), error = function(e) NULL)
  if (!is.null(m) && ncol(m) > 0) return(m)
  m <- tryCatch(SeuratObject::GetAssayData(seu_obj, assay = assay, slot = layer_name), error = function(e) NULL)
  if (is.null(m) || ncol(m) == 0) stop("Could not retrieve expression layer '", layer_name, "' from assay '", assay, "'.")
  m
}

pc1_module_score <- function(expr_mat, genes) {
  g <- intersect(genes, rownames(expr_mat))
  if (length(g) < 2) return(rep(NA_real_, ncol(expr_mat)))
  x <- as.matrix(expr_mat[g, , drop = FALSE])
  x <- t(scale(t(x), center = TRUE, scale = TRUE))
  x <- x[apply(x, 1, function(v) all(is.finite(v))), , drop = FALSE]
  if (nrow(x) < 2) return(rep(NA_real_, ncol(expr_mat)))

  score <- tryCatch({
    if (requireNamespace("irlba", quietly = TRUE)) {
      s <- irlba::irlba(x, nv = 1, nu = 1)
      as.numeric(s$v[, 1] * s$d[1])
    } else {
      as.numeric(stats::prcomp(t(x), center = FALSE, scale. = FALSE)$x[, 1])
    }
  }, error = function(e) rep(NA_real_, ncol(expr_mat)))

  if (all(is.na(score))) return(score)
  m <- colMeans(x, na.rm = TRUE)
  if (is.finite(stats::cor(score, m, use = "pairwise.complete.obs")) &&
      stats::cor(score, m, use = "pairwise.complete.obs") < 0) {
    score <- -score
  }
  score
}

bootstrap_cluster_ci <- function(df, boot_n = 500, seed = 42) {
  set.seed(seed)
  gr <- unique(df[, c("donor", "cell_type", "condition")])
  out <- vector("list", nrow(gr))
  for (i in seq_len(nrow(gr))) {
    s <- df[df$donor == gr$donor[i] & df$cell_type == gr$cell_type[i] & df$condition == gr$condition[i], , drop = FALSE]
    if (nrow(s) == 0) next
    b <- replicate(boot_n, mean(sample(s$MISI, size = nrow(s), replace = TRUE), na.rm = TRUE))
    out[[i]] <- data.frame(
      donor = gr$donor[i],
      cell_type = gr$cell_type[i],
      condition = gr$condition[i],
      n_cells = nrow(s),
      MISI_mean = mean(s$MISI, na.rm = TRUE),
      ci_low = stats::quantile(b, 0.025, na.rm = TRUE),
      ci_high = stats::quantile(b, 0.975, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(out)
}

moran_multiradii <- function(df, value_col = "MISI", radii = c(50, 100, 200)) {
  if (!requireNamespace("spdep", quietly = TRUE)) return(data.frame())
  xy <- as.matrix(df[, c("x", "y")])
  out <- list()
  for (r in radii) {
    nb <- tryCatch(spdep::dnearneigh(xy, 0, r, longlat = FALSE), error = function(e) NULL)
    if (is.null(nb)) next
    lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
    lm <- spdep::localmoran(df[[value_col]], lw, zero.policy = TRUE)
    out[[as.character(r)]] <- data.frame(
      cell = df$cell,
      radius = r,
      Ii = lm[, 1],
      Z = lm[, 4],
      p = lm[, 5],
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(out)
}

obj_candidates <- c(
  file.path(CFG$dirs$objects, "seu_with_ea_flt1_proxies.qs"),
  file.path(CFG$dirs$objects, "seu_with_nk_flags.qs"),
  file.path(CFG$dirs$objects, "seu_with_scores.qs"),
  file.path(CFG$dirs$objects, "seu_core_architecture.qs"),
  file.path(CFG$dirs$objects, "seu_core_figs.qs"),
  file.path(CFG$dirs$objects, "seu_qc_checked.qs"),
  file.path(CFG$dirs$objects, "seu_clean.qs"),
  file.path(CFG$dirs$objects, "seu_with_architecture_transfer_lite.qs"),
  file.path(CFG$dirs$objects, "seu_with_architecture_transfer.qs"),
  file.path(CFG$dirs$objects, "seu_with_integration_checks.qs")
)
obj_candidates <- obj_candidates[file.exists(obj_candidates)]
if (length(obj_candidates) == 0) stop("No Seurat object found for script 16_compute_misi")

load_seurat_qs <- function(path) {
  attempts <- list(
    function() qs::qread(path, nthreads = 1, use_alt_rep = TRUE),
    function() qs::qread(path, use_alt_rep = TRUE),
    function() qs::qread(path, nthreads = 1),
    function() qs::qread(path)
  )
  for (fn in attempts) {
    obj <- tryCatch(fn(), error = function(e) NULL)
    if (!is.null(obj)) return(obj)
  }
  NULL
}

seu <- NULL
obj_path <- NULL
load_errors <- list()
for (cand in obj_candidates) {
  s <- load_seurat_qs(cand)
  if (!is.null(s) && inherits(s, "Seurat")) {
    seu <- s
    obj_path <- cand
    break
  }
  load_errors[[cand]] <- if (is.null(s)) "qread_failed_or_not_seurat" else paste0("loaded_class=", class(s)[1])
}
if (is.null(seu)) {
  err_txt <- paste(paste(names(load_errors), unlist(load_errors), sep = " => "), collapse = "; ")
  stop("Could not load any candidate Seurat object. Details: ", err_txt)
}
log_msg("Using object: ", obj_path, log_file = log_file)

if (!"RNA" %in% Assays(seu)) stop("RNA assay is required for MISI script.")
DefaultAssay(seu) <- "RNA"

if (!CFG$cols$cell_type %in% colnames(seu@meta.data)) {
  if ("celltype_refined" %in% colnames(seu@meta.data)) {
    seu[[CFG$cols$cell_type]] <- seu$celltype_refined
  } else if ("celltype_author" %in% colnames(seu@meta.data)) {
    seu[[CFG$cols$cell_type]] <- seu$celltype_author
  } else {
    stop("No cell-type column found.")
  }
}
if (!CFG$cols$infection %in% colnames(seu@meta.data) || !CFG$cols$hpi %in% colnames(seu@meta.data)) {
  stop("infection and hpi metadata columns are required.")
}
if (!"condition" %in% colnames(seu@meta.data)) {
  seu$condition <- paste0(seu[[CFG$cols$infection]][, 1], "_", seu[[CFG$cols$hpi]][, 1])
}
if (!CFG$cols$donor_id %in% colnames(seu@meta.data)) {
  seu[[CFG$cols$donor_id]] <- "donor_unknown"
}

expr <- get_expr_data(seu, assay = "RNA", layer_name = "data")

modules <- list(
  Smet = c("PLD1", "GDPD1", "GDPD5", "ETNK1", "FAAH", "PCYT2"),
  Simm = c("TLR2", "TLR4", "TNF", "IL6", "IL1B", "NFKB1", "NFKB2", "IL18", "CSF3"),
  Sbar = c("CDH1", "CDH5", "CLDN1", "CLDN3", "OCLN", "TJP1", "EPCAM"),
  Smmp = c("MMP1", "MMP2", "MMP9", "MMP14"),
  TIMP = c("TIMP1", "TIMP2", "TIMP3"),
  Sdef = c("DEFB1", "DEFB4A", "CAMP", "C3", "C5AR1", "LCN2"),
  Shome = c("GALNT1", "GALNT2", "GALNT3", "C1GALT1", "C1GALT1C1", "MGAT5", "GCNT1"),
  Hypoxia = c("HIF1A", "VEGFA", "LDHA", "PGK1", "ENO1"),
  Svaso = c("NOS3", "EDN1", "VCAM1", "ICAM1", "SELE"),
  Swnt = c("MYC", "CCND1", "AXIN2", "FOS", "JUN"),
  Slyso = c("LAMP1", "LAMP2", "MAP1LC3B", "SQSTM1")
)

presence <- data.frame(
  module = names(modules),
  n_present = vapply(modules, function(g) length(intersect(g, rownames(expr))), integer(1)),
  stringsAsFactors = FALSE
)
safe_write_csv(presence, file.path(CFG$dirs$tables, "misi_module_gene_presence.csv"))

module_scores <- lapply(modules, function(g) pc1_module_score(expr, g))
module_scores_df <- as.data.frame(module_scores)
module_scores_df$cell <- colnames(expr)
module_scores_df <- module_scores_df[, c("cell", names(modules))]
safe_write_csv(module_scores_df, file.path(CFG$dirs$tables, "misi_per_cell_module_scores.csv"))

z <- as.data.frame(scale(module_scores_df[, names(modules), drop = FALSE]))
z$Z_Sbar <- -z$Sbar
z$Z_Sdef <- -z$Sdef
z$Z_Smmp_composite <- z$Smmp - z$TIMP

X <- as.matrix(data.frame(
  Z_Smet = z$Smet,
  Z_Simm = z$Simm,
  Z_Sbar = z$Z_Sbar,
  Z_Smmp = z$Z_Smmp_composite,
  Z_Sdef = z$Z_Sdef,
  Z_Shome = z$Shome
))

misi_naive <- rowMeans(X, na.rm = TRUE)
weights <- rep(1 / ncol(X), ncol(X)); names(weights) <- colnames(X)
fit_status <- "equal_weights_fallback"

donor <- as.character(seu[[CFG$cols$donor_id]][, 1])
ctype <- as.character(seu[[CFG$cols$cell_type]][, 1])
cond <- as.character(seu$condition)

ratio <- NULL
if ("FLT1_PGF_ratio" %in% colnames(seu@meta.data)) {
  ratio <- as.numeric(seu$FLT1_PGF_ratio)
} else if (all(c("FLT1", "PGF") %in% rownames(expr))) {
  troph_idx <- grepl("SCT|STB|VCT|vCTB|EVT", ctype, ignore.case = TRUE)
  if (sum(troph_idx) > 20) {
    ratio <- rep(NA_real_, ncol(seu))
    ratio[troph_idx] <- as.numeric(log1p(expr["FLT1", troph_idx]) - log1p(expr["PGF", troph_idx]))
  }
}

if (!is.null(ratio) && requireNamespace("glmnet", quietly = TRUE)) {
  donor_df <- data.frame(donor = donor, ratio = ratio, X)
  donor_df <- donor_df[is.finite(donor_df$ratio), , drop = FALSE]
  donor_agg <- donor_df %>%
    group_by(donor) %>%
    summarise(ratio = mean(ratio, na.rm = TRUE), across(starts_with("Z_"), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

  if (nrow(donor_agg) >= 6) {
    xd <- as.matrix(donor_agg[, grepl("^Z_", colnames(donor_agg)), drop = FALSE])
    yd <- donor_agg$ratio

    cv <- tryCatch(glmnet::cv.glmnet(
      x = xd,
      y = yd,
      family = "gaussian",
      alpha = 1,
      nfolds = min(5, nrow(donor_agg)),
      lower.limits = 0
    ), error = function(e) NULL)

    if (!is.null(cv)) {
      coef_vec <- as.numeric(stats::coef(cv, s = "lambda.1se"))[-1]
      if (sum(coef_vec, na.rm = TRUE) > 0) {
        weights <- coef_vec / sum(coef_vec)
        names(weights) <- colnames(xd)
        fit_status <- "glmnet_gaussian_donor_ratio"
      }
      sens <- data.frame(
        feature = rownames(as.matrix(stats::coef(cv, s = "lambda.min"))),
        coef_lambda_min = as.numeric(as.matrix(stats::coef(cv, s = "lambda.min"))),
        coef_lambda_1se = as.numeric(as.matrix(stats::coef(cv, s = "lambda.1se"))),
        stringsAsFactors = FALSE
      )
      safe_write_csv(sens, file.path(CFG$dirs$tables, "glmnet_weight_sensitivity.csv"))
    }
  }
}

seu$MISI <- as.numeric(X %*% weights)
seu$MISI_naive <- misi_naive
seu$MISI_hotspot_flag <- as.integer(seu$MISI >= stats::quantile(seu$MISI, 0.90, na.rm = TRUE))

per_cell <- data.frame(
  cell = colnames(seu),
  donor = donor,
  cell_type = ctype,
  condition = cond,
  MISI = seu$MISI,
  MISI_naive = seu$MISI_naive,
  MISI_hotspot_flag = seu$MISI_hotspot_flag,
  X,
  stringsAsFactors = FALSE
)
safe_write_csv(per_cell, file.path(CFG$dirs$tables, "misi_per_cell.csv"))

boot <- bootstrap_cluster_ci(per_cell[, c("donor", "cell_type", "condition", "MISI")], boot_n = 500, seed = 42)
safe_write_csv(boot, file.path(CFG$dirs$tables, "bootstrap_cluster_CIs.csv"))

loo <- lapply(names(modules), function(mn) {
  g <- intersect(modules[[mn]], rownames(expr))
  if (length(g) < 3) return(NULL)
  full <- pc1_module_score(expr, g)
  corrs <- sapply(g, function(drop_g) {
    alt <- pc1_module_score(expr, setdiff(g, drop_g))
    suppressWarnings(stats::cor(full, alt, use = "pairwise.complete.obs"))
  })
  data.frame(module = mn, dropped_gene = g, corr_with_full = as.numeric(corrs), stringsAsFactors = FALSE)
})
safe_write_csv(dplyr::bind_rows(loo), file.path(CFG$dirs$tables, "gene_leaveoneout_correlations.csv"))

coord_cols <- if (all(c("x", "y") %in% colnames(seu@meta.data))) c("x", "y") else if (all(c("X", "Y") %in% colnames(seu@meta.data))) c("X", "Y") else NULL
if (!is.null(coord_cols) && requireNamespace("spdep", quietly = TRUE)) {
  spdf <- data.frame(cell = colnames(seu), x = as.numeric(seu@meta.data[[coord_cols[1]]]), y = as.numeric(seu@meta.data[[coord_cols[2]]]), MISI = seu$MISI, stringsAsFactors = FALSE)
  spdf <- spdf[is.finite(spdf$x) & is.finite(spdf$y) & is.finite(spdf$MISI), , drop = FALSE]
  if (nrow(spdf) > 20) {
    mor <- moran_multiradii(spdf, value_col = "MISI", radii = c(50, 100, 200))
    safe_write_csv(mor, file.path(CFG$dirs$tables, "spatial_moran_multiradii.csv"))
    if (nrow(mor) > 0) {
      thr <- mor %>% group_by(radius) %>% summarise(t = quantile(Ii, 0.90, na.rm = TRUE), .groups = "drop")
      hot <- mor %>% left_join(thr, by = "radius") %>% filter(Ii >= t & p < 0.05)
      safe_write_csv(hot, file.path(CFG$dirs$tables, "toxic_blast_zone_cells.csv"))
    }
  }
}

weights_tbl <- data.frame(component = names(weights), weight = as.numeric(weights), fit_status = fit_status, stringsAsFactors = FALSE)
safe_write_csv(weights_tbl, file.path(CFG$dirs$tables, "misi_weights.csv"))

sum_tbl <- per_cell %>%
  group_by(cell_type, condition) %>%
  summarise(n_cells = dplyr::n(), MISI_mean = mean(MISI, na.rm = TRUE), MISI_median = median(MISI, na.rm = TRUE), hotspot_pct = mean(MISI_hotspot_flag, na.rm = TRUE), .groups = "drop")
safe_write_csv(sum_tbl, file.path(CFG$dirs$tables, "misi_summary_by_celltype_condition.csv"))

p <- ggplot(per_cell, aes(x = cell_type, y = MISI, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "MISI by cell type and condition", x = "Cell type", y = "MISI")
save_plot(file.path(CFG$dirs$figures, "FigMISI_violin_by_celltype_infection.pdf"), p, w = 14, h = 7)

qs::qsave(seu, file.path(CFG$dirs$objects, "seu_with_misi.qs"), preset = "high")

write_legend(
  "FigMISI", "Mechanistic Infection Susceptibility Index (MISI)",
  hypothesis = "A composite vulnerability score integrating metabolic, immune, barrier, MMP/TIMP, defense, and glyco-homeostasis signals identifies infection-sensitive placental niches.",
  methods = "PC1 module scoring (with signed vulnerability orientation), optional donor-level non-negative glmnet weighting against independent FLT1/PGF ratio, bootstrap CIs, leave-one-out, and optional multi-radius local Moran analysis for spatial data.",
  readout = "Higher MISI and hotspot fraction indicate stronger vulnerability signatures.",
  interpretation_template = "Confirm inferential claims with donor-aware tests and compare weighted vs naive MISI agreement.",
  outfile = file.path(CFG$dirs$legends, "FigMISI_legend.md")
)

status <- data.frame(
  object_used = obj_path,
  glmnet_fit_status = fit_status,
  spatial_coords_detected = !is.null(coord_cols),
  radii_tested = paste(c(50, 100, 200), collapse = ";"),
  stringsAsFactors = FALSE
)
safe_write_csv(status, file.path(CFG$dirs$tables, "misi_run_status.csv"))


append_methods_draft(
  script_name = "16_compute_misi.R",
  hypothesis = "A multicomponent spatial vulnerability index better captures pathogenic injury risk than any single marker threshold.",
  method_details = list(
    "Compute PC1 module scores for metabolic, immune, barrier, MMP/TIMP, defense, and glyco-homeostasis modules",
    "Construct MISI from signed z-scored components with optional donor-level non-negative glmnet weighting",
    "Bootstrap confidence intervals by donor/cell-type/condition",
    "Optional local Moran's I across radii to localize toxic blast zones"
  ),
  rationale = "Single-gene cutoffs are unstable across donors and cell states, whereas a composite score integrates convergent biological axes (metabolic stress, barrier failure, inflammatory activation, and matrix remodeling) to estimate local vulnerability with higher robustness.",
  citations = c(
    "Hao Y, Hao S, Andersen-Nissen E, et al. Integrated analysis of multimodal single-cell data. Cell. 2021;184(13):3573-3587.e29. doi:10.1016/j.cell.2021.04.048",
    "Palla G, Spitzer H, Klein M, et al. Squidpy: a scalable framework for spatial omics analysis. Nature Methods. 2022;19:171-178. doi:10.1038/s41592-021-01358-2"
  )
)

log_msg("16_compute_misi (v2) done.", log_file = log_file)
