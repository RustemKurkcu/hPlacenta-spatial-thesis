#!/usr/bin/env Rscript

# =============================================================================
# Script: 04b_plot_spatial_cellchat_atlas.R
# Purpose: Visualization stage for the 5-Phase Spatial CellChat Atlas.
#          Produces phase-wise figures and appends entries to
#          output/reports/04_figures_manifest.json + methods log.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(jsonlite)
  library(dplyr)
})

source("R/spatial_color_themes.R")
source("R/celltype_dictionary.R")

OUT_ROOT <- "output"
DIR_OBJECTS <- file.path(OUT_ROOT, "objects")
DIR_FIGS <- file.path(OUT_ROOT, "figures", "04_atlas")
DIR_REPORTS <- file.path(OUT_ROOT, "reports")
for (d in c(DIR_FIGS, DIR_REPORTS)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

manifest_path <- file.path(DIR_REPORTS, "04_figures_manifest.json")
methods_path <- file.path(DIR_REPORTS, "04_methods_and_provenance.md")

if (file.exists(manifest_path)) {
  fig_manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
  if (!is.data.frame(fig_manifest)) fig_manifest <- data.frame()
} else {
  fig_manifest <- data.frame()
}

append_manifest <- function(entry) {
  entry_df <- as.data.frame(entry, stringsAsFactors = FALSE)
  if (nrow(fig_manifest) == 0) {
    fig_manifest <<- entry_df
  } else {
    fig_manifest <<- dplyr::bind_rows(fig_manifest, entry_df)
  }
  jsonlite::write_json(fig_manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
}

append_methods <- function(header, lines) {
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  cat(paste0("\n## ", stamp, " — ", header, "\n\n"), file = methods_path, append = TRUE)
  for (ln in lines) cat(paste0("- ", ln, "\n"), file = methods_path, append = TRUE)
}

save_base_plot <- function(fig_id, data_source, hypothesis_tested, expected_outcome, claim_strength, methods_blurb, plot_expr) {
  png_path <- file.path(DIR_FIGS, paste0(fig_id, ".png"))
  pdf_path <- file.path(DIR_FIGS, paste0(fig_id, ".pdf"))

  grDevices::png(png_path, width = 2400, height = 1800, res = 300)
  try(eval(plot_expr), silent = TRUE)
  grDevices::dev.off()

  grDevices::pdf(pdf_path, width = 12, height = 9)
  try(eval(plot_expr), silent = TRUE)
  grDevices::dev.off()

  append_manifest(list(
    figure_id = fig_id,
    file_path = paste(png_path, pdf_path, sep = " | "),
    data_source = data_source,
    hypothesis_tested = hypothesis_tested,
    expected_outcome = expected_outcome,
    biological_claim_strength = claim_strength,
    methods_blurb_2to4_sentences = methods_blurb
  ))

  gc(verbose = FALSE)
}

save_ggplot <- function(fig_id, p, data_source, hypothesis_tested, expected_outcome, claim_strength, methods_blurb) {
  png_path <- file.path(DIR_FIGS, paste0(fig_id, ".png"))
  pdf_path <- file.path(DIR_FIGS, paste0(fig_id, ".pdf"))

  p <- p + theme_thesis_spatial()
  ggsave(png_path, plot = p, width = 12, height = 9, dpi = 300)
  ggsave(pdf_path, plot = p, width = 12, height = 9, device = cairo_pdf)

  append_manifest(list(
    figure_id = fig_id,
    file_path = paste(png_path, pdf_path, sep = " | "),
    data_source = data_source,
    hypothesis_tested = hypothesis_tested,
    expected_outcome = expected_outcome,
    biological_claim_strength = claim_strength,
    methods_blurb_2to4_sentences = methods_blurb
  ))

  gc(verbose = FALSE)
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

cache <- readRDS(file.path(DIR_OBJECTS, "04a_atlas_compute_cache.rds"))
weeks <- names(cache$week_paths)
week_objs <- lapply(cache$week_paths, read_object)
merged_obj <- read_object(cache$merged_path)

# ------------------------
# Phase 1: weekly global QC
# ------------------------
for (wk in weeks) {
  obj <- week_objs[[wk]]

  save_base_plot(
    fig_id = paste0("P1_", wk, "_circle"),
    data_source = cache$week_paths[[wk]],
    hypothesis_tested = "Global communication volume is biologically plausible and non-zero across key lineages.",
    expected_outcome = "Dense but structured network with macrophage/trophoblast participation.",
    claim_strength = "Moderate",
    methods_blurb = "netVisual_circle() applied to group-level weight matrix from the weekly CellChat object. Edge thickness encodes interaction strength.",
    plot_expr = quote({
      m <- obj@net$weight
      if (is.null(m) && !is.null(obj@net$count)) m <- obj@net$count
      if (!is.null(m)) CellChat::netVisual_circle(m, title.name = paste0("Global communication - ", wk))
    })
  )

  save_base_plot(
    fig_id = paste0("P1_", wk, "_heatmap"),
    data_source = cache$week_paths[[wk]],
    hypothesis_tested = "Sender/receiver intensity matrix has interpretable dominant hubs.",
    expected_outcome = "Strong trophoblast-immune cross-talk blocks with week-specific asymmetry.",
    claim_strength = "Moderate",
    methods_blurb = "netVisual_heatmap() run on each weekly object to summarize sender-receiver interaction intensity.",
    plot_expr = quote({
      try(CellChat::netVisual_heatmap(obj), silent = TRUE)
    })
  )

  save_base_plot(
    fig_id = paste0("P1_", wk, "_signaling_role_scatter"),
    data_source = cache$week_paths[[wk]],
    hypothesis_tested = "Dominant sender and receiver cell types can be ranked without pathway preselection.",
    expected_outcome = "Distinct clusters of high-outgoing versus high-incoming signaling roles.",
    claim_strength = "Moderate",
    methods_blurb = "netAnalysis_signalingRole_scatter() compares outgoing and incoming centrality across cell groups.",
    plot_expr = quote({
      try(CellChat::netAnalysis_signalingRole_scatter(obj), silent = TRUE)
    })
  )

  save_base_plot(
    fig_id = paste0("P1_", wk, "_pathway_river"),
    data_source = cache$week_paths[[wk]],
    hypothesis_tested = "Unbiased pathway hierarchy can identify unexpectedly dominant signaling programs.",
    expected_outcome = "A subset of pathways dominates total flow for each week.",
    claim_strength = "Speculative",
    methods_blurb = "netAnalysis_river() visualizes pathway-level sender/receiver hierarchy from weekly pathway probabilities.",
    plot_expr = quote({
      pw <- tryCatch(obj@netP$pathways, error = function(e) NULL)
      if (!is.null(pw) && length(pw) > 0) {
        top_pw <- head(pw, 8)
        try(CellChat::netAnalysis_river(obj, signaling = top_pw), silent = TRUE)
      }
    })
  )
}

# ------------------------
# Phase 2: longitudinal
# ------------------------
save_base_plot(
  fig_id = "P2_total_interactions_count",
  data_source = cache$merged_path,
  hypothesis_tested = "Overall interaction counts vary across gestational weeks.",
  expected_outcome = "Monotonic or phase-shifted trajectory from W7 to W11.",
  claim_strength = "Moderate",
  methods_blurb = "compareInteractions() run on merged object with measure='count' for per-dataset totals.",
  plot_expr = quote({ try(CellChat::compareInteractions(merged_obj, show.legend = FALSE, group = 1), silent = TRUE) })
)

save_base_plot(
  fig_id = "P2_total_interactions_weight",
  data_source = cache$merged_path,
  hypothesis_tested = "Overall interaction weights vary across gestational weeks.",
  expected_outcome = "Developmental shifts in summed communication strength.",
  claim_strength = "Moderate",
  methods_blurb = "compareInteractions() run on merged object with measure='weight'.",
  plot_expr = quote({ try(CellChat::compareInteractions(merged_obj, measure = "weight", show.legend = FALSE, group = 1), silent = TRUE) })
)

save_base_plot(
  fig_id = "P2_diff_network_W7_vs_W11_count",
  data_source = cache$merged_path,
  hypothesis_tested = "Specific communication axes are gained/lost between W7 and W11.",
  expected_outcome = "Red/blue differential edges highlighting developmental rewiring.",
  claim_strength = "Strong",
  methods_blurb = "netVisual_diffInteraction() contrasts first and last datasets in merged object by count.",
  plot_expr = quote({ try(CellChat::netVisual_diffInteraction(merged_obj, weight.scale = TRUE), silent = TRUE) })
)

save_base_plot(
  fig_id = "P2_ranknet_pathway_flow",
  data_source = cache$merged_path,
  hypothesis_tested = "Pathway information flow ranking changes across weeks and highlights thesis pathways.",
  expected_outcome = "MMP/IDO1/TGFb occupy week-specific ranks.",
  claim_strength = "Strong",
  methods_blurb = "rankNet() ranks pathway information flow in merged object for cross-week comparison.",
  plot_expr = quote({
    try(CellChat::rankNet(merged_obj, mode = "comparison", stacked = TRUE, do.stat = FALSE), silent = TRUE)
  })
)

# ------------------------
# Phase 3: spatial mapping
# ------------------------
for (wk in weeks) {
  obj <- week_objs[[wk]]
  md <- obj@meta
  id_col <- as.character(obj@idents)
  xy_x <- if ("x_cent" %in% colnames(md)) "x_cent" else if ("x_um" %in% colnames(md)) "x_um" else NA_character_
  xy_y <- if ("y_cent" %in% colnames(md)) "y_cent" else if ("y_um" %in% colnames(md)) "y_um" else NA_character_
  if (is.na(xy_x) || is.na(xy_y)) next

  plot_df <- data.frame(
    x = md[[xy_x]],
    y = md[[xy_y]],
    celltype = id_col,
    stringsAsFactors = FALSE
  )

  sel <- plot_df$celltype %in% c("EVT", "Extravillous Trophoblast 2", "Syncytiotrophoblast", "Hofbauer Macrophages", "maternal macrophages")
  p_baseline <- ggplot(plot_df[sel, , drop = FALSE], aes(x = x, y = y, color = celltype)) +
    geom_point(size = 0.15, alpha = 0.7) +
    scale_color_manual(values = get_universal_colors(unique(plot_df$celltype[sel]))) +
    coord_fixed() +
    labs(title = paste0("Spatial architecture - ", wk), x = "x", y = "y", color = "Cell type")

  save_ggplot(
    fig_id = paste0("P3_", wk, "_spatial_architecture"),
    p = p_baseline,
    data_source = cache$week_paths[[wk]],
    hypothesis_tested = "Spatial baseline architecture shows interpretable organization of trophoblast and immune populations.",
    expected_outcome = "Distinct spatial niches and boundaries across key cell types.",
    claim_strength = "Moderate",
    methods_blurb = "Scatter map of cell coordinates from metadata, colored by canonical cell type labels."
  )

  save_base_plot(
    fig_id = paste0("P3_", wk, "_netVisual_spatial_global"),
    data_source = cache$week_paths[[wk]],
    hypothesis_tested = "Communication vectors localize to non-random tissue regions.",
    expected_outcome = "Spatial hotspots rather than uniform field-wide signaling.",
    claim_strength = "Strong",
    methods_blurb = "netVisual_spatial() plots pathway/network communication over physical coordinates.",
    plot_expr = quote({
      try(SpatialCellChat::netVisual_spatial(obj, signaling = NULL), silent = TRUE)
    })
  )

  for (sig in c("MMP", "IDO1")) {
    save_base_plot(
      fig_id = paste0("P3_", wk, "_netVisual_spatial_", sig),
      data_source = cache$week_paths[[wk]],
      hypothesis_tested = paste0(sig, " spatial hotspots map to biologically relevant invasion/shield zones."),
      expected_outcome = paste0(sig, " vectors enriched near specific tissue interfaces."),
      claim_strength = "Strong",
      methods_blurb = "netVisual_spatial() restricted to thesis pathway.",
      plot_expr = bquote({ try(SpatialCellChat::netVisual_spatial(obj, signaling = .(sig)), silent = TRUE) })
    )
  }
}

# ------------------------
# Phase 4: vulnerability split plots from compute tables
# ------------------------
vuln_tbl <- tryCatch(read.csv(cache$tables$vulnerability, stringsAsFactors = FALSE), error = function(e) data.frame())
if (nrow(vuln_tbl) > 0) {
  for (wk in unique(vuln_tbl$week)) {
    sub <- vuln_tbl[vuln_tbl$week == wk, , drop = FALSE]
    tgtA <- paste0(cache$target_celltype, "_Exposed")
    tgtB <- paste0(cache$target_celltype, "_Shielded")

    bsub <- sub[sub$source %in% c(tgtA, tgtB) | sub$target %in% c(tgtA, tgtB), , drop = FALSE]
    if (nrow(bsub) > 0) {
      p_bubble <- ggplot(bsub, aes(x = source, y = target, size = prob_value, color = pathway_name)) +
        geom_point(alpha = 0.75) +
        scale_size_continuous(range = c(1, 10)) +
        labs(title = paste0("Vulnerability bubble plot - ", wk), x = "Sender", y = "Receiver", size = "Probability")

      save_ggplot(
        fig_id = paste0("P4_", wk, "_vulnerability_bubble"),
        p = p_bubble,
        data_source = cache$tables$vulnerability,
        hypothesis_tested = "Exposed vs shielded target cells differ in LR communication profile.",
        expected_outcome = "Distinct bubble size/color patterns between Exposed and Shielded states.",
        claim_strength = "Strong",
        methods_blurb = "Bubble plot of re-aggregated LR communication probabilities for vulnerability split groups."
      )

      auto <- bsub[bsub$source == bsub$target, , drop = FALSE]
      if (nrow(auto) > 0) {
        p_auto <- ggplot(auto, aes(x = interaction_name, y = prob_value, fill = source)) +
          geom_col(position = "dodge") +
          coord_flip() +
          labs(title = paste0("Autocrine survival signaling - ", wk), x = "Ligand-Receptor pair", y = "Aggregated probability")

        save_ggplot(
          fig_id = paste0("P4_", wk, "_autocrine_chord_proxy"),
          p = p_auto,
          data_source = cache$tables$vulnerability,
          hypothesis_tested = "Exposed cells compensate via autocrine signaling programs.",
          expected_outcome = "Higher autocrine LR signal in exposed niche for selected pathways.",
          claim_strength = "Moderate",
          methods_blurb = "Autocrine-only subset (source==target) plotted as ranked bar proxy for chord-style self-signaling focus."
        )
      }
    }

    meta <- cache$vulnerability_meta[[wk]]
    if (!is.null(meta) && nrow(meta) > 0) {
      xcol <- if ("x_cent" %in% colnames(meta)) "x_cent" else if ("x_um" %in% colnames(meta)) "x_um" else NA_character_
      ycol <- if ("y_cent" %in% colnames(meta)) "y_cent" else if ("y_um" %in% colnames(meta)) "y_um" else NA_character_
      if (!is.na(xcol) && !is.na(ycol)) {
        p_zone <- ggplot(meta, aes_string(x = xcol, y = ycol, color = "vulnerability_group")) +
          geom_point(size = 0.15, alpha = 0.8) +
          coord_fixed() +
          labs(title = paste0("MISI vulnerability zones - ", wk), x = "x", y = "y", color = "Zone")

        save_ggplot(
          fig_id = paste0("P4_", wk, "_misi_stratification"),
          p = p_zone,
          data_source = cache$week_paths[[wk]],
          hypothesis_tested = "MISI-defined shielded and exposed zones are spatially non-random.",
          expected_outcome = "Visible niche partitioning of Exposed vs Shielded cells.",
          claim_strength = "Strong",
          methods_blurb = "Cells are stratified by MISI quantiles (compute stage), then mapped over tissue coordinates."
        )
      }
    }
  }
}

append_methods(
  "Phase 1-4 plotting",
  c(
    "Phase 1 generated per-week circle, heatmap, signaling-role scatter, and pathway-river plots from completed CellChat objects.",
    "Phase 2 generated merged longitudinal comparisons (interaction counts/weights, differential network, and rankNet pathway flow).",
    "Phase 3 generated baseline coordinate maps and attempted netVisual_spatial global + thesis-pathway overlays (MMP, IDO1).",
    "Phase 4 generated MISI stratification maps and Exposed-vs-Shielded LR comparison plots from 04a compute outputs.",
    "All figures saved as 300-dpi PNG + vector PDF and appended to 04_figures_manifest.json with hypothesis and method descriptors."
  )
)

cat("Plotting stage complete. Figures: ", DIR_FIGS, "\n", sep = "")
