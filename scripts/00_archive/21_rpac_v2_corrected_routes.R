# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  21_rpac_v2_corrected_routes.R                                              ║
# ║  Re-run rPAC with HIF-2α–corrected route definitions (9 routes)             ║
# ║                                                                             ║
# ║  AIM:   Aim 3 — Mechanistic Pathway Activity Profiling                      ║
# ║  THESIS CHAPTER:  Results, Section 3.3                                      ║
# ║                                                                             ║
# ║  FIGURES PRODUCED:                                                          ║
# ║    Fig21a  — FeaturePlot per route (×9)                [Aim 3, Fig 21a]     ║
# ║    Fig21b  — ARS barplot (cohort-level summary)        [Aim 3, Fig 21b]     ║
# ║    Fig21c  — New routes combined panel (×3)            [Aim 3, Fig 21c]     ║
# ║    Fig21d  — Violin plots by condition (×3 new routes) [Aim 3, Fig 21d]     ║
# ║    Fig21e  — Routes × Conditions heatmap               [Aim 3, Fig 21e]     ║
# ║    Fig21f  — Inter-route correlation heatmap            [Aim 3, Fig 21f]     ║
# ║                                                                             ║
# ║  TABLES PRODUCED:                                                           ║
# ║    Tab21a  — Per-cell rPAC scores (cells × routes)     [Aim 3, Tab 21a]     ║
# ║    Tab21b  — Per-cell p-values (cells × routes)        [Aim 3, Tab 21b]     ║
# ║    Tab21c  — Route summary metrics (ARS, PS, n_sig)    [Aim 3, Tab 21c]     ║
# ║                                                                             ║
# ║  KEY CHANGES FROM SCRIPT 19:                                                ║
# ║    • Uses rpac_example_routes_v2() (HIF-2α correction)                      ║
# ║    • 9 routes instead of 6 (3 new: EPAS1 anti-angiogenic, NF-κB→EPAS1      ║
# ║      bridge, EPAS1 Notch/endothelial)                                       ║
# ║    • Route 3 corrected: FLT1 removed from HIF-1α, SLC2A1/PDK1 added        ║
# ║    • All plots saved as both PDF and PNG                                     ║
# ║                                                                             ║
# ║  INPUT:                                                                     ║
# ║    - outputs/objects/sandbox/seu_misi_extremes.qs  (from 17b)               ║
# ║                                                                             ║
# ║  OUTPUT DIRS:                                                               ║
# ║    - outputs/figures/sandbox/21_rpac_v2/                                    ║
# ║    - outputs/tables/sandbox/21_rpac_v2/                                     ║
# ║    - outputs/objects/sandbox/21_rpac_v2/                                    ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

cat("
╔══════════════════════════════════════════════════════╗
║  21 · rPAC v2 — HIF-2α Corrected Route Analysis     ║
╚══════════════════════════════════════════════════════╝
")

# ── 0. Environment ──────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
})

# ── 1. Source pipeline infrastructure ───────────────────────────────────────
source("R/config.R")
source("R/helpers_io.R")
source("R/helpers_plot.R")
source("R/helpers_scores.R")
source("R/helpers_methods.R")
source("R/rpac_core.R")
source("R/rpac_routes_v2.R")

# ── 2. Directory setup (matches pipeline convention) ────────────────────────
obj_dir <- file.path(CFG$dirs$objects, "sandbox", "21_rpac_v2")
fig_dir <- file.path(CFG$dirs$figures, "sandbox", "21_rpac_v2")
tab_dir <- file.path(CFG$dirs$tables,  "sandbox", "21_rpac_v2")
ensure_dir(obj_dir)
ensure_dir(fig_dir)
ensure_dir(tab_dir)

log_file <- file.path(CFG$dirs$logs, "pipeline.log")
log_msg("sandbox/21_rpac_v2_corrected_routes starting.", log_file = log_file)
log_msg("  Object dir: ", obj_dir)
log_msg("  Figure dir: ", fig_dir)
log_msg("  Table dir:  ", tab_dir)

# ── Helper: save both PDF and PNG ──────────────────────────────────────────
save_dual <- function(path_pdf, p, w = 10, h = 7) {
  save_plot(path_pdf, p, w = w, h = h, dpi = 300)
  path_png <- sub("\\.pdf$", ".png", path_pdf)
  ggsave(filename = path_png, plot = p, width = w, height = h,
         units = "in", dpi = 300, bg = "white", limitsize = FALSE)
  log_msg("  Saved: ", basename(path_pdf), " + ", basename(path_png))
  invisible(c(path_pdf, path_png))
}

# ── 3. Load data ────────────────────────────────────────────────────────────
seu_path <- file.path(CFG$dirs$objects, "sandbox", "seu_misi_extremes.qs")
if (!file.exists(seu_path)) {
  stop("seu_misi_extremes.qs not found at: ", seu_path,
       "\nRun 17b_misi_stratified_spotlight.R first.")
}
log_msg("Loading MISI extremes object ...")
seu <- qread(seu_path)
log_msg("  Cells: ", ncol(seu), "  Genes: ", nrow(seu))

# ── 4. Load corrected routes ────────────────────────────────────────────────
routes_v2 <- rpac_example_routes_v2()
log_msg("Loaded ", length(routes_v2), " routes (v2, HIF-2α corrected)")

# Print route summary
route_summary <- rpac_route_summary_v2()
cat("\n── Route Summary ──\n")
print(route_summary[, c("route_id", "primary_tf", "route_type", "n_nodes")])
cat("\n")

# ── 5. Check gene availability for new routes ──────────────────────────────
available_genes <- rownames(seu)
new_route_ids <- c("EPAS1_ER_antiangiogenic", "NFKB_EPAS1_bridge", "EPAS1_ER_notch")
for (rid in new_route_ids) {
  r <- routes_v2[[rid]]
  rg <- route_genes(r)
  present <- rg[rg %in% available_genes]
  missing <- rg[!rg %in% available_genes]
  log_msg("  Route ", rid, ": ", length(present), "/", length(rg),
          " genes present",
          if (length(missing) > 0) paste0(" [missing: ", paste(missing, collapse = ", "), "]") else "")
}

# ── 6. Run rPAC scoring ────────────────────────────────────────────────────
log_msg("Running rPAC scoring with control_condition = 'UI' ...")
rpac_out <- compute_rpac_scores(
  seu               = seu,
  routes            = routes_v2,
  control_condition = "UI",
  condition_col     = CFG$cols$infection,
  c_max             = 1.5,
  u_min_and         = 0.5,
  u_min_or          = 0.2,
  n_null            = 1000L,
  p_threshold       = 0.05,
  seed              = 42L,
  verbose           = TRUE
)

scores  <- rpac_out$scores
pvals   <- rpac_out$pvalues
summary <- rpac_out$summary
log_msg("rPAC scoring complete.")
log_msg("  Score matrix: ", nrow(scores), " cells × ", ncol(scores), " routes")

# ── 7. Save tables ──────────────────────────────────────────────────────────
write.csv(scores,  file.path(tab_dir, "21_rpac_v2_scores.csv"),  row.names = TRUE)
write.csv(pvals,   file.path(tab_dir, "21_rpac_v2_pvalues.csv"), row.names = TRUE)
write.csv(summary, file.path(tab_dir, "21_rpac_v2_summary.csv"), row.names = FALSE)
log_msg("  Saved: 3 CSVs → ", tab_dir)

# ── 8. Embed scores in Seurat ──────────────────────────────────────────────
for (rid in colnames(scores)) {
  col_name <- paste0("rPAC_", rid)
  seu[[col_name]] <- scores[colnames(seu), rid]
}
log_msg("  Embedded ", ncol(scores), " rPAC score columns into Seurat metadata.")

# ── 9. Detect UMAP reduction name ──────────────────────────────────────────
umap_reduction <- if ("umap" %in% names(seu@reductions)) {
  "umap"
} else if ("ref.umap" %in% names(seu@reductions)) {
  "ref.umap"
} else {
  names(seu@reductions)[1]
}
log_msg("  Using UMAP reduction: ", umap_reduction)

# =============================================================================
# 10. FIGURE 21a — FeaturePlots (one per route)
#     [Aim 3, Fig 21a] — Per-cell rPAC activity on UMAP
# =============================================================================
log_msg("Generating Fig21a: FeaturePlots ...")
plot_list <- list()

for (i in seq_along(routes_v2)) {
  rt  <- routes_v2[[i]]
  rid <- rt$route_id
  col_name <- paste0("rPAC_", rid)

  p <- FeaturePlot(seu, features = col_name, reduction = umap_reduction,
                   cols = c("blue", "white", "red"),
                   order = TRUE, pt.size = 0.3, raster = TRUE) +
    ggtitle(paste0(rid, " (", rt$route_type, ")"),
            subtitle = rt$description) +
    theme(plot.title    = element_text(face = "bold", size = 10),
          plot.subtitle = element_text(size = 7, color = "grey40"))

  plot_list[[rid]] <- p

  save_dual(
    file.path(fig_dir, paste0("Fig21a_featureplot_", rid, ".pdf")),
    p, w = 7, h = 6
  )
}

# Combined panel (all 9 routes)
n_routes <- length(plot_list)
n_cols   <- min(3, n_routes)
n_rows   <- ceiling(n_routes / n_cols)

p_combined_all <- wrap_plots(plot_list, ncol = n_cols) +
  plot_annotation(
    title    = "rPAC v2 Route Activity Scores — All 9 Routes (HIF-2α Corrected)",
    subtitle = paste0("Control: UI | Cmax=1.5 | ", ncol(seu), " cells | ", umap_reduction),
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

save_dual(
  file.path(fig_dir, "Fig21a_featureplot_combined_all9.pdf"),
  p_combined_all, w = 7 * n_cols, h = 6 * n_rows
)
log_msg("  Fig21a complete: ", length(plot_list), " FeaturePlots + combined panel.")

# =============================================================================
# 11. FIGURE 21b — ARS Barplot (cohort-level summary)
#     [Aim 3, Fig 21b] — Average Route Scores with PS annotation
# =============================================================================
log_msg("Generating Fig21b: ARS barplot ...")

summary$route_id <- factor(summary$route_id, levels = summary$route_id[order(summary$ARS)])
summary$is_new <- summary$route_id %in% c("EPAS1_ER_antiangiogenic", "NFKB_EPAS1_bridge", "EPAS1_ER_notch")

p_ars <- ggplot(summary, aes(x = route_id, y = ARS, fill = is_new)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0("PS=", round(PS, 2))),
            hjust = ifelse(summary$ARS >= 0, -0.1, 1.1), size = 3) +
  coord_flip() +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick"),
                    labels = c("FALSE" = "Original (v1)", "TRUE" = "New (v2)"),
                    name = "Route origin") +
  labs(title = "rPAC v2 — Average Route Score (ARS)",
       subtitle = paste0("9 routes | ", ncol(seu), " MISI-extreme cells | PS = proportion significant"),
       x = NULL, y = "Average Route Score (ARS)") +
  theme_minimal(base_size = 12) +
  theme(legend.position   = "bottom",
        panel.grid.major.y = element_blank())

save_dual(
  file.path(fig_dir, "Fig21b_ars_barplot_v2.pdf"),
  p_ars, w = 10, h = 6
)
log_msg("  Fig21b complete.")

# =============================================================================
# 12. FIGURE 21c — New Routes Combined Panel (3 EPAS1/bridge routes)
#     [Aim 3, Fig 21c] — Focus on HIF-2α–corrected routes
# =============================================================================
log_msg("Generating Fig21c: New routes combined panel ...")

new_plots <- plot_list[new_route_ids]
if (length(new_plots) == 3) {
  p_new_combined <- wrap_plots(new_plots, ncol = 3) +
    plot_annotation(
      title    = "New HIF-2α Routes — rPAC v2",
      subtitle = "EPAS1 anti-angiogenic | NF-κB→EPAS1 bridge | EPAS1 Notch/endothelial",
      theme = theme(
        plot.title    = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey40")
      )
    )

  save_dual(
    file.path(fig_dir, "Fig21c_new_routes_combined_v2.pdf"),
    p_new_combined, w = 21, h = 6
  )
  log_msg("  Fig21c complete.")
} else {
  log_msg("  WARNING: Could not create Fig21c — expected 3 new routes, got ", length(new_plots))
}

# =============================================================================
# 13. FIGURE 21d — Violin Plots by Condition (new routes only)
#     [Aim 3, Fig 21d] — Distribution of new route scores across conditions
# =============================================================================
log_msg("Generating Fig21d: Violin plots by condition ...")

cond_col <- CFG$cols$infection
if (cond_col %in% colnames(seu@meta.data)) {
  violin_list <- list()
  for (rid in new_route_ids) {
    col_name <- paste0("rPAC_", rid)
    p_violin <- VlnPlot(seu, features = col_name, group.by = cond_col,
                        pt.size = 0, cols = brewer.pal(4, "Set2")) +
      geom_hline(yintercept = 0, lty = 2, color = "grey50") +
      labs(title = rid, y = "rPAC Score") +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold", size = 10))

    violin_list[[rid]] <- p_violin

    save_dual(
      file.path(fig_dir, paste0("Fig21d_violin_", rid, ".pdf")),
      p_violin, w = 7, h = 5
    )
  }

  # Combined violin panel
  p_violin_combined <- wrap_plots(violin_list, ncol = 3) +
    plot_annotation(
      title = "rPAC v2 — New Route Scores by Infection Condition",
      theme = theme(plot.title = element_text(size = 13, face = "bold"))
    )
  save_dual(
    file.path(fig_dir, "Fig21d_violins_combined_new_routes.pdf"),
    p_violin_combined, w = 21, h = 6
  )
  log_msg("  Fig21d complete.")
} else {
  log_msg("  WARNING: Condition column '", cond_col, "' not found — skipping violins.")
}

# =============================================================================
# 14. FIGURE 21e — Routes × Conditions Heatmap
#     [Aim 3, Fig 21e] — Mean rPAC score per route per condition
# =============================================================================
log_msg("Generating Fig21e: Routes × Conditions heatmap ...")

if (cond_col %in% colnames(seu@meta.data)) {
  score_cols <- paste0("rPAC_", sapply(routes_v2, `[[`, "route_id"))
  md <- seu@meta.data[, c(cond_col, score_cols), drop = FALSE]
  md_long <- md %>%
    pivot_longer(cols = all_of(score_cols), names_to = "route", values_to = "score") %>%
    mutate(route = gsub("^rPAC_", "", route))

  heat_df <- md_long %>%
    group_by(!!sym(cond_col), route) %>%
    summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = route, values_from = mean_score)

  heat_mat <- as.matrix(heat_df[, -1])
  rownames(heat_mat) <- heat_df[[1]]

  heat_long <- heat_df %>%
    pivot_longer(-1, names_to = "route", values_to = "mean_score")
  colnames(heat_long)[1] <- "condition"

  p_heat <- ggplot(heat_long, aes(x = route, y = condition, fill = mean_score)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = round(mean_score, 3)), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         name = "Mean\nrPAC") +
    labs(title = "rPAC v2 — Mean Route Activity by Condition",
         subtitle = "Blue = suppressed, Red = activated relative to UI control",
         x = "Route", y = "Condition") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          panel.grid = element_blank())

  save_dual(
    file.path(fig_dir, "Fig21e_heatmap_routes_conditions_v2.pdf"),
    p_heat, w = 12, h = 5
  )
  log_msg("  Fig21e complete.")
} else {
  log_msg("  WARNING: Skipping heatmap — no condition column.")
}

# =============================================================================
# 15. FIGURE 21f — Inter-route Correlation Heatmap
#     [Aim 3, Fig 21f] — Pairwise Pearson correlation of route scores
# =============================================================================
log_msg("Generating Fig21f: Inter-route correlation heatmap ...")

cor_mat <- cor(scores, use = "pairwise.complete.obs")

cor_long <- as.data.frame(as.table(cor_mat))
colnames(cor_long) <- c("Route1", "Route2", "r")

p_cor <- ggplot(cor_long, aes(x = Route1, y = Route2, fill = r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(r, 2)), size = 2.8) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  labs(title = "rPAC v2 — Inter-Route Correlation",
       subtitle = "Pairwise Pearson correlation across all scored cells") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        panel.grid = element_blank())

save_dual(
  file.path(fig_dir, "Fig21f_correlation_heatmap_v2.pdf"),
  p_cor, w = 9, h = 8
)
log_msg("  Fig21f complete.")

# =============================================================================
# 16. SAVE UPDATED SEURAT OBJECT
# =============================================================================
qsave(seu, file.path(obj_dir, "seu_misi_extremes_rpac_v2.qs"), preset = "high")
log_msg("  Saved: seu_misi_extremes_rpac_v2.qs (with rPAC v2 metadata)")

# =============================================================================
# 17. FIGURE LEGEND — Detailed, thesis-ready
# =============================================================================
write_legend(
  fig_id  = "Fig21",
  title   = "rPAC v2 Route Activity Scores — HIF-2α Corrected (MISI Extreme Cells)",
  hypothesis = paste(
    "The Fn Ethanolamine → MegL Toxic Switch hypothesis predicts that specific",
    "TF-centric signaling routes are differentially activated between MISI-high",
    "(vulnerable) and MISI-low (resilient) cells. The HIF-2α correction",
    "(Sasagawa et al. 2018, 2021; Colson et al. 2023) establishes that EPAS1",
    "(HIF-2α), not HIF1A (HIF-1α), drives FLT1/sFLT1 upregulation in",
    "placental trophoblasts. Three new EPAS1-based routes test the",
    "infection-to-preeclampsia bridge hypothesis: (1) EPAS1 anti-angiogenic",
    "(HIF-2α/ARNT → FLT1/PGF), (2) NF-κB→EPAS1 bridge (RELA → EPAS1 →",
    "FLT1/PGF), and (3) EPAS1 Notch/endothelial (EPAS1 → DLL4/ANGPT2/NOTCH1)."
  ),
  methods = paste(
    "The rPAC algorithm (Joshi et al., Methods 2022) was applied to MISI-extreme",
    "cells (Q1 vs Q3 of Target_MISI_Score, from script 17b). Nine TF-centric",
    "routes (v2) were defined: 6 original routes with FLT1 removed from HIF-1α",
    "(Route 3 now targets SLC2A1/PDK1 glycolysis), plus 3 new EPAS1-based routes.",
    "For each cell, gene-level log2 fold-changes were computed relative to",
    "uninfected (UI) control cells (Eq 1, α=1). Node expected values propagated",
    "from the primary TF (Eq 2), bundle nodes evaluated with AND/OR thresholds",
    "(Eq 3, Umin_AND=0.5, Umin_OR=0.2), contributions clamped at Cmax=1.5 (Eq 4),",
    "and activity scores computed as mean weighted node contributions (Eq 5).",
    "Statistical significance assessed via 1000 null samples from N(0,1) (Eq 6).",
    "Summary metrics PS and ARS computed per route (Eqs 7–8). All plots saved as",
    "both PDF (vector) and PNG (300 dpi raster)."
  ),
  readout = paste(
    "Fig21a: FeaturePlots show per-cell rPAC scores on UMAP using a diverging",
    "blue–white–red palette (blue = route suppressed, white = neutral, red =",
    "activated). Fig21b: ARS barplot shows cohort-level average activation per",
    "route with PS annotated. New routes highlighted in red, original in blue.",
    "Fig21c: Combined panel of the 3 new HIF-2α routes for direct comparison.",
    "Fig21d: Violin plots of new route scores split by infection condition",
    "(dashed line = 0). Fig21e: Heatmap of mean rPAC score per route per",
    "condition — tests whether infection type modulates route activation.",
    "Fig21f: Pearson correlation heatmap reveals co-regulation structure",
    "(e.g., NF-κB and EPAS1 routes should positively correlate if the bridge",
    "hypothesis holds). Tab21a–c: Full score matrices and summary statistics."
  ),
  interpretation_template = paste(
    "If EPAS1_ER_antiangiogenic shows elevated ARS in Fn-infected cells, this",
    "supports the hypothesis that Fn infection drives HIF-2α–mediated sFLT1",
    "upregulation, a hallmark of preeclampsia. Positive correlation between",
    "NFKB_SR_innate and NFKB_EPAS1_bridge would confirm the NF-κB→EPAS1→FLT1",
    "super-route, linking innate immune activation to anti-angiogenic crisis.",
    "EPAS1_ER_notch activation alongside EPAS1_ER_antiangiogenic would indicate",
    "a broader endothelial disruption program. Negative EA_nutrient_axis scores",
    "in infected conditions would indicate EA liberation exceeding sequestration,",
    "consistent with increased EA availability fueling Fn colonization. The",
    "HIF-1α route (now metabolic-only: SLC2A1/PDK1/VEGFA) should decouple from",
    "FLT1 upregulation, confirming the dual-HIF model."
  )
)
log_msg("  Saved: Fig21 legend → outputs/legends/")

# =============================================================================
# 18. METHODS DRAFT — Thesis-ready paragraph
# =============================================================================
if (exists("append_methods_draft", mode = "function")) {
  append_methods_draft(
    script_name    = "sandbox/21_rpac_v2_corrected_routes.R",
    hypothesis     = paste(
      "Route-based pathway analysis with HIF-2α correction reveals that EPAS1,",
      "not HIF1A, drives the anti-angiogenic FLT1/sFLT1 axis in infection-vulnerable",
      "placental cells, and that an NF-κB→EPAS1→FLT1 super-route links innate",
      "immune activation to preeclampsia-like vascular disruption."
    ),
    method_details = list(
      "Algorithm: rPAC (Joshi et al., Methods 2022) adapted for single-cell data",
      "Routes: 9 hand-curated TF-centric routes (v2, HIF-2α corrected)",
      "  - Original 6 with FLT1 removed from HIF-1α (Route 3 → SLC2A1/PDK1/VEGFA)",
      "  - Route 7 NEW: EPAS1_ER_antiangiogenic (EPAS1/ARNT → FLT1/PGF)",
      "  - Route 8 NEW: NFKB_EPAS1_bridge (RELA → EPAS1 → FLT1/PGF)",
      "  - Route 9 NEW: EPAS1_ER_notch (EPAS1 → DLL4/ANGPT2/NOTCH1)",
      "Control: mean expression of UI (uninfected) cells as baseline",
      paste0("Hyperparameters: Cmax=1.5, Umin_AND=0.5, Umin_OR=0.2, n_null=1000, seed=42"),
      paste0("Cells scored: ", ncol(seu), " (MISI extreme subset from 17b)")
    ),
    rationale      = paste(
      "The original rPAC routes (v1, script 19) assigned FLT1 to the HIF-1α effector",
      "route. However, Sasagawa et al. (2018, Mol. Cell. Endocrinol.) demonstrated",
      "that EPAS1 (HIF-2α), not HIF1A (HIF-1α), directly transactivates the FLT1",
      "promoter in human trophoblasts via an HIF-2α/ARNT complex (confirmed by",
      "Sasagawa et al. 2021, Hum. Cell). In vivo validation by Colson et al. (2023,",
      "Hypertension) showed that the oral HIF-2α inhibitor PT2385 prevented sFLT1",
      "elevation and hypertension in the RUPP rat model. This correction strengthens",
      "the thesis by establishing a dual-HIF model: HIF-1α drives metabolic",
      "reprogramming (glycolysis, VEGFA), while HIF-2α drives the anti-angiogenic",
      "crisis (FLT1/sFLT1 ↑, PGF ↓). The NF-κB→EPAS1 bridge route tests whether",
      "pathogen-driven inflammation directly activates this preeclampsia pathway."
    ),
    citations      = list(
      "Joshi et al., Methods 198 (2022): 76-87. doi:10.1016/j.ymeth.2021.10.002",
      "Sasagawa et al., Mol. Cell. Endocrinol. 472 (2018): 68-75. doi:10.1016/j.mce.2017.12.003",
      "Sasagawa et al., Hum. Cell 34 (2021): 2-9. doi:10.1007/s13577-020-00432-5",
      "Colson et al., Hypertension 80 (2023): e26-e37. doi:10.1161/HYPERTENSIONAHA.122.20149",
      "Skuli et al., Blood 119 (2012): 2304-2314. doi:10.1182/blood-2011-09-380410"
    )
  )
  log_msg("  Saved: methods draft → outputs/reports/THESIS_METHODS_DRAFT.md")
}

# =============================================================================
# 19. COMPLETION
# =============================================================================
log_msg("sandbox/21_rpac_v2_corrected_routes done.", log_file = log_file)

cat("
╔══════════════════════════════════════════════════════════════════════════╗
║  Script 21 complete!                                                    ║
║                                                                         ║
║  Aim 3 — Mechanistic Pathway Activity Profiling (HIF-2α Corrected)     ║
║                                                                         ║
║  Figures (PDF + PNG):                                                   ║
║    Fig21a — FeaturePlot per route (×9) + combined panel                 ║
║    Fig21b — ARS barplot                                                 ║
║    Fig21c — New routes combined panel                                   ║
║    Fig21d — Violin plots by condition (×3 new routes) + combined        ║
║    Fig21e — Routes × Conditions heatmap                                 ║
║    Fig21f — Inter-route correlation heatmap                             ║
║                                                                         ║
║  Tables:                                                                ║
║    Tab21a — Per-cell scores  (21_rpac_v2_scores.csv)                    ║
║    Tab21b — Per-cell p-values (21_rpac_v2_pvalues.csv)                  ║
║    Tab21c — Route summary    (21_rpac_v2_summary.csv)                   ║
║                                                                         ║
║  Object: seu_misi_extremes_rpac_v2.qs                                   ║
║  Legend: outputs/legends/Fig21_legend.md                                 ║
║  Methods: outputs/reports/THESIS_METHODS_DRAFT.md (appended)            ║
╚══════════════════════════════════════════════════════════════════════════╝
")