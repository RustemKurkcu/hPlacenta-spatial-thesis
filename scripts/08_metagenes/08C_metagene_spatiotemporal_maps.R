# ======================================================================
# scripts/08_metagenes/08C_metagene_spatiotemporal_maps.R
# Apply metagene modules to Slide-tags and STARmap; plot in space/time.
#
# PURPOSE:
#   Map discovered metagene modules onto spatial data to identify
#   spatiotemporal patterns of biological programs (remodeling,
#   tolerance, metabolism, etc.)
#
# METHOD:
#   1. Load metagene modules from 08B
#   2. Score modules in Slide-tags (spatial)
#   3. Score modules in STARmap (spatial, using imputed assay)
#   4. Generate spatial maps for top modules
#   5. Generate timecourse plots by week and cell type
#   6. Identify "permissiveness windows"
#
# BIOLOGICAL RATIONALE:
#   Metagene modules represent coordinated biological programs.
#   Mapping them spatially reveals:
#   - "Remodeling Highway" zones (high MMP/ECM module)
#   - "Immune privilege" zones (high tolerance module)
#   - "BAIT signal" zones (high ethanolamine module)
#   - Temporal dynamics (week 8-10 vulnerability window)
#
# OUTPUT:
#   - Spatial maps of top metagene modules
#   - Timecourse plots (week x cell type)
#   - Module variance tables
#   - Permissiveness window identification
#
# USAGE:
#   source("scripts/08_metagenes/08C_metagene_spatiotemporal_maps.R")
#   (Requires 08B to be run first)
#
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

source("config/config.R")
source("scripts/R/utils.R")

# Load enhanced utilities if available
if (file.exists("scripts/R/utils_enhanced.R")) {
  source("scripts/R/utils_enhanced.R")
}

ensure_dir(DIR_TABLES); ensure_dir(DIR_FIGURES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "08C_metagene_spatiotemporal_maps.log")

# Create metagenes subdirectories
dir.create(file.path(DIR_TABLES, "metagenes"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(DIR_FIGURES, "metagenes"), showWarnings = FALSE, recursive = TRUE)

log_msg("Starting metagene spatiotemporal mapping...", logfile)

# ---- Load module gene lists ----
mod_path <- file.path(DIR_TABLES, "metagenes", "metagene_modules_gene_lists.csv")
if (!file.exists(mod_path)) {
  stop("Missing module list: run 08B first: ", mod_path)
}

log_msg(paste0("Loading modules from: ", mod_path), logfile)
mods <- read.csv(mod_path) %>% as_tibble()
modules <- split(mods$gene, mods$module)

log_msg(sprintf("Loaded %d modules", length(modules)), logfile)

# ---- Load mapped objects ----
sl_path <- file.path(DIR_OBJS, "slidetags_mapped_to_multiome.rds")
st_path <- file.path(DIR_OBJS, "starmap_mapped_to_multiome.rds")

if (!file.exists(sl_path)) stop("Missing: ", sl_path)
if (!file.exists(st_path)) stop("Missing: ", st_path)

log_msg("Loading Slide-tags...", logfile)
sl <- readRDS(sl_path)
sl <- ensure_week_column(sl, COL_WEEK_CANDIDATES)
sl <- ensure_spatial_coords(sl, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)

log_msg("Loading STARmap...", logfile)
st <- readRDS(st_path)
st <- ensure_week_column(st, COL_WEEK_CANDIDATES)
st <- ensure_spatial_coords(st, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)

# ---- Score modules on Slide-tags ----
log_msg("Scoring modules on Slide-tags...", logfile)

DefaultAssay(sl) <- if ("RNA" %in% names(sl@assays)) "RNA" else DefaultAssay(sl)
if (!("data" %in% Layers(sl[[DefaultAssay(sl)]]))) {
  sl <- NormalizeData(sl, verbose=FALSE)
}

# Use safe scoring if available
if (exists("add_modules_from_list_enhanced")) {
  sl <- add_modules_from_list_enhanced(sl, modules, prefix = "MG_", 
                                       assay = DefaultAssay(sl), 
                                       seed = SEED, verbose = FALSE)
} else {
  sl <- AddModuleScore(sl, features = modules, name="MG_", 
                       assay=DefaultAssay(sl), search=FALSE)
}

log_msg("  Module scoring complete for Slide-tags", logfile)

# ---- Score modules on STARmap (use imputed if available) ----
log_msg("Scoring modules on STARmap...", logfile)

# Use imputed assay for better gene coverage
assay_st <- if ("imputed" %in% names(st@assays)) "imputed" else "RNA_raw"
log_msg(paste0("  Using assay: ", assay_st), logfile)

DefaultAssay(st) <- assay_st
if (assay_st == "RNA_raw" && !("data" %in% Layers(st[["RNA_raw"]]))) {
  st <- NormalizeData(st, verbose=FALSE)
}

# Use safe scoring if available
if (exists("add_modules_from_list_enhanced")) {
  st <- add_modules_from_list_enhanced(st, modules, prefix = "MG_", 
                                       assay = assay_st, 
                                       seed = SEED, verbose = FALSE)
} else {
  st <- AddModuleScore(st, features = modules, name="MG_", 
                       assay=assay_st, search=FALSE)
}

log_msg("  Module scoring complete for STARmap", logfile)

# ---- Identify top modules by variance ----
mg_cols_sl <- grep("^MG_", colnames(sl@meta.data), value=TRUE)

var_df <- tibble(
  module = mg_cols_sl,
  variance = sapply(mg_cols_sl, function(cn) var(sl@meta.data[[cn]], na.rm=TRUE))
) %>%
  arrange(desc(variance))

write.csv(var_df, 
          file.path(DIR_TABLES, "metagenes", "metagenes_slidetags_module_variances.csv"), 
          row.names=FALSE)

top_mods <- head(var_df$module, 6)
log_msg(sprintf("Top 6 modules by variance: %s", paste(top_mods, collapse=", ")), logfile)

# ---- Spatial plots for Slide-tags ----
log_msg("Creating spatial maps for top modules...", logfile)

md_sl <- sl@meta.data %>% tibble::rownames_to_column("cell")

for (m in top_mods) {
  log_msg(paste0("  Plotting: ", m), logfile)
  
  p <- scatter_spatial(md_sl, x = "spatial_x_use", y = "spatial_y_use", 
                       color = m, 
                       title = paste0("Slide-tags Spatial: ", m),
                       continuous = TRUE) +
    scale_color_viridis_c(option = "magma") +
    theme(legend.position = "right")
  
  save_plot(p, 
            file.path(DIR_FIGURES, "metagenes", paste0("slidetags_spatial_", m, ".png")), 
            w=10, h=8)
}

# ---- Spatial plots for STARmap (top 3 modules) ----
log_msg("Creating spatial maps for STARmap (top 3 modules)...", logfile)

md_st <- st@meta.data %>% tibble::rownames_to_column("cell")
top_mods_st <- head(top_mods, 3)

for (m in top_mods_st) {
  if (m %in% colnames(md_st)) {
    log_msg(paste0("  Plotting: ", m), logfile)
    
    p <- scatter_spatial(md_st, x = "spatial_x_use", y = "spatial_y_use",
                         color = m,
                         title = paste0("STARmap Spatial: ", m),
                         continuous = TRUE) +
      scale_color_viridis_c(option = "magma") +
      theme(legend.position = "right")
    
    save_plot(p,
              file.path(DIR_FIGURES, "metagenes", paste0("starmap_spatial_", m, ".png")),
              w=10, h=8)
  }
}

# ---- Timecourse summary for Slide-tags ----
log_msg("Creating timecourse summaries...", logfile)

ct_col <- if (COL_PRED_CELLTYPE %in% colnames(sl@meta.data)) {
  COL_PRED_CELLTYPE
} else {
  "cell_type_spatial"
}

tc <- bind_rows(lapply(top_mods, function(m) {
  sl@meta.data %>%
    mutate(
      module = m, 
      score = .data[[m]], 
      celltype = .data[[ct_col]]
    ) %>%
    group_by(week, celltype, module) %>%
    summarize(
      mean_score = mean(score, na.rm=TRUE),
      median_score = median(score, na.rm=TRUE),
      sd_score = sd(score, na.rm=TRUE),
      n = n(), 
      .groups="drop"
    )
}))

write.csv(tc, 
          file.path(DIR_TABLES, "metagenes", "slidetags_metagene_timecourse_by_week_celltype.csv"), 
          row.names=FALSE)

# ---- Timecourse plot ----
log_msg("Creating timecourse plot...", logfile)

p_tc <- tc %>%
  filter(module %in% head(top_mods, 4)) %>%
  ggplot(aes(x = week, y = mean_score, color = celltype, group = celltype)) +
  geom_line(alpha = 0.7, size = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~module, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Metagene Module Dynamics Over Gestational Time",
    subtitle = "Top 4 modules by variance",
    x = "Gestational Week",
    y = "Mean Module Score",
    color = "Cell Type"
  )

save_plot(p_tc,
          file.path(DIR_FIGURES, "metagenes", "metagenes_timecourse_top4.png"),
          w=12, h=10)

# ---- Identify "permissiveness window" ----
log_msg("Identifying permissiveness windows...", logfile)

# Look for modules that peak in weeks 8-10
window_analysis <- tc %>%
  filter(week >= 8, week <= 10) %>%
  group_by(module, celltype) %>%
  summarize(
    mean_score_window = mean(mean_score, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_score_window))

write.csv(window_analysis,
          file.path(DIR_TABLES, "metagenes", "permissiveness_window_week8-10.csv"),
          row.names=FALSE)

log_msg("Top module-celltype combinations in week 8-10 window:", logfile)
for (i in 1:min(10, nrow(window_analysis))) {
  log_msg(sprintf("  %d. %s in %s (score = %.3f)",
                  i, window_analysis$module[i], 
                  window_analysis$celltype[i],
                  window_analysis$mean_score_window[i]), logfile)
}

log_msg("08C complete.", logfile)

cat("\n")
cat(strrep("=", 70), "\n")
cat("METAGENE SPATIOTEMPORAL MAPPING COMPLETE\n")
cat(strrep("=", 70), "\n")
cat(sprintf("Spatial maps: %s/metagenes/slidetags_spatial_MG_*.png\n", DIR_FIGURES))
cat(sprintf("              %s/metagenes/starmap_spatial_MG_*.png\n", DIR_FIGURES))
cat(sprintf("Timecourse: %s/metagenes/slidetags_metagene_timecourse_by_week_celltype.csv\n", DIR_TABLES))
cat(sprintf("            %s/metagenes/metagenes_timecourse_top4.png\n", DIR_FIGURES))
cat(sprintf("Window analysis: %s/metagenes/permissiveness_window_week8-10.csv\n", DIR_TABLES))
cat(strrep("=", 70), "\n\n")

cat("BIOLOGICAL INTERPRETATION:\n")
cat("1. High-variance modules represent coordinated biological programs\n")
cat("2. Spatial clustering of modules indicates functional zones\n")
cat("3. Temporal peaks in weeks 8-10 suggest vulnerability windows\n")
cat("4. Cell type-specific module expression reveals susceptibility patterns\n\n")

cat("GRANT/THESIS APPLICATIONS:\n")
cat("- Use spatial maps to show 'Remodeling Highway' zones\n")
cat("- Use timecourse to demonstrate 'permissiveness window'\n")
cat("- Connect high-scoring zones to infection susceptibility\n")
cat("- Validate with known biological programs (ECM, immune, metabolic)\n\n")

# ======================================================================
# END OF SCRIPT
# ======================================================================