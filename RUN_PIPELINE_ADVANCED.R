# Run the full “cutting-edge” pipeline.
#
# This run adds optional modules (CellChat, NicheNet, metagene discovery,
# spatial neighborhood DE) that may require extra packages and runtime.
# The script will *skip* modules whose dependencies are not installed.

source("config/config.R")

message("\n=== RUN_PIPELINE_ADVANCED ===")
message("Project root: ", getwd())
message("Seurat version: ", as.character(utils::packageVersion("Seurat")))

# --- Core pipeline (same as minimal) ---
source("scripts/02_preprocess/02A_preprocess_multiome_reference.R")
source("scripts/03_mapping/03A_map_slidetags_to_multiome.R")
source("scripts/03_mapping/03B_map_starmap_to_multiome.R")
source("scripts/03_mapping/03C_harmonize_celltype_labels.R")
source("scripts/04_timecourse/04A_gene_of_interest_timecourse.R")
source("scripts/04_timecourse/04B_immune_subsets_refinement.R")
source("scripts/05_spatial/05A_spatial_overview_plots.R")
source("scripts/05_spatial/05B_neighborhood_enrichment.R")
source("scripts/05_spatial/05C_permissiveness_score_maps.R")
source("scripts/06_cell_communication/06B_simple_LR_scoring.R")
source("scripts/07_export/07A_export_shareable_outputs.R")

# --- Advanced / optional modules ---

if (file.exists("scripts/08_metagene/08A_housekeeping_qc.R")) {
  source("scripts/08_metagene/08A_housekeeping_qc.R")
}

if (file.exists("scripts/08_metagene/08B_metagene_nmf.R")) {
  source("scripts/08_metagene/08B_metagene_nmf.R")
}

if (file.exists("scripts/05_spatial/05D_neighborhood_DE.R")) {
  source("scripts/05_spatial/05D_neighborhood_DE.R")
}

if (file.exists("scripts/06_cell_communication/06C_cellchat_optional.R")) {
  source("scripts/06_cell_communication/06C_cellchat_optional.R")
}

message("\n=== ADVANCED RUN COMPLETE ===")