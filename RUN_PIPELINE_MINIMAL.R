# Run the minimal, stable pipeline end-to-end.
# This script is intentionally conservative: it avoids heavyweight dependencies
# and prefers deterministic outputs.

source("config/config.R")

message("\n=== RUN_PIPELINE_MINIMAL ===")
message("Project root: ", getwd())
message("Seurat version: ", as.character(utils::packageVersion("Seurat")))
message("SeuratObject version: ", as.character(utils::packageVersion("SeuratObject")))

source("scripts/02_preprocess/02A_preprocess_multiome_reference.R")

source("scripts/03_mapping/03A_map_slidetags_to_multiome.R")
source("scripts/03_mapping/03B_map_starmap_to_multiome.R")
source("scripts/03_mapping/03C_harmonize_celltype_labels.R")

source("scripts/04_timecourse/04A_gene_of_interest_timecourse.R")

source("scripts/05_spatial/05A_spatial_overview_plots.R")
source("scripts/05_spatial/05B_neighborhood_enrichment.R")
source("scripts/05_spatial/05C_permissiveness_score_maps.R")

source("scripts/06_cell_communication/06B_simple_LR_scoring.R")

source("scripts/07_export/07A_export_shareable_outputs.R")

message("=== DONE (MINIMAL) ===\n")
