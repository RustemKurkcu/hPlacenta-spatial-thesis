# Run order

Recommended order (from project root):

1. Preprocess reference:
   Rscript scripts/02_preprocess/02A_preprocess_multiome_reference.R

2. Map Slide-tags:
   Rscript scripts/03_mapping/03A_map_slidetags_to_multiome.R

3. Map STARmap:
   Rscript scripts/03_mapping/03B_map_starmap_to_multiome.R

4. Harmonize labels (preserve author + predicted; add conservative/refined/fallback):
   Rscript scripts/03_mapping/03C_harmonize_celltype_labels.R

4. Timecourse (genes + module scores):
   Rscript scripts/04_timecourse/04A_gene_of_interest_timecourse.R

5. Immune subtype refinement (macrophage/NK; lightweight score-based):
   Rscript scripts/04_timecourse/04B_immune_subsets_refinement.R

Optional: Greenbaum-style gene coordination score:
   Rscript scripts/04_timecourse/04C_gene_coordination_score.R

5. Spatial overview:
   Rscript scripts/05_spatial/05A_spatial_overview_plots.R

6. Neighborhood enrichment:
   Rscript scripts/05_spatial/05B_neighborhood_enrichment.R

7. Permissiveness score maps (Remodeling + Tolerance - NK + EA):
   Rscript scripts/05_spatial/05C_permissiveness_score_maps.R

7. Cell–cell interaction scoring:
   - Lightweight (recommended): Rscript scripts/06_cell_communication/06B_simple_LR_scoring.R
   - Optional heavy (CellChat):  Rscript scripts/06_cell_communication/06A_cellchat_spatial_constrained.R
     (requires RUN_OPTIONAL_HEAVY <- TRUE in config/config.R)

8. Export a shareable package (no large Seurat objects):
   Rscript scripts/07_export/07A_export_shareable_outputs.R
