# 05F and 05H description

This file describes the outputs and interpretation guidelines for:
- scripts/05_spatial/05C_permissiveness_score_maps.R (05C)
- scripts/05_spatial/05H_spatial_permissiveness_panels_global.R (05H)

### Purpose
05C: compute dataset-relative permissiveness per cell: permissiveness = z_MMP + z_Tolerance - z_NK + z_Ethanolamine (z computed within each dataset).
05H: pooled global calibration across datasets to allow cross-sample/time comparisons: z_global computed across pooled cells.

### Figure types created by 05H (per-sample)
1. 3-panel global: (A) image/background, (B) KDE of top 10% permissive cells, (C) global-permissiveness overlay (viridis).
2. 3-panel white->red: (A) image, (B) same KDE, (C) white->red overlay emphasizing positive permissiveness.
3. Cell-type highlighted variant: (A) image, (B) same KDE, (C) overlay with highlighted biologically-important cell types (filled circles) on top of permissiveness.
4. Protected/diff panels: (A) top-density, (B) bottom-density (protected), (C) hotspot - protected difference map.

### Interpretation guidance (for figure legends / methods)
- Hotspot defined as KDE of the top 10% permissive cells (dataset-pooled percentile when using permissiveness_global).
- Protected regions: KDE of the bottom 10% permissive cells.
- Difference map = hotspot_density - protected_density (positive = permissive, negative = protected).
- Include a note whether permissiveness is dataset-relative (05C) or globally-calibrated (05H).

### Files produced
- output/figures/05_spatial/<dataset>/<sample>/ : PNG/PDF figures and caption .txt files
- output/tables/05_spatial/ : permissiveness_global_allcells.csv, summary_table_with_highlights.csv, PROCESS_LOG.txt, per-sample hotspot/protected CSVs and meta.json files

### Methods notes
See the Methods bullets saved to scripts/05_spatial/METHODS_05C_05H.md for full reproducible details.
