# 01 Methods and Data Provenance Log

This file is auto-appended by `01_preprocess_harmony_embeddings.R`.
Each run records the exact modalities/files used for each week.

## Run: 2026-04-15 10:58:34 EDT
- Script: 01_preprocess_harmony_embeddings v1.0.4
- Week: W7
- Modalities used:
  - Expression matrix: data/raw/zenodo_spatial/STARmap-ISS_sample_W7_imputed_expression.csv
  - Spatial coordinates (spots metadata): data/raw/zenodo_spatial/STARmap-ISS_sample_W7_spots_metadata.csv
  - Cell metadata (optional): data/raw/zenodo_spatial/STARmap-ISS_sample_W7_cell_metadata.csv
- Method notes:
  - Expression + spots metadata were sourced from split raw directories (Broad + Zenodo).
  - Spatial coordinates are harmonized to x_um/y_um before Seurat construction.
  - SCTransform + Harmony embeddings are generated downstream in this script.

