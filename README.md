# Spatial Architecture of the Human Placenta: Infection & Non-Infection Pipelines

![R Version](https://img.shields.io/badge/R-4.3%2B-blue)
![Seurat](https://img.shields.io/badge/Seurat-v5-success)
![Status](https://img.shields.io/badge/Status-Active_Development-orange)

## 📌 Overview
This repository contains the code and spatial analysis pipelines for my thesis, focusing on the spatial architecture of the human placenta. It merges two previously separate analytical tracks—the "non-infection" baseline pipeline and the "infection-informed" pipeline—to model how localized pathobiont activity (such as *Fusobacterium nucleatum*) alters placental microenvironments.

### Key Features
* **Memory-Safe Spatial Data Loading:** Optimized for large `.qs` and `.rds` spatial matrices.
* **MISI v2 Scoring:** Spatial mapping of the Microenvironment Infection Susceptibility Index (MISI) using UCell.
* **Toxic Blast Zone Detection:** Multi-radius Moran's I sweeps to identify spatially autocorrelated hotspots of pathology.
* **Spatial CellChat:** Ligand-receptor interaction analysis weighted by adjacency and spot deconvolution.

## 📂 Repository Structure
```text
├── data/                  # Raw and processed spatial data (Note: Managed via Git LFS)
├── R/                     # Helper functions (I/O, spatial checks, MISI math)
├── scripts/               # Sequential pipeline scripts
│   ├── 00_config.R        # Central parameters and radii configurations
│   ├── 01_load_make_seurat.R 
│   ├── 05_spatial_contract_check.R
│   ├── 16_compute_misi_v2_spatial.R
│   └── 17_spatial_cellchat.R
├── outputs/               # Generated figures, hotspot polygons (.gpkg), and tables
└── README.md              # Project documentation

## Repository sync note

Date: 2026-04-09 16:26:37
Repository synced from the current local tree at: E:\hPlacenta-architecture

Rules used:
- Files <= 100 MiB added with regular Git
- Files > 100 MiB and <= 500 MiB added with Git LFS
- Excluded: .tmp.drivedownload, .tmp.driveupload, .Rproj.user, desktop.ini, Thumbs.db, *.tmp, *.lock

This note was added automatically during the sync.

## Repository sync note

Date: 2026-04-14 18:46:27
Repository synced from the current local tree at: E:\hPlacenta-architecture

Rules used:
- Files <= 100 MiB added with regular Git
- Files > 100 MiB and <= 500 MiB added with Git LFS
- Included if within size rules: .qs, .rds, .RData, .rda
- Excluded: .tmp.drivedownload, .tmp.driveupload, .Rproj.user, desktop.ini, Thumbs.db, *.tmp, *.lock

Remote main backup branch: backup/origin-main-20260414_184548
Local backup branch: backup/local-main-20260414_184548

## Repository sync note

Date: 2026-04-14 18:59:25
Repository synced from the current local tree at: E:\hPlacenta-architecture

Rules used:
- Files <= 100 MiB added with regular Git
- Files > 100 MiB and <= 500 MiB added with Git LFS
- Files > 500 MiB excluded
- Included if within size rules: .qs, .rds, .RData, .rda
- Excluded: .tmp.drivedownload, .tmp.driveupload, .Rproj.user, desktop.ini, Thumbs.db, *.tmp, *.lock

## Repository sync note

Date: 2026-04-14 21:09:16
Repository synced from the current local tree at: E:\hPlacenta-architecture

Rules used:
- Files <= 100 MiB added with regular Git
- Files > 100 MiB and <= 500 MiB added with Git LFS
- Files > 500 MiB excluded
- Included if within size rules: .qs, .rds, .RData, .rda
- Excluded: .tmp.drivedownload, .tmp.driveupload, .Rproj.user, desktop.ini, Thumbs.db, *.tmp, *.lock

## Repository sync note

Date: 2026-04-15 09:30:10
Repository synced from the current local tree at: E:\hPlacenta-architecture

Rules used:
- Files <= 100 MiB added with regular Git
- Files > 100 MiB and <= 500 MiB added with Git LFS
- Files > 500 MiB excluded
- Included if within size rules: .qs, .rds, .RData, .rda
- Excluded: .tmp.drivedownload, .tmp.driveupload, .Rproj.user, desktop.ini, Thumbs.db, *.tmp, *.lock

## Repository sync note

Date: 2026-04-15 09:32:29
Repository synced from the current local tree at: E:\hPlacenta-architecture

Rules used:
- Files <= 100 MiB added with regular Git
- Files > 100 MiB and <= 500 MiB added with Git LFS
- Files > 500 MiB excluded
- Included if within size rules: .qs, .rds, .RData, .rda
- Excluded: .tmp.drivedownload, .tmp.driveupload, .Rproj.user, desktop.ini, Thumbs.db, *.tmp, *.lock

## Repository sync note

Date: 2026-04-15 10:00:18
Repository synced from the current local tree at: E:\hPlacenta-architecture

Rules used:
- Files <= 100 MiB added with regular Git
- Files > 100 MiB and <= 500 MiB added with Git LFS
- Files > 500 MiB excluded
- Included if within size rules: .qs, .rds, .RData, .rda
- Excluded: .tmp.drivedownload, .tmp.driveupload, .Rproj.user, desktop.ini, Thumbs.db, *.tmp, *.lock

## Repository sync note

Date: 2026-04-15 10:48:02
Repository synced from the current local tree at: E:\hPlacenta-architecture

Rules used:
- Files <= 100 MiB added with regular Git
- Files > 100 MiB and <= 500 MiB added with Git LFS
- Files > 500 MiB excluded
- Included if within size rules: .qs, .rds, .RData, .rda
- Excluded: .tmp.drivedownload, .tmp.driveupload, .Rproj.user, desktop.ini, Thumbs.db, *.tmp, *.lock

## Repository sync note

Date: 2026-04-15 11:12:18
Repository synced from the current local tree at: E:\hPlacenta-architecture

Rules used:
- Files <= 100 MiB added with regular Git
- Files > 100 MiB and <= 500 MiB added with Git LFS
- Files > 500 MiB excluded
- Included if within size rules: .qs, .rds, .RData, .rda
- Excluded: .tmp.drivedownload, .tmp.driveupload, .Rproj.user, desktop.ini, Thumbs.db, *.tmp, *.lock
