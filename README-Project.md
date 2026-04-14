# Placenta Vicious Cycle Analysis Pipeline

**Project:** Spatiotemporal Dynamics of Placental Nutritional Immunity  
**Hypothesis:** *F. nucleatum* exploits placental immune privilege and metabolic gradients to establish a self-reinforcing "vicious cycle" driving preeclampsia  
**Version:** Final — Seurat v5.4.0 Tested  
**Date:** 2026-02-15

---

## 📚 Documentation

| Document | Description |
|----------|-------------|
| **[docs/ANALYSIS_REPORT.md](docs/ANALYSIS_REPORT.md)** | Comprehensive methods, hypotheses, and rationale for every analytical step |
| **[docs/SOURCE_MANIFEST.md](docs/SOURCE_MANIFEST.md)** | Complete academic bibliography with 30 cited references |
| [docs/00_COUNTS_AND_PLATFORMS.md](docs/00_COUNTS_AND_PLATFORMS.md) | Data types and platform details |
| [docs/01_ANALYSIS_GUARDRAILS.md](docs/01_ANALYSIS_GUARDRAILS.md) | Best practices and caveats |
| [docs/02_FIGURE_PLAN_AND_METHODS_SNIPPETS.md](docs/02_FIGURE_PLAN_AND_METHODS_SNIPPETS.md) | Figure generation guide |
| [docs/03_RUN_ORDER.md](docs/03_RUN_ORDER.md) | Script execution order |
| [docs/05_METHODS_WRITEUP_TEMPLATE.md](docs/05_METHODS_WRITEUP_TEMPLATE.md) | Methods section template for manuscripts |

---

## 🔧 All Known Issues Fixed

| Issue | Root Cause | Fix |
|-------|-----------|-----|
| `FeaturePlot` crash: "non-numeric argument" | `umap_1`/`umap_2` metadata columns clash with UMAP reduction | Auto-rename to `meta_umap_1`/`meta_umap_2` |
| `Layer 'data' is empty` warning | Seurat v5 split layers not joined | `safe_join_layers()` before all processing |
| `Insufficient data values to produce 24 bins` | Fixed nbin=24 too large for small datasets | Adaptive nbin (3–24) based on cell count |
| `'Layers' is not exported` | Seurat v4 code in v5 environment | `has_data_layer()` with try/catch |
| `could not find function "%R%"` | Non-standard operator | Replaced with `strrep()` |
| Wrong assay for STARmap | Active assay RNA_raw has only 1004 genes | `select_starmap_assay()` prefers imputed |

---

## 📁 Directory Structure

```
├── README.md                              # This file
├── config/
│   ├── config.R                           # Paths, parameters, column names
│   ├── config_01_gene_signatures.R        # 957-line gene signature database
│   ├── gene_sets.R                        # Hypothesis-driven gene sets
│   ├── ligand_receptor_pairs.csv          # Curated L-R pairs
│   └── celltype_refinement_map.csv        # Label refinement rules
├── scripts/
│   ├── R/
│   │   └── utils.R                        # ✨ Unified utilities (all fixes)
│   ├── 02_preprocess/
│   │   └── 02A_preprocess_multiome_reference.R
│   ├── 03_mapping/
│   │   ├── 03A_map_slidetags_to_multiome.R    # ✨ Fixed metadata clash
│   │   ├── 03B_map_starmap_to_multiome.R      # ✨ Fixed split layers
│   │   └── 03C_harmonize_celltype_labels.R    # ✨ Fixed adaptive scoring
│   ├── 04_timecourse/
│   │   ├── 04A_gene_of_interest_timecourse.R
│   │   ├── 04A_gene_of_interest_timecourse_ENHANCED.R
│   │   ├── 04B_immune_subsets_refinement.R
│   │   └── 04C_gene_coordination_score.R
│   ├── 05_spatial/
│   │   ├── 05A_spatial_overview_plots.R
│   │   ├── 05B_neighborhood_enrichment.R
│   │   └── 05C_permissiveness_score_maps.R
│   ├── 06_cell_communication/
│   │   ├── 06A_cellchat_spatial_constrained.R
│   │   └── 06B_simple_LR_scoring.R
│   ├── 07_export/
│   │   └── 07A_export_shareable_outputs.R
│   └── 08_metagenes/
│       ├── 08A_housekeeping_diagnostics.R
│       ├── 08B_metagene_module_discovery.R
│       └── 08C_metagene_spatiotemporal_maps.R
├── docs/
│   ├── ANALYSIS_REPORT.md                 # ✨ Comprehensive methods & hypotheses
│   ├── SOURCE_MANIFEST.md                 # ✨ Academic bibliography (30 refs)
│   └── ... (additional documentation)
├── data/
│   ├── raw/                               # Original data files
│   └── processed/                         # Seurat objects (.rds)
└── output/
    ├── objects/                            # Pipeline-generated objects
    ├── figures/                            # All plots
    ├── tables/                             # Summary tables
    └── logs/                              # Execution logs
```

---

## 🚀 Quick Start

### 1. Setup

```r
# Set working directory
setwd("path/to/FINAL_PIPELINE")

# Verify Seurat version
packageVersion("Seurat")  # Should be 5.x
```

### 2. Update Data Paths

Edit `config/config.R`:
```r
PATH_MULTIOME_RDS  <- "data/processed/multiome_rna_seurat.rds"
PATH_SLIDETAGS_RDS <- "data/processed/slidetags_mapped_to_multiome.rds"
PATH_STARMAP_RDS   <- "data/processed/starmap_spatial_raw_plus_imputed_seurat.rds"
```

### 3. Run Pipeline

```r
# Step 1: Preprocess reference
source("scripts/02_preprocess/02A_preprocess_multiome_reference.R")

# Step 2: Map spatial datasets
source("scripts/03_mapping/03A_map_slidetags_to_multiome.R")
source("scripts/03_mapping/03B_map_starmap_to_multiome.R")
source("scripts/03_mapping/03C_harmonize_celltype_labels.R")

# Step 3: Downstream analyses
source("scripts/04_timecourse/04A_gene_of_interest_timecourse.R")
source("scripts/05_spatial/05A_spatial_overview_plots.R")
source("scripts/05_spatial/05B_neighborhood_enrichment.R")
source("scripts/06_cell_communication/06B_simple_LR_scoring.R")

# Step 4: Metagene analysis
source("scripts/08_metagenes/08A_housekeeping_diagnostics.R")
source("scripts/08_metagenes/08B_metagene_module_discovery.R")

# Step 5: Export
source("scripts/07_export/07A_export_shareable_outputs.R")
```

---

## 🔬 The "Vicious Cycle" Hypothesis

```
    ┌─────────────────────────────────────────────┐
    │         F. nucleatum Colonization            │
    │  (Fap2 → GALNT1 on EVT = "Address Label")   │
    └──────────────────┬──────────────────────────┘
                       │
                       ▼
    ┌─────────────────────────────────────────────┐
    │         Nutrient Exploitation                │
    │  (Ethanolamine = "Dinner Bell")              │
    │  PLD1 → EA release → bacterial growth        │
    └──────────────────┬──────────────────────────┘
                       │
                       ▼
    ┌─────────────────────────────────────────────┐
    │         Barrier Breakdown                    │
    │  (MMP2/9/14 = "Remodeling Highway")          │
    │  ECM degradation → bacterial spread          │
    └──────────────────┬──────────────────────────┘
                       │
                       ▼
    ┌─────────────────────────────────────────────┐
    │         Immune Dysregulation                 │
    │  (HLA-G/IDO1 collapse → inflammation)        │
    │  NK activation, macrophage polarization      │
    └──────────────────┬──────────────────────────┘
                       │
                       ▼
    ┌─────────────────────────────────────────────┐
    │         Toxic Amplification                  │
    │  (H2S, ammonia → tissue damage)              │
    │  Endothelial dysfunction → preeclampsia      │
    └──────────────────┬──────────────────────────┘
                       │
                       └──────── CYCLE REPEATS ──→ ↑
```

---

## 📊 Key Outputs

| Output | Description |
|--------|-------------|
| `multiome_reference_processed.rds` | Normalized, clustered reference with UMAP/tSNE |
| `slidetags_mapped.rds` | Slide-tags with cell types, UMAP, spatial coords |
| `starmap_mapped.rds` | STARmap with cell types, UMAP, spatial coords |
| `*_harmonized.rds` | Final harmonized cell type labels + lineage scores |
| `*_celltype_harmonization_summary.csv` | Cell type assignment breakdown |
| `output/figures/*.png` | All QC and analysis plots |

---

## 📖 Citation

If you use this pipeline, please cite:

1. Greenbaum et al. (2023). *Nature*, 619, 801–810. [Primary data]
2. Hao et al. (2024). *Nature Biotechnology*, 42, 293–304. [Seurat v5]
3. See `docs/SOURCE_MANIFEST.md` for the complete bibliography.

---

## 🆘 Troubleshooting

| Symptom | Solution |
|---------|----------|
| "non-numeric argument to mathematical function" | Fixed: metadata clash auto-renamed |
| "Layer 'data' is empty" | Fixed: JoinLayers called automatically |
| "Insufficient data values to produce N bins" | Fixed: adaptive nbin parameter |
| Script fails on STARmap | Check that imputed assay exists; see logs |
| Gene set scores are all NA | Check gene availability with `check_gene_availability()` |

**Check logs:** `output/logs/*.log` for detailed execution traces.

---

*Pipeline ready for production use. All scripts tested with Seurat v5.4.0.*