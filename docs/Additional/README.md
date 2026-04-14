<h1>Placenta Vicious Cycle Analysis Pipeline</h1><p><strong>Project:</strong> Spatiotemporal Dynamics of Placental Nutritional Immunity<br><strong>Hypothesis:</strong> <em>F. nucleatum</em> exploits placental immune privilege and metabolic gradients to establish a self-reinforcing "vicious cycle" driving preeclampsia<br><strong>Version:</strong> Final — Seurat v5.4.0 Tested<br><strong>Date:</strong> 2026-02-15</p><hr><h2>📚 Documentation</h2><table class="e-rte-table"> <thead> <tr> <th>Document</th> <th>Description</th> </tr> </thead> <tbody><tr> <td><strong><a href="docs/ANALYSIS_REPORT.md">docs/ANALYSIS_REPORT.md</a></strong></td> <td>Comprehensive methods, hypotheses, and rationale for every analytical step</td> </tr> <tr> <td><strong><a href="docs/SOURCE_MANIFEST.md">docs/SOURCE_MANIFEST.md</a></strong></td> <td>Complete academic bibliography with 30 cited references</td> </tr> <tr> <td><a href="docs/00_COUNTS_AND_PLATFORMS.md">docs/00_COUNTS_AND_PLATFORMS.md</a></td> <td>Data types and platform details</td> </tr> <tr> <td><a href="docs/01_ANALYSIS_GUARDRAILS.md">docs/01_ANALYSIS_GUARDRAILS.md</a></td> <td>Best practices and caveats</td> </tr> <tr> <td><a href="docs/02_FIGURE_PLAN_AND_METHODS_SNIPPETS.md">docs/02_FIGURE_PLAN_AND_METHODS_SNIPPETS.md</a></td> <td>Figure generation guide</td> </tr> <tr> <td><a href="docs/03_RUN_ORDER.md">docs/03_RUN_ORDER.md</a></td> <td>Script execution order</td> </tr> <tr> <td><a href="docs/05_METHODS_WRITEUP_TEMPLATE.md">docs/05_METHODS_WRITEUP_TEMPLATE.md</a></td> <td>Methods section template for manuscripts</td> </tr> </tbody></table><hr><h2>🔧 All Known Issues Fixed</h2><table class="e-rte-table"> <thead> <tr> <th>Issue</th> <th>Root Cause</th> <th>Fix</th> </tr> </thead> <tbody><tr> <td><code>FeaturePlot</code> crash: "non-numeric argument"</td> <td><code>umap_1</code>/<code>umap_2</code> metadata columns clash with UMAP reduction</td> <td>Auto-rename to <code>meta_umap_1</code>/<code>meta_umap_2</code></td> </tr> <tr> <td><code>Layer 'data' is empty</code> warning</td> <td>Seurat v5 split layers not joined</td> <td><code>safe_join_layers()</code> before all processing</td> </tr> <tr> <td><code>Insufficient data values to produce 24 bins</code></td> <td>Fixed nbin=24 too large for small datasets</td> <td>Adaptive nbin (3–24) based on cell count</td> </tr> <tr> <td><code>'Layers' is not exported</code></td> <td>Seurat v4 code in v5 environment</td> <td><code>has_data_layer()</code> with try/catch</td> </tr> <tr> <td><code>could not find function "%R%"</code></td> <td>Non-standard operator</td> <td>Replaced with <code>strrep()</code></td> </tr> <tr> <td>Wrong assay for STARmap</td> <td>Active assay RNA_raw has only 1004 genes</td> <td><code>select_starmap_assay()</code> prefers imputed</td> </tr> </tbody></table><hr><h2>📁 Directory Structure</h2><pre><code>├── README.md                              # This file
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
│   ├── ANALYSIS_REPORT.md                 # ✨ Comprehensive methods &amp; hypotheses
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
</code></pre><hr><h2>🚀 Quick Start</h2><h3>1. Setup</h3><pre><code class="language-r"># Set working directory
setwd("path/to/FINAL_PIPELINE")

# Verify Seurat version
packageVersion("Seurat")  # Should be 5.x
</code></pre><h3>2. Update Data Paths</h3><p>Edit <code>config/config.R</code>:</p><pre><code class="language-r">PATH_MULTIOME_RDS  &lt;- "data/processed/multiome_rna_seurat.rds"
PATH_SLIDETAGS_RDS &lt;- "data/processed/slidetags_mapped_to_multiome.rds"
PATH_STARMAP_RDS   &lt;- "data/processed/starmap_spatial_raw_plus_imputed_seurat.rds"
</code></pre><h3>3. Run Pipeline</h3><pre><code class="language-r"># Step 1: Preprocess reference
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
</code></pre><hr><h2>🔬 The "Vicious Cycle" Hypothesis</h2><pre><code>    ┌─────────────────────────────────────────────┐
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
</code></pre><hr><h2>📊 Key Outputs</h2><table class="e-rte-table"> <thead> <tr> <th>Output</th> <th>Description</th> </tr> </thead> <tbody><tr> <td><code>multiome_reference_processed.rds</code></td> <td>Normalized, clustered reference with UMAP/tSNE</td> </tr> <tr> <td><code>slidetags_mapped.rds</code></td> <td>Slide-tags with cell types, UMAP, spatial coords</td> </tr> <tr> <td><code>starmap_mapped.rds</code></td> <td>STARmap with cell types, UMAP, spatial coords</td> </tr> <tr> <td><code>*_harmonized.rds</code></td> <td>Final harmonized cell type labels + lineage scores</td> </tr> <tr> <td><code>*_celltype_harmonization_summary.csv</code></td> <td>Cell type assignment breakdown</td> </tr> <tr> <td><code>output/figures/*.png</code></td> <td>All QC and analysis plots</td> </tr> </tbody></table><hr><h2>📖 Citation</h2><p>If you use this pipeline, please cite:</p><ol> <li>Greenbaum et al. (2023). <em>Nature</em>, 619, 801–810. [Primary data]</li> <li>Hao et al. (2024). <em>Nature Biotechnology</em>, 42, 293–304. [Seurat v5]</li> <li>See <code>docs/SOURCE_MANIFEST.md</code> for the complete bibliography.</li> </ol><hr><h2>🆘 Troubleshooting</h2><table class="e-rte-table"> <thead> <tr> <th>Symptom</th> <th>Solution</th> </tr> </thead> <tbody><tr> <td>"non-numeric argument to mathematical function"</td> <td>Fixed: metadata clash auto-renamed</td> </tr> <tr> <td>"Layer 'data' is empty"</td> <td>Fixed: JoinLayers called automatically</td> </tr> <tr> <td>"Insufficient data values to produce N bins"</td> <td>Fixed: adaptive nbin parameter</td> </tr> <tr> <td>Script fails on STARmap</td> <td>Check that imputed assay exists; see logs</td> </tr> <tr> <td>Gene set scores are all NA</td> <td>Check gene availability with <code>check_gene_availability()</code></td> </tr> </tbody></table><p><strong>Check logs:</strong> <code>output/logs/*.log</code> for detailed execution traces.</p><hr><p><em>Pipeline ready for production use. All scripts tested with Seurat v5.4.0.</em></p>