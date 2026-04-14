# Placenta “Vicious Cycle” Spatial Pipeline (Seurat v5)

This repository contains analysis code to integrate:

* a **multiome-derived RNA reference** (Seurat object), and
* two spatial modalities:
  * **Slide-tags** (spatial RNA)
  * **STARmap-ISS** (measured panel + optional imputation)

The biological goal is to quantify **spatiotemporal vulnerability** at the
maternal–fetal interface (remodeling zones + immune tolerance niches) and relate that
to infection susceptibility.

## What you get

Two runnable pipelines:

1. **Minimal (stable) pipeline** — bug-fixes + robust plotting + robust module scoring.
2. **Advanced (cutting-edge) pipeline** — adds microenvironment discovery, immune
   subtyping, and optional spatially constrained cell–cell communication.

## Key docs

* `docs/ANALYSIS_REPORT.md` — hypotheses, methods, and how to interpret outputs.
* `docs/SOURCE_MANIFEST.bib` — bibliography (BibTeX, easy to extend).
* `docs/SOURCE_MANIFEST.md` — the same bibliography in human-readable form.

When writing a thesis/grant/paper, cite sources from `docs/SOURCE_MANIFEST.bib`.

## Seurat v5 note (layers)

Seurat v5 stores assay matrices as **layers**. Functions like `Layers()` and
`JoinLayers()` live in **SeuratObject**, not the **Seurat** namespace. This pipeline
therefore calls `SeuratObject::Layers()` and `SeuratObject::JoinLayers()` internally.

## Quick start

1. Open `config/config.R` and set the three input paths:
   * `PATH_MULTIOME_RDS`
   * `PATH_SLIDETAGS_RDS` (or `PATH_SLIDETAGS_RAW`)
   * `PATH_STARMAP_RDS`
2. Run the minimal pipeline:

```r
source("RUN_PIPELINE_MINIMAL.R")
```

3. If you want the expanded analysis (recommended once minimal is clean):

```r
source("RUN_PIPELINE_ADVANCED.R")
```

Outputs are written to:

* `output/objects/` (RDS)
* `output/figures/` (PNG)
* `output/tables/` (CSV)
* `output/logs/` (LOG)

## Reproducibility

* Set `SEED` in `config/config.R`.
* Keep a copy of your session info (`sessionInfo()`) with each major run.
* @@ -40,25 +40,205 @@ therefore calls `SeuratObject::Layers()` and `SeuratObject::JoinLayers()` intern
   * `PATH_SLIDETAGS_RDS` (or `PATH_SLIDETAGS_RAW`)
   * `PATH_STARMAP_RDS`
2. Run the minimal pipeline:

```r
source("RUN_PIPELINE_MINIMAL.R")
```

3. If you want the expanded analysis (recommended once minimal is clean):

```r
source("RUN_PIPELINE_ADVANCED.R")
```

Outputs are written to:

* `output/objects/` (RDS)
* `output/figures/` (PNG)
* `output/tables/` (CSV)
* `output/logs/` (LOG)

## Reproducibility

* Set `SEED` in `config/config.R`.
* Keep a copy of your session info (`sessionInfo()`) with each major run.


## Spatial v2 (publication-focused) updates

The 05_spatial series now includes a stronger **enhanced v2** methodology:

- **K-sensitivity robustness panel** in neighborhood and LR analyses (`k = 8, 15, 25, 40`) to avoid over-claims from single-k results.
- **Week-stratified neighborhood rewiring** statistics for niche composition change over gestation.
- **Spatial-density stratified neighborhood enrichment** (quartiles) to reduce artifacts from local crowding differences.
- **Spatially constrained LR null model** (receiver-label shuffling within week and density bins) to better calibrate significance.
- **Pseudo-bulk niche DE by donor/week** to complement marker-style single-cell DE with stronger inferential structure.

### New/updated 05 outputs

- `scripts/05_spatial/05B_neighborhood_enrichment.R`
  - `*_neighbor_enrichment_long.csv`
  - `*_neighbor_enrichment_k_robustness.csv`
  - `*_neighbor_enrichment_trends.csv`
  - `*_neighbor_enrichment_trend_shortlist.csv`
  - `*_neighbor_enrichment_effect_summary.csv`
  - Heatmaps in both palettes:
    - white→red (`*_neighbor_z_heatmap_*`)
    - blue→white→red (`*_neighbor_z_heatmap_*_bwr.png`)
- `scripts/05_spatial/05D_neighborhood_DE.R`
  - `*_niche_week_composition.csv`
  - `*_niche_fraction_trends.csv` (which niches increase/decrease with week)
  - `*_niche_week_rewiring_jsd.csv`
  - `*_niche_center_celltype_composition.csv` and `*_niche_center_celltype_top5.csv`
  - `*_niche_geneset_scores.csv` and `*_niche_week_geneset_scores.csv`
  - `*_niche_pseudobulk_de.csv`
- `scripts/05_spatial/05E_spatial_lr_proximity.R`
  - `*_spatial_lr_edges.csv`
  - `*_spatial_lr_summary.csv`
  - `*_spatial_lr_k_robustness.csv`
- `scripts/05_spatial/05F_interaction_adjacency_followup.R`
  - `*_adjacency_followup_summary.csv`
  - Spatial follow-up figures for shortlisted interaction pairs

### Interpreting 05B neighbor heatmaps and tables

- Axes are both cell types by design:
  - y-axis = **center (focal) cell label**
  - x-axis = **neighbor label**
- Each tile is a z-score comparing observed center→neighbor adjacency vs a label-shuffled null.
- Positive (red) = more frequent than random; negative (blue in BWR plots) = less frequent than random.
- Use these companion tables for quantitative interpretation:
  - `*_neighbor_enrichment_long.csv` (observed, expected, z, fold-enrichment, FDR)
  - `*_neighbor_enrichment_effect_summary.csv` (median effects, fraction significant)
  - `*_neighbor_enrichment_trends.csv` (increasing/decreasing over week)
  - `*_neighbor_enrichment_trend_shortlist.csv` (compact shortlist of changing interactions)

### Time-change and adjacency follow-up workflow

1. Run `05B` to estimate pairwise interactions by week and k.
2. Use `*_neighbor_enrichment_trend_shortlist.csv` to identify interactions increasing/decreasing over time.
3. Run `05F` to compare, for each shortlisted pair:
   - center cells adjacent to target neighbor type vs not adjacent,
   - fraction adjacent and closeness summary per week,
   - spatial plots highlighting adjacency status and closeness.

### Why this matters biologically

These changes better support publication-quality claims that placental permissive niches are:

1. **Spatially structured** (non-random adjacency),
2. **Temporally dynamic** (week-specific rewiring), and
3. **Mechanistically plausible** (LR support beyond constrained null expectations),

instead of being driven by a single neighborhood parameter choice or density confounding.


### Quick interpretation guide (new)

If you are new to the spatial outputs, use:

- `docs/05_SPATIAL_INTERPRETATION_CHEATSHEET.md`

It provides a one-page, non-technical explanation of what each 05 module does, how to read each figure/table, and common interpretation pitfalls.

## Fast post-hoc figure refresh (no heavy reruns)

If 04A/04C already ran and you want clearer figure panels without recomputing heavy steps:

- From existing 04A enhanced summary table:
  - `Rscript scripts/04_timecourse/04A_timecourse_posthoc_significance_plots.R`
  - Produces **all-trends** and **significant-only** plots plus trend tables.
- From existing 04C coordination table:
  - `Rscript scripts/04_timecourse/04C_gene_coordination_posthoc_plots.R`
  - Produces coordination dotplots (all + significant-only), top-hit barplots,
    and updated stats table with BH-adjusted significance.

This is useful when 04C takes a long time and you only want improved visualization and interpretation outputs.

## Recommended additional graphs (high value) + how to interpret

The pipeline now supports three publication-friendly figure families for 04A/04C:

1. **Timecourse panels (all vs significant-only)**
   - Built by:
     - `scripts/04_timecourse/04A_gene_of_interest_timecourse_ENHANCED.R` (during full run), and
     - `scripts/04_timecourse/04A_timecourse_posthoc_significance_plots.R` (fast refresh from CSV).
   - Interpretation:
     - **All** view = complete biological context.
     - **Significant-only** view = concise headline trends (BH FDR<0.05).
     - Line type encodes trend significance; use these for manuscript summary figures.

2. **Coordination dotplot (all + significant-only)**
   - Built by:
     - `scripts/04_timecourse/04C_gene_coordination_posthoc_plots.R`.
   - Interpretation:
     - Dot fill = `coordination_z` (higher = tighter temporal coordination vs null).
     - Dot size = geneset size (`n_genes`) to contextualize confidence/stability.
     - Significant-only panel helps focus on robust pathway/celltype coordination signals.

3. **Top coordination hits barplot**
   - Built by:
     - `scripts/04_timecourse/04C_gene_coordination_posthoc_plots.R`.
   - Interpretation:
     - Ranks strongest coordination findings (by FDR then absolute z).
     - Useful as a compact “main findings” panel and figure legend entry.

### Method notes for interpretation

- **Timecourse significance**: Spearman correlation across weeks per dataset × celltype,
  BH-adjusted across tested celltypes/metrics.
- **Coordination score**: compares observed geneset COM dispersion (Gini) against random null sets;
  higher `coordination_z` indicates stronger temporal coordination than expected by chance.
- **Reporting best practice**:
  - show both all-data and significant-only figures,
  - report `n_genes` and week coverage,
  - avoid over-interpreting datasets with only 2 weeks (trend tests are underpowered).



## Biological Markers and Literature Sources

To improve biological interpretability and manuscript/grant traceability, the pipeline includes high-resolution marker signatures in `config/gene_sets.R` for trophoblast, macrophage, NK, stromal, and vascular subtyping.

### Trophoblast substates (Ounadjela 2024; Cui 2022)
- `Pan_Trophoblast`: `EGFR`
- `vCTB_Core`: `PAGE4, PEG10, MKI67, TOP2A`
- `vCTB_Substates`: `TBL1X, SMAGP, IFI6, LRP5`
- `EVT_Core`: `HLA-G, CCNE1, NOTUM`
- `EVT_Substates`: `UTRN, HAPLN3, LY6E, AOC1, PAPPA2`
- `STB_Core`: `CYP19A1, ERVFRD-1, CGA, GDF15, ENDOU`

### Macrophage subsets (Greenbaum 2023; Hoo 2024; Vento-Tormo 2018)
- `Mac_Pan`: `CD14, CD68, SPP1`
- `Mac_M2_Tolerogenic`: `CD163, MRC1, HAVCR2, LGALS9, CD274, HIF1A, VEGFA`
- `Mac_M1_Inflammatory`: `HLA-DRA, CD40, CD80, IL1B, IL6, CXCL8`
- `Mac_Fetal_HBC`: `LYVE1, ADAMTS17, FPR2, MSR1`
- `Mac_Maternal_PAMM1`: `CD74, LYZ`
- `Mac_Spatial_Decidual`: `CD209, ITGAX`

### NK subsets (Vento-Tormo 2018)
- `dNK_Pan`: `NCAM1, CD9, ITGA1`
- `dNK1`: `LILRB1, ENTPD1`
- `dNK2`: `ANXA1, ITGB2`
- `dNK3`: `CD160, ITGAE`

### Stromal and vascular subsets (Ounadjela 2024)
- `FIB_Pan`: `COL3A1, COL6A2, VIM`
- `FIB_Fetal`: `PDGFRB, AGTR1, PDGFRA, CXCL14`
- `FIB_Maternal`: `ALDH1A2, FAM155A`
- `Endo_Vascular`: `PECAM1, KDR, CDH5, VWF`

### Why these signatures are used
These signatures support biologically specific spatial interpretation beyond broad lineage labels, including:
- distinguishing fetal Hofbauer-like vs maternal PAMM-like macrophage programs,
- mapping dNK1/dNK2/dNK3-like immune states,
- separating tolerogenic vs inflammatory myeloid niches,
- resolving trophoblast substate geography and transitions,
- improving grant/paper methods reproducibility with explicit literature-linked marker definitions.

### Primary literature
- Ounadjela et al., *Nature Medicine* (2024)
- Cui et al., *Bioengineering & Translational Medicine* (2022)
- Greenbaum et al., *Nature* (2023)
- Vento-Tormo et al., *Nature* (2018)
- Hoo et al., *Cell Host & Microbe* (2024)
