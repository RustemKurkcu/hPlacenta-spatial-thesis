# Placenta “Vicious Cycle” Project — Analysis Report (v1)

**Project goal:** use **single-cell reference** + **spatial multiomics** (Slide-tags + STARmap) to
quantify *when*, *where*, and *why* the maternal–fetal interface becomes permissive to
pathogen invasion (Fn/Listeria framing), and to connect that permissive window to
remodeling programs, immune tolerance, and metabolic axes.

This document does two things:

1) **Interprets the outputs you already generated** (tables/figures in `output/`)
2) Defines a **two-track pipeline** (minimal fixed vs cutting-edge extended) that you can
re-run reproducibly without breaking Seurat v5.

> **Sources:** see `docs/SOURCE_MANIFEST.md` (BibTeX in `docs/SOURCE_MANIFEST.bib`).

---

## 0. Datasets and what each can/cannot support

### 0.1 Multiome RNA reference (10x multiome RNA)
**Role:** high-coverage *taxonomy + marker discovery*; provides anchors for label transfer.

**Strengths:** genome-wide expression; strong for differential expression and program scoring.

**Limitations:** if only the RNA modality is present, TF motif/gene activity inference is limited
unless ATAC is also available.

### 0.2 Slide-tags (spatial RNA; UMI-like)
**Role:** spatial distribution across subjects/weeks; good for spatial neighborhoods and
cell-type co-localization.

**Strengths:** broad gene coverage; supports module scoring, spatial autocorrelation, and
program-by-celltype analyses.

**Limitations:** sparse per-spot/cell depending on preprocessing; careful normalization required.

### 0.3 STARmap-ISS (panel-based spatial, with optional imputed expression)
**Role:** high-quality spatial cellular geometry and local microenvironment structure.

**Two assays matter (and we treat them differently):**

* `RNA_raw` (measured panel): use for **mapping / anchors / DE on panel genes**.
* `imputed` (genome-wide imputation): use for **program scoring and hypothesis generation**
  (not for formal DE unless the method’s assumptions are explicitly justified).

This split is essential for defensible Methods.

---

## 1. Core biological hypotheses (what we test)

### H1 — A “Remodeling Highway” creates a permissive spatial corridor
**Hypothesis:** remodeling programs (ECM/MMP, EVT invasion/remodeling) form a spatially
structured corridor; invasion risk is highest where remodeling and tolerance programs co-occur.

*Operational tests*
* Program scores (ECM/MMP, EVT invasion) are **spatially clustered** and peak during a
  developmental window.
* Cell-type adjacency (EVT ↔ fibroblast ↔ myeloid) is **enriched** beyond a randomized null.

### H2 — Immune tolerance niches coincide with reduced cytotoxic surveillance
**Hypothesis:** tolerance/immune checkpoint programs align with altered NK / myeloid states,
increasing permissiveness.

*Operational tests*
* NK cytotoxic module scores are **locally depleted** where tolerance programs peak.
* Ligand–receptor axes (e.g., checkpoint pairs highlighted in placenta atlases) are elevated in
  neighborhoods that score high for permissiveness.

### H3 — Metabolic “bait” and “nutritional immunity” program shifts are temporally structured
**Hypothesis:** ethanolamine/PLD1-related programs fluctuate across the same developmental
window, consistent with “bait” (baseline) and “defense” (downregulation upon infection).

*Operational tests*
* PLD1/EA-related program scores show time-varying patterns by trophoblast subtype.

---

## 2. What your current outputs already show (quick read)

This section is based on the current tables in `output/tables/`.

### 2.1 Neighborhood enrichment (Slide-tags)
Files: `SlideTags_week_*_neighbor_z.csv`

**Week 8:** strong positive adjacency within the EVT lineage and between EVT and
EVT-progenitor-type labels (high positive Z). This is consistent with a structured trophoblast
remodeling niche.

**Week 11:** Hofbauer ↔ Endothelial adjacency is among the strongest signals, suggesting a
shift toward immune–vascular microenvironments later in the timeline.

> Interpretation: you are already seeing *time-dependent rewiring* of “who sits next to whom”,
> which is exactly what you need for the “spatiotemporal window” story.

### 2.2 Neighborhood enrichment (STARmap)
Files: `STARmap_week_*_neighbor_z.csv`

The panel-based nature and lower gene dimensionality makes these Z-scores **more sensitive**
to null model details (cell counts, k, permutation count). Your early outputs show very large
negative values for some pairs, which usually means one of:

* very sparse pair counts (near zero observed)
* too few permutations / unstable variance in the null
* a distance threshold (or k) that is too small or too large for the STARmap cell density

**Fix:** in the minimal pipeline we clamp k and ensure robust permutation counts; in the
advanced pipeline we move to a *Greenbaum-style neighborhood vector* approach.

### 2.3 Timecourse gene/program summaries
File: `timecourse_gene_module_summaries.csv`

This table is the backbone of the thesis/grant “timeline” figures:

* It gives **per-week × per-celltype** mean/median for genes and program scores.
* Slide-tags already shows **week-dependent trends in PLD1** at the dataset level.
* STARmap raw panel cannot report PLD1 if it is not in the panel (expected); this motivates
  using the `imputed` assay for program scoring.

---

## 3. Two pipelines you asked for

### Pipeline A — “Minimal + Correct” (stable, thesis-safe)
Goal: eliminate errors, enforce Seurat v5 correctness, and regenerate the core deliverables.

Key properties:
* Strict assay rules (`RNA_raw` for STARmap mapping; `imputed` only for scoring/plots)
* Robust module scoring that cannot fail on small panels (adaptive `nbin`/`ctrl`)
* UMAP/tSNE generated for **reference + each query** (for global vs local structure)

Produces:
* Mapping QC plots (UMAP + tSNE)
* Spatial celltype maps for each subject/week
* Neighborhood enrichment Z matrices by week
* Timecourse tables for marker genes + programs
* Harmonized celltype tables that **never overwrite author labels**

### Pipeline B — “Cutting-edge Extended” (more ambitious)
Goal: add microenvironment discovery + spatially constrained signaling + metagene programs.

Adds:
1) **Neighborhood vectors → niche clustering** (in the style of recent placenta spatial atlases)
2) **Spatially constrained ligand–receptor inference** (CellChat if available; fallback simple LR)
3) **Immune subtyping** (NK and macrophage state scores; conservative assignment)
4) **Metagene program discovery** (NMF/WGCNA optional) + housekeeping QC diagnostics

Produces:
* Niche maps (neighborhood cluster IDs in space)
* “Permissiveness surface” + top/bottom neighborhood DE comparisons
* Prioritized candidate signaling axes per niche and week

---

## 4. Methods (what we do, why we do it, what the null is)

### 4.1 Normalization
* **LogNormalize** is used when counts are not UMI-like (or when the data layer is missing).
* **SCTransform** is used when counts are UMI-like and SCT is desired for variance
  stabilization.

**Why:** anchor finding and program scoring depend on the scale of the expression matrix.
SCTransform is a regularized NB regression method designed to stabilize variance across genes
[(Hafemeister2019SCTransform)].

### 4.2 Dimensionality reduction
* **PCA** as the linear basis for neighbor graph construction.
* **UMAP** for global manifold visualization [(McInnes2018UMAP)].
* **t-SNE** (optional) for local neighborhood visualization [(VanDerMaaten2008tSNE)].

**Hypothesis framing:** embeddings are visualization tools; statistical testing should be done in
expression/program space, not in 2D.

### 4.3 Label transfer / mapping
We use Seurat’s anchor-based transfer from the multiome RNA reference to each spatial query
[(Stuart2019SeuratIntegration)]. If you later add multi-modal reference integration (RNA+ATAC),
Seurat’s WNN framework is the recommended pattern [(Hao2021WNN)].

**Null:** if no shared structure exists between query and reference, transfer confidence scores
should be low and assignments unstable.

### 4.4 Conservative label harmonization (your “don’t erase author labels” rule)
We keep:
* `celltype_author` (never overwritten)
* `predicted.id` + confidence

And create:
* `celltype_final_conservative`: author if present; else predicted if high confidence; else unknown.
* `celltype_lineage_rescue`: for unknowns only, assign broad lineage based on marker scoring.

### 4.5 Module scoring (robust to gene panel limits)
We use Seurat’s AddModuleScore concept, but adapt `nbin` and `ctrl` based on assay size and
unique mean-expression values, so it won’t fail on panel data [(Satija2015Seurat)].

**Null:** genes are compared to control features of similar mean expression.

### 4.6 Spatial neighborhood enrichment
Minimal pipeline: KNN adjacency in (x,y), compute observed pair counts, compare to permuted
labels to get Z-scores (the same statistical framing used in placenta spatial atlases that compare
observed adjacency to randomized nulls) [(Greenbaum2024SpatialTimelineMFI)].

Advanced pipeline: build neighborhood composition vectors (cell-type frequencies around each
cell), cluster neighborhoods to define niches, then track niche dynamics across week, following
the “spatial niche / neighborhood composition” strategy used in recent placenta spatial work
[(Greenbaum2024SpatialTimelineMFI; Ounadjela2024PlacentaMultiomic)].

---

## 5. Deliverables checklist (what you will get after re-run)

### Tables
* `*_counts_by_week_*`
* `timecourse_gene_module_summaries.csv`
* `*_neighbor_z.csv` per week and dataset
* `*_celltype_harmonization_summary.csv`
* (Advanced) `*_niche_composition.csv`, `*_niche_LR_rankings.csv`

### Figures
* multiome QC + UMAP + tSNE
* Slide-tags mapping QC (UMAP+tSNE) + spatial maps (per subject)
* STARmap mapping QC + spatial maps
* neighborhood enrichment heatmaps
* permissiveness score maps (advanced)

---

## 6. Known limitations and how we handle them

1) **STARmap imputation:** use imputed assay for scoring only; interpret as hypothesis, not
   definitive DE.
2) **Batch/subject effects:** do not pool blindly; always stratify outputs by `subject`/`week`.
3) **Spatial resolution mismatch:** Slide-tags vs STARmap have different spot/cell geometry;
   neighborhood parameters (k, radius) must be dataset-specific.

---

## 7. What to run

*Minimal pipeline runner:* `RUN_PIPELINE_MINIMAL.R`

*Advanced pipeline runner:* `RUN_PIPELINE_ADVANCED.R`

---

## 8. Next “advisor-ready” figure (suggested)

**3-panel figure:**

1) Week 8–11 permissiveness surface (Slide-tags)
2) Neighborhood enrichment heatmap per week (EVT/fibroblast/myeloid focus)
3) Candidate checkpoint / tolerizing LR axes by niche (advanced; CellChat or fallback LR)

This directly mirrors the spatial-tolerance niche logic described in recent placenta atlases
(Greenbaum2024SpatialTimelineMFI; VentoTormo2018MFIAtlas) and connects to your
invasion-window hypothesis. If you run the optional communication layers, keep the
interpretation aligned with the underlying statistical model (e.g., curated LR scoring or
CellChat/CellPhoneDB/NicheNet-style inference) (Jin2021CellChat; Efremova2020CellPhoneDB;
Browaeys2020NicheNet).
