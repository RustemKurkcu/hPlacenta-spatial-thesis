# Counts and platforms: STARmap vs Slide-tags vs 10x Multiome

This document explains **what a “count” means** in each platform you’re using, why they differ,
and how to do analyses that stay scientifically defensible.

## 1) Big picture

You are integrating three *different measurement technologies*:

- **10x Multiome (RNA)**: droplet-based sequencing with **UMIs** (molecule counts).
- **Slide-tags spatial RNA**: sequencing with **UMIs**, plus spatial coordinates (cell-resolved).
- **STARmap-ISS**: imaging/in situ sequencing with a **targeted gene panel** (here ~1,001 genes),
  and an additional **imputed transcriptome matrix** that is *not directly measured*.

Because the measurement process differs, the “counts” differ in meaning, noise model, and normalization.

## 2) What STARmap “counts” are

### 2.1 STARmap-ISS raw panel counts (REAL MEASUREMENT)
**What it is:**
- STARmap-ISS is an in situ method: RNA molecules are detected in tissue by targeted chemistry and read out by imaging/sequencing-in-place.
- For your dataset, the “raw” matrix is a **1,001-gene panel** measured per cell.

**How a number becomes a count (conceptual):**
1) Tissue is stained/processed with a STARmap-ISS chemistry
2) Imaging produces fluorescent signals per sequencing cycle
3) Base calls / reads are decoded into “transcripts”
4) Detected transcripts are assigned to segmented cells
5) Each detected transcript increments a gene’s count for that cell

**Important consequence:** These are *not UMIs*. They are “detected molecule events” from imaging.
They can behave like counts, but the error profile differs (optical noise, segmentation errors, per-gene detection efficiency).

### 2.2 STARmap imputed transcriptome (NOT DIRECTLY MEASURED)
You also have an “imputed_expression” matrix (genome-wide genes).
This is created computationally by using similarities between:
- STARmap panel expression
- a sequencing-based reference (here, the multiome/snRNA data)

It is useful for:
- visualization (heatmaps, UMAP coloring)
- hypothesis generation (candidate pathways)
- ligand–receptor exploration (with distance constraints)

It is **not appropriate** for:
- reporting “raw expression levels” as if they were measured
- differential expression claims that depend on count-based statistics

## 3) What Slide-tags “counts” are

Slide-tags is sequencing-based. It yields UMI counts similar to scRNA/snRNA-seq,
but also provides *cell-resolved spatial coordinates* by using slide-attached tags.

**Key consequence:** Slide-tags counts are much more comparable to Multiome RNA counts than STARmap counts are.

## 4) What 10x Multiome RNA “counts” are

10x Multiome RNA uses UMIs → each UMI ideally corresponds to one captured transcript molecule.
Counts are sparse and affected by capture efficiency and ambient RNA, but are “classic” single-cell UMI counts.

## 5) Practical rules for cross-platform analysis

### Rule A — never compare raw counts across platforms
Instead:
- normalize within each dataset
- use **log-normalized expression**, **module scores**, or **z-scored per-dataset values**
- compare **directions of effect** (up/down over time) rather than absolute levels

### Rule B — treat imputed STARmap as a different data type
Label it everywhere as “imputed” and keep it out of any claim that depends on measurement.

### Rule C — do per-cell-type, per-timepoint summaries
For each dataset:
- define a harmonized cell-type label
- summarize mean/median expression per (week × celltype)
- do statistics across donors where possible (pseudobulk)

### Rule D — ratios are OK, but do them carefully
If you want gene ratios (e.g., PLD1 / MMP9 or EA metabolism / MMP module):
- use log space: `log1p(A) - log1p(B)` to avoid division by zeros
- compute per cell type + per sample; compare distributions across weeks
- interpret as relative shifts, not absolute molecule ratios

## 6) Recommended outputs (for grant/thesis figures)

1) **Data modality overview** (schematic + summary table)
2) **Cell type maps**:
   - Multiome UMAP
   - Slide-tags spatial map colored by predicted cell type
   - STARmap spatial map colored by predicted cell type
3) **Temporal window** (weeks 7–11):
   - cell-type abundance vs week
   - MMP “remodeling highway” score vs week (EVT/FIB)
4) **Spatial microenvironment**:
   - neighborhood enrichment (EVT ↔ FIB, EVT ↔ endo, EVT ↔ myeloid)
   - distance-to-STB gradients of key modules
5) **Communication**:
   - CellChat (distance constrained) interactions relevant to EVT invasion / immune tolerance

## 7) References (for your thesis/grant)
- Wang et al. STARmap (in situ sequencing / imaging-based spatial transcriptomics)
- Rodriques et al. Slide-tags (UMI sequencing + spatial tags)
- Ounadjela et al. / Greenbaum et al. placenta spatial multiomic landscape (Slide-tags + STARmap-ISS + imputation)
- Seurat label transfer / integration methodology (Stuart et al. 2019)
