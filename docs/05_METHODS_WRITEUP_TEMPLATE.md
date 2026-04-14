# Methods writeup template (thesis / grant)

Use this as a starting point for a formal Methods section. Replace bracketed text.

---

## Data sources

### Multiome reference
- [Describe sample source, gestational ages, library preparation, QC thresholds]
- Object: `multiome_rna_seurat.rds`

### Slide-tags spatial RNA-seq
- [Describe Slide-tags chemistry, tissue handling, spatial coordinate generation]
- Object: `slidetags_mapped_to_multiome.rds` (or `slidetags_harmonized.rds`)

### STARmap-ISS
- [Describe STARmap panel, imaging pipeline, segmentation, count extraction]
- Object: `starmap_spatial_raw_plus_imputed_seurat.rds`

---

## Preprocessing and normalization

### Multiome reference
We infer whether the RNA `counts` layer is UMI-like (integer dominated).  
- If UMI-like: we apply SCTransform.
- Otherwise: LogNormalize + FindVariableFeatures + ScaleData.

Dimensional reduction (PCA/UMAP) and clustering are run on the normalized assay (SCT or RNA).

Script: `scripts/02_preprocess/02A_preprocess_multiome_reference.R`

### Spatial datasets
- Slide-tags: SCT if present; otherwise, LogNormalize on RNA.
- STARmap: mapping uses `RNA_raw` (measured panel), with JoinLayers() followed by LogNormalize.

---

## Reference mapping / label transfer

We map spatial cells to the multiome taxonomy using Seurat anchor-based integration:
- FindTransferAnchors(reference, query)
- TransferData(refdata = celltype_ref)

Script: `scripts/03_mapping/03A_map_slidetags_to_multiome.R`  
Script: `scripts/03_mapping/03B_map_starmap_to_multiome.R`

**Optional**: MapQuery projection into the reference UMAP is disabled by default due to strict assay/normalization matching requirements.

---

## Cell-type label harmonization

We define:
- `celltype_author` (original labels, if provided)
- `predicted.id` + `prediction.score.max` (Seurat label transfer)

We then derive:
- `celltype_final_conservative`
- `celltype_final_refined` (optional refinement of broad author labels if a mapping table is provided)
- `celltype_lineage_rescue` (lineage-level rescue for unknown cells based on marker scoring)

Script: `scripts/03_mapping/03C_harmonize_celltype_labels.R`

---

## Module scoring

We compute hypothesis-driven module scores (Seurat AddModuleScore), intersecting each gene set with available features.

Key modules include:
- ECM/MMP remodeling (“Remodeling Highway”)
- Immune tolerance checkpoints (“Immune privilege”)
- NK cytotoxicity (“Cytolytic pressure”)
- Ethanolamine metabolism (“Nutrient gradient”)

Gene sets are stored in: `config/gene_sets.R`  
Timecourse script: `scripts/04_timecourse/04A_gene_of_interest_timecourse.R`

---

## Immune subtyping

We perform a lightweight, score-based refinement of immune populations:
- NK cytotoxic-high vs other NK
- inflammatory-high macrophage vs ISG-high macrophage vs other

Script: `scripts/04_timecourse/04B_immune_subsets_refinement.R`

---

## Spatial neighborhood enrichment (proximity analysis)

We quantify whether cell types co-localize beyond chance:
1. Construct a k-nearest-neighbor graph in (x,y) space.
2. Count directed edges between cell-type labels.
3. Permute labels to build a null distribution.
4. Report z-score enrichment for each label pair.

Script: `scripts/05_spatial/05B_neighborhood_enrichment.R`

---

## Composite permissiveness score

We compute a per-cell “permissiveness” score:
z(MMP/ECM) + z(tolerance) – z(NK cytotoxic) + z(ethanolamine)

We visualize permissiveness spatially and test which cell types dominate the high-permissiveness tail.

Script: `scripts/05_spatial/05C_permissiveness_score_maps.R`

---

## Cell–cell interaction hypotheses (lightweight LR scoring)

For curated ligand–receptor pairs, we compute:
score(sender, receiver) = mean(ligand in sender) × mean(receptor in receiver)

Script: `scripts/06_cell_communication/06B_simple_LR_scoring.R`

Optional full modeling:
- CellChat (requires additional dependencies)
- `scripts/06_cell_communication/06A_cellchat_spatial_constrained.R`

---

## Statistics and visualization
- [Specify tests used for differential abundance / score comparisons]
- Figures are produced using ggplot2 and Seurat plotting functions.
- All random processes use a fixed seed (`SEED` in config/config.R).
