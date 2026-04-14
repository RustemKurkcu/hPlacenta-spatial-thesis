# Analysis guardrails (what is valid vs invalid)

## 1) Label every matrix by measurement type

- **UMI counts**: Multiome RNA, Slide-tags
- **Imaging counts**: STARmap raw panel
- **Imputed**: STARmap imputed transcriptome

Never describe imputed values as “measured”.

## 2) Differential expression (DE) rules

✅ OK:
- DE within Multiome RNA
- DE within Slide-tags RNA
- DE within STARmap raw *panel* (panel genes only), using non-parametric or appropriate models

⚠️ Risky:
- DE on STARmap imputed (treat as hypothesis generation only)

❌ Not OK:
- Comparing “PLD1 counts” between STARmap and Multiome as absolute values

## 3) Cross-dataset validation pattern

A defensible validation looks like:

1) Detect a trend in Multiome (cell-type specific temporal change)
2) Show the same **direction** of trend in Slide-tags
3) If gene is in STARmap panel, show spatial localization of the same program
4) If gene is *not* in panel, show STARmap-imputed *only as supportive visualization*.

## 4) Spatial analyses

- Treat coordinates as the ground truth for local microenvironment
- Use multiple definitions of neighborhood:
  - KNN (fixed k)
  - radius graph (fixed physical distance, if units are microns)

For key results, confirm that findings are robust to:
- changing k (e.g., 10, 20, 30)
- changing radius (e.g., 25µm, 50µm, 100µm)

## 5) Reporting checklist (for each figure)

- dataset name(s)
- week(s)
- assay used (RNA counts vs STARmap raw vs imputed)
- normalization
- number of cells and donors
- statistical test and multiple testing method
- if spatial: neighborhood definition (k or radius)
