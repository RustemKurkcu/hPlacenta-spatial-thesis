# STARmap Pipeline Alignment Report (Author Method vs Thesis Pipeline)

**Date:** 2026-04-15  
**Scope:** Align `scripts/01_active_pipeline/01_preprocess_harmony_embeddings.R` with methods described by Jian Shu Lab resources and the Nature Medicine placenta paper.

## 1) Source-of-truth references consulted

1. Jian Shu Lab repository: https://github.com/jian-shu-lab/hPlacenta-architecture  
2. Nature Medicine paper: https://www.nature.com/articles/s41591-024-03073-9  
3. Paper public mirror summary (data/code availability excerpt): https://pmc.ncbi.nlm.nih.gov/articles/PMC12660148/

## 2) What the original study reports doing

From the paper + repository structure:

- Technologies integrated: snRNA-seq + snATAC-seq + Slide-tags + STARmap-ISS + STARmap-ISH.
- STARmap-ISS used targeted in situ profiling and downstream whole-transcriptome imputation.
- Early gestation STARmap samples include W7, W8-2, W9, W11.
- Harmony is explicitly used for batch-effect correction in PCA space.
- Code repository modules map to this structure:
  - `STARmap_ISS_data_processing`
  - `STARmap_ISS_Imputation`
  - `Slide-tags`
  - `multiome-analysis`

## 3) Why “No MT genes detected” can be expected here

Your current input files are STARmap **imputed** matrices with a fixed gene panel/derived imputation layer, not raw 10x whole-transcriptome UMI matrices. In this context, absence of `^MT-` features is plausible and not automatically a QC failure.

## 4) Current pipeline status and fixes applied

### A. Parsing and ingestion robustness
- Handles split roots:
  - `data/raw/Broad_SCP2601human-placenta-architecture`
  - `data/raw/zenodo_spatial`
- Delimiter and header-offset detection added.
- Fallback parser logic added.
- Troubleshooting bundle emits compact artifacts per week.

### B. SCTransform memory/future crash fix
- Added `future` controls to avoid 500 MiB globals crash:
  - configurable `future.globals.maxSize` (GiB)
  - default `future::plan("sequential")`
- Added resume checkpoints:
  - per-week post-QC checkpoints
  - merged post-QC checkpoint
  - post-SCTransform checkpoint

This enables restarting from where runs fail, instead of reprocessing all weeks.

## 5) Recommended best-practice run order on your machine

1. Run Script 01 once with checkpoints ON.  
2. If interrupted, rerun Script 01 (it resumes from checkpoints).  
3. After successful SCT, continue PCA/Harmony/UMAP/tSNE in same script.  
4. Keep `future::plan("sequential")` unless you explicitly tune parallel workers and memory.

## 6) Data provenance to keep for thesis defense

For every run:
- Keep `output/reports/01_methods_and_provenance.md`.
- Keep troubleshooting folder snapshots for failures.
- Record final cell counts per week after QC and after merge.
- Record exact script version and commit hash.

## 7) Open methodological cautions

- SCTransform on very large, dense imputed matrices is computationally heavy and may not always be the most biologically conservative choice; sensitivity checks with alternate normalization for specific downstream tasks are recommended.
- “Global + local embeddings” (UMAP + tSNE) is reasonable for exploratory structure, but conclusions should be cross-validated with spatial neighborhood evidence.

