# Greenbaum / maternal–fetal interface logic: how we operationalize it here

This note bridges **Greenbaum-style “spatiotemporal immune remodeling”** concepts to the concrete analyses in this pipeline.

## What the Greenbaum paper contributes conceptually

### 1) Immune composition changes with gestational age (GA)
A central message of the spatiotemporal atlas is that **GA and tissue context reshape the immune milieu** at the maternal–fetal interface.

### 2) Macrophage and NK diversity is structured (not “one bucket”)
The atlas highlights **substructure inside macrophage and NK compartments**, e.g.:

- **Macrophage subsets** (e.g., “Mac1” and “Mac2”) can be separated by markers such as **DC-SIGN**, **CD11c**, and **HLA-DR**.
- **NK subsets** (e.g., “NK1–NK4”) can be separated by markers such as **CD57**, **CD11c**, and **CD8**.

### 3) A tolerogenic program emerges and is time-regulated
A key immune “tone” described is **tolerogenic**, supported by expression of markers like **CD206**, **CD163** and inhibitory checkpoints like **TIM-3 / Galectin-9**, plus metabolic immunoregulators like **IDO-1**.

**Why this matters for your project:** these are the ingredients of the “immune privileged zone” exploited during the **Remodeling Highway** window.

---

## How we map these ideas into our pipeline

### A) Immune subtyping (lightweight, score-based)
Script: `scripts/04_timecourse/04B_immune_subsets_refinement.R`

We keep your original labels untouched and add:
- `immune_class` (NK vs Macrophage)
- `immune_subtype` (e.g., cytotoxic-high NK; inflammatory-high vs ISG-high macrophage)

This is designed to work even when STARmap has limited genes, and it produces tables for week-by-week immune shifts.

### B) Spatial proximity enrichment (who neighbors whom)
Script: `scripts/05_spatial/05B_neighborhood_enrichment.R`

We compute KNN adjacency in physical (x,y) space and test label-pair enrichment versus a permutation null.

This operationalizes:
- “Remodeling Highway” neighborhoods (EVT/fibroblast/vascular adjacency)
- “Immune privileged” neighborhoods (tolerogenic macrophage / low cytotoxic NK)

### C) “Permissiveness score” maps (one number per cell)
Script: `scripts/05_spatial/05C_permissiveness_score_maps.R`

Composite score:
- + ECM remodeling
- + immune tolerance checkpoints
- – NK cytotoxicity
- + ethanolamine metabolism

This generates:
- Spatial maps of permissiveness per week
- Composition of the top 10% permissive cells

### D) Immune signaling as cell–cell interaction hypotheses
Script: `scripts/06_cell_communication/06B_simple_LR_scoring.R`

Lightweight ligand–receptor scoring produces interpretable “grant-ready” interactions:
- HLA-G → LILRB1/2 (EVT tolerance)
- LGALS9 → HAVCR2 (TIM-3)
- TGFB1 → TGFBR1/2

---

## Where to plug in deeper immunology later

If you want to go beyond score-based immune subtyping:
- Subcluster immune cells and annotate with canonical markers
- Run CellChat (optional; `RUN_OPTIONAL_HEAVY <- TRUE`)
- Stratify proximity enrichment by immune subtype, not only broad cell type
