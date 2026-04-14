# Source Manifest (Bibliography)

This project keeps citations in **two synchronized formats**:

1. `docs/SOURCE_MANIFEST.bib` — BibTeX (best for manuscripts, grants, Overleaf)
2. `docs/SOURCE_MANIFEST.md` — human-readable list (this file)

**How to add a new source**

* Add a BibTeX entry to `docs/SOURCE_MANIFEST.bib` (DOI preferred).
* Add a short bullet here with the same citation key.
* If a source justifies a specific pipeline choice, mention the step/script.

---

## Core single-cell / Seurat methods

* **Stuart2019SeuratIntegration** — Stuart et al. (2019) *Cell*. Seurat anchors + label transfer.
* **Hao2021SeuratMultimodal** — Hao et al. (2021) *Cell*. Multimodal integration; reference framing.
* **Hafemeister2019SCTransform** — Hafemeister & Satija (2019) *Genome Biology*. SCTransform.
* **McInnes2018UMAP** — McInnes et al. (2018) arXiv. UMAP.
* **vanDerMaaten2008tSNE** — van der Maaten & Hinton (2008) JMLR. t-SNE.

## Cell programs / scoring

* **Satija2015Seurat** — Satija et al. (2015) *Nature Biotechnology*. Seurat + early “module score” concept.

## Cell–cell communication

* **Jin2021CellChat** — Jin et al. (2021) *Nature Communications*. CellChat framework.
* **Chen2025CellChatProtocol** — Chen et al. (2025) *Nature Protocols*. Practical CellChat workflow.
* **Browaeys2020NicheNet** — Browaeys et al. (2020) *Nature Methods*. NicheNet ligand→target.
* **Efremova2020CellPhoneDB** — Efremova et al. (2020) *Nature Protocols*. CellPhoneDB.

## Placenta / maternal–fetal interface atlases (your biological ground truth)

* **Greenbaum2024SpatialTimelineMFI** — Greenbaum et al. (2024) *Nature*. Spatially resolved timeline; immune+trophoblast niches.
* **Ounadjela2024SpatialMultiomePlacenta** — Ounadjela et al. (2024) *Nature Medicine*. Spatial multiomic placenta at molecular resolution.
* **VentoTormo2018MFIAtlas** — Vento-Tormo et al. (2018) *Nature*. Single-cell atlas of early maternal–fetal interface.
