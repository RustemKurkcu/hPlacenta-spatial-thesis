# ======================================================================
# config/gene_sets.R
# Curated (and hypothesis-driven) host gene sets used for module scoring.
#
# Notes
# - These sets are meant to operationalize the "Vicious Cycle" blueprint.
# - Not all genes will be present in all assays/platforms (e.g., STARmap panel).
# - The scoring code will automatically intersect with available features.
#
# Citations / rationale:
# - ECM remodeling / MMP activity is central to EVT invasion and barrier remodeling.
# - Immune tolerance checkpoint signals help define immune-privileged zones.
# - Ethanolamine (EA) metabolism genes operationalize the "nutrient gradient" hypothesis.
# ======================================================================

GENESETS <- list(
  # -----------------------
  # Barrier remodeling ("Remodeling Highway")
  # -----------------------
  MMP_ECM_Remodeling = c(
    "MMP2","MMP9","MMP14","ADAMTS4","ADAMTS5","CTSK","CTSL","CTSB",
    "TIMP1","TIMP2","SERPINE1","PLAUR"
  ),
  ECM_Structure = c(
    "COL1A1","COL1A2","COL3A1","COL5A1","COL6A1","COL6A2","FN1","DCN","LUM"
  ),

  # -----------------------
  # Immune tolerance ("Immune Privilege")
  # -----------------------
  Immune_Tolerance = c(
    "HLA-G","CD274","PDCD1LG2","IDO1","TGFB1","IL10","LGALS9","HAVCR2",
    "VSIR","LILRB1","LILRB2"
  ),
  Cytotoxic_NK = c(
    "NKG7","GNLY","PRF1","GZMB","GZMA","GZMH","CTSW","XCL1","XCL2",
    "KLRD1","KLRK1","NCR1","FCGR3A","TYROBP","TRAC"
  ),
  NK_Expanded = c(
    "NKG7","GNLY","PRF1","GZMB","GZMA","GZMH","CTSW","XCL1","XCL2",
    "KLRD1","KLRK1","KLRB1","KLRC1","KLRC2","NCR1","NCR3","FCGR3A","TYROBP",
    "TRAC","TRBC1","TRBC2","CD247","IL7R",
    "B3GAT1","ITGA1","ITGAX","NCAM1","NOS2","HAVCR2","LGALS9","IDO1","TIGIT",
    "KIR2DL1","KIR3DL1"
  ),
  dNK_Residency = c(
    "KLRD1","KLRB1","NCR1","ITGAE","CXCR4","XCL1","XCL2","GZMK","CD69"
  ),
  NK_Core = c(
    "NKG7","PRF1","GZMB","GNLY","KLRK1"  ),
  Myeloid_Inflammation = c(
    "IL1B","TNF","NFKBIA","CXCL8","S100A8","S100A9","CCL2","CCL3","CCL4"
  ),
  Interferon_Stimulated = c(
    "ISG15","IFIT1","IFIT2","IFIT3","MX1","OAS1","OAS2","OAS3","IFI6","IFI27"
  ),

  # -----------------------
  # Nutrient gradient / ethanolamine axis ("Dinner bell" / "Nutritional trap")
  # -----------------------
  Ethanolamine_Metabolism = c(
    "PLD1","ETNK1","ETNK2","PCYT2","SELENOI","CEPT1","PLA2G6"
  ),

  # -----------------------
  # Vascular / endothelial stress (proxy for PE-like vascular toxicity)
  # -----------------------
  Endothelial_Activation = c(
    "VCAM1","ICAM1","SELE","ENG","EDN1","NOS3","KDR","FLT1","PGF","VEGFA"
  ),
  Hypoxia_Response = c(
    "HIF1A","VEGFA","SLC2A1","LDHA","CA9","BNIP3","EGLN1"
  ),

  # -----------------------
  # Tropism ("Lock vs Label")
  # -----------------------
  Entry_Receptors = c(
    "CDH1","EPCAM","GALNT1"
  )
)

# Convenience: a smaller "core" set for quick demo runs
GENESETS_CORE <- list(
  MMP_ECM_Remodeling = GENESETS$MMP_ECM_Remodeling,
  Immune_Tolerance   = GENESETS$Immune_Tolerance,
  Cytotoxic_NK       = GENESETS$Cytotoxic_NK,
  Ethanolamine_Metabolism = GENESETS$Ethanolamine_Metabolism
)


# -----------------------
# Major lineage marker sets (for conservative label rescue / sanity checks)
# -----------------------
LINEAGE_MARKERS <- list(
  EVT = c("HLA-G","ITGA1","MMP2","MMP9","TACSTD2","KRT7"),
  vCTB = c("EPCAM","KRT7","KRT8","KRT18","TACSTD2"),
  STB = c("CSH1","CSH2","CYP19A1","PSG1","PSG3"),
  Endothelial = c("PECAM1","VWF","KDR","ESAM","RAMP2"),
  Fibroblast = c("COL1A1","COL1A2","DCN","LUM","COL3A1"),
  Macrophage = c("LST1","C1QA","C1QB","C1QC","APOE","TYROBP"),
  NK = c("NKG7","GNLY","PRF1","GZMB","FCER1G")
)
