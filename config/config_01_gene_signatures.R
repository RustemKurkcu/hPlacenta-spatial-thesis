# ==============================================================================
# config/01_gene_signatures.R
# Comprehensive Gene Signature Database for Placenta "Vicious Cycle" Project
# ==============================================================================
#
# PURPOSE:
#   Centralized repository of all gene signatures used in the analysis,
#   with literature citations and biological rationale.
#
# ORGANIZATION:
#   1. Bacterial Adhesion Targets (Fap2, FadA)
#   2. Metabolic Pathways (Ethanolamine, PLD1)
#   3. ECM Remodeling (MMP, Collagen)
#   4. Immune Signatures (NK, Macrophage, T cell)
#   5. Toxic Switch Pathways (H2S, Ammonia)
#   6. Barrier Function (STB, Tight junctions)
#   7. Cell Type Markers (Comprehensive)
#   8. Signaling Pathways (for CellChat)
#
# USAGE:
#   source("config/01_gene_signatures.R")
#   genes <- GENE_SIGNATURES$fap2_targets$all
#
# ==============================================================================

# ==============================================================================
# 1. BACTERIAL ADHESION TARGETS
# ==============================================================================

# Fap2 Targets - Gal-GalNAc Glycosylation Enzymes
# Reference: Abed et al. 2016 (Cell Host & Microbe)
#            Parhi et al. 2020 (EMBO Reports)
# 
# BIOLOGICAL RATIONALE:
#   Fap2 is the primary adhesin of F. nucleatum in placental tissue.
#   It recognizes Gal-GalNAc glycans on cell surfaces, particularly
#   on EVT cells which express high levels of GALNT1 (the "address label").
#   This creates tissue-specific tropism.

FAP2_TARGETS <- list(
  # Core GALNT family - GalNAc transferases
  galnt_family = c(
    "GALNT1",   # Primary placental target - "address label"
    "GALNT2",   # Ubiquitous GalNAc-T
    "GALNT3",   # Expressed in trophoblasts
    "GALNT7",   # Mucin-type O-glycosylation
    "GALNT10",  # Tissue-specific expression
    "GALNT14"   # Placental expression
  ),
  
  # Core 1/2 synthesis enzymes
  core_synthesis = c(
    "C1GALT1",    # Core 1 synthase (T-synthase)
    "C1GALT1C1",  # Core 1 synthase chaperone (Cosmc)
    "GCNT1",      # Core 2 GlcNAc-T
    "GCNT3"       # Core 2/4 GlcNAc-T
  ),
  
  # Sialylation enzymes
  sialylation = c(
    "ST6GALNAC1",  # Sialyltransferase
    "ST6GALNAC2"   # Sialyltransferase
  ),
  
  # Beta-3-GlcNAc transferases
  beta3_transferases = c(
    "B3GNT6"  # Core 3 synthase
  ),
  
  # All Fap2 targets combined
  all = c(
    "GALNT1", "GALNT2", "GALNT3", "GALNT7", "GALNT10", "GALNT14",
    "C1GALT1", "C1GALT1C1", "GCNT1", "GCNT3",
    "ST6GALNAC1", "ST6GALNAC2", "B3GNT6"
  )
)

# FadA Targets - E-Cadherin Adhesion Complex
# Reference: Rubinstein et al. 2013 (Cell Host & Microbe)
#
# BIOLOGICAL RATIONALE:
#   FadA is a secondary adhesin that binds E-cadherin (CDH1) and
#   modulates β-catenin signaling. This is shared with Listeria
#   (making Listeria a good infection proxy). FadA binding disrupts
#   cell-cell junctions and activates oncogenic pathways.

FADA_TARGETS <- list(
  # Core E-cadherin complex
  core_complex = c(
    "CDH1",     # E-cadherin - primary FadA receptor
    "CTNNB1",   # Beta-catenin - signaling mediator
    "CTNNA1",   # Alpha-catenin - actin linker
    "CTNND1",   # Delta-catenin (p120) - stability
    "JUP"       # Junction plakoglobin (gamma-catenin)
  ),
  
  # Additional junction components
  junction_associated = c(
    "CDH2",     # N-cadherin (neural)
    "CDH3",     # P-cadherin (placental)
    "CTNNB1",   # Beta-catenin
    "CTNNA2",   # Alpha-catenin 2
    "CTNNA3"    # Alpha-catenin 3
  ),
  
  # All FadA targets
  all = c(
    "CDH1", "CTNNB1", "CTNNA1", "CTNND1", "JUP",
    "CDH2", "CDH3", "CTNNA2", "CTNNA3"
  )
)

# ==============================================================================
# 2. METABOLIC PATHWAYS - ETHANOLAMINE (BAIT SIGNAL)
# ==============================================================================

# Ethanolamine Metabolism - The "BAIT" Signal
# Reference: Garsin 2010 (Nature Rev Microbiol)
#            Thiennimitr et al. 2011 (PNAS)
#
# BIOLOGICAL RATIONALE:
#   Ethanolamine (EA) is released during membrane turnover and cell death.
#   F. nucleatum uses EA as a nutrient source and chemoattractant via the
#   Eut operon. The placenta produces EA via PLD1 to fuel trophoblast
#   invasion, inadvertently creating a "dinner bell" for bacteria.
#   
#   THE VICIOUS CYCLE:
#   1. Placenta produces EA (BAIT) → attracts Fn
#   2. Host detects infection → downregulates PLD1 (DEFENSE)
#   3. Fn starves → activates toxic switch (ACCIDENT)
#   4. H2S/ammonia → vascular collapse (COLLAPSE)

ETHANOLAMINE_PATHWAY <- list(
  # Ethanolamine kinases - first step
  kinases = c(
    "ETNK1",  # Ethanolamine kinase 1
    "ETNK2"   # Ethanolamine kinase 2
  ),
  
  # CDP-ethanolamine synthesis
  cdp_synthesis = c(
    "PCYT2"   # CTP:phosphoethanolamine cytidylyltransferase
  ),
  
  # Phosphatidylethanolamine synthesis
  pe_synthesis = c(
    "SELENOI",  # Ethanolamine phosphotransferase (EPT)
    "EPT1",     # Ethanolamine phosphotransferase 1
    "CEPT1"     # Choline/ethanolamine phosphotransferase
  ),
  
  # PE modification
  pe_modification = c(
    "PEMT"    # PE N-methyltransferase
  ),
  
  # KEY BAIT MARKER - Phospholipase D1
  key_marker = c(
    "PLD1"    # *** PRIMARY BAIT SIGNAL ***
  ),
  
  # All ethanolamine pathway genes
  all = c(
    "ETNK1", "ETNK2", "PCYT2", "SELENOI", "EPT1", 
    "CEPT1", "PEMT", "PLD1"
  )
)

# Phospholipase D Family - EA Production
PLD_FAMILY <- list(
  phospholipases = c(
    "PLD1",   # Primary EA producer - KEY BAIT
    "PLD2",   # Secondary EA producer
    "PLD3",   # Endosomal PLD
    "PLD4",   # Phagocytic PLD
    "PLD5",   # Mitochondrial PLD
    "PLD6"    # MitoPLD
  ),
  
  key_bait = "PLD1"  # *** MOST IMPORTANT ***
)

# ==============================================================================
# 3. ECM REMODELING - "REMODELING HIGHWAY"
# ==============================================================================

# Matrix Metalloproteinases - Physical Barrier Breach
# Reference: Greenbaum et al. 2023 (Nature)
#
# BIOLOGICAL RATIONALE:
#   MMPs digest extracellular matrix to allow trophoblast invasion.
#   This creates "physical holes" in the barrier that Fn exploits.
#   MMP2 and MMP9 are the primary gelatinases active in placenta.
#   The "Remodeling Highway" is the spatial zone of high MMP activity.

MMP_REMODELING <- list(
  # Primary gelatinases - KEY HIGHWAY MARKERS
  gelatinases = c(
    "MMP2",   # Gelatinase A - constitutive
    "MMP9"    # Gelatinase B - inducible
  ),
  
  # Additional MMPs in placenta
  collagenases = c(
    "MMP1",   # Collagenase-1
    "MMP8",   # Collagenase-2
    "MMP13"   # Collagenase-3
  ),
  
  stromelysins = c(
    "MMP3",   # Stromelysin-1
    "MMP10",  # Stromelysin-2
    "MMP11"   # Stromelysin-3
  ),
  
  membrane_type = c(
    "MMP14",  # MT1-MMP - cell surface
    "MMP15",  # MT2-MMP
    "MMP16",  # MT3-MMP
    "MMP17"   # MT4-MMP
  ),
  
  # Tissue inhibitors - BARRIER PROTECTION
  inhibitors = c(
    "TIMP1",  # TIMP-1
    "TIMP2",  # TIMP-2
    "TIMP3",  # TIMP-3
    "TIMP4"   # TIMP-4
  ),
  
  # All MMP pathway genes
  all = c(
    "MMP1", "MMP2", "MMP3", "MMP7", "MMP8", "MMP9", "MMP10", 
    "MMP11", "MMP13", "MMP14", "MMP15", "MMP16", "MMP17",
    "TIMP1", "TIMP2", "TIMP3", "TIMP4"
  )
)

# Extracellular Matrix Components
ECM_COMPONENTS <- list(
  # Collagens - structural integrity
  collagens = c(
    "COL1A1", "COL1A2",  # Type I collagen
    "COL3A1",            # Type III collagen
    "COL4A1", "COL4A2",  # Type IV collagen (basement membrane)
    "COL5A1", "COL5A2",  # Type V collagen
    "COL6A1", "COL6A2", "COL6A3"  # Type VI collagen
  ),
  
  # Laminins - basement membrane
  laminins = c(
    "LAMA1", "LAMA2", "LAMA3", "LAMA4", "LAMA5",
    "LAMB1", "LAMB2", "LAMB3",
    "LAMC1", "LAMC2", "LAMC3"
  ),
  
  # Fibronectin
  fibronectin = c(
    "FN1"   # Fibronectin
  ),
  
  # Proteoglycans
  proteoglycans = c(
    "HSPG2",  # Perlecan
    "DCN",    # Decorin
    "BGN",    # Biglycan
    "LUM",    # Lumican
    "VCAN"    # Versican
  ),
  
  # All ECM components
  all = c(
    "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", 
    "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3",
    "LAMA1", "LAMA2", "LAMA3", "LAMA4", "LAMA5",
    "LAMB1", "LAMB2", "LAMB3", "LAMC1", "LAMC2", "LAMC3",
    "FN1", "HSPG2", "DCN", "BGN", "LUM", "VCAN"
  )
)

# ==============================================================================
# 4. IMMUNE SIGNATURES - COMPREHENSIVE
# ==============================================================================

# Hofbauer Cells - Fetal Macrophages (TROJAN_HORSE)
# Reference: Thomas et al. 2021 (Immunity)
#
# BIOLOGICAL RATIONALE:
#   Hofbauer cells are fetal-derived macrophages that reside in placental
#   villi. They are M2-polarized and have high phagocytic capacity,
#   making them ideal intracellular reservoirs for Fn (TROJAN_HORSE).
#   They can harbor bacteria and disseminate infection.

HOFBAUER_CELLS <- list(
  # Core macrophage markers
  core_macrophage = c(
    "CD68",    # Pan-macrophage
    "CD163",   # M2 macrophage / hemoglobin scavenger
    "CD14",    # LPS receptor
    "CSF1R",   # M-CSF receptor
    "FCGR1A"   # Fc gamma receptor I (CD64)
  ),
  
  # Hofbauer-specific markers
  hofbauer_specific = c(
    "FOLR2",   # Folate receptor beta - HIGHLY SPECIFIC
    "LYVE1",   # Lymphatic vessel endothelial HA receptor
    "F13A1",   # Coagulation factor XIII A
    "MAF"      # Transcription factor
  ),
  
  # M2 polarization markers
  m2_polarization = c(
    "MRC1",    # Mannose receptor (CD206)
    "CD209",   # DC-SIGN
    "MSR1",    # Scavenger receptor A
    "MARCO",   # Macrophage receptor with collagenous structure
    "MERTK",   # MER tyrosine kinase
    "AXL",     # AXL receptor tyrosine kinase
    "TIMD4"    # T-cell immunoglobulin mucin-4
  ),
  
  # Phagocytosis and intracellular survival
  phagocytosis = c(
    "CD36",    # Scavenger receptor
    "SCARB1",  # Scavenger receptor B1
    "FCGR2A",  # Fc gamma receptor IIA
    "FCGR2B",  # Fc gamma receptor IIB
    "FCGR3A"   # Fc gamma receptor IIIA (CD16)
  ),
  
  # All Hofbauer markers
  all = c(
    "CD68", "CD163", "CD14", "CSF1R", "FCGR1A",
    "FOLR2", "LYVE1", "F13A1", "MAF",
    "MRC1", "CD209", "MSR1", "MARCO", "MERTK", "AXL", "TIMD4",
    "CD36", "SCARB1", "FCGR2A", "FCGR2B", "FCGR3A"
  )
)

# Maternal Macrophages - Decidual Macrophages
MATERNAL_MACROPHAGES <- list(
  # Maternal-specific markers
  maternal_specific = c(
    "S100A8",  # Calgranulin A
    "S100A9",  # Calgranulin B
    "FCN1",    # Ficolin-1
    "VCAN",    # Versican
    "CD86"     # Costimulatory molecule
  ),
  
  # Pro-inflammatory markers
  proinflammatory = c(
    "IL1B",    # Interleukin-1 beta
    "TNF",     # Tumor necrosis factor
    "IL6",     # Interleukin-6
    "CXCL8",   # IL-8
    "CCL2"     # MCP-1
  ),
  
  # All maternal macrophage markers
  all = c(
    "S100A8", "S100A9", "FCN1", "VCAN", "CD86",
    "IL1B", "TNF", "IL6", "CXCL8", "CCL2"
  )
)

# NK Cells - Decidual NK Cells
# Reference: Vento-Tormo et al. 2018 (Nature)
#
# BIOLOGICAL RATIONALE:
#   Decidual NK (dNK) cells are specialized NK cells that regulate
#   trophoblast invasion and vascular remodeling. They are less cytotoxic
#   than peripheral NK cells and create an "immune privileged" zone.
#   Fn may exploit this tolerance.

NK_CELLS <- list(
  # Core NK markers
  core_nk = c(
    "NCAM1",   # CD56 - NK cell marker
    "NKG7",    # Natural killer granule protein 7
    "GNLY"     # Granulysin
  ),
  
  # Cytotoxic markers
  cytotoxic = c(
    "GZMB",    # Granzyme B
    "GZMA",    # Granzyme A
    "PRF1",    # Perforin
    "GZMK",    # Granzyme K
    "GZMH"     # Granzyme H
  ),
  
  # NK receptors
  receptors = c(
    "KLRD1",   # CD94
    "KLRF1",   # NKp80
    "KLRB1",   # CD161
    "KLRC1",   # NKG2A
    "KLRC2",   # NKG2C
    "KLRK1",   # NKG2D
    "NCR1",    # NKp46
    "NCR3"     # NKp30
  ),
  
  # Decidual NK-specific
  decidual_specific = c(
    "CD9",     # Tetraspanin
    "ITGAE",   # CD103 - integrin alpha E
    "KIR2DL4"  # Killer cell Ig-like receptor
  ),
  
  # All NK markers
  all = c(
    "NCAM1", "NKG7", "GNLY",
    "GZMB", "GZMA", "PRF1", "GZMK", "GZMH",
    "KLRD1", "KLRF1", "KLRB1", "KLRC1", "KLRC2", "KLRK1", "NCR1", "NCR3",
    "CD9", "ITGAE", "KIR2DL4"
  )
)

# T Cells
T_CELLS <- list(
  # Core T cell markers
  core_t = c(
    "CD3D",    # CD3 delta
    "CD3E",    # CD3 epsilon
    "CD3G"     # CD3 gamma
  ),
  
  # CD4+ T cells
  cd4 = c(
    "CD4"      # CD4 molecule
  ),
  
  # CD8+ T cells
  cd8 = c(
    "CD8A",    # CD8 alpha
    "CD8B"     # CD8 beta
  ),
  
  # Regulatory T cells
  treg = c(
    "FOXP3",   # Forkhead box P3
    "IL2RA",   # CD25
    "CTLA4",   # CTLA-4
    "IKZF2"    # Helios
  ),
  
  # All T cell markers
  all = c(
    "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B",
    "FOXP3", "IL2RA", "CTLA4", "IKZF2"
  )
)

# Dendritic Cells
DENDRITIC_CELLS <- list(
  # Core DC markers
  core_dc = c(
    "ITGAX",   # CD11c
    "CD1C",    # BDCA-1
    "CLEC9A",  # DNGR-1
    "CLEC10A"  # CD301
  ),
  
  # Plasmacytoid DC
  pdc = c(
    "CLEC4C",  # BDCA-2
    "IL3RA",   # CD123
    "TCF4"     # E2-2
  ),
  
  # All DC markers
  all = c(
    "ITGAX", "CD1C", "CLEC9A", "CLEC10A",
    "CLEC4C", "IL3RA", "TCF4"
  )
)

# ==============================================================================
# 5. IMMUNE TOLERANCE / PRIVILEGE
# ==============================================================================

# Immune Tolerance Markers
# Reference: Vento-Tormo et al. 2018 (Nature)
#
# BIOLOGICAL RATIONALE:
#   The placenta creates an "immune privileged" zone to prevent
#   maternal rejection of the semi-allogeneic fetus. This involves
#   expression of non-classical MHC (HLA-G), checkpoint molecules
#   (PD-L1), and immunosuppressive cytokines (TGF-β, IL-10).
#   Fn may exploit this tolerance to evade immune clearance.

IMMUNE_TOLERANCE <- list(
  # Non-classical MHC
  hla_g = c(
    "HLA-G"    # Non-classical MHC class I - KEY TOLERANCE MARKER
  ),
  
  # Checkpoint molecules
  checkpoints = c(
    "CD274",    # PD-L1
    "PDCD1LG2", # PD-L2
    "CTLA4",    # CTLA-4
    "LAG3",     # LAG-3
    "HAVCR2"    # TIM-3
  ),
  
  # Immunosuppressive enzymes
  enzymes = c(
    "IDO1",     # Indoleamine 2,3-dioxygenase 1
    "IDO2",     # Indoleamine 2,3-dioxygenase 2
    "ARG1",     # Arginase 1
    "ARG2"      # Arginase 2
  ),
  
  # Immunosuppressive cytokines
  cytokines = c(
    "TGFB1",    # TGF-beta 1
    "TGFB2",    # TGF-beta 2
    "TGFB3",    # TGF-beta 3
    "IL10",     # Interleukin-10
    "IL35"      # Interleukin-35 (EBI3 + IL12A)
  ),
  
  # All tolerance markers
  all = c(
    "HLA-G", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2",
    "IDO1", "IDO2", "ARG1", "ARG2",
    "TGFB1", "TGFB2", "TGFB3", "IL10"
  )
)

# ==============================================================================
# 6. TOXIC SWITCH PATHWAYS
# ==============================================================================

# Hydrogen Sulfide (H2S) Metabolism
# Reference: Attene-Ramos et al. 2006
#            Szabo et al. 2014 (Nature Rev Drug Discov)
#
# BIOLOGICAL RATIONALE:
#   When starved of ethanolamine, Fn activates MegL (methionine gamma-lyase)
#   which produces H2S from methionine. H2S is a potent vascular toxin that
#   causes endothelial dysfunction and hypertension - hallmarks of preeclampsia.
#   This is the "TOXIC SWITCH" - the metabolic accident.

H2S_METABOLISM <- list(
  # H2S production enzymes
  production = c(
    "CBS",     # Cystathionine beta-synthase
    "CTH",     # Cystathionine gamma-lyase
    "MPST"     # Mercaptopyruvate sulfurtransferase
  ),
  
  # H2S detoxification
  detoxification = c(
    "SQOR",    # Sulfide quinone oxidoreductase
    "ETHE1",   # Persulfide dioxygenase
    "TST"      # Thiosulfate sulfurtransferase
  ),
  
  # Vascular targets of H2S
  vascular_targets = c(
    "NOS3",    # eNOS - endothelial nitric oxide synthase
    "EDN1",    # Endothelin-1
    "VEGFA",   # VEGF-A
    "FLT1",    # VEGFR1 (sFlt-1 in preeclampsia)
    "KDR"      # VEGFR2
  ),
  
  # All H2S pathway genes
  all = c(
    "CBS", "CTH", "MPST", "SQOR", "ETHE1", "TST",
    "NOS3", "EDN1", "VEGFA", "FLT1", "KDR"
  )
)

# Ammonia Metabolism
AMMONIA_METABOLISM <- list(
  # Ammonia production
  production = c(
    "GLUL",    # Glutamine synthetase
    "GLUD1",   # Glutamate dehydrogenase 1
    "GLUD2"    # Glutamate dehydrogenase 2
  ),
  
  # Urea cycle - ammonia detoxification
  urea_cycle = c(
    "CPS1",    # Carbamoyl phosphate synthetase 1
    "OTC",     # Ornithine transcarbamylase
    "ASS1",    # Argininosuccinate synthase
    "ASL",     # Argininosuccinate lyase
    "ARG1"     # Arginase 1
  ),
  
  # All ammonia pathway genes
  all = c(
    "GLUL", "GLUD1", "GLUD2",
    "CPS1", "OTC", "ASS1", "ASL", "ARG1"
  )
)

# ==============================================================================
# 7. BARRIER FUNCTION - PROTECTIVE
# ==============================================================================

# Syncytiotrophoblast (STB) - Primary Barrier
# Reference: Turco & Moffett 2019 (Development)
#
# BIOLOGICAL RATIONALE:
#   The STB is a multinucleated syncytium that forms the primary
#   barrier between maternal and fetal blood. It is highly resistant
#   to infection due to lack of intercellular junctions and
#   continuous membrane. STB differentiation is protective.

SYNCYTIOTROPHOBLAST <- list(
  # Syncytialization transcription factors
  transcription_factors = c(
    "GCM1",    # Glial cells missing 1 - MASTER REGULATOR
    "TFAP2A",  # AP-2 alpha
    "TFAP2C",  # AP-2 gamma
    "GATA2",   # GATA binding protein 2
    "GATA3"    # GATA binding protein 3
  ),
  
  # Syncytins - fusogenic proteins
  syncytins = c(
    "ERVW-1",   # Syncytin-1 (HERV-W env)
    "ERVFRD-1"  # Syncytin-2 (HERV-FRD env)
  ),
  
  # Hormones - STB function markers
  hormones = c(
    "CGB",     # hCG beta (CGB3, CGB5, CGB7, CGB8)
    "CGB3",    # hCG beta 3
    "CGB5",    # hCG beta 5
    "CGB7",    # hCG beta 7
    "CGB8",    # hCG beta 8
    "CSH1",    # Placental lactogen (hPL)
    "CSH2",    # Placental lactogen 2
    "CSHL1"    # Chorionic somatomammotropin hormone-like 1
  ),
  
  # Pregnancy-specific glycoproteins
  glycoproteins = c(
    "PSG1", "PSG2", "PSG3", "PSG4", "PSG5",
    "PSG6", "PSG7", "PSG8", "PSG9"
  ),
  
  # All STB markers
  all = c(
    "GCM1", "TFAP2A", "TFAP2C", "GATA2", "GATA3",
    "ERVW-1", "ERVFRD-1",
    "CGB", "CGB3", "CGB5", "CGB7", "CGB8", "CSH1", "CSH2", "CSHL1",
    "PSG1", "PSG2", "PSG3", "PSG4", "PSG5", "PSG6", "PSG7", "PSG8", "PSG9"
  )
)

# Tight Junctions - Barrier Integrity
TIGHT_JUNCTIONS <- list(
  # Claudins
  claudins = c(
    "CLDN1", "CLDN3", "CLDN4", "CLDN5", "CLDN7"
  ),
  
  # Occludin
  occludin = c(
    "OCLN"
  ),
  
  # Zonula occludens
  zonula_occludens = c(
    "TJP1",    # ZO-1
    "TJP2",    # ZO-2
    "TJP3"     # ZO-3
  ),
  
  # All tight junction genes
  all = c(
    "CLDN1", "CLDN3", "CLDN4", "CLDN5", "CLDN7",
    "OCLN", "TJP1", "TJP2", "TJP3"
  )
)

# ==============================================================================
# 8. CELL TYPE MARKERS - COMPREHENSIVE
# ==============================================================================

# Extravillous Trophoblast (EVT) - PRIMARY VICTIM
EVT_MARKERS <- list(
  # Core EVT markers
  core = c(
    "HLA-G",   # Non-classical MHC - HIGHLY SPECIFIC
    "ITGA5",   # Integrin alpha 5
    "ITGA1",   # Integrin alpha 1
    "MMP2",    # Matrix metalloproteinase 2
    "MMP9"     # Matrix metalloproteinase 9
  ),
  
  # Invasion markers
  invasion = c(
    "PAPPA",   # Pregnancy-associated plasma protein A
    "PAPPA2",  # PAPPA2
    "ADAM12",  # ADAM metallopeptidase 12
    "PLAC8"    # Placenta-specific 8
  ),
  
  # All EVT markers
  all = c(
    "HLA-G", "ITGA5", "ITGA1", "MMP2", "MMP9",
    "PAPPA", "PAPPA2", "ADAM12", "PLAC8"
  )
)

# Villous Cytotrophoblast (vCTB)
VCTB_MARKERS <- list(
  # Core vCTB markers
  core = c(
    "TEAD4",   # TEA domain transcription factor 4
    "TP63",    # Tumor protein p63
    "ITGA6",   # Integrin alpha 6
    "ITGB4",   # Integrin beta 4
    "KRT7"     # Keratin 7
  ),
  
  # Proliferation markers
  proliferation = c(
    "MKI67",   # Ki-67
    "PCNA",    # Proliferating cell nuclear antigen
    "TOP2A"    # Topoisomerase 2A
  ),
  
  # All vCTB markers
  all = c(
    "TEAD4", "TP63", "ITGA6", "ITGB4", "KRT7",
    "MKI67", "PCNA", "TOP2A"
  )
)

# Fibroblasts (FIB2 - BAIT cells)
FIBROBLAST_MARKERS <- list(
  # Core fibroblast markers
  core = c(
    "COL1A1",  # Collagen type I alpha 1
    "COL1A2",  # Collagen type I alpha 2
    "COL3A1",  # Collagen type III alpha 1
    "DCN",     # Decorin
    "LUM"      # Lumican
  ),
  
  # FIB2-specific (from Ounadjela)
  fib2_specific = c(
    "PDGFRA",  # PDGF receptor alpha
    "PDGFRB",  # PDGF receptor beta
    "THY1",    # CD90
    "VIM"      # Vimentin
  ),
  
  # All fibroblast markers
  all = c(
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM",
    "PDGFRA", "PDGFRB", "THY1", "VIM"
  )
)

# Endothelial Cells
ENDOTHELIAL_MARKERS <- list(
  # Core endothelial markers
  core = c(
    "PECAM1",  # CD31
    "VWF",     # von Willebrand factor
    "CDH5",    # VE-cadherin
    "KDR",     # VEGFR2
    "FLT1"     # VEGFR1
  ),
  
  # All endothelial markers
  all = c(
    "PECAM1", "VWF", "CDH5", "KDR", "FLT1"
  )
)

# ==============================================================================
# 9. SIGNALING PATHWAYS (for CellChat)
# ==============================================================================

# Key signaling pathways relevant to placental infection
SIGNALING_PATHWAYS <- list(
  # Growth factors
  growth_factors = c("VEGF", "HGF", "EGF", "FGF", "PDGF", "IGF"),
  
  # Morphogens
  morphogens = c("BMP", "WNT", "NOTCH", "SHH"),
  
  # Immune/inflammatory
  immune = c("TNF", "IL1", "IL6", "IL10", "TGFB", "IFNG"),
  
  # Chemokines
  chemokines = c("CCL", "CXCL", "CX3CL1"),
  
  # Cell adhesion
  adhesion = c("CDH", "NCAM", "ICAM", "VCAM"),
  
  # ECM-receptor
  ecm_receptor = c("ITGA", "ITGB", "CD44", "SDC")
)

# ==============================================================================
# 10. COMBINED SIGNATURES FOR MISI SCORING
# ==============================================================================

# Multi-dimensional Infection Susceptibility Index (MISI) Components
MISI_COMPONENTS <- list(
  fap2_targets = FAP2_TARGETS$all,
  fada_targets = FADA_TARGETS$all,
  ethanolamine_axis = ETHANOLAMINE_PATHWAY$all,
  intracellular_survival = HOFBAUER_CELLS$phagocytosis,
  h2s_sensitivity = H2S_METABOLISM$all,
  barrier_function = SYNCYTIOTROPHOBLAST$all  # Negative weight
)

# MISI Weights (sum to 1.0, barrier is negative)
MISI_WEIGHTS <- c(
  fap2_targets = 0.30,           # Primary adhesion
  fada_targets = 0.15,           # Secondary adhesion
  ethanolamine_axis = 0.20,      # BAIT signal
  intracellular_survival = 0.20, # TROJAN_HORSE
  h2s_sensitivity = 0.05,        # Collateral damage
  barrier_function = -0.10       # Protective (negative)
)

# ==============================================================================
# 11. EXPORT ALL SIGNATURES AS A SINGLE LIST
# ==============================================================================

GENE_SIGNATURES <- list(
  # Bacterial targets
  fap2_targets = FAP2_TARGETS,
  fada_targets = FADA_TARGETS,
  
  # Metabolic pathways
  ethanolamine = ETHANOLAMINE_PATHWAY,
  pld_family = PLD_FAMILY,
  
  # ECM remodeling
  mmp_remodeling = MMP_REMODELING,
  ecm_components = ECM_COMPONENTS,
  
  # Immune cells
  hofbauer_cells = HOFBAUER_CELLS,
  maternal_macrophages = MATERNAL_MACROPHAGES,
  nk_cells = NK_CELLS,
  t_cells = T_CELLS,
  dendritic_cells = DENDRITIC_CELLS,
  
  # Immune tolerance
  immune_tolerance = IMMUNE_TOLERANCE,
  
  # Toxic switch
  h2s_metabolism = H2S_METABOLISM,
  ammonia_metabolism = AMMONIA_METABOLISM,
  
  # Barrier function
  syncytiotrophoblast = SYNCYTIOTROPHOBLAST,
  tight_junctions = TIGHT_JUNCTIONS,
  
  # Cell type markers
  evt_markers = EVT_MARKERS,
  vctb_markers = VCTB_MARKERS,
  fibroblast_markers = FIBROBLAST_MARKERS,
  endothelial_markers = ENDOTHELIAL_MARKERS,
  
  # Signaling pathways
  signaling_pathways = SIGNALING_PATHWAYS,
  
  # MISI components
  misi_components = MISI_COMPONENTS,
  misi_weights = MISI_WEIGHTS
)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Get all genes from a signature
#' @param signature_name Name of the signature (e.g., "fap2_targets")
#' @return Character vector of gene symbols
get_signature_genes <- function(signature_name) {
  sig <- GENE_SIGNATURES[[signature_name]]
  if (is.null(sig)) {
    stop(sprintf("Signature '%s' not found", signature_name))
  }
  if ("all" %in% names(sig)) {
    return(sig$all)
  } else {
    return(unlist(sig, use.names = FALSE))
  }
}

#' Get citation for a signature
#' @param signature_name Name of the signature
#' @return Character string with citation
get_signature_citation <- function(signature_name) {
  citations <- list(
    fap2_targets = "Abed et al. 2016 (Cell Host & Microbe); Parhi et al. 2020 (EMBO Reports)",
    fada_targets = "Rubinstein et al. 2013 (Cell Host & Microbe)",
    ethanolamine = "Garsin 2010 (Nature Rev Microbiol); Thiennimitr et al. 2011 (PNAS)",
    mmp_remodeling = "Greenbaum et al. 2023 (Nature)",
    hofbauer_cells = "Thomas et al. 2021 (Immunity)",
    nk_cells = "Vento-Tormo et al. 2018 (Nature)",
    immune_tolerance = "Vento-Tormo et al. 2018 (Nature)",
    h2s_metabolism = "Attene-Ramos et al. 2006; Szabo et al. 2014 (Nature Rev Drug Discov)",
    syncytiotrophoblast = "Turco & Moffett 2019 (Development)"
  )
  
  citation <- citations[[signature_name]]
  if (is.null(citation)) {
    return("See 00_CITATIONS_DATABASE.md")
  }
  return(citation)
}

#' Print summary of all signatures
print_signature_summary <- function() {
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("GENE SIGNATURE DATABASE SUMMARY\n")
  cat(strrep("=", 70), "\n\n")
  
  for (sig_name in names(GENE_SIGNATURES)) {
    sig <- GENE_SIGNATURES[[sig_name]]
    n_genes <- length(get_signature_genes(sig_name))
    citation <- get_signature_citation(sig_name)
    
    cat(sprintf("%-30s: %3d genes\n", sig_name, n_genes))
    cat(sprintf("  Citation: %s\n\n", citation))
  }
  
  cat(strrep("=", 70), "\n")
  cat(sprintf("Total signatures: %d\n", length(GENE_SIGNATURES)))
  cat(strrep("=", 70), "\n\n")
}

# Print summary when sourced
if (interactive()) {
  print_signature_summary()
}

cat("✓ Gene signatures loaded successfully\n")
cat("  Use GENE_SIGNATURES$<name> to access signatures\n")
cat("  Use get_signature_genes('<name>') to get gene list\n")
cat("  Use get_signature_citation('<name>') for citations\n\n")

# ==============================================================================
# END OF GENE SIGNATURES
# ==============================================================================