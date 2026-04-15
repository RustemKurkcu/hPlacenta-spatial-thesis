# Master cell-type naming and lineage-grouped colors for CellChat and ggplot use.

rename_with_true_names <- function(x) {
  mapped <- sapply(as.character(x), function(v) {
    switch(v,
           # Immune
           "Hofbauer.cells" = "Hofbauer Macrophages",
           "HBC" = "HBC (Minor)",
           "HBC_p" = "Prolif. Hofbauer Macrophages",
           "PAMM1" = "PAMM1 Macrophages",

           # Fibroblasts
           "FIB1" = "Fibroblast 1",
           "FIB2" = "Fibroblast 2",
           "F" = "Stromal Fibroblast",
           "F_p" = "Prolif. Fibroblast",
           "F_sm" = "Smooth Muscle Fibroblast",

           # Vascular
           "Endothelial" = "Maternal Endothelial",
           "Endo_f" = "Fetal Endothelial",
           "PV" = "Perivascular Cells",

           # Trophoblasts (villous/syncytial)
           "vCTB" = "Villous Cytotrophoblast",
           "VCT_fusing" = "Fusing Cytotrophoblast",
           "VCT_p" = "Prolif. Cytotrophoblast",
           "VCT_CCC" = "Column Cytotrophoblast",
           "VCT" = "VCT",
           "STB" = "Syncytiotrophoblast",
           "STB.progenitor" = "STB Progenitor",

           # Trophoblasts (extravillous)
           "EVT_1" = "Extravillous Trophoblast 1",
           "EVT_2" = "Extravillous Trophoblast 2",
           "iEVT" = "Interstitial EVT",
           "EVT.progenitor" = "EVT Progenitor",
           "EVT" = "EVT",

           # Other
           "Erythroblasts" = "Erythroblasts",
           v)
  }, USE.NAMES = FALSE)
  as.character(mapped)
}

# -------------------------------------------------------------------------
# Universal lineage-grouped palette
# -------------------------------------------------------------------------
universal_colors <- c(
  # Immune
  "Hofbauer Macrophages" = "#E31A1C",
  "HBC (Minor)" = "#FF7F00",
  "Prolif. Hofbauer Macrophages" = "#FB9A99",
  "PAMM1 Macrophages" = "#FDBF6F",

  # Fibroblasts
  "Fibroblast 1" = "#33A02C",
  "Fibroblast 2" = "#B2DF8A",
  "Stromal Fibroblast" = "#006D2C",
  "Prolif. Fibroblast" = "#A6DBA0",
  "Smooth Muscle Fibroblast" = "#B3DE69",

  # Vascular
  "Maternal Endothelial" = "#1F78B4",
  "Fetal Endothelial" = "#A6CEE3",
  "Perivascular Cells" = "#08519C",

  # Villous trophoblasts
  "Villous Cytotrophoblast" = "#6A3D9A",
  "Fusing Cytotrophoblast" = "#CAB2D6",
  "Prolif. Cytotrophoblast" = "#9E9AC8",
  "Column Cytotrophoblast" = "#807DBA",
  "VCT" = "#DADAEB",

  # Syncytial
  "Syncytiotrophoblast" = "#FFD700",
  "STB Progenitor" = "#E6D5B8",

  # Extravillous
  "Extravillous Trophoblast 1" = "#D01C8B",
  "Extravillous Trophoblast 2" = "#F1B6DA",
  "Interstitial EVT" = "#F781BF",
  "EVT Progenitor" = "#DF65B0",
  "EVT" = "#C51B7D",

  # Other
  "Erythroblasts" = "#969696",
  "Myeloid (Unknown)" = "#BDBDBD",
  "Unknown Cell Type" = "#D9D9D9"
)

master_celltype_order <- names(universal_colors)

get_universal_colors <- function(cell_names) {
  cols <- universal_colors[as.character(cell_names)]
  missing <- is.na(cols)
  if (any(missing)) cols[missing] <- "#BDBDBD"
  names(cols) <- as.character(cell_names)
  cols
}

pick_celltype_source_column <- function(md) {
  candidates <- c(
    "cell_label_display",
    "celltype_corrected",
    "celltype_original",
    "predicted.celltype",
    "celltype_final_refined",
    "celltype_final_conservative",
    "celltype_author",
    "predicted.id",
    "seurat_clusters"
  )
  hit <- candidates[candidates %in% colnames(md)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

add_true_celltype_metadata <- function(seu_obj, source_col = NULL) {
  md <- seu_obj@meta.data
  src <- if (!is.null(source_col)) source_col else pick_celltype_source_column(md)
  if (is.na(src) || !(src %in% colnames(md))) {
    md$celltype_original <- "Unknown Cell Type"
    md$celltype_true_name <- "Unknown Cell Type"
    md$cell_label_display <- "Unknown Cell Type"
  } else {
    orig <- as.character(md[[src]])
    true <- rename_with_true_names(orig)
    true[is.na(true) | !nzchar(true)] <- "Unknown Cell Type"
    md$celltype_original <- orig
    md$celltype_true_name <- true
    md$cell_label_display <- true
  }
  seu_obj@meta.data <- md
  seu_obj
}

augment_seurat_with_predicted_celltypes <- function(
    seu_obj,
    roots = c("data/raw/zenodo_spatial", "data/raw/Broad_SCP2601human-placenta-architecture")) {
  md <- seu_obj@meta.data
  files <- unique(unlist(lapply(roots, function(r) {
    if (!dir.exists(r)) return(character(0))
    list.files(r, pattern = "_cell_metadata\\.csv$", recursive = TRUE, full.names = TRUE)
  })))
  if (length(files) == 0) return(seu_obj)

  read_one <- function(path) {
    x <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE), error = function(...) NULL)
    if (is.null(x) || !("cell_id" %in% colnames(x))) return(NULL)
    keep <- intersect(c("cell_id", "predicted.celltype", "predicted.celltype.score"), colnames(x))
    if (!("predicted.celltype" %in% keep)) return(NULL)
    x <- x[, keep, drop = FALSE]
    x$cell_id <- as.character(x$cell_id)
    x
  }
  tbls <- Filter(Negate(is.null), lapply(files, read_one))
  if (length(tbls) == 0) return(seu_obj)
  pred <- unique(do.call(rbind, tbls))
  if (nrow(pred) == 0) return(seu_obj)

  row_id <- if ("cell_id" %in% colnames(md)) as.character(md$cell_id) else rownames(md)
  idx <- match(row_id, pred$cell_id)

  if (!("predicted.celltype" %in% colnames(md))) md$predicted.celltype <- NA_character_
  fill_ct <- is.na(md$predicted.celltype) | !nzchar(as.character(md$predicted.celltype))
  md$predicted.celltype[fill_ct] <- as.character(pred$predicted.celltype[idx[fill_ct]])

  if ("predicted.celltype.score" %in% colnames(pred)) {
    if (!("predicted.celltype.score" %in% colnames(md))) md$predicted.celltype.score <- NA_real_
    fill_sc <- is.na(md$predicted.celltype.score)
    md$predicted.celltype.score[fill_sc] <- suppressWarnings(as.numeric(pred$predicted.celltype.score[idx[fill_sc]]))
  }

  seu_obj@meta.data <- md
  seu_obj
}
