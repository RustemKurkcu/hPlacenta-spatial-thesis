# scripts/05_spatial/05I_DE_hotspot_vs_background.R
# Differential expression: hotspot (top X% permissiveness) vs background cells
# Run from repo root: source("scripts/05_spatial/05I_DE_hotspot_vs_background.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
})

repo_root <- normalizePath(".")
seurat_dir <- file.path(repo_root, "output", "objects")
tables_dir <- file.path(repo_root, "output", "tables", "05_spatial")
out_dir <- file.path(repo_root, "output", "tables", "05_spatial", "DE_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters (tweak as desired)
min.pct <- 0.10          # feature present fraction to consider
logfc.threshold <- 0.25  # log fold-change threshold
adj_p_cutoff <- 0.05     # for reporting / top genes
n_top_report <- 200      # how many top genes to save & plot
test.use <- "wilcox"     # "wilcox" or "MAST" (if MAST installed and preferable)
latent.vars <- c("nCount_RNA")  # covariates in FindMarkers (if present in metadata)

# Helper listing of samples (scans per-sample folders)
sample_dirs <- list.dirs(tables_dir, recursive = FALSE, full.names = TRUE)
if (length(sample_dirs) == 0) stop("No per-dataset sample folders found under ", tables_dir)

log_msg <- function(...) cat(sprintf(...), "\n")

for (ds_dir in sample_dirs) {
  ds_name <- basename(ds_dir)
  sample_dirs2 <- list.dirs(ds_dir, recursive = FALSE, full.names = TRUE)
  if (length(sample_dirs2) == 0) {
    log_msg("No sample folders under dataset dir: %s (skip)", ds_dir)
    next
  }
  for (sdir in sample_dirs2) {
    sname <- basename(sdir)
    log_msg("Processing dataset=%s sample=%s", ds_name, sname)
    
    # hotspot file (as produced by 05H)
    hotspot_file <- list.files(sdir, pattern = paste0(ds_name, ".*hotspot_cells_top.*pct\\.csv$"), full.names = TRUE)
    if (length(hotspot_file) == 0) {
      log_msg("  No hotspot CSV found in %s - skipping", sdir); next
    }
    hotspot_file <- hotspot_file[1]
    
    # read hotspot barcodes
    hotspots <- tryCatch(read_csv(hotspot_file, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(hotspots) || nrow(hotspots) == 0) {
      log_msg("  Hotspot file empty or unreadable: %s - skip", hotspot_file); next
    }
    
    # Determine Seurat RDS for this dataset (prefer harmonized/_with_permissiveness)
    # try common names
    possible_rds <- c(
      file.path(seurat_dir, paste0(ds_name, "_with_permissiveness.rds")),
      file.path(seurat_dir, paste0(ds_name, "_harmonized.rds")),
      file.path(seurat_dir, paste0(ds_name, ".rds"))
    )
    rds_exists <- possible_rds[file.exists(possible_rds)]
    if (length(rds_exists) == 0) {
      # fallback: look through all rds and attempt to find ds_name in filename
      all_rds <- list.files(seurat_dir, pattern = "\\.rds$", full.names = TRUE)
      rds_match <- all_rds[grepl(ds_name, tolower(basename(all_rds)))]
      if (length(rds_match) >= 1) rds_exists <- rds_match[1]
    }
    if (length(rds_exists) == 0) {
      log_msg("  No Seurat RDS found for dataset %s - skipping", ds_name); next
    }
    seurat_path <- rds_exists[1]
    log_msg("  Loading Seurat: %s", seurat_path)
    so <- tryCatch(readRDS(seurat_path), error = function(e) NULL)
    if (is.null(so) || !inherits(so, "Seurat")) {
      log_msg("  Failed to read Seurat object: %s", seurat_path); next
    }
    
    # Ensure cell names matching: hotspot CSV must include a column with barcodes
    # Many of our 05H hotspot CSVs will have a column "cell" from rownames. Try to detect column containing barcodes.
    hotspot_cols <- colnames(hotspots)
    barcode_col <- hotspot_cols[which(tolower(hotspot_cols) %in% c("cell","barcode","cell_id","cellname","cell_id"))]
    if (length(barcode_col) == 0) {
      # try first col
      barcode_col <- hotspot_cols[1]
      log_msg("  Warning: assuming first column %s contains cell barcodes", barcode_col)
    } else {
      barcode_col <- barcode_col[1]
    }
    hot_cells <- as.character(hotspots[[barcode_col]])
    # Some barcodes may have been written with rownames, ensure proper format
    hot_cells <- hot_cells[hot_cells != "" & !is.na(hot_cells)]
    
    # Check which hotspots are present in Seurat
    in_so <- hot_cells %in% colnames(so)
    if (sum(in_so) == 0) {
      log_msg("  None of hotspot barcodes found in Seurat object (0/%d). Trying to match by prefix...", length(hot_cells))
      # try partial matching by prefix if identifiers lost
      matched <- hot_cells[sapply(hot_cells, function(x) any(startsWith(colnames(so), x)))]
      if (length(matched) == 0) { log_msg("  No matches found -> skip"); next }
      hot_cells <- matched
      log_msg("  Found %d matched hotspots via prefix.", length(hot_cells))
    } else if (sum(in_so) < length(hot_cells)) {
      log_msg("  Warning: only %d/%d hotspot barcodes found in Seurat object; continuing with matched subset.", sum(in_so), length(hot_cells))
      hot_cells <- hot_cells[in_so]
    }
    
    # Build labels in Seurat: hotspot vs background
    # Make a temporary metadata column "hotspot_label" with "hotspot" and "background"
    so$hotspot_label <- "background"
    common_cells <- intersect(colnames(so), hot_cells)
    so$hotspot_label[common_cells] <- "hotspot"
    
    # Subset to cells on the image / with permissiveness_global (optional)
    # If 'permissiveness' or 'permissiveness_global' are present in metadata, you can filter NA
    meta <- so@meta.data
    keep_cells <- rownames(meta) # default: all
    if ("permissiveness_global" %in% colnames(meta)) {
      keep_cells <- rownames(meta)[is.finite(meta$permissiveness_global)]
    } else if ("permissiveness" %in% colnames(meta)) {
      keep_cells <- rownames(meta)[is.finite(meta$permissiveness)]
    }
    if (length(keep_cells) == 0) { log_msg("  No finite permissiveness values present -> skip"); next }
    so_sub <- subset(so, cells = keep_cells)
    
    # Pre-filter features: keep genes expressed in >= min.pct fraction in either group
    # Run DE
    de_out <- tryCatch({
      # set identity to hotspot_label
      Idents(so_sub) <- so_sub$hotspot_label
      # ensure two groups exist
      if (!all(c("hotspot", "background") %in% levels(Idents(so_sub)))) {
        log_msg("  Not both hotspot/background present in subset levels: %s", paste(levels(Idents(so_sub)), collapse=", "))
      }
      # run FindMarkers controlling for latent vars if present
      latent_present <- latent.vars[latent.vars %in% colnames(so_sub@meta.data)]
      if (length(latent_present) == 0) latent_present <- NULL
      FindMarkers(
        object = so_sub,
        ident.1 = "hotspot",
        ident.2 = "background",
        min.pct = min.pct,
        logfc.threshold = logfc.threshold,
        test.use = test.use,
        latent.vars = latent_present,
        only.pos = FALSE
      )
    }, error = function(e) { log_msg("  FindMarkers failed: %s", e$message); NULL })
    
    if (is.null(de_out) || nrow(de_out) == 0) {
      log_msg("  No DE results for %s/%s", ds_name, sname); next
    }
    
    # tidy: add gene name column, p_val_adj column may already exist
    de_out <- de_out %>% tibble::rownames_to_column("gene")
    if (!("p_val_adj" %in% colnames(de_out)) && ("p_val" %in% colnames(de_out))) {
      de_out$p_val_adj <- p.adjust(de_out$p_val, method = "BH")
    }
    # rank and save top genes
    de_out <- de_out %>% arrange(p_val_adj, desc(avg_log2FC))
    out_prefix <- file.path(out_dir, ds_name, sname)
    dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
    out_csv <- paste0(out_prefix, "_DE_hotspot_vs_background.csv")
    write_csv(de_out, out_csv)
    log_msg("  Wrote DE CSV: %s (%d genes)", out_csv, nrow(de_out))
    
    # Volcano plot of top genes
    v <- de_out %>% mutate(logp = -log10(p_val_adj + 1e-300)) %>% head(n_top_report)
    plt <- ggplot(de_out, aes(x = avg_log2FC, y = -log10(p_val_adj + 1e-300))) +
      geom_point(alpha = 0.4, size = 0.6) +
      theme_classic() +
      xlab("avg_log2FC (hotspot vs background)") + ylab("-log10(adj p-value)") +
      geom_vline(xintercept = c(-logfc.threshold, logfc.threshold), linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = -log10(adj_p_cutoff), linetype = "dashed", color = "grey50")
    # label top genes
    top_label <- de_out %>% filter(!is.na(p_val_adj)) %>% arrange(p_val_adj) %>% slice_head(n = 8)
    plt <- plt + geom_text_repel(data = top_label, aes(label = gene), size = 3)
    
    ggsave(filename = paste0(out_prefix, "_volcano.png"), plot = plt, width = 6, height = 4, dpi = 200)
    log_msg("  Wrote volcano: %s", paste0(out_prefix, "_volcano.png"))
  } # samples
} # datasets

log_msg("DE analysis finished. Outputs in: ", out_dir)

