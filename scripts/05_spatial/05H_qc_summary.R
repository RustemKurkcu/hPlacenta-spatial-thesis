# scripts/05_spatial/05H_qc_summary.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

dir_tab <- file.path("output", "tables", "05_spatial")
if (!dir.exists(dir_tab)) stop("Directory missing: ", dir_tab)

# 1) NK coverage: find *_NK_module_coverage_qc.csv
nk_qc_files <- list.files(dir_tab, pattern = "_NK_module_coverage_qc\\.csv$", recursive=TRUE, full.names=TRUE)
nk_qc <- bind_rows(lapply(nk_qc_files, function(f) {
  df <- tryCatch(read_csv(f, show_col_types = FALSE), error=function(e) NULL)
  if (is.null(df)) return(NULL)
  df$source_file <- f
  df
}))

if (nrow(nk_qc) > 0) {
  nk_summary <- nk_qc %>% select(dataset = dataset, assay, nk_gene_set_used, nk_found_genes = nk_found_genes, nk_missing_genes = nk_missing_genes)
  write_csv(nk_summary, file.path(dir_tab, "QC_NK_coverage_summary.csv"))
  print("Wrote QC_NK_coverage_summary.csv")
} else {
  message("No NK QC files found.")
}

# 2) Global permissiveness NA check
pg_file <- file.path(dir_tab, "permissiveness_global_allcells.csv")
if (file.exists(pg_file)) {
  all_md <- read_csv(pg_file, show_col_types = FALSE)
  per_dataset_na <- all_md %>% group_by(dataset) %>% summarize(
    n_total = n(),
    n_finite = sum(is.finite(permissiveness_global)),
    frac_finite = n_finite / n_total
  )
  write_csv(per_dataset_na, file.path(dir_tab, "QC_permissiveness_global_finiteness.csv"))
  print("Wrote QC_permissiveness_global_finiteness.csv")
} else {
  message("No permissiveness_global_allcells.csv found.")
}

# 3) Hotspot/protected counts per sample (scan subfolders for *_hotspot_cells_*.csv and *_protected_cells_*.csv)
hot_files <- list.files(dir_tab, pattern = "_hotspot_cells_top.*pct\\.csv$", recursive = TRUE, full.names = TRUE)
prot_files <- list.files(dir_tab, pattern = "_protected_cells_bottom.*pct\\.csv$", recursive = TRUE, full.names = TRUE)

read_counts <- function(files, tag) {
  if (length(files) == 0) return(tibble())
  bind_rows(lapply(files, function(f) {
    df <- tryCatch(read_csv(f, show_col_types = FALSE), error = function(e) tibble())
    n <- nrow(df)
    # parse dataset/sample from path
    parts <- strsplit(f, "/")[[1]]
    # try to find dataset and sample (assumes DIR_TAB/dataset/sample/filename)
    i_tab <- which(parts == "05_spatial")
    ds <- NA; sp <- NA
    if (length(i_tab)>=1 && length(parts) >= i_tab + 2) {
      ds <- parts[i_tab + 1]
      sp <- parts[i_tab + 2]
    }
    tibble(dataset = ds, sample = sp, file = f, tag = tag, n = n)
  }))
}

hot_counts <- read_counts(hot_files, "hotspot")
prot_counts <- read_counts(prot_files, "protected")
counts <- bind_rows(hot_counts, prot_counts) %>% tidyr::pivot_wider(names_from = tag, values_from = n, values_fill = 0)
write_csv(counts, file.path(dir_tab, "QC_hotspot_protected_counts_per_sample.csv"))
print("Wrote QC_hotspot_protected_counts_per_sample.csv")