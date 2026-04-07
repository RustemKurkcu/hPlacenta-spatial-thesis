# ======================================================================
# scripts/07_export/07A_export_shareable_outputs.R
# Collects key outputs into a small, shareable folder and zips it.
#
# This is useful for sending to advisors without including huge Seurat objects.
# ======================================================================

suppressPackageStartupMessages({
  library(fs)
})

source("config/config.R")
source("scripts/R/utils.R")

ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "07A_export_shareable_outputs.log")

share_dir <- file.path("output", "shareable_package")
if (dir.exists(share_dir)) fs::dir_delete(share_dir)
fs::dir_create(share_dir)

# Copy tables and figures
if (dir.exists(DIR_TABLES)) fs::dir_copy(DIR_TABLES, file.path(share_dir, "tables"))
if (dir.exists(DIR_FIGURES)) fs::dir_copy(DIR_FIGURES, file.path(share_dir, "figures"))

# Copy docs (methods templates, etc.)
if (dir.exists(DIR_DOCS)) fs::dir_copy(DIR_DOCS, file.path(share_dir, "docs"))

# Include README + config for reproducibility
fs::file_copy("README.md", file.path(share_dir, "README.md"), overwrite = TRUE)
fs::dir_copy("config", file.path(share_dir, "config"), overwrite = TRUE)

# Zip it
zip_path <- file.path("output", "shareable_package.zip")
if (file.exists(zip_path)) file.remove(zip_path)
utils::zip(zipfile = zip_path, files = share_dir, flags = "-r9Xq")

log_msg(paste0("Created shareable package: ", zip_path), logfile)
