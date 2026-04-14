# ==============================================================================
# MULTI-VOLUME RESEARCH ARCHIVER & METADATA EXTRACTOR
# ==============================================================================
# Purpose: 
# 1. Scan and inspect Seurat objects (Dimensions, Layers, Reductions).
# 2. Split files into multiple ZIP volumes (Max 500MB per zip).
# 3. Generate a single Master Manifest tracking where every file ended up.
# ==============================================================================

# 1. SETUP
# ------------------------------------------------------------------------------
# Update this path to your project folder
work_dir <- "D:/Dr.Das_Placenta_SpatialAnalysis/Placenta_ViciousCycle_Pipeline_v2/COMPLETE_FIXED_PIPELINE"
setwd(work_dir)

# Configuration
FILE_SIZE_LIMIT_MB <- 100      # Skip individual files bigger than this
CHUNK_LIMIT_MB     <- 200      # Max size per Zip file (Leaving 5MB buffer for safety)

MANIFEST_FILE <- "MASTER_PROJECT_MANIFEST.txt"

# Install/Load required packages
if (!require("zip")) install.packages("zip")
if (!require("Seurat")) install.packages("Seurat")
library(zip)
library(Seurat)

# Constants for Byte conversion
BYTES_PER_MB <- 1024 * 1024
FILE_LIMIT_BYTES <- FILE_SIZE_LIMIT_MB * BYTES_PER_MB
CHUNK_LIMIT_BYTES <- CHUNK_LIMIT_MB * BYTES_PER_MB

# 2. HELPER FUNCTION: INSPECT SEURAT OBJECTS
# ------------------------------------------------------------------------------
get_seurat_info <- function(file_path) {
  tryCatch({
    obj <- readRDS(file_path)
    
    lines <- c(
      paste0("  [OBJECT]: ", basename(file_path)),
      paste0("  [CLASS]: ", class(obj)[1]),
      paste0("  [DIMENSIONS]: ", nrow(obj), " features x ", ncol(obj), " cells"),
      paste0("  [ACTIVE ASSAY]: ", DefaultAssay(obj)),
      paste0("  [LAYERS]: ", paste(Layers(obj), collapse = ", "))
    )
    
    reds <- Reductions(obj)
    if (length(reds) > 0) {
      lines <- c(lines, paste0("  [REDUCTIONS CALCULATED]: ", paste(reds, collapse = ", ")))
      for (r in reds) {
        coords <- Embeddings(obj, reduction = r)
        preview <- head(coords, 10)
        lines <- c(lines, paste0("    -> PREVIEW: First 10 rows of ", r, " coordinates:"))
        capture_preview <- capture.output(print(preview))
        lines <- c(lines, paste0("       ", capture_preview))
      }
    } else {
      lines <- c(lines, "  [REDUCTIONS]: None found")
    }
    
    lines <- c(lines, paste0("  [META DATA COLUMNS]: ", paste(colnames(obj@meta.data), collapse = ", ")))
    separator <- paste0("  ", paste(rep("-", 50), collapse = ""))
    lines <- c(lines, separator) 
    return(lines)
    
  }, error = function(e) {
    return(paste0("  [ERROR READING OBJECT]: ", e$message))
  })
}

# 3. SCAN AND ORGANIZE FILES
# ------------------------------------------------------------------------------
message("Scanning directory: ", work_dir)

# Get all files
all_files <- list.files(work_dir, recursive = TRUE, full.names = FALSE)
# Exclude self (script outputs)
all_files <- all_files[!grepl("Project_Part_|MASTER_PROJECT_MANIFEST", all_files)]

full_paths <- file.path(work_dir, all_files)
file_info <- file.info(full_paths)

# Lists for Manifest
manifest_rows <- c()
seurat_details <- c()
files_to_zip <- list() # Logic: List of lists. files_to_zip[[1]] = files for zip 1

# Variables for Chunking Logic
current_chunk_idx <- 1
current_chunk_size <- 0
current_chunk_files <- c()

message("Analyzing ", length(all_files), " files for packing...")

for (i in seq_along(all_files)) {
  fname <- all_files[i]
  fsize <- file_info$size[i]
  fsize_mb <- fsize / BYTES_PER_MB
  
  status_msg <- ""
  location_msg <- ""
  
  # --- LOGIC: IS FILE TOO BIG TO ZIP AT ALL? ---
  if (fsize > FILE_LIMIT_BYTES) {
    status_msg <- "[SKIPPED - TOO LARGE]"
    location_msg <- "Not Archived"
    
  } else {
    # --- LOGIC: BIN PACKING (SPLIT INTO CHUNKS) ---
    
    # Check if adding this file exceeds the 500MB chunk limit
    if ((current_chunk_size + fsize) > CHUNK_LIMIT_BYTES) {
      # Save current batch
      files_to_zip[[current_chunk_idx]] <- current_chunk_files
      
      # Start new batch
      current_chunk_idx <- current_chunk_idx + 1
      current_chunk_size <- 0
      current_chunk_files <- c()
    }
    
    # Add file to current batch
    current_chunk_files <- c(current_chunk_files, fname)
    current_chunk_size <- current_chunk_size + fsize
    
    status_msg <- "[ZIPPED]"
    location_msg <- paste0("Project_Part_", current_chunk_idx, ".zip")
  }
  
  # --- MANIFEST ENTRY ---
  entry <- sprintf("%-20s | %-20s | %8.2f MB | %s", status_msg, location_msg, fsize_mb, fname)
  manifest_rows <- c(manifest_rows, entry)
  
  # --- SEURAT INSPECTION ---
  if (grepl("\\.rds$", fname, ignore.case = TRUE)) {
    message("  Inspecting Seurat object: ", fname)
    obj_info <- get_seurat_info(full_paths[i])
    seurat_details <- c(seurat_details, obj_info)
  }
}

# Save the final partial batch
if (length(current_chunk_files) > 0) {
  files_to_zip[[current_chunk_idx]] <- current_chunk_files
}

# 4. GENERATE MANIFEST
# ------------------------------------------------------------------------------
header <- c(
  "================================================================================",
  " MASTER PROJECT MANIFEST",
  paste(" Generated on:", Sys.time()),
  paste(" Root Directory:", work_dir),
  "================================================================================",
  sprintf("%-20s | %-20s | %-11s | %s", "STATUS", "LOCATION", "SIZE", "FILENAME"),
  paste(rep("-", 100), collapse = "")
)

section_header <- c(
  "", "",
  "================================================================================",
  " DETAILED SEURAT OBJECT INSPECTION",
  "================================================================================",
  ""
)

final_text <- c(header, manifest_rows, section_header, seurat_details)
writeLines(final_text, MANIFEST_FILE)
message("Manifest created: ", MANIFEST_FILE)

# 5. CREATE ZIP VOLUMES
# ------------------------------------------------------------------------------
message("----------------------------------------------------------------")
message("Starting compression of ", length(files_to_zip), " volume(s)...")

for (idx in seq_along(files_to_zip)) {
  zip_name <- paste0("Project_Part_", idx, ".zip")
  files_in_batch <- files_to_zip[[idx]]
  
  # Add the manifest to the FIRST zip only (so it's easy to find)
  if (idx == 1) {
    files_in_batch <- c(files_in_batch, MANIFEST_FILE)
  }
  
  message("Creating ", zip_name, " containing ", length(files_in_batch), " files...")
  
  zip::zip(
    zipfile = zip_name,
    files = files_in_batch,
    mode = "mirror"
  )
}

message("----------------------------------------------------------------")
message("SUCCESS! All valid files have been split into ", length(files_to_zip), " zip volumes.")