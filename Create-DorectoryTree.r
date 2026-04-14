print_dir_tree_to_file <- function(root = ".", output_file = "directory_tree.txt") {
  if (!dir.exists(root)) {
    stop("Directory does not exist: ", root)
  }
  
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  output_path <- file.path(root, output_file)
  
  tree_lines <- character()
  
  add_line <- function(x) {
    tree_lines <<- c(tree_lines, x)
  }
  
  print_branch <- function(path, prefix = "") {
    items <- list.files(
      path,
      full.names = TRUE,
      all.files = FALSE,
      no.. = TRUE
    )
    
    if (length(items) == 0) return(invisible(NULL))
    
    info <- file.info(items)
    dirs  <- items[info$isdir %in% TRUE]
    files <- items[!(info$isdir %in% TRUE)]
    
    # folders first, then files
    items <- c(sort(dirs), sort(files))
    
    for (i in seq_along(items)) {
      item <- items[i]
      is_last <- i == length(items)
      connector <- if (is_last) "└── " else "├── "
      
      add_line(paste0(prefix, connector, basename(item)))
      
      if (isTRUE(file.info(item)$isdir)) {
        new_prefix <- paste0(prefix, if (is_last) "    " else "│   ")
        print_branch(item, new_prefix)
      }
    }
  }
  
  add_line(basename(root))
  print_branch(root)
  
  writeLines(tree_lines, con = output_path, useBytes = TRUE)
  
  message("Directory tree saved to: ", output_path)
  
  invisible(tree_lines)
}

#############
# print_dir_tree_to_file("E:/hPlacenta-architecture")
#######