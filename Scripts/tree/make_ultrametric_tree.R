#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript make_ultrametric_tree.R <input_tree> <output_tree> [root_height]\n")
  cat("\nArguments:\n")
  cat("  input_tree   : Path to input Newick tree file (from IQ-TREE)\n")
  cat("  output_tree  : Path to output ultrametric tree file\n")
  cat("  root_height  : Root height for scaled tree (default: 1.0)\n")
  cat("\nExamples:\n")
  cat("  # Basic usage with root height = 1\n")
  cat("  Rscript make_ultrametric_tree.R tree.nwk tree_ultrametric.nwk\n")
  cat("\n  # Custom root height\n")
  cat("  Rscript make_ultrametric_tree.R tree.nwk tree_ultrametric.nwk 2.0\n")
  quit(status = 1)
}

input_tree_file <- args[[1]]
output_tree_file <- args[[2]]
root_height <- ifelse(length(args) >= 3, as.numeric(args[[3]]), 1.0)

# Check if input file exists
if (!file.exists(input_tree_file)) {
  stop("Error: Input tree file not found: ", input_tree_file)
}

message(paste(rep("=", 60), collapse = ""))
message("Making Tree Ultrametric")
message(paste(rep("=", 60), collapse = ""))
message("Input tree: ", input_tree_file)
message("Output tree: ", output_tree_file)
message("Target root height: ", root_height)
message(paste(rep("=", 60), collapse = ""))

# Load tree
message("\nLoading tree...")
tree <- read.tree(input_tree_file)

# Check if tree is already rooted
if (!is.rooted(tree)) {
  message("\nTree is not rooted. Attempting to root at midpoint...")
  tree <- midpoint.root(tree)
  message("Tree rooted at midpoint")
}

message("\nTree information:")
message("  Number of tips: ", length(tree$tip.label))
message("  Number of nodes: ", tree$Nnode)
message("  Root height (original): ", max(node.depth.edgelength(tree)))

# Check if tree is already ultrametric
is_ultrametric <- is.ultrametric(tree, tol = 1e-6)
message("  Is ultrametric: ", is_ultrametric)

if (is_ultrametric) {
  message("\nTree is already ultrametric. Scaling to root height = ", root_height, "...")
  current_root_height <- max(node.depth.edgelength(tree))
  if (current_root_height > 0) {
    scaling_factor <- root_height / current_root_height
    tree$edge.length <- tree$edge.length * scaling_factor
    message("  Scaling factor: ", scaling_factor)
  } else {
    warning("  Warning: Current root height is 0, cannot scale")
  }
} else {
  message("\nMaking tree ultrametric...")
  message("  Using method hierarchy: chronos() -> force.ultrametric(nnls) -> force.ultrametric(extend)")
  
  # Method 1: chronos() - Maximum likelihood under relaxed molecular clock (MOST ACCURATE)
  # This is the most statistically rigorous method for dating trees
  message("\n  Attempting Method 1: chronos() (maximum likelihood, relaxed clock)...")
  success <- FALSE
  tryCatch({
    # Try relaxed clock model first (most flexible)
    tree_ultra <- chronos(tree, lambda = 1, model = "relaxed", quiet = TRUE)
    if (is.ultrametric(tree_ultra, tol = 1e-6)) {
      tree <- tree_ultra
      message("  ✓ Successfully made tree ultrametric using chronos() with relaxed clock model")
      success <- TRUE
    } else {
      message("  chronos(relaxed) did not produce ultrametric tree, trying strict clock...")
      # Try strict clock model
      tree_ultra <- chronos(tree, lambda = 1, model = "strict", quiet = TRUE)
      if (is.ultrametric(tree_ultra, tol = 1e-6)) {
        tree <- tree_ultra
        message("  ✓ Successfully made tree ultrametric using chronos() with strict clock model")
        success <- TRUE
      }
    }
  }, error = function(e) {
    message("  ✗ chronos() failed: ", e$message)
  })
  
  # Method 2: force.ultrametric with NNLS (Non-Negative Least Squares) - better than extend
  if (!success) {
    message("\n  Attempting Method 2: force.ultrametric() with NNLS method...")
    tryCatch({
      tree_ultra <- force.ultrametric(tree, method = "nnls")
      if (is.ultrametric(tree_ultra, tol = 1e-6)) {
        tree <- tree_ultra
        message("  ✓ Successfully made tree ultrametric using force.ultrametric(nnls)")
        success <- TRUE
      }
    }, error = function(e) {
      message("  ✗ force.ultrametric(nnls) failed: ", e$message)
    })
  }
  
  # Method 3: force.ultrametric with extend (last resort - least accurate)
  if (!success) {
    message("\n  Attempting Method 3: force.ultrametric() with extend method (last resort)...")
    message("  WARNING: This method is less accurate and should only be used if other methods fail")
    tryCatch({
      tree_ultra <- force.ultrametric(tree, method = "extend")
      if (is.ultrametric(tree_ultra, tol = 1e-6)) {
        tree <- tree_ultra
        message("  ✓ Successfully made tree ultrametric using force.ultrametric(extend)")
        success <- TRUE
      } else {
        stop("force.ultrametric(extend) did not produce ultrametric tree")
      }
    }, error = function(e) {
      stop("All methods failed. Last error: ", e$message)
    })
  }
  
  if (!success) {
    stop("Failed to make tree ultrametric with any method")
  }
  
  # Scale to desired root height
  current_root_height <- max(node.depth.edgelength(tree))
  if (current_root_height > 0) {
    scaling_factor <- root_height / current_root_height
    tree$edge.length <- tree$edge.length * scaling_factor
    message("\nScaling tree to root height = ", root_height, "...")
    message("  Original root height: ", current_root_height)
    message("  Scaling factor: ", scaling_factor)
  } else {
    warning("  Warning: Root height is 0 after ultrametric conversion, cannot scale")
  }
}

# Verify final tree
final_root_height <- max(node.depth.edgelength(tree))
message("\nFinal tree information:")
message("  Root height: ", final_root_height)
message("  Is ultrametric: ", is.ultrametric(tree, tol = 1e-6))

# Check for zero-length branches
zero_branches <- sum(tree$edge.length == 0)
if (zero_branches > 0) {
  message("  Warning: Tree has ", zero_branches, " zero-length branches")
  message("  (These may cause issues in some analyses)")
}

# Save tree
message("\nSaving ultrametric tree to: ", output_tree_file)
write.tree(tree, file = output_tree_file)

# Verify saved tree
tree_check <- read.tree(output_tree_file)
message("  Tree saved successfully")
message("  Verification - Root height: ", max(node.depth.edgelength(tree_check)))
message("  Verification - Is ultrametric: ", is.ultrametric(tree_check, tol = 1e-6))

message("\n", paste(rep("=", 60), collapse = ""))
message("Done!")
message(paste(rep("=", 60), collapse = ""))

