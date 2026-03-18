#!/usr/bin/env Rscript

# Universal script to create 3D phylomorphospace plots from any ultrametric tree and PCA scores
# Usage: Rscript create_3d_phylomorphospace_universal.R --tree <tree_file> --pca <pca_file> --output <output_html> [options]
#
# Required arguments:
#   --tree, -t    : Path to ultrametric tree file (Newick format)
#   --pca, -p     : Path to PCA scores CSV file (must contain PC1 and PC2 columns)
#   --output, -o  : Path to output HTML file
#
# Optional arguments:
#   --name, -n    : Dataset name for plot title (default: basename of output file)
#   --id-col, -i  : Name of ID column in PCA file (default: auto-detect from 'clean_id', 'seq_id', or first column)
#   --metadata, -m: Comma-separated list of metadata columns to include (e.g., "family,genus,species")
#   --help, -h    : Show this help message
#
# CLADE COLOR/SHAPE CONFIGURATION:
# =================================
# To customize colors and shapes for different clades, edit the section below.
# If left empty (NULL), all points will be displayed in grey.
#
# Format: clade_name = c("color", "shape")
# Available plotly shapes: "circle", "square", "diamond", "triangle-up", "triangle-down",
#                          "star", "cross", "x", "pentagon", "hexagon", "octagon", etc.
# Available colors: Any valid color name (e.g., "red", "blue", "purple", "#FF5733") or hex code

# Classification functions (harmonized with phylomorphospace_colored.R for primates)
get_genus <- function(seq_id) {
  if (grepl("_", seq_id)) {
    return(strsplit(seq_id, "_")[[1]][1])
  } else if (grepl("\\|", seq_id)) {
    return(strsplit(seq_id, "\\|")[[1]][1])
  } else {
    return(seq_id)
  }
}

# Build a genus -> family lookup for this dataset.
# (Used to reproduce the exact legend/shapes from the updated 2D phylomorphospace.)
genus_to_family <- c(
  # Cetacea
  "Balaenoptera" = "Cetacea",
  "Delphinapterus" = "Cetacea",
  "Delphinus" = "Cetacea",
  "Eschrichtius" = "Cetacea",
  "Eubalaena" = "Cetacea",
  "Globicephala" = "Cetacea",
  "Kogia" = "Cetacea",
  "Lagenorhynchus" = "Cetacea",
  "Orcinus" = "Cetacea",
  "Lipotes" = "Cetacea",
  "Monodon" = "Cetacea",
  "Neophocaena" = "Cetacea",
  "Phocoena" = "Cetacea",
  "Physeter" = "Cetacea",
  "Pseudorca" = "Cetacea",
  "Sousa" = "Cetacea",
  "Tursiops" = "Cetacea",

  # Bovidae
  "Bison" = "Bovidae",
  "Bos" = "Bovidae",
  "Bubalus" = "Bovidae",
  "Budorcas" = "Bovidae",
  "Capra" = "Bovidae",
  "Capricornis" = "Bovidae",
  "Oryx" = "Bovidae",
  "Ovibos" = "Bovidae",
  "Ovis" = "Bovidae",
  "Pantholops" = "Bovidae",

  # Camelidae
  "Camelus" = "Camelidae",
  "Vicugna" = "Camelidae",

  # Cervidae
  "Cervus" = "Cervidae",
  "Dama" = "Cervidae",
  "Muntiacus" = "Cervidae",
  "Odocoileus" = "Cervidae",
  "Rangifer" = "Cervidae",

  # Hippopotamidae
  "Hippopotamus" = "Hippopotamidae",

  # Hominidae
  "Homo" = "Hominidae",
  "HomoIPR003599IG" = "Hominidae",

  # Moschidae
  "Moschus" = "Moschidae",

  # Suidae
  "Phacochoerus" = "Suidae",
  "Sus" = "Suidae"
)

classify_family <- function(seq_id) {
  genus <- get_genus(seq_id)
  fam <- unname(genus_to_family[genus])
  if (length(fam) == 0 || is.na(fam)) return("Unknown")
  return(as.character(fam))
}

# Species groups matching the 3-color palette in the updated 2D phylomorphospace.
classify_species_group <- function(seq_id) {
  fam <- classify_family(seq_id)
  if (fam == "Cetacea") return("Cetacean")
  if (fam == "Hominidae") return("Other")
  if (fam %in% c("Bovidae", "Camelidae", "Cervidae", "Hippopotamidae", "Moschidae", "Suidae")) {
    return("Terrestrial_Artiodactyl")
  }
  return("Other")
}

# Map R pch shapes to plotly shapes
# R pch: 21=circle, 22=square, 23=diamond, 24=triangle-up, 25=triangle-down
#        0=open square, 1=open circle, 2=open triangle up, 3=plus, 4=X, 5=diamond, 6=triangle down
# Plotly symbols: circle, square, diamond, triangle, triangle-down, cross, x, star, etc.
# Note: plotly uses "triangle" not "triangle-up"
# For open shapes, we'll use filled equivalents in plotly (plotly doesn't have "open" variants)
pch_to_plotly <- function(pch) {
  shapes <- c("21" = "circle", "22" = "square", "23" = "diamond", 
              "24" = "triangle", "25" = "triangle-down",
              "0" = "square",    # Open square -> square (plotly doesn't have open shapes)
              "1" = "circle",    # Open circle -> circle
              "2" = "triangle",  # Open triangle -> triangle
              "3" = "cross",     # Plus -> cross
              "4" = "x",         # X -> x
              "5" = "diamond",   # Diamond -> diamond
              "6" = "triangle-down",  # Triangle down -> triangle-down
              "8" = "star")      # Asterisk -> star
  if (as.character(pch) %in% names(shapes)) {
    return(shapes[as.character(pch)])
  }
  return("circle")
}

# Family shapes (from phylomorphospace_colored.R)
# R pch mapping used by base graphics:
#   21 circle, 22 square, 23 diamond, 24 triangle-up, 25 triangle-down, 8 asterisk/star
family_shapes_pch <- c(
  "Bovidae" = 24,          # triangle up
  "Camelidae" = 21,       # filled circle
  "Cervidae" = 23,        # filled diamond
  "Cetacea" = 21,         # filled circle
  "Hippopotamidae" = 8,   # asterisk
  "Hominidae" = 21,       # filled circle
  "Moschidae" = 25,       # triangle down
  "Suidae" = 22           # square
)

# Species group colors matching the target 2D phylomorphospace legend.
group_colors <- c(
  "Cetacean" = "#0072B2",                  # Blue
  "Terrestrial_Artiodactyl" = "#E69F00", # Orange
  "Other" = "#F0E442"                     # Yellow
)

# Ensure pandoc is in PATH (for conda installations)
conda_bin <- Sys.getenv("CONDA_PREFIX")
if (conda_bin != "") {
  conda_bin_path <- file.path(conda_bin, "bin")
  if (dir.exists(conda_bin_path)) {
    current_path <- Sys.getenv("PATH")
    Sys.setenv(PATH = paste(conda_bin_path, current_path, sep = ":"))
  }
}
# Also try common conda locations
common_conda_paths <- c(
  file.path(Sys.getenv("HOME"), "miniforge3", "bin"),
  file.path(Sys.getenv("HOME"), "anaconda3", "bin"),
  file.path(Sys.getenv("HOME"), "conda", "bin")
)
for (conda_path in common_conda_paths) {
  if (dir.exists(conda_path) && !grepl(conda_path, Sys.getenv("PATH"))) {
    current_path <- Sys.getenv("PATH")
    Sys.setenv(PATH = paste(conda_path, current_path, sep = ":"))
    break
  }
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(phytools)
  library(geiger)
  library(castor)
  library(geometry)
  library(plotly)
  library(htmlwidgets)
})

# Helper function to convert phylomorphospace to data frame
phylomorphospace_to_df <- function(phylomorphospace_obj, tree) {
  # Extract coordinates from phylomorphospace object
  # phylomorphospace returns a list with $xx (x-coords) and $yy (y-coords) for nodes
  # We need to extract edge information
  
  # Get node coordinates
  node_coords <- data.frame(
    node = 1:(length(tree$tip.label) + tree$Nnode),
    x = c(phylomorphospace_obj$xx[1:length(tree$tip.label)], 
          phylomorphospace_obj$xx[(length(tree$tip.label)+1):length(phylomorphospace_obj$xx)]),
    y = c(phylomorphospace_obj$yy[1:length(tree$tip.label)], 
          phylomorphospace_obj$yy[(length(tree$tip.label)+1):length(phylomorphospace_obj$yy)])
  )
  
  # Create edge data frame
  edge_df <- data.frame(
    nodestart = tree$edge[, 1],
    nodestop = tree$edge[, 2]
  )
  
  # Add coordinates
  edge_df <- merge(edge_df, node_coords, by.x = "nodestart", by.y = "node", all.x = TRUE)
  edge_df$xstart <- edge_df$x
  edge_df$ystart <- edge_df$y
  edge_df <- edge_df[, c("nodestart", "nodestop", "xstart", "ystart")]
  
  edge_df <- merge(edge_df, node_coords, by.x = "nodestop", by.y = "node", all.x = TRUE)
  edge_df$xstop <- edge_df$x
  edge_df$ystop <- edge_df$y
  edge_df <- edge_df[, c("nodestart", "nodestop", "xstart", "ystart", "xstop", "ystop")]
  
  # Add node type and seq_id
  bound <- length(tree$tip.label)
  edge_df$nodeType <- ifelse(edge_df$nodestop <= bound, "Tip", "Internal")
  edge_df$seq_id <- ifelse(edge_df$nodestop <= bound, 
                           tree$tip.label[edge_df$nodestop], 
                           "NA")
  
  return(edge_df)
}

# Fix for morphospaceCreator to handle missing heights
morphospaceCreator_fixed <- function(phylomorphospace_object, tree) {
  df <- phylomorphospace_to_df(phylomorphospace_object, tree)
  
  # Use nodeHeights from phytools - this gives heights for each edge
  # nodeHeights returns a matrix with 2 columns: [start_height, end_height] for each edge
  node_heights <- nodeHeights(tree)
  max_height <- max(node_heights)
  
  # Initialize zstart and zstop
  df$zstart <- NA_real_
  df$zstop <- NA_real_
  
  # Match each edge in df to the corresponding edge in the tree
  for (i in 1:nrow(df)) {
    # Find the corresponding edge in the tree
    edge_match <- which(tree$edge[,1] == df$nodestart[i] & tree$edge[,2] == df$nodestop[i])
    if (length(edge_match) > 0) {
      # Get original heights from nodeHeights
      # nodeHeights gives distance FROM root:
      # - Column 1: distance from root to start node (closer to root = smaller value, root = 0)
      # - Column 2: distance from root to end node (closer to tips = larger value, tips = max_height)
      orig_start <- node_heights[edge_match[1], 1]  # Distance from root to start node
      orig_stop <- node_heights[edge_match[1], 2]    # Distance from root to end node
      
      # We want: root (distance 0) at z=1 (top), tips (distance max_height) at z=0 (bottom)
      # So: z = 1 - (distance_from_root / max_height)
      if (max_height > 0) {
        df$zstart[i] <- 1 - (orig_start / max_height)  # Root (0) -> 1, internal nodes -> intermediate
        df$zstop[i] <- 1 - (orig_stop / max_height)    # Tips (max_height) -> 0, internal nodes -> intermediate
      } else {
        df$zstart[i] <- 1
        df$zstop[i] <- 0
      }
    } else {
      # Fallback: if edge not found, try to calculate from node depths
      warning("Edge not found for row ", i, ". Using fallback calculation.")
      df$zstart[i] <- 1
      df$zstop[i] <- 0
    }
  }
  
  df2 <- df
  
  return(df2)
}

# Function to parse command-line arguments
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Simple argument parser
  tree_file <- NULL
  pca_file <- NULL
  output_file <- NULL
  dataset_name <- NULL
  id_col <- NULL
  metadata_cols <- NULL
  
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    
    if (arg %in% c("--tree", "-t")) {
      if (i + 1 <= length(args)) {
        tree_file <- args[i + 1]
        i <- i + 2
      } else {
        stop("--tree requires a file path")
      }
    } else if (arg %in% c("--pca", "-p")) {
      if (i + 1 <= length(args)) {
        pca_file <- args[i + 1]
        i <- i + 2
      } else {
        stop("--pca requires a file path")
      }
    } else if (arg %in% c("--output", "-o")) {
      if (i + 1 <= length(args)) {
        output_file <- args[i + 1]
        i <- i + 2
      } else {
        stop("--output requires a file path")
      }
    } else if (arg %in% c("--name", "-n")) {
      if (i + 1 <= length(args)) {
        dataset_name <- args[i + 1]
        i <- i + 2
      } else {
        stop("--name requires a name")
      }
    } else if (arg %in% c("--id-col", "-i")) {
      if (i + 1 <= length(args)) {
        id_col <- args[i + 1]
        i <- i + 2
      } else {
        stop("--id-col requires a column name")
      }
    } else if (arg %in% c("--metadata", "-m")) {
      if (i + 1 <= length(args)) {
        metadata_cols <- strsplit(args[i + 1], ",")[[1]]
        metadata_cols <- trimws(metadata_cols)
        i <- i + 2
      } else {
        stop("--metadata requires a comma-separated list")
      }
    } else if (arg %in% c("--help", "-h")) {
      cat("Usage: Rscript create_3d_phylomorphospace_universal.R --tree <tree_file> --pca <pca_file> --output <output_html> [options]\n\n")
      cat("Required arguments:\n")
      cat("  --tree, -t    : Path to ultrametric tree file (Newick format)\n")
      cat("  --pca, -p     : Path to PCA scores CSV file (must contain PC1 and PC2 columns)\n")
      cat("  --output, -o  : Path to output HTML file\n\n")
      cat("Optional arguments:\n")
      cat("  --name, -n    : Dataset name for plot title (default: basename of output file)\n")
      cat("  --id-col, -i  : Name of ID column in PCA file (default: auto-detect)\n")
      cat("  --metadata, -m: Comma-separated list of metadata columns (e.g., \"family,genus,species\")\n")
      cat("  --help, -h    : Show this help message\n")
      quit(status = 0)
    } else {
      stop("Unknown argument: ", arg, "\nUse --help for usage information")
    }
  }
  
  # Check required arguments
  if (is.null(tree_file)) {
    stop("--tree is required. Use --help for usage information")
  }
  if (is.null(pca_file)) {
    stop("--pca is required. Use --help for usage information")
  }
  if (is.null(output_file)) {
    stop("--output is required. Use --help for usage information")
  }
  
  # Set default dataset name if not provided
  if (is.null(dataset_name)) {
    dataset_name <- tools::file_path_sans_ext(basename(output_file))
  }
  
  return(list(
    tree_file = tree_file,
    pca_file = pca_file,
    output_file = output_file,
    dataset_name = dataset_name,
    id_col = id_col,
    metadata_cols = metadata_cols
  ))
}

# Function to auto-detect ID column
detect_id_column <- function(pca_data) {
  # Try common ID column names
  id_candidates <- c("clean_id", "seq_id", "id", "ID", "sample_id", "sample", "taxon", "species", "original_sequence_id", "row_id")
  
  for (candidate in id_candidates) {
    if (candidate %in% colnames(pca_data)) {
      return(candidate)
    }
  }
  
  # If first column is empty or named "X", use row names instead
  first_col <- colnames(pca_data)[1]
  if (first_col == "" || first_col == "X") {
    return(NULL)  # Signal to use row names
  }
  
  # If none found, use first column that's not PC1, PC2, pPC1, pPC2, or any other PC column
  pc_cols <- grep("^PC[0-9]+$|^pPC[0-9]+$", colnames(pca_data), value = TRUE)
  non_pc_cols <- setdiff(colnames(pca_data), pc_cols)
  if (length(non_pc_cols) > 0) {
    return(non_pc_cols[1])
  }
  
  # Last resort: use row names
  return(NULL)
}

# Function to create 3D phylomorphospace plot
create_3d_phylomorphospace <- function(pca_scores_file, tree_file, output_html, dataset_name, id_col = NULL, metadata_cols = NULL) {
  
  message("\n", paste(rep("=", 60), collapse = ""))
  message("Creating 3D Phylomorphospace for: ", dataset_name)
  message(paste(rep("=", 60), collapse = ""))
  
  # Read PCA scores
  message("\nReading PCA scores from: ", pca_scores_file)
  if (!file.exists(pca_scores_file)) {
    stop("PCA scores file not found: ", pca_scores_file)
  }
  # Try reading with row names first (if first column is empty/contains IDs)
  pca_data <- tryCatch({
    read.csv(pca_scores_file, stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
  }, error = function(e) {
    # Fallback: read without row names
    read.csv(pca_scores_file, stringsAsFactors = FALSE)
  })
  
  # If we have row names, add them as a column for matching
  if (length(rownames(pca_data)) > 0 && !all(rownames(pca_data) == as.character(1:nrow(pca_data)))) {
    pca_data$row_id <- rownames(pca_data)
  }
  
  # Check if we have the required columns (try PC1/PC2 or pPC1/pPC2)
  has_pc1_pc2 <- "PC1" %in% colnames(pca_data) && "PC2" %in% colnames(pca_data)
  has_ppc1_ppc2 <- "pPC1" %in% colnames(pca_data) && "pPC2" %in% colnames(pca_data)
  
  if (!has_pc1_pc2 && !has_ppc1_ppc2) {
    stop("PCA scores file must contain PC1/PC2 or pPC1/pPC2 columns")
  }
  
  # Use pPC columns if available, otherwise PC columns
  pc1_col <- if (has_ppc1_ppc2) "pPC1" else "PC1"
  pc2_col <- if (has_ppc1_ppc2) "pPC2" else "PC2"
  pc1_label <- if (has_ppc1_ppc2) "pPC1" else "PC1"
  pc2_label <- if (has_ppc1_ppc2) "pPC2" else "PC2"
  message("Using columns: ", pc1_col, " and ", pc2_col)
  
  # Auto-detect ID column if not provided
  if (is.null(id_col)) {
    # First check if we have row_id (from row names)
    if ("row_id" %in% colnames(pca_data)) {
      id_col <- "row_id"
      message("Using row names as sequence IDs")
    } else {
    id_col <- detect_id_column(pca_data)
      if (is.null(id_col)) {
        # Use row names as IDs
        pca_data$row_id <- rownames(pca_data)
        id_col <- "row_id"
        message("Using row names as sequence IDs")
      } else {
    message("Auto-detected ID column: ", id_col)
      }
    }
  } else {
    if (!id_col %in% colnames(pca_data)) {
      stop("Specified ID column '", id_col, "' not found in PCA data")
    }
    message("Using ID column: ", id_col)
  }
  
  # Handle duplicate IDs - keep only first occurrence for phylomorphospace
  if (any(duplicated(pca_data[[id_col]]))) {
    message("Warning: Found duplicate IDs. Keeping first occurrence of each...")
    dup_ids <- duplicated(pca_data[[id_col]])
    pca_data <- pca_data[!dup_ids, ]
    message("  Removed ", sum(dup_ids), " duplicate entries")
  }
  
  # Filter out rows with NA or infinite values
  valid_rows <- !is.na(pca_data[[pc1_col]]) & !is.na(pca_data[[pc2_col]]) & 
                is.finite(pca_data[[pc1_col]]) & is.finite(pca_data[[pc2_col]])
  if (sum(!valid_rows) > 0) {
    message("Warning: Removing ", sum(!valid_rows), " rows with NA/Inf values")
    pca_data <- pca_data[valid_rows, ]
  }
  
  # Prepare data matrix with PC1/PC2 or pPC1/pPC2
  pc_matrix <- pca_data[, c(pc1_col, pc2_col)]
  colnames(pc_matrix) <- c("PC1", "PC2")  # Standardize column names for phylomorphospace
  rownames(pc_matrix) <- pca_data[[id_col]]
  
  # Read tree
  message("Reading tree from: ", tree_file)
  if (!file.exists(tree_file)) {
    stop("Tree file not found: ", tree_file)
  }
  tree <- read.tree(tree_file)
  
  # Check if tree is ultrametric
  if (!is.ultrametric(tree)) {
    message("Warning: Tree is not ultrametric. Attempting to force ultrametric...")
    tree <- force.ultrametric(tree, method = "extend")
  }
  
  # Normalize sequence IDs for matching (matches phylomorphospace_colored.R)
  normalize_id <- function(seq_id) {
    # Replace spaces with underscores
    normalized <- gsub(" ", "_", seq_id)
    # Remove pipe separators (|) and everything after them (e.g., version numbers)
    normalized <- gsub("\\|.*", "", normalized)
    # Remove dots from version numbers (e.g., XP_042112342.1 -> XP_0421123421)
    # Pattern: number followed by dot followed by number at end of word or before underscore/end
    normalized <- gsub("([0-9])\\.([0-9]+)([^0-9]|$)", "\\1\\2\\3", normalized)
    # Remove IPR003599IG_3c suffix and range (e.g., IPR003599IG_3c26-128) from PCA row names
    # This matches what tree tip labels have after removing pipes
    normalized <- gsub("IPR003599IG_3c[0-9]+-[0-9]+$", "", normalized)
    # Normalize keywords followed by colons/underscores: PREDICTED, LOW_QUALITY_PROTEIN, TPA
    normalized <- gsub("PREDICTED[_:]+", "PREDICTED_", normalized)
    normalized <- gsub("LOW_QUALITY_PROTEIN[_:]+", "LOW_QUALITY_PROTEIN_", normalized)
    normalized <- gsub("TPA[_:]+", "TPA_", normalized)
    # Normalize partial annotations: __partial or ,_partial -> _partial
    normalized <- gsub("_{1,2}partial|,_partial", "_partial", normalized)
    # Remove commas (they appear in some names like "group_C,_group_5" -> "group_C_group_5")
    normalized <- gsub(",", "", normalized)
    # Remove any remaining double underscores
    normalized <- gsub("_{2,}", "_", normalized)
    return(normalized)
  }
  
  # Normalize IDs before matching
  normalized_tree_tips <- sapply(tree$tip.label, normalize_id)
  normalized_pca_ids <- sapply(rownames(pc_matrix), normalize_id)
  
  # Debug: show some examples
  message("\nID matching debug:")
  message("  Sample tree tip labels (first 3): ", paste(head(tree$tip.label, 3), collapse = ", "))
  message("  Sample normalized tree tips (first 3): ", paste(head(normalized_tree_tips, 3), collapse = ", "))
  message("  Sample PCA row names (first 3): ", paste(head(rownames(pc_matrix), 3), collapse = ", "))
  message("  Sample normalized PCA IDs (first 3): ", paste(head(normalized_pca_ids, 3), collapse = ", "))
  
  # Check name overlap and prune tree to match PCA data
  test <- pc_matrix$PC1
  names(test) <- normalized_pca_ids  # Use normalized IDs for matching
  
  # Find common tips using normalized IDs
  common_tips_normalized <- intersect(normalized_tree_tips, normalized_pca_ids)
  message("  Common tips found: ", length(common_tips_normalized))
  if (length(common_tips_normalized) == 0) {
    # Show some examples of what we're trying to match
    message("\n  First 5 normalized tree tips: ", paste(head(normalized_tree_tips, 5), collapse = ", "))
    message("  First 5 normalized PCA IDs: ", paste(head(normalized_pca_ids, 5), collapse = ", "))
    stop("No matching tips between tree and PCA data. Check sequence IDs.")
  }
  
  # Create mapping from normalized IDs back to original IDs
  tree_normalized_to_original <- setNames(tree$tip.label, normalized_tree_tips)
  pca_normalized_to_original <- setNames(rownames(pc_matrix), normalized_pca_ids)
  
  # Get original IDs for the common tips
  common_tips_original <- tree_normalized_to_original[common_tips_normalized]
  common_pca_ids_original <- pca_normalized_to_original[common_tips_normalized]
  
  message("\nName matching (with normalized IDs):")
  message("  Tree tips: ", length(tree$tip.label))
  message("  PCA sequences: ", length(rownames(pc_matrix)))
  message("  Common (normalized): ", length(common_tips_normalized))
  
  if (length(common_tips_normalized) < length(tree$tip.label)) {
    message("Pruning tree to match PCA data...")
    tree <- keep.tip(tree, common_tips_original)
  }
  
  if (length(common_tips_normalized) < length(rownames(pc_matrix))) {
    message("Filtering PCA data to match tree...")
    # Filter using original IDs
    pc_matrix <- pc_matrix[common_pca_ids_original, , drop = FALSE]
    pca_data <- pca_data[pca_data[[id_col]] %in% common_pca_ids_original, ]
    # Update rownames to match tree tip labels (original IDs)
    rownames(pc_matrix) <- common_pca_ids_original
  } else {
    # Reorder matrix to match tree tip order
    # Map normalized tree tips back to original PCA IDs
    tree_to_pca_map <- setNames(common_pca_ids_original, common_tips_normalized)
    tree_tips_normalized <- sapply(tree$tip.label, normalize_id)
    matching_pca_ids <- tree_to_pca_map[tree_tips_normalized]
    pc_matrix <- pc_matrix[matching_pca_ids, , drop = FALSE]
    rownames(pc_matrix) <- tree$tip.label  # Use tree tip labels as rownames
  }
  
  # Final check for valid data
  if (any(!is.finite(pc_matrix$PC1)) || any(!is.finite(pc_matrix$PC2))) {
    stop("Invalid PC values detected after filtering. Check PCA scores.")
  }
  
  message("\nFinal dataset:")
  message("  Taxa: ", nrow(pc_matrix))
  message("  PC1 range: [", round(min(pc_matrix$PC1), 3), ", ", round(max(pc_matrix$PC1), 3), "]")
  message("  PC2 range: [", round(min(pc_matrix$PC2), 3), ", ", round(max(pc_matrix$PC2), 3), "]")
  
  # Calculate ancestral states manually
  message("\nCalculating ancestral states...")
  # Handle zero-length branches
  if (any(tree$edge.length == 0)) {
    message("  Warning: Tree has ", sum(tree$edge.length == 0), " zero-length branches")
    message("  Adding small value to zero-length branches...")
    tree$edge.length[tree$edge.length == 0] <- min(tree$edge.length[tree$edge.length > 0], na.rm = TRUE) * 0.001
  }
  
  # Calculate ancestral states for PC1 and PC2
  anc_pc1 <- fastAnc(tree, pc_matrix[, 1])
  anc_pc2 <- fastAnc(tree, pc_matrix[, 2])
  
  # Create edge data frame with coordinates
  message("Processing 3D coordinates...")
  edge_df <- data.frame(
    nodestart = tree$edge[, 1],
    nodestop = tree$edge[, 2]
  )
  
  # Get node coordinates
  n_tips <- length(tree$tip.label)
  node_x <- c(pc_matrix[, 1], anc_pc1)
  node_y <- c(pc_matrix[, 2], anc_pc2)
  
  # Add coordinates to edges
  edge_df$xstart <- node_x[edge_df$nodestart]
  edge_df$ystart <- node_y[edge_df$nodestart]
  edge_df$xstop <- node_x[edge_df$nodestop]
  edge_df$ystop <- node_y[edge_df$nodestop]
  
  # Add node type and seq_id
  edge_df$nodeType <- ifelse(edge_df$nodestop <= n_tips, "Tip", "Internal")
  edge_df$seq_id <- ifelse(edge_df$nodestop <= n_tips, 
                           tree$tip.label[edge_df$nodestop], 
                           "NA")
  
  ggplot_phylomorphospace_object <- edge_df
  
  # Calculate z-axis (time/branch length) from root to tips
  node_heights <- nodeHeights(tree)
  max_height <- max(node_heights)
  
  # Map node heights to z coordinates (root at z=1, tips at z=0)
  node_z <- numeric(length(node_x))
  for (i in 1:n_tips) {
    # Tips: z = 0
    node_z[i] <- 0
  }
  for (i in (n_tips+1):length(node_x)) {
    # Internal nodes: find height from root
    node_id <- i
    # Find edges ending at this node
    edge_idx <- which(tree$edge[, 2] == node_id)
    if (length(edge_idx) > 0) {
      # Get height from nodeHeights (distance from root)
      height_from_root <- node_heights[edge_idx[1], 2]
      node_z[i] <- 1 - (height_from_root / max_height)
    } else {
      node_z[i] <- 1  # Root
    }
  }
  
  # Add z coordinates to edges
  ggplot_phylomorphospace_object$zstart <- node_z[edge_df$nodestart]
  ggplot_phylomorphospace_object$zstop <- node_z[edge_df$nodestop]
  
  max_z <- max(c(ggplot_phylomorphospace_object$zstart, ggplot_phylomorphospace_object$zstop), na.rm = TRUE)
  message("Tree height range: root at z=", round(max_z, 4), ", tips at z=0")
  
  # Add metadata if available
  # The seq_id in phylomorphospace object comes from tree tip labels
  # We need to match it with pca_data using normalized IDs
  # Create normalized seq_id for joining
  ggplot_phylomorphospace_object$seq_id_normalized <- sapply(ggplot_phylomorphospace_object$seq_id, normalize_id)
  pca_data$seq_id_normalized <- sapply(pca_data[[id_col]], normalize_id)
  
  # Determine which metadata columns to use
  if (!is.null(metadata_cols)) {
    # Use specified metadata columns
    available_metadata <- intersect(metadata_cols, colnames(pca_data))
    if (length(available_metadata) > 0) {
      message("Adding metadata columns: ", paste(available_metadata, collapse = ", "))
      metadata <- pca_data %>% select(all_of(c("seq_id_normalized", available_metadata)))
      ggplot_phylomorphospace_object <- ggplot_phylomorphospace_object %>%
        left_join(metadata, by = "seq_id_normalized")
      
      # Create clade column from first metadata column if clade doesn't exist
      if (!"clade" %in% colnames(ggplot_phylomorphospace_object)) {
        ggplot_phylomorphospace_object$clade <- ggplot_phylomorphospace_object[[available_metadata[1]]]
      }
    } else {
      message("Warning: None of the specified metadata columns found in data")
    }
  } else {
    # Auto-detect common metadata columns
    common_metadata <- c("family", "genus", "species", "clade", "group", "taxon")
    available_metadata <- intersect(common_metadata, colnames(pca_data))
    if (length(available_metadata) > 0) {
      message("Auto-detected metadata columns: ", paste(available_metadata, collapse = ", "))
      metadata <- pca_data %>% select(all_of(c("seq_id_normalized", available_metadata)))
      ggplot_phylomorphospace_object <- ggplot_phylomorphospace_object %>%
        left_join(metadata, by = "seq_id_normalized")
      
      # Create clade column from first metadata column if clade doesn't exist
      if (!"clade" %in% colnames(ggplot_phylomorphospace_object)) {
        ggplot_phylomorphospace_object$clade <- ggplot_phylomorphospace_object[[available_metadata[1]]]
      }
    }
  }
  
  # Ensure we have clade, family, and genus columns (fill with Unknown if missing)
  if (!"clade" %in% colnames(ggplot_phylomorphospace_object)) {
    ggplot_phylomorphospace_object$clade <- "Unknown"
  }
  if (!"family" %in% colnames(ggplot_phylomorphospace_object)) {
    ggplot_phylomorphospace_object$family <- "Unknown"
  }
  if (!"genus" %in% colnames(ggplot_phylomorphospace_object)) {
    ggplot_phylomorphospace_object$genus <- "Unknown"
  }
  
  # Prepare data for plotting
  df2 <- ggplot_phylomorphospace_object
  segments <- data.frame(
    x = c(rbind(df2$xstart, df2$xstop, NA)),
    y = c(rbind(df2$ystart, df2$ystop, NA)),
    z = c(rbind(df2$zstart, df2$zstop, NA))
  )
  tips <- subset(df2, nodeType == "Tip")
  
  # Report how many tips we have
  message("Tips in phylomorphospace object: ", nrow(tips))
  message("Expected tips (from tree): ", length(tree$tip.label))
  
  # Classify families and species groups for tips (primate groups)
  message("\nClassifying families and species groups...")
  tips$family <- sapply(tips$seq_id, classify_family)
  tips$species_group <- sapply(tips$seq_id, classify_species_group)
  
  # Colors by species group using Okabe-Ito palette
  tips$color <- group_colors[tips$species_group]
  tips$color[is.na(tips$color)] <- "gray50"
  
  # Shapes: derive plotly symbols from family-specific pch
  tips$plotly_symbol <- sapply(tips$family, function(fam) pch_to_plotly(family_shapes_pch[fam]))
  tips$plotly_symbol[is.na(tips$plotly_symbol)] <- "circle"
  
  # Get unique families for legend (order matches the 2D legend)
  family_order <- c("Bovidae", "Camelidae", "Cervidae", "Cetacea", "Hippopotamidae", "Hominidae", "Moschidae", "Suidae")
  unique_families <- intersect(family_order, unique(tips$family))
  unique_families <- unique_families[!is.na(unique_families) & unique_families != "Unknown"]
  
  message("Family counts:")
  for (fam in unique_families) {
    count <- sum(tips$family == fam, na.rm = TRUE)
    message("  ", fam, ": ", count)
  }
  
  message("\nSpecies group counts:")
  group_counts <- table(tips$species_group)
  for (group in names(group_counts)) {
    message("  ", group, ": ", group_counts[group])
  }
  
  message("\nCreating 3D plot...")
  
  # Calculate fixed axis ranges from ALL data (not just visible)
  # This ensures axes don't change when toggling groups
  all_x <- c(df2$xstart, df2$xstop, tips$xstop)
  all_y <- c(df2$ystart, df2$ystop, tips$ystop)
  all_z <- c(df2$zstart, df2$zstop, tips$zstop)
  
  # Remove NA values
  all_x <- all_x[!is.na(all_x) & is.finite(all_x)]
  all_y <- all_y[!is.na(all_y) & is.finite(all_y)]
  all_z <- all_z[!is.na(all_z) & is.finite(all_z)]
  
  # Calculate ranges with padding
  x_range <- range(all_x)
  y_range <- range(all_y)
  z_range <- range(all_z)
  
  # Add 5% padding to each axis
  x_padding <- diff(x_range) * 0.05
  y_padding <- diff(y_range) * 0.05
  z_padding <- diff(z_range) * 0.05
  
  x_axis_range <- c(x_range[1] - x_padding, x_range[2] + x_padding)
  y_axis_range <- c(y_range[1] - y_padding, y_range[2] + y_padding)
  z_axis_range <- c(z_range[1] - z_padding, z_range[2] + z_padding)
  
  message("Fixed X-axis range: [", round(x_axis_range[1], 2), ", ", round(x_axis_range[2], 2), "]")
  message("Fixed Y-axis range: [", round(y_axis_range[1], 2), ", ", round(y_axis_range[2], 2), "]")
  message("Fixed Z-axis range: [", round(z_axis_range[1], 4), ", ", round(z_axis_range[2], 4), "]")
  
  # Organize tips by family for legend - each family is its own toggleable item
  # Create base plot
  p <- plot_ly()
  
  # Add branches - make them toggleable
  p <- p %>%
    add_trace(
      data = segments,
      x = ~x,
      y = ~y,
      z = ~z,
      type = "scatter3d",
      mode = "lines",
      inherit = FALSE,
      color = I("lightblue"),
      line = list(width = 2),
      showlegend = TRUE,
      name = "Tree Branches",
      legendgroup = "Branches"
    )
  
  # Count tips per family for reporting
  total_tips_plotted <- 0
  
  # Add tips by family (one trace per family, toggleable in legend)
  for (fam_name in unique_families) {
    fam_tips <- subset(tips, family == fam_name)
    if (nrow(fam_tips) == 0) next
    
    # Color for this group (already set in tips$color)
    fam_color <- unique(fam_tips$color)[1]
    
    shape_pch <- family_shapes_pch[fam_name]
    plotly_symbol <- pch_to_plotly(shape_pch)
    if (is.na(plotly_symbol)) plotly_symbol <- "circle"
    
    message("  Plotting ", nrow(fam_tips), " tips for family: ", fam_name, " (color: ", fam_color, ", symbol: ", plotly_symbol, ")")
    total_tips_plotted <- total_tips_plotted + nrow(fam_tips)
    
    marker_list <- list(
      size = 8,
      color = fam_color,
      line = list(width = 1.5, color = "black"),
      symbol = plotly_symbol
    )
    
    p <- p %>%
      add_trace(
        data = fam_tips,
        x = ~xstop,
        y = ~ystop,
        z = ~zstop,
        type = "scatter3d",
        mode = "markers",
        name = fam_name,
        marker = marker_list,
        legendgroup = fam_name,
        text = ~paste(
          "seq_id:", seq_id,
          "<br>Group:", species_group
        ),
        hoverinfo = "text"
      )
  }
  
  message("Total tips plotted: ", total_tips_plotted, " / ", nrow(tips))
  
  # Add layout with interactive legend and FIXED axis ranges
  # Fixed ranges prevent axes from changing when toggling groups
  p <- p %>%
    layout(
      scene = list(
        xaxis = list(
          title = pc1_label,
          range = x_axis_range,  # Fixed range - won't auto-adjust
          autorange = FALSE       # Disable auto-ranging
        ),
        yaxis = list(
          title = pc2_label,
          range = y_axis_range,  # Fixed range - won't auto-adjust
          autorange = FALSE       # Disable auto-ranging
        ),
        zaxis = list(
          title = "Time (Branch Length)",
          range = z_axis_range,  # Fixed range - won't auto-adjust
          autorange = FALSE       # Disable auto-ranging
        )
      ),
      title = paste("3D Phylomorphospace:", dataset_name),
      legend = list(
        title = list(text = "<b>Click to toggle visibility</b><br>Double-click to isolate"),
        itemsizing = "constant",
        itemclick = "toggleothers",  # Click toggles other items
        itemdoubleclick = "toggle"   # Double-click isolates this item
      )
    )
  
  # Convex hulls removed per user request
  
  # Save as HTML
  message("\nSaving plot to: ", output_html)
  # Create output directory if needed
  output_dir <- dirname(output_html)
  if (output_dir != "." && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  # Try self-contained first, fall back to non-self-contained if pandoc not available
  tryCatch({
    htmlwidgets::saveWidget(p, output_html, selfcontained = TRUE)
  }, error = function(e) {
    message("Warning: Could not save as self-contained HTML (pandoc may be missing). Saving as non-self-contained...")
    htmlwidgets::saveWidget(p, output_html, selfcontained = FALSE)
  })
  message("Done! Plot saved to: ", output_html)
  
  return(p)
}

# Main execution
main <- function() {
  # Parse command-line arguments
  args <- parse_args()
  
  # Create the 3D phylomorphospace plot
  create_3d_phylomorphospace(
    pca_scores_file = args$pca_file,
    tree_file = args$tree_file,
    output_html = args$output_file,
    dataset_name = args$dataset_name,
    id_col = args$id_col,
    metadata_cols = args$metadata_cols
  )
  
  message("\n", paste(rep("=", 80), collapse = ""))
  message("All done!")
  message(paste(rep("=", 80), collapse = ""))
}

# Run if executed as script
if (!interactive()) {
  main()
}

