#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(dplyr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("Usage: Rscript phylomorphospace_from_pca.R <pca_scores_csv> <tree_file> <output_prefix> [PC_x] [PC_y]\n")
  cat("\nArguments:\n")
  cat("  pca_scores_csv  : Path to CSV file with PCA scores (from structural_pca.py)\n")
  cat("  tree_file       : Path to Newick tree file (from IQ-TREE or similar)\n")
  cat("  output_prefix   : Prefix for output files (e.g., 'phylomorphospace')\n")
  cat("  PC_x            : Principal component for X-axis (default: PC1)\n")
  cat("  PC_y            : Principal component for Y-axis (default: PC2)\n")
  cat("\nExamples:\n")
  cat("  # Basic usage with default PC1 vs PC2\n")
  cat("  Rscript phylomorphospace_from_pca.R pca_scores.csv tree.nwk output\n")
  cat("\n  # Specify different PCs\n")
  cat("  Rscript phylomorphospace_from_pca.R pca_scores.csv tree.nwk output PC1 PC3\n")
  cat("\n  # Use full paths\n")
  cat("  Rscript phylomorphospace_from_pca.R /path/to/scores.csv /path/to/tree.nwk /path/to/output\n")
  quit(status = 1)
}

# Get arguments
pca_scores_file <- args[[1]]
tree_file <- args[[2]]
out_prefix <- args[[3]]
pc_x <- ifelse(length(args) >= 4, args[[4]], "PC1")
pc_y <- ifelse(length(args) >= 5, args[[5]], "PC2")

# Check if files exist and normalize paths
if (!file.exists(pca_scores_file)) {
  stop("Error: PCA scores file not found: ", pca_scores_file)
}
if (!file.exists(tree_file)) {
  stop("Error: Tree file not found: ", tree_file)
}

pca_scores_file <- normalizePath(pca_scores_file, mustWork = TRUE)
tree_file <- normalizePath(tree_file, mustWork = TRUE)

# Create output directory if needed
output_dir <- dirname(out_prefix)
if (output_dir != "." && !dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("Created output directory: ", output_dir)
}

message(paste(rep("=", 60), collapse = ""))
message("Phylomorphospace Analysis")
message(paste(rep("=", 60), collapse = ""))
message("PCA scores file: ", pca_scores_file)
message("Tree file: ", tree_file)
message("Output prefix: ", out_prefix)
message("X-axis (PC): ", pc_x)
message("Y-axis (PC): ", pc_y)
message(paste(rep("=", 60), collapse = ""))

# 1) Load PCA scores
message("\nLoading PCA scores...")
pca_scores <- read.csv(pca_scores_file, stringsAsFactors = FALSE)

# Check if required PCs exist
if (!pc_x %in% colnames(pca_scores)) {
  stop("PC '", pc_x, "' not found in PCA scores. Available: ", 
       paste(colnames(pca_scores)[grep("^PC", colnames(pca_scores))], collapse = ", "))
}
if (!pc_y %in% colnames(pca_scores)) {
  stop("PC '", pc_y, "' not found in PCA scores. Available: ", 
       paste(colnames(pca_scores)[grep("^PC", colnames(pca_scores))], collapse = ", "))
}

# Determine which column to use for matching with tree
# IQ-TREE tree uses original_sequence_id format (with pipes)
if ("original_sequence_id" %in% colnames(pca_scores)) {
  pca_scores$seq_id <- pca_scores$original_sequence_id
  message("Using 'original_sequence_id' column to match with tree")
} else if ("clean_id" %in% colnames(pca_scores)) {
  pca_scores$seq_id <- pca_scores$clean_id
  message("Using 'clean_id' column to match with tree")
  message("  WARNING: Tree may use original_sequence_id format - check for mismatches!")
} else {
  # Use first column as fallback
  pca_scores$seq_id <- pca_scores[, 1]
  message("Using first column '", colnames(pca_scores)[1], "' to match with tree")
}

message("Total sequences in PCA data: ", nrow(pca_scores))
message("  Sample sequence IDs: ", paste(head(pca_scores$seq_id, 3), collapse = ", "))

# 2) Load phylogenetic tree
message("\nLoading phylogenetic tree...")
tree <- read.tree(tree_file)
tree$tip.label <- gsub(" ", "_", tree$tip.label)
message("Tree loaded: ", length(tree$tip.label), " tips")

# 3) Match sequences between tree and PCA data
common_seqs <- intersect(tree$tip.label, pca_scores$seq_id)
message("\nMatching sequences...")
message("  Sequences in tree: ", length(tree$tip.label))
message("  Sequences in PCA: ", nrow(pca_scores))
message("  Common sequences: ", length(common_seqs))

if (length(common_seqs) < 3) {
  # Show some examples of mismatches for debugging
  tree_only <- setdiff(tree$tip.label, pca_scores$seq_id)
  pca_only <- setdiff(pca_scores$seq_id, tree$tip.label)
  if (length(tree_only) > 0) {
    message("\n  Sample sequences in tree but not in PCA:")
    message("    ", paste(head(tree_only, 3), collapse = ", "))
  }
  if (length(pca_only) > 0) {
    message("\n  Sample sequences in PCA but not in tree:")
    message("    ", paste(head(pca_only, 3), collapse = ", "))
  }
  stop("Too few common sequences (", length(common_seqs), ") between tree and PCA data. Need at least 3.")
}

# Show mismatch info if significant
if (length(common_seqs) < min(length(tree$tip.label), nrow(pca_scores)) * 0.9) {
  tree_only <- setdiff(tree$tip.label, pca_scores$seq_id)
  pca_only <- setdiff(pca_scores$seq_id, tree$tip.label)
  if (length(tree_only) > 0) {
    message("\n  Warning: ", length(tree_only), " sequences in tree but not in PCA")
    message("    Sample: ", paste(head(tree_only, 2), collapse = ", "))
  }
  if (length(pca_only) > 0) {
    message("  Warning: ", length(pca_only), " sequences in PCA but not in tree")
    message("    Sample: ", paste(head(pca_only, 2), collapse = ", "))
  }
}

# Filter and reorder data
pca_filtered <- pca_scores %>%
  filter(seq_id %in% common_seqs) %>%
  arrange(match(seq_id, tree$tip.label))

# Create matrix for phylomorphospace (in tree tip order)
pc_matrix <- pca_filtered %>%
  select(all_of(c(pc_x, pc_y))) %>%
  as.matrix()
rownames(pc_matrix) <- pca_filtered$seq_id

# Check for invalid values (NaN, Inf, NA)
valid_rows <- complete.cases(pc_matrix) & 
               is.finite(pc_matrix[, 1]) & 
               is.finite(pc_matrix[, 2])
pc_matrix <- pc_matrix[valid_rows, , drop = FALSE]
valid_seqs <- rownames(pc_matrix)

if (nrow(pc_matrix) < 3) {
  stop("Too few sequences with valid PC values (", nrow(pc_matrix), "). Need at least 3.")
}

if (nrow(pc_matrix) < length(common_seqs)) {
  message("  Warning: ", length(common_seqs) - nrow(pc_matrix), 
          " sequences removed due to invalid PC values")
}

# Prune tree to sequences with valid PC values
tree_pruned <- keep.tip(tree, valid_seqs)

# Reorder matrix to match tree tip order exactly
pc_matrix <- pc_matrix[tree_pruned$tip.label, , drop = FALSE]

# Final check for valid ranges
if (any(!is.finite(pc_matrix))) {
  stop("Invalid PC values detected after filtering. Check PCA scores.")
}

message("\nFinal dataset:")
message("  Taxa: ", nrow(pc_matrix))
message("  Dimensions: ", ncol(pc_matrix))
message("  PC", substr(pc_x, 3, 10), " range: [", round(min(pc_matrix[, 1]), 3), ", ", round(max(pc_matrix[, 1]), 3), "]")
message("  PC", substr(pc_y, 3, 10), " range: [", round(min(pc_matrix[, 2]), 3), ", ", round(max(pc_matrix[, 2]), 3), "]")

# 4) Create phylomorphospace plot
message("\nCreating phylomorphospace plot...")

# Get variance explained if available
pca_all <- read.csv(pca_scores_file, stringsAsFactors = FALSE)
variance_file <- sub("_pca_scores\\.csv$", "_explained_variance.csv", pca_scores_file)
if (file.exists(variance_file)) {
  variance_df <- read.csv(variance_file)
  pc_x_idx <- which(variance_df$PC == pc_x)
  pc_y_idx <- which(variance_df$PC == pc_y)
  var_x <- if (length(pc_x_idx) > 0) variance_df$Explained_Variance_Ratio[pc_x_idx] * 100 else NA
  var_y <- if (length(pc_y_idx) > 0) variance_df$Explained_Variance_Ratio[pc_y_idx] * 100 else NA
} else {
  var_x <- NA
  var_y <- NA
}

# Create plot
pdf(paste0(out_prefix, "_phylomorphospace_", pc_x, "_", pc_y, ".pdf"), 
    width = 10, height = 10)
par(mar = c(5, 5, 2, 2))

xlab <- if (!is.na(var_x)) {
  paste0(pc_x, " (", round(var_x, 1), "% variance)")
} else {
  pc_x
}

ylab <- if (!is.na(var_y)) {
  paste0(pc_y, " (", round(var_y, 1), "% variance)")
} else {
  pc_y
}

# Get tree info
n_tips <- length(tree_pruned$tip.label)
n_nodes <- tree_pruned$Nnode

# Handle zero-length branches (fastAnc can have issues with these)
# Add small value to zero-length branches to avoid NaN in ancestral reconstruction
if (any(tree_pruned$edge.length == 0)) {
  message("  Warning: Tree has ", sum(tree_pruned$edge.length == 0), " zero-length branches")
  message("  Adding small value to zero-length branches for ancestral state calculation...")
  tree_pruned$edge.length[tree_pruned$edge.length == 0] <- min(tree_pruned$edge.length[tree_pruned$edge.length > 0]) * 0.001
}

# Calculate ancestral states for internal nodes
# Method: Maximum Likelihood (ML) with Brownian Motion model
# Note: Maximum Parsimony (MP) is designed for discrete characters, not continuous traits.
#       For continuous PCA scores, ML with Brownian Motion is the standard approach.
#       This assumes traits evolve under a random walk (Brownian motion) model.
message("  Calculating ancestral states for internal nodes...")
message("  Method: Maximum Likelihood with Brownian Motion model")
message("  (MP is for discrete characters; ML/BM is standard for continuous traits)")
anc_x <- fastAnc(tree_pruned, pc_matrix[, 1])
anc_y <- fastAnc(tree_pruned, pc_matrix[, 2])

# Check for NaN values and replace with interpolated values if needed
if (any(is.nan(anc_x)) || any(is.nan(anc_y))) {
  message("  Warning: Some ancestral states are NaN, attempting to fix...")
  # For NaN values, use the mean of child nodes
  for (i in 1:length(anc_x)) {
    if (is.nan(anc_x[i])) {
      node_id <- as.integer(names(anc_x)[i])
      # Find children of this node
      children <- tree_pruned$edge[tree_pruned$edge[, 1] == node_id, 2]
      if (length(children) > 0) {
        # Get trait values of children and average
        child_vals <- sapply(children, function(cid) {
          if (cid <= n_tips) {
            tip_name <- tree_pruned$tip.label[cid]
            return(pc_matrix[tip_name, 1])
          } else {
            cid_str <- as.character(cid)
            if (cid_str %in% names(anc_x) && !is.nan(anc_x[cid_str])) {
              return(anc_x[cid_str])
            }
          }
          return(NA)
        })
        child_vals <- child_vals[!is.na(child_vals)]
        if (length(child_vals) > 0) {
          anc_x[i] <- mean(child_vals)
        }
      }
    }
    if (is.nan(anc_y[i])) {
      node_id <- as.integer(names(anc_y)[i])
      children <- tree_pruned$edge[tree_pruned$edge[, 1] == node_id, 2]
      if (length(children) > 0) {
        child_vals <- sapply(children, function(cid) {
          if (cid <= n_tips) {
            tip_name <- tree_pruned$tip.label[cid]
            return(pc_matrix[tip_name, 2])
          } else {
            cid_str <- as.character(cid)
            if (cid_str %in% names(anc_y) && !is.nan(anc_y[cid_str])) {
              return(anc_y[cid_str])
            }
          }
          return(NA)
        })
        child_vals <- child_vals[!is.na(child_vals)]
        if (length(child_vals) > 0) {
          anc_y[i] <- mean(child_vals)
        }
      }
    }
  }
}

# Create a mapping of all node IDs to their trait values

# Create a complete mapping: node_id -> c(x, y)
node_traits <- list()

# Add tips
for (i in 1:n_tips) {
  tip_name <- tree_pruned$tip.label[i]
  node_traits[[as.character(i)]] <- c(pc_matrix[tip_name, 1], pc_matrix[tip_name, 2])
}

# Add internal nodes - fastAnc returns named vector with node IDs as names
for (node_id_str in names(anc_x)) {
  node_id <- as.integer(node_id_str)
  if (!is.na(node_id) && node_id > n_tips) {
    node_traits[[node_id_str]] <- c(anc_x[node_id_str], anc_y[node_id_str])
  }
}

# Function to get trait values for any node
get_node_trait <- function(node_id) {
  node_id_str <- as.character(node_id)
  if (node_id_str %in% names(node_traits)) {
    return(node_traits[[node_id_str]])
  } else {
    # Fallback: if it's a tip, try direct lookup
    if (node_id <= n_tips) {
      tip_name <- tree_pruned$tip.label[node_id]
      if (tip_name %in% rownames(pc_matrix)) {
        return(c(pc_matrix[tip_name, 1], pc_matrix[tip_name, 2]))
      }
    }
    return(c(NA, NA))
  }
}

# Calculate axis limits including all nodes (filter invalid values)
all_x <- c(pc_matrix[, 1], anc_x)
all_y <- c(pc_matrix[, 2], anc_y)
all_x <- all_x[is.finite(all_x)]
all_y <- all_y[is.finite(all_y)]

if (length(all_x) == 0 || length(all_y) == 0) {
  stop("Error: No valid values for axis limits after ancestral state calculation")
}

x_range <- range(all_x)
y_range <- range(all_y)
x_padding <- diff(x_range) * 0.1
y_padding <- diff(y_range) * 0.1
xlim <- c(x_range[1] - x_padding, x_range[2] + x_padding)
ylim <- c(y_range[1] - y_padding, y_range[2] + y_padding)

# Create empty plot
plot(1, 1, type = "n", 
     xlim = xlim, ylim = ylim,
     xlab = xlab, ylab = ylab,
     main = paste("Phylomorphospace:", pc_x, "vs", pc_y))

# Draw tree branches (edges) in morphospace
edge <- tree_pruned$edge

# Debug: check node mapping
message("  Node mapping summary:")
message("    Tips mapped: ", sum(1:n_tips %in% as.integer(names(node_traits))), " / ", n_tips)
message("    Internal nodes mapped: ", length(node_traits) - n_tips, " / ", n_nodes)
message("    Total edges in tree: ", nrow(edge))
message("    Sample node IDs in tree edges (parents): ", paste(head(unique(edge[, 1]), 5), collapse = ", "))
message("    Sample node IDs in tree edges (children): ", paste(head(unique(edge[, 2]), 5), collapse = ", "))
message("    Node IDs in node_traits: ", paste(head(names(node_traits), 5), collapse = ", "))

edges_drawn <- 0
edges_skipped <- 0

for (i in 1:nrow(edge)) {
  parent <- edge[i, 1]
  child <- edge[i, 2]
  
  parent_trait <- get_node_trait(parent)
  child_trait <- get_node_trait(child)
  
  # Draw line connecting parent and child (only if both are valid)
  if (all(is.finite(parent_trait)) && all(is.finite(child_trait))) {
    lines(c(parent_trait[1], child_trait[1]), 
          c(parent_trait[2], child_trait[2]),
          col = "gray60", lwd = 0.8)
    edges_drawn <- edges_drawn + 1
  } else {
    edges_skipped <- edges_skipped + 1
  }
}

message("  Edges drawn: ", edges_drawn, " / ", nrow(edge))
if (edges_skipped > 0) {
  message("  Edges skipped (invalid traits): ", edges_skipped)
}

# Plot internal nodes (ancestors) as small gray points
if (n_nodes > 0 && length(anc_x) > 0 && all(is.finite(anc_x)) && all(is.finite(anc_y))) {
  points(anc_x, anc_y, 
         pch = 21, 
         bg = "lightgray", 
         cex = 0.6,
         col = "gray40",
         lwd = 0.5)
}

# Plot tips (species) as blue points
points(pc_matrix[, 1], pc_matrix[, 2], 
       pch = 21, 
       bg = "steelblue", 
       cex = 1.0,
       col = "darkblue",
       lwd = 1.0)

dev.off()

message("Saved phylomorphospace plot: ", paste0(out_prefix, "_phylomorphospace_", pc_x, "_", pc_y, ".pdf"))

# 5) Save matched data
write.csv(pc_matrix, 
          file = paste0(out_prefix, "_matched_pca_scores.csv"), 
          quote = FALSE)
message("Saved matched PCA scores: ", paste0(out_prefix, "_matched_pca_scores.csv"))

# 6) Create additional plots (PC1 vs PC2, PC2 vs PC3 if available)
if (pc_x == "PC1" && pc_y == "PC2" && "PC3" %in% colnames(pca_scores)) {
  message("\nCreating additional plots...")
  
  # PC2 vs PC3
  pca_filtered_23 <- pca_scores %>%
    filter(seq_id %in% valid_seqs) %>%
    arrange(match(seq_id, tree_pruned$tip.label))
  
  pc_matrix_23 <- pca_filtered_23 %>%
    select(PC2, PC3) %>%
    as.matrix()
  rownames(pc_matrix_23) <- pca_filtered_23$seq_id
  
  # Check for invalid values
  valid_rows_23 <- complete.cases(pc_matrix_23) & 
                   is.finite(pc_matrix_23[, 1]) & 
                   is.finite(pc_matrix_23[, 2])
  pc_matrix_23 <- pc_matrix_23[valid_rows_23, , drop = FALSE]
  valid_seqs_23 <- rownames(pc_matrix_23)
  tree_pruned_23 <- keep.tip(tree_pruned, valid_seqs_23)
  pc_matrix_23 <- pc_matrix_23[tree_pruned_23$tip.label, , drop = FALSE]
  
  var_2 <- if (file.exists(variance_file)) {
    variance_df <- read.csv(variance_file)
    idx <- which(variance_df$PC == "PC2")
    if (length(idx) > 0) variance_df$Explained_Variance_Ratio[idx] * 100 else NA
  } else NA
  
  var_3 <- if (file.exists(variance_file)) {
    variance_df <- read.csv(variance_file)
    idx <- which(variance_df$PC == "PC3")
    if (length(idx) > 0) variance_df$Explained_Variance_Ratio[idx] * 100 else NA
  } else NA
  
  pdf(paste0(out_prefix, "_phylomorphospace_PC2_PC3.pdf"), width = 10, height = 10)
  par(mar = c(5, 5, 2, 2))
  
  xlab_23 <- if (!is.na(var_2)) {
    paste0("PC2 (", round(var_2, 1), "% variance)")
  } else {
    "PC2"
  }
  
  ylab_23 <- if (!is.na(var_3)) {
    paste0("PC3 (", round(var_3, 1), "% variance)")
  } else {
    "PC3"
  }
  
  # Get tree info
  n_tips_23 <- length(tree_pruned_23$tip.label)
  n_nodes_23 <- tree_pruned_23$Nnode
  
  # Handle zero-length branches
  if (any(tree_pruned_23$edge.length == 0)) {
    tree_pruned_23$edge.length[tree_pruned_23$edge.length == 0] <- min(tree_pruned_23$edge.length[tree_pruned_23$edge.length > 0]) * 0.001
  }
  
  # Calculate ancestral states for PC2 vs PC3
  # Method: Maximum Likelihood (ML) with Brownian Motion model
  message("  Calculating ancestral states for PC2 vs PC3...")
  message("  Method: Maximum Likelihood with Brownian Motion model")
  anc_x_23 <- fastAnc(tree_pruned_23, pc_matrix_23[, 1])
  anc_y_23 <- fastAnc(tree_pruned_23, pc_matrix_23[, 2])
  
  # Check for NaN values and fix
  if (any(is.nan(anc_x_23)) || any(is.nan(anc_y_23))) {
    message("  Warning: Some ancestral states are NaN, attempting to fix...")
    for (i in 1:length(anc_x_23)) {
      if (is.nan(anc_x_23[i])) {
        node_id <- as.integer(names(anc_x_23)[i])
        children <- tree_pruned_23$edge[tree_pruned_23$edge[, 1] == node_id, 2]
        if (length(children) > 0) {
          child_vals <- sapply(children, function(cid) {
            if (cid <= n_tips_23) {
              tip_name <- tree_pruned_23$tip.label[cid]
              return(pc_matrix_23[tip_name, 1])
            } else {
              cid_str <- as.character(cid)
              if (cid_str %in% names(anc_x_23) && !is.nan(anc_x_23[cid_str])) {
                return(anc_x_23[cid_str])
              }
            }
            return(NA)
          })
          child_vals <- child_vals[!is.na(child_vals)]
          if (length(child_vals) > 0) {
            anc_x_23[i] <- mean(child_vals)
          }
        }
      }
      if (is.nan(anc_y_23[i])) {
        node_id <- as.integer(names(anc_y_23)[i])
        children <- tree_pruned_23$edge[tree_pruned_23$edge[, 1] == node_id, 2]
        if (length(children) > 0) {
          child_vals <- sapply(children, function(cid) {
            if (cid <= n_tips_23) {
              tip_name <- tree_pruned_23$tip.label[cid]
              return(pc_matrix_23[tip_name, 2])
            } else {
              cid_str <- as.character(cid)
              if (cid_str %in% names(anc_y_23) && !is.nan(anc_y_23[cid_str])) {
                return(anc_y_23[cid_str])
              }
            }
            return(NA)
          })
          child_vals <- child_vals[!is.na(child_vals)]
          if (length(child_vals) > 0) {
            anc_y_23[i] <- mean(child_vals)
          }
        }
      }
    }
  }
  
  # Create mapping of all node IDs to their trait values
  
  node_traits_23 <- list()
  
  # Add tips
  for (i in 1:n_tips_23) {
    tip_name <- tree_pruned_23$tip.label[i]
    node_traits_23[[as.character(i)]] <- c(pc_matrix_23[tip_name, 1], pc_matrix_23[tip_name, 2])
  }
  
  # Add internal nodes
  for (node_id_str in names(anc_x_23)) {
    node_id <- as.integer(node_id_str)
    if (!is.na(node_id) && node_id > n_tips_23) {
      node_traits_23[[node_id_str]] <- c(anc_x_23[node_id_str], anc_y_23[node_id_str])
    }
  }
  
  # Function to get trait values for any node
  get_node_trait_23 <- function(node_id) {
    node_id_str <- as.character(node_id)
    if (node_id_str %in% names(node_traits_23)) {
      return(node_traits_23[[node_id_str]])
    } else {
      if (node_id <= n_tips_23) {
        tip_name <- tree_pruned_23$tip.label[node_id]
        if (tip_name %in% rownames(pc_matrix_23)) {
          return(c(pc_matrix_23[tip_name, 1], pc_matrix_23[tip_name, 2]))
        }
      }
      return(c(NA, NA))
    }
  }
  
  # Calculate axis limits (filter invalid values)
  all_x_23 <- c(pc_matrix_23[, 1], anc_x_23)
  all_y_23 <- c(pc_matrix_23[, 2], anc_y_23)
  all_x_23 <- all_x_23[is.finite(all_x_23)]
  all_y_23 <- all_y_23[is.finite(all_y_23)]
  
  if (length(all_x_23) == 0 || length(all_y_23) == 0) {
    message("  Warning: No valid values for PC2 vs PC3 plot, skipping...")
    dev.off()
  } else {
  
  x_range_23 <- range(all_x_23)
  y_range_23 <- range(all_y_23)
  x_padding_23 <- diff(x_range_23) * 0.1
  y_padding_23 <- diff(y_range_23) * 0.1
  xlim_23 <- c(x_range_23[1] - x_padding_23, x_range_23[2] + x_padding_23)
  ylim_23 <- c(y_range_23[1] - y_padding_23, y_range_23[2] + y_padding_23)
  
  # Create plot
  plot(1, 1, type = "n", 
       xlim = xlim_23, ylim = ylim_23,
       xlab = xlab_23, ylab = ylab_23,
       main = "Phylomorphospace: PC2 vs PC3")
  
  # Draw tree branches
  edge_23 <- tree_pruned_23$edge
  edges_drawn_23 <- 0
  edges_skipped_23 <- 0
  
  for (i in 1:nrow(edge_23)) {
    parent <- edge_23[i, 1]
    child <- edge_23[i, 2]
    
    parent_trait <- get_node_trait_23(parent)
    child_trait <- get_node_trait_23(child)
    
    if (all(is.finite(parent_trait)) && all(is.finite(child_trait))) {
      lines(c(parent_trait[1], child_trait[1]), 
            c(parent_trait[2], child_trait[2]),
            col = "gray60", lwd = 0.8)
      edges_drawn_23 <- edges_drawn_23 + 1
    } else {
      edges_skipped_23 <- edges_skipped_23 + 1
    }
  }
  
  message("  Edges drawn (PC2 vs PC3): ", edges_drawn_23, " / ", nrow(edge_23))
  if (edges_skipped_23 > 0) {
    message("  Edges skipped (invalid traits): ", edges_skipped_23)
  }
  
  # Plot internal nodes (ancestors) as small gray points
  if (n_nodes_23 > 0 && length(anc_x_23) > 0 && all(is.finite(anc_x_23)) && all(is.finite(anc_y_23))) {
    points(anc_x_23, anc_y_23, 
           pch = 21, 
           bg = "lightgray", 
           cex = 0.6,
           col = "gray40",
           lwd = 0.5)
  }
  
  # Plot tips (species) as blue points
  points(pc_matrix_23[, 1], pc_matrix_23[, 2], 
         pch = 21, 
         bg = "steelblue", 
         cex = 1.0,
         col = "darkblue",
         lwd = 1.0)
  
  dev.off()
  
  message("Saved PC2 vs PC3 plot: ", paste0(out_prefix, "_phylomorphospace_PC2_PC3.pdf"))
  }
}

message("\nDone! Output files:")
message("  - ", paste0(out_prefix, "_phylomorphospace_", pc_x, "_", pc_y, ".pdf"))
message("  - ", paste0(out_prefix, "_matched_pca_scores.csv"))



