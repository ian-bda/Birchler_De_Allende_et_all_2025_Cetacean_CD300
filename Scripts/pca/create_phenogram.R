#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(dplyr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("Usage: Rscript create_phenogram.R <pca_scores_csv> <tree_file> <output_prefix> [PC_components]\n")
  cat("\nArguments:\n")
  cat("  pca_scores_csv  : Path to CSV file with PCA scores\n")
  cat("  tree_file       : Path to ultrametric Newick tree file\n")
  cat("  output_prefix   : Prefix for output files\n")
  cat("  PC_components   : Comma-separated list of PCs to plot (default: PC1,PC2)\n")
  cat("\nThis script creates phenograms showing trait evolution along the phylogeny.\n")
  quit(status = 1)
}

# Get arguments
pca_scores_file <- args[[1]]
tree_file <- args[[2]]
out_prefix <- args[[3]]
pc_components <- ifelse(length(args) >= 4, args[[4]], "PC1,PC2")

# Parse PC components
pc_list <- unlist(strsplit(pc_components, ","))
pc_list <- trimws(pc_list)

# Check if files exist
if (!file.exists(pca_scores_file)) {
  stop("Error: PCA scores file not found: ", pca_scores_file)
}
if (!file.exists(tree_file)) {
  stop("Error: Tree file not found: ", tree_file)
}

message(paste(rep("=", 60), collapse = ""))
message("Creating Phenograms for PCA Components")
message(paste(rep("=", 60), collapse = ""))
message("PCA scores file: ", pca_scores_file)
message("Tree file: ", tree_file)
message("Output prefix: ", out_prefix)
message("PC components: ", paste(pc_list, collapse = ", "))
message(paste(rep("=", 60), collapse = ""))

# 1) Load PCA scores
message("\n1. Loading PCA scores...")
pca_scores <- read.csv(pca_scores_file, stringsAsFactors = FALSE)

# Determine ID column
if ("original_sequence_id" %in% colnames(pca_scores)) {
  pca_scores$seq_id <- pca_scores$original_sequence_id
  message("Using 'original_sequence_id' column")
} else if ("clean_id" %in% colnames(pca_scores)) {
  pca_scores$seq_id <- pca_scores$clean_id
  message("Using 'clean_id' column")
} else {
  pca_scores$seq_id <- pca_scores[, 1]
  message("Using first column as sequence ID")
}

# Check if required PCs exist
missing_pcs <- setdiff(pc_list, colnames(pca_scores))
if (length(missing_pcs) > 0) {
  stop("PC components not found: ", paste(missing_pcs, collapse = ", "))
}

# 2) Load phylogenetic tree
message("\n2. Loading phylogenetic tree...")
tree <- read.tree(tree_file)
tree$tip.label <- gsub(" ", "_", tree$tip.label)
message("Tree loaded: ", length(tree$tip.label), " tips")

# Check if tree is ultrametric
if (!is.ultrametric(tree, tol = 1e-6)) {
  warning("WARNING: Tree may not be ultrametric. Phenogram may be inaccurate.")
}

# 3) Match sequences
common_seqs <- intersect(tree$tip.label, pca_scores$seq_id)
message("\n3. Matching sequences...")
message("  Sequences in tree: ", length(tree$tip.label))
message("  Sequences in PCA: ", nrow(pca_scores))
message("  Common sequences: ", length(common_seqs))

if (length(common_seqs) < 10) {
  stop("Too few common sequences (", length(common_seqs), "). Need at least 10.")
}

# Filter and prepare data
pca_filtered <- pca_scores %>%
  filter(seq_id %in% common_seqs) %>%
  arrange(match(seq_id, tree$tip.label))

# Extract PC scores matrix
pc_matrix <- pca_filtered %>%
  select(all_of(pc_list)) %>%
  as.matrix()
rownames(pc_matrix) <- pca_filtered$seq_id

# Prune tree to common sequences
tree_pruned <- keep.tip(tree, common_seqs)

# Reorder matrix to match tree
pc_matrix <- pc_matrix[tree_pruned$tip.label, , drop = FALSE]

# Handle zero-length branches
if (any(tree_pruned$edge.length == 0)) {
  message("  Warning: Tree has ", sum(tree_pruned$edge.length == 0), " zero-length branches")
  message("  Adding small value to zero-length branches...")
  min_edge <- min(tree_pruned$edge.length[tree_pruned$edge.length > 0])
  tree_pruned$edge.length[tree_pruned$edge.length == 0] <- min_edge * 0.001
}

message("\n4. Final dataset:")
message("  Taxa: ", nrow(pc_matrix))
message("  PC dimensions: ", ncol(pc_matrix))

# 5) Create phenograms for each PC
message("\n", paste(rep("=", 60), collapse = ""))
message("Creating Phenograms")
message(paste(rep("=", 60), collapse = ""))

# Create output PDF with larger dimensions for better readability
pdf(paste0(out_prefix, "_phenograms.pdf"), width = 16, height = 10 * length(pc_list))

par(mfrow = c(length(pc_list), 1), mar = c(5, 4, 4, 2) + 0.1)

for (pc_idx in 1:length(pc_list)) {
  pc_name <- pc_list[pc_idx]
  trait <- pc_matrix[, pc_name]
  names(trait) <- rownames(pc_matrix)
  
  message("\nCreating phenogram for ", pc_name, "...")
  
  # Create phenogram
  # phenogram plots tree with trait values on y-axis and time on x-axis
  tryCatch({
    phenogram(tree_pruned, trait,
              xlab = "Time (from root)",
              ylab = paste(pc_name, "value"),
              main = paste("Phenogram:", pc_name),
              spread.labels = TRUE,
              spread.cost = c(1, 0.5),  # Better label spreading (higher cost = more spreading)
              fsize = 0.35,  # Smaller font to reduce overlap
              colors = "steelblue",
              offset = 0.8)  # Offset labels further from tips
    
    message("  ✓ Phenogram created for ", pc_name)
  }, error = function(e) {
    message("  ✗ Error creating phenogram for ", pc_name, ": ", e$message)
    # Create a simple plot as fallback
    plot(1, 1, type = "n", 
         xlab = "Time", ylab = pc_name,
         main = paste("Phenogram:", pc_name, "(error occurred)"))
    text(0.5, 0.5, paste("Error:", e$message), col = "red")
  })
}

dev.off()
message("\nSaved: ", paste0(out_prefix, "_phenograms.pdf"))

# 6) Create individual phenograms with better formatting
message("\nCreating individual high-quality phenograms...")

for (pc_idx in 1:length(pc_list)) {
  pc_name <- pc_list[pc_idx]
  trait <- pc_matrix[, pc_name]
  names(trait) <- rownames(pc_matrix)
  
  # Use larger plot dimensions for better label readability
  pdf(paste0(out_prefix, "_phenogram_", pc_name, ".pdf"), width = 18, height = 12)
  
  tryCatch({
    # Create phenogram with better formatting
    # Use smaller font and better label spreading to avoid overlap
    phenogram(tree_pruned, trait,
              xlab = "Time (from root)",
              ylab = paste(pc_name, "value"),
              main = paste("Phenogram:", pc_name, "\nTrait Evolution Along Phylogeny"),
              spread.labels = TRUE,
              spread.cost = c(1, 0.5),  # Better label spreading (higher cost = more spreading)
              fsize = 0.35,  # Smaller font size to reduce overlap
              colors = "steelblue",
              lwd = 1.5,
              offset = 0.8)  # Offset labels further from tips
    
    # Add summary statistics
    trait_mean <- mean(trait)
    trait_sd <- sd(trait)
    trait_range <- range(trait)
    
    legend("topright",
           legend = c(
             paste("Mean:", round(trait_mean, 3)),
             paste("SD:", round(trait_sd, 3)),
             paste("Range: [", round(trait_range[1], 3), ", ", round(trait_range[2], 3), "]", sep = "")
           ),
           bty = "n",
           cex = 0.9)
    
    message("  ✓ Individual phenogram created for ", pc_name)
  }, error = function(e) {
    message("  ✗ Error creating individual phenogram for ", pc_name, ": ", e$message)
    plot(1, 1, type = "n", 
         xlab = "Time", ylab = pc_name,
         main = paste("Phenogram:", pc_name, "(error occurred)"))
    text(0.5, 0.5, paste("Error:", e$message), col = "red")
  })
  
  dev.off()
  message("  Saved: ", paste0(out_prefix, "_phenogram_", pc_name, ".pdf"))
}

# 7) Create a combined phenogram showing both PC1 and PC2
if (length(pc_list) >= 2 && "PC1" %in% pc_list && "PC2" %in% pc_list) {
  message("\nCreating combined phenogram (PC1 and PC2)...")
  
  # Use larger plot dimensions for better label readability
  pdf(paste0(out_prefix, "_phenogram_combined_PC1_PC2.pdf"), width = 20, height = 14)
  
  par(mfrow = c(2, 1), mar = c(5, 4, 3, 2) + 0.1)
  
  # PC1
  trait_pc1 <- pc_matrix[, "PC1"]
  names(trait_pc1) <- rownames(pc_matrix)
  
  tryCatch({
    phenogram(tree_pruned, trait_pc1,
              xlab = "Time (from root)",
              ylab = "PC1 value",
              main = "Phenogram: PC1",
              spread.labels = TRUE,
              spread.cost = c(1, 0.5),  # Better label spreading
              fsize = 0.35,  # Smaller font to reduce overlap
              colors = "steelblue",
              lwd = 1.5,
              offset = 0.8)  # Offset labels further from tips
  }, error = function(e) {
    plot(1, 1, type = "n", main = "PC1 (error occurred)")
    text(0.5, 0.5, paste("Error:", e$message), col = "red")
  })
  
  # PC2
  trait_pc2 <- pc_matrix[, "PC2"]
  names(trait_pc2) <- rownames(pc_matrix)
  
  tryCatch({
    phenogram(tree_pruned, trait_pc2,
              xlab = "Time (from root)",
              ylab = "PC2 value",
              main = "Phenogram: PC2",
              spread.labels = TRUE,
              spread.cost = c(1, 0.5),  # Better label spreading
              fsize = 0.35,  # Smaller font to reduce overlap
              colors = "darkorange",
              lwd = 1.5,
              offset = 0.8)  # Offset labels further from tips
  }, error = function(e) {
    plot(1, 1, type = "n", main = "PC2 (error occurred)")
    text(0.5, 0.5, paste("Error:", e$message), col = "red")
  })
  
  dev.off()
  message("Saved: ", paste0(out_prefix, "_phenogram_combined_PC1_PC2.pdf"))
}

message("\n", paste(rep("=", 60), collapse = ""))
message("Analysis Complete!")
message(paste(rep("=", 60), collapse = ""))
message("\nPhenograms show trait evolution along the phylogeny:")
message("  - X-axis: Time from root (using ultrametric tree)")
message("  - Y-axis: Trait value (PC score)")
message("  - Each branch shows how the trait changes over time")
message("  - Tips show current trait values for each species")

