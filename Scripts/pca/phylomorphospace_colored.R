#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(dplyr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("Usage: Rscript phylomorphospace_colored.R <pca_scores_csv> <tree_file> <output_prefix> [PC_x] [PC_y]\n")
  cat("\nArguments:\n")
  cat("  pca_scores_csv  : Path to CSV file with PCA scores (from structural_pca.py)\n")
  cat("  tree_file       : Path to Newick tree file (from IQ-TREE or similar)\n")
  cat("  output_prefix   : Prefix for output files (e.g., 'phylomorphospace')\n")
  cat("  PC_x            : Principal component for X-axis (default: PC1)\n")
  cat("  PC_y            : Principal component for Y-axis (default: PC2)\n")
  cat("\nThis script colors tips by species group:\n")
  cat("  - Cetacean: Blue\n")
  cat("  - Terrestrial Artiodactyl: Orange/Brown\n")
  cat("  - Other (Human): Green\n")
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
message("Phylomorphospace Analysis (Colored by Species Group)")
message(paste(rep("=", 60), collapse = ""))
message("PCA scores file: ", pca_scores_file)
message("Tree file: ", tree_file)
message("Output prefix: ", out_prefix)
message("X-axis (PC): ", pc_x)
message("Y-axis (PC): ", pc_y)
message(paste(rep("=", 60), collapse = ""))

# Function to classify family from sequence ID
classify_family <- function(seq_id) {
  # Extract genus (first part before underscore or pipe)
  # Try underscore first, then pipe separator
  if (grepl("_", seq_id)) {
    genus <- strsplit(seq_id, "_")[[1]][1]
  } else if (grepl("\\|", seq_id)) {
    genus <- strsplit(seq_id, "\\|")[[1]][1]
  } else {
    genus <- seq_id  # Fallback: use whole string
  }
  
  # Family classifications based on Artiodactyla taxonomy
  # Camelidae (Tylopoda)
  if (genus %in% c("Camelus", "Vicugna")) {
    return("Camelidae")
  }
  
  # Suidae (pigs)
  if (genus %in% c("Sus", "Phacochoerus")) {
    return("Suidae")
  }
  
  # Cervidae (deer)
  if (genus %in% c("Cervus", "Dama", "Muntiacus", "Odocoileus", "Rangifer")) {
    return("Cervidae")
  }
  
  # Bovidae (cattle, sheep, antelopes)
  if (genus %in% c("Bison", "Bos", "Bubalus", "Budorcas", "Capra", 
                    "Capricornis", "Oryx", "Ovibos", "Ovis", "Pantholops")) {
    return("Bovidae")
  }
  
  # Moschidae (musk deer)
  if (genus == "Moschus") {
    return("Moschidae")
  }
  
  # Hippopotamidae (hippos)
  if (genus == "Hippopotamus") {
    return("Hippopotamidae")
  }
  
  # Cetaceans (whales and dolphins)
  if (genus %in% c("Balaenoptera", "Delphinapterus", "Delphinus", "Eschrichtius", 
                   "Eubalaena", "Globicephala", "Kogia", "Lagenorhynchus", 
                   "Lipotes", "Monodon", "Neophocaena", "Orcinus", "Phocoena", 
                   "Physeter", "Pseudorca", "Tursiops", "Sousa", "Mesoplodon")) {
    return("Cetacea")
  }
  
  # Hominidae (humans)
  if (genus == "Homo") {
    return("Hominidae")
  }
  
  # Unknown/Other - classify as Hominidae (these are all hominidae)
  return("Hominidae")
}

# Function to classify species group from sequence ID (for colors)
classify_species_group <- function(seq_id) {
  # Extract genus (first part before underscore or pipe)
  # Try underscore first, then pipe separator
  if (grepl("_", seq_id)) {
    genus <- strsplit(seq_id, "_")[[1]][1]
  } else if (grepl("\\|", seq_id)) {
    genus <- strsplit(seq_id, "\\|")[[1]][1]
  } else {
    genus <- seq_id  # Fallback: use whole string
  }
  
  # Cetacean genera
  cetacean_genera <- c(
    "Balaenoptera", "Delphinapterus", "Delphinus", "Eschrichtius", 
    "Eubalaena", "Globicephala", "Kogia", "Lagenorhynchus", 
    "Lipotes", "Monodon", "Neophocaena", "Orcinus", "Phocoena", 
    "Physeter", "Pseudorca", "Tursiops", "Sousa", "Mesoplodon"
  )
  
  # Terrestrial Artiodactyl genera (including Homo)
  artiodactyl_genera <- c(
    "Bison", "Bos", "Bubalus", "Budorcas", "Camelus", "Capra", 
    "Capricornis", "Cervus", "Dama", "Hippopotamus", "Moschus", 
    "Muntiacus", "Odocoileus", "Oryx", "Ovibos", "Ovis", 
    "Pantholops", "Phacochoerus", "Rangifer", "Sus", "Vicugna", "Homo"
  )
  
  if (genus %in% cetacean_genera) {
    return("Cetacean")
  } else if (genus %in% artiodactyl_genera) {
    return("Terrestrial_Artiodactyl")
  } else {
    return("Other")
  }
}

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

# Function to normalize sequence IDs for matching
normalize_seq_id <- function(seq_id) {
  # Replace spaces with underscores
  normalized <- gsub(" ", "_", seq_id)
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

# Determine which column to use for matching with tree
if ("original_sequence_id" %in% colnames(pca_scores)) {
  pca_scores$seq_id <- pca_scores$original_sequence_id
  message("Using 'original_sequence_id' column to match with tree")
} else if ("clean_id" %in% colnames(pca_scores)) {
  pca_scores$seq_id <- pca_scores$clean_id
  message("Using 'clean_id' column to match with tree")
} else {
  pca_scores$seq_id <- pca_scores[, 1]
  message("Using first column '", colnames(pca_scores)[1], "' to match with tree")
}

# Normalize PCA sequence IDs for matching
pca_scores$seq_id_normalized <- normalize_seq_id(pca_scores$seq_id)

# Classify species groups (for colors)
message("\nClassifying species groups...")
pca_scores$species_group <- sapply(pca_scores$seq_id, classify_species_group)
group_counts <- table(pca_scores$species_group)
message("Species group counts:")
for (group in names(group_counts)) {
  message("  ", group, ": ", group_counts[group])
}

# Classify families (for shapes)
message("\nClassifying families...")
pca_scores$family <- sapply(pca_scores$seq_id, classify_family)
family_counts <- table(pca_scores$family)
message("Family counts:")
for (fam in names(family_counts)) {
  message("  ", fam, ": ", family_counts[fam])
}

message("Total sequences in PCA data: ", nrow(pca_scores))

# 2) Load phylogenetic tree
message("\nLoading phylogenetic tree...")
tree <- read.tree(tree_file)
tree$tip.label <- gsub(" ", "_", tree$tip.label)
# Normalize tree tip labels
tree$tip.label <- normalize_seq_id(tree$tip.label)
message("Tree loaded: ", length(tree$tip.label), " tips")

# 3) Match sequences between tree and PCA data (using normalized IDs)
common_seqs_normalized <- intersect(tree$tip.label, pca_scores$seq_id_normalized)
message("\nMatching sequences (with normalized IDs)...")
message("  Sequences in tree: ", length(tree$tip.label))
message("  Sequences in PCA: ", nrow(pca_scores))
message("  Common sequences (normalized): ", length(common_seqs_normalized))

# Create mapping from normalized IDs to original IDs
tree_normalized_to_original <- setNames(tree$tip.label, tree$tip.label)  # Tree IDs are already normalized
pca_normalized_to_original <- setNames(pca_scores$seq_id, pca_scores$seq_id_normalized)

# Get original IDs for matching
common_seqs <- common_seqs_normalized  # Use normalized IDs for matching

if (length(common_seqs) < 3) {
  stop("Too few common sequences (", length(common_seqs), ") between tree and PCA data. Need at least 3.")
}

# Filter and reorder data using normalized IDs
pca_filtered <- pca_scores %>%
  filter(seq_id_normalized %in% common_seqs) %>%
  arrange(match(seq_id_normalized, tree$tip.label))

# Create matrix for phylomorphospace (in tree tip order)
# Use normalized IDs for rownames to match tree tip labels
pc_matrix <- pca_filtered %>%
  select(all_of(c(pc_x, pc_y))) %>%
  as.matrix()
rownames(pc_matrix) <- pca_filtered$seq_id_normalized

# Get species groups for tips (in tree tip order) - for colors
tip_groups <- pca_filtered$species_group
names(tip_groups) <- pca_filtered$seq_id_normalized

# Get families for tips (in tree tip order) - for shapes
tip_families <- pca_filtered$family
names(tip_families) <- pca_filtered$seq_id_normalized

# Check for invalid values
valid_rows <- complete.cases(pc_matrix) & 
               is.finite(pc_matrix[, 1]) & 
               is.finite(pc_matrix[, 2])
pc_matrix <- pc_matrix[valid_rows, , drop = FALSE]
valid_seqs <- rownames(pc_matrix)
tip_groups <- tip_groups[valid_seqs]
tip_families <- tip_families[valid_seqs]

if (nrow(pc_matrix) < 3) {
  stop("Too few sequences with valid PC values (", nrow(pc_matrix), "). Need at least 3.")
}

# Prune tree to sequences with valid PC values
tree_pruned <- keep.tip(tree, valid_seqs)

# Reorder matrix to match tree tip order exactly
pc_matrix <- pc_matrix[tree_pruned$tip.label, , drop = FALSE]
tip_groups <- tip_groups[tree_pruned$tip.label]
tip_families <- tip_families[tree_pruned$tip.label]

# Final check
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

# Define colors for species groups (with transparency for better visibility of overlapping points)
# Using alpha = 0.7 (70% opacity) so overlapping points appear darker
group_colors <- c(
  "Cetacean" = adjustcolor("steelblue", alpha.f = 0.7),
  "Terrestrial_Artiodactyl" = adjustcolor("darkorange", alpha.f = 0.7),
  "Other" = adjustcolor("gray50", alpha.f = 0.7)
)

# Define point shapes (pch) for families
# pch values: 21=filled circle, 22=filled square, 23=filled diamond, 24=filled triangle up, 25=filled triangle down
#             0=open square, 1=open circle, 2=open triangle up, 3=plus, 4=X, 5=diamond, 6=triangle down
family_shapes <- c(
  "Camelidae" = 21,        # Filled circle
  "Suidae" = 22,          # Filled square
  "Cervidae" = 23,        # Filled diamond
  "Bovidae" = 24,         # Filled triangle up
  "Moschidae" = 25,       # Filled triangle down
  "Hippopotamidae" = 21,  # Filled circle (different color from Cetacea)
  "Cetacea" = 21,         # Filled circle (different color from Hippopotamidae)
  "Hominidae" = 21,       # Filled circle (yellow color)
  "Unknown" = 21          # Filled circle (same as Hominidae - these are all hominidae)
)

# Create color vector for tips
tip_colors <- group_colors[tip_groups]
tip_colors[is.na(tip_colors)] <- adjustcolor("gray50", alpha.f = 0.7)  # Fallback for unclassified

# Override color for Hominidae to yellow
hominidae_mask <- tip_families == "Hominidae"
if (any(hominidae_mask)) {
  tip_colors[hominidae_mask] <- "yellow"
}

# Create shape vector for tips
tip_shapes <- family_shapes[tip_families]
tip_shapes[is.na(tip_shapes)] <- 21  # Default to circle if family not found

# Create plot
pdf(paste0(out_prefix, "_phylomorphospace_colored_", pc_x, "_", pc_y, ".pdf"), 
    width = 12, height = 10)
par(mar = c(5, 5, 3, 2))

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

# Handle zero-length branches
if (any(tree_pruned$edge.length == 0)) {
  message("  Warning: Tree has ", sum(tree_pruned$edge.length == 0), " zero-length branches")
  message("  Adding small value to zero-length branches for ancestral state calculation...")
  tree_pruned$edge.length[tree_pruned$edge.length == 0] <- min(tree_pruned$edge.length[tree_pruned$edge.length > 0]) * 0.001
}

# Calculate ancestral states
message("  Calculating ancestral states for internal nodes...")
message("  Method: Maximum Likelihood with Brownian Motion model")
anc_x <- fastAnc(tree_pruned, pc_matrix[, 1])
anc_y <- fastAnc(tree_pruned, pc_matrix[, 2])

# Check for NaN values and fix
if (any(is.nan(anc_x)) || any(is.nan(anc_y))) {
  message("  Warning: Some ancestral states are NaN, attempting to fix...")
  for (i in 1:length(anc_x)) {
    if (is.nan(anc_x[i])) {
      node_id <- as.integer(names(anc_x)[i])
      children <- tree_pruned$edge[tree_pruned$edge[, 1] == node_id, 2]
      if (length(children) > 0) {
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

# Calculate axis limits
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
     main = paste("Phylomorphospace (Colored by Species Group):", pc_x, "vs", pc_y))

# Draw tree branches (edges) in morphospace
edge <- tree_pruned$edge

# Create node trait mapping
node_traits <- list()
for (i in 1:n_tips) {
  tip_name <- tree_pruned$tip.label[i]
  node_traits[[as.character(i)]] <- c(pc_matrix[tip_name, 1], pc_matrix[tip_name, 2])
}
for (node_id_str in names(anc_x)) {
  node_id <- as.integer(node_id_str)
  if (!is.na(node_id) && node_id > n_tips) {
    node_traits[[node_id_str]] <- c(anc_x[node_id_str], anc_y[node_id_str])
  }
}

get_node_trait <- function(node_id) {
  node_id_str <- as.character(node_id)
  if (node_id_str %in% names(node_traits)) {
    return(node_traits[[node_id_str]])
  } else {
    if (node_id <= n_tips) {
      tip_name <- tree_pruned$tip.label[node_id]
      if (tip_name %in% rownames(pc_matrix)) {
        return(c(pc_matrix[tip_name, 1], pc_matrix[tip_name, 2]))
      }
    }
    return(c(NA, NA))
  }
}

edges_drawn <- 0
edges_skipped <- 0

for (i in 1:nrow(edge)) {
  parent <- edge[i, 1]
  child <- edge[i, 2]
  
  parent_trait <- get_node_trait(parent)
  child_trait <- get_node_trait(child)
  
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

# Plot tips (species) colored by group and shaped by family
# Plot each family separately to use different shapes
unique_families <- unique(tip_families)
for (fam in unique_families) {
  if (!is.na(fam) && fam != "Unknown") {
    family_mask <- tip_families == fam
    if (sum(family_mask) > 0) {
      # Get the color for this family's group
      # Special case: Hominidae gets yellow color
      if (fam == "Hominidae") {
        family_color <- "yellow"
      } else {
        family_group <- unique(tip_groups[family_mask])[1]
        family_color <- group_colors[family_group]
        if (is.na(family_color)) {
          family_color <- adjustcolor("gray50", alpha.f = 0.7)
        }
      }
      
      # Determine if shape is filled (21-25) or open (0-6)
      shape_pch <- family_shapes[fam]
      is_filled <- shape_pch >= 21 && shape_pch <= 25
      
      if (is_filled) {
        # Filled shapes use bg for fill and col for border
        points(pc_matrix[family_mask, 1], pc_matrix[family_mask, 2], 
               pch = shape_pch, 
               bg = family_color, 
               cex = 1.4,
               col = adjustcolor("black", alpha.f = 0.8),
               lwd = 1.2)
      } else {
        # Open shapes use col for the shape itself
        points(pc_matrix[family_mask, 1], pc_matrix[family_mask, 2], 
               pch = shape_pch, 
               col = family_color, 
               cex = 1.4,
               lwd = 1.5)
      }
    }
  }
}

# Add legend with both colors and shapes
# First, create a combined legend showing families with their shapes and colors
legend_families <- unique(tip_families)
legend_families <- legend_families[!is.na(legend_families) & legend_families != "Unknown"]
legend_families <- sort(legend_families)

# Create legend entries
legend_labels <- character()
legend_cols <- character()
legend_pch <- integer()
legend_bg <- character()

for (fam in legend_families) {
  # Get the group color for this family
  # Special case: Hominidae gets yellow color
  if (fam == "Hominidae") {
    fam_color <- "yellow"
  } else {
    fam_group <- unique(tip_groups[tip_families == fam])[1]
    fam_color <- if (fam_group == "Cetacean") "steelblue" else if (fam_group == "Terrestrial_Artiodactyl") "darkorange" else "gray50"
  }
  
  shape_pch <- family_shapes[fam]
  is_filled <- shape_pch >= 21 && shape_pch <= 25
  
  legend_labels <- c(legend_labels, fam)
  legend_pch <- c(legend_pch, shape_pch)
  
  if (is_filled) {
    # Filled shapes: use bg for fill, black for border
    legend_cols <- c(legend_cols, "black")
    legend_bg <- c(legend_bg, fam_color)
  } else {
    # Open shapes: use col for the shape itself, no bg
    legend_cols <- c(legend_cols, fam_color)
    legend_bg <- c(legend_bg, NA)
  }
}

# Add legend
legend("topright", 
       legend = legend_labels,
       pch = legend_pch,
       pt.bg = legend_bg,
       col = legend_cols,
       pt.cex = 1.2,
       cex = 0.8,
       title = "Family",
       bg = adjustcolor("white", alpha.f = 0.9),
       box.lwd = 1)

dev.off()

message("Saved colored phylomorphospace plot: ", paste0(out_prefix, "_phylomorphospace_colored_", pc_x, "_", pc_y, ".pdf"))

# Save matched data with species groups and families
# Map normalized IDs back to original IDs for output
original_seq_ids <- pca_filtered$seq_id[match(rownames(pc_matrix), pca_filtered$seq_id_normalized)]
output_df <- data.frame(
  sequence_id = original_seq_ids,  # Use original IDs in output
  sequence_id_normalized = rownames(pc_matrix),  # Also include normalized IDs for reference
  species_group = tip_groups,
  family = tip_families,
  pc_matrix
)
colnames(output_df)[5:6] <- c(pc_x, pc_y)

write.csv(output_df, 
          file = paste0(out_prefix, "_matched_pca_scores_colored.csv"), 
          quote = FALSE, row.names = FALSE)
message("Saved matched PCA scores with species groups: ", paste0(out_prefix, "_matched_pca_scores_colored.csv"))

message("\nDone! Output files:")
message("  - ", paste0(out_prefix, "_phylomorphospace_colored_", pc_x, "_", pc_y, ".pdf"))
message("  - ", paste0(out_prefix, "_matched_pca_scores_colored.csv"))

