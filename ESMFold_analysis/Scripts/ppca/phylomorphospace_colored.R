#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(dplyr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Default values (pPCA by default)
# Assume script is run from project root directory
default_pca_file <- "ppca/ppca_results_csv/ppca_scores.csv"
default_tree_file <- "Tree_for_Analysis/Arteodactyl_ig_hits_no_PIGR_single_domain_iqtree_ultrametric_reltime.nwk"
default_output_prefix <- "ppca/ppca_results_figures/phylomorphospace_ppca"

# Show usage if --help is requested
if (length(args) >= 1 && (args[[1]] == "--help" || args[[1]] == "-h")) {
  cat("Usage: Rscript phylomorphospace_colored.R [pca_scores_csv] [tree_file] [output_prefix] [PC_x] [PC_y]\n")
  cat("\nArguments (all optional, defaults to pPCA):\n")
  cat("  pca_scores_csv  : Path to CSV file with PCA/pPCA scores (default: ppca/ppca_results_csv/ppca_scores.csv)\n")
  cat("  tree_file       : Path to Newick tree file (default: Tree_for_Analysis/...ultrametric_reltime.nwk)\n")
  cat("  output_prefix   : Prefix for output files (default: ppca/ppca_results_figures/phylomorphospace_ppca)\n")
  cat("  PC_x            : Principal component for X-axis (default: pPC1)\n")
  cat("  PC_y            : Principal component for Y-axis (default: pPC2)\n")
  cat("\nThis script colors tips by species group:\n")
  cat("  - Cetacean: Blue\n")
  cat("  - Terrestrial Artiodactyl: Orange/Brown\n")
  cat("  - Other (Human): Green\n")
  quit(status = 0)
}

# Get arguments with defaults
pca_scores_file <- ifelse(length(args) >= 1 && args[[1]] != "", args[[1]], default_pca_file)
tree_file <- ifelse(length(args) >= 2 && args[[2]] != "", args[[2]], default_tree_file)
out_prefix <- ifelse(length(args) >= 3 && args[[3]] != "", args[[3]], default_output_prefix)
pc_x <- ifelse(length(args) >= 4 && args[[4]] != "", args[[4]], "pPC1")
pc_y <- ifelse(length(args) >= 5 && args[[5]] != "", args[[5]], "pPC2")

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

# Function to extract genus from sequence ID
get_genus <- function(seq_id) {
  if (grepl("_", seq_id)) {
    strsplit(seq_id, "_")[[1]][1]
  } else if (grepl("\\|", seq_id)) {
    strsplit(seq_id, "\\|")[[1]][1]
  } else {
    seq_id
  }
}

# Build a genus -> family lookup for this dataset.
# (Used to reproduce the exact legend/shapes from the target PDF.)
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

# Classify family from sequence ID using genus parsing.
classify_family <- function(seq_id) {
  genus <- get_genus(seq_id)
  fam <- unname(genus_to_family[genus])
  if (length(fam) == 0 || is.na(fam)) return("Unknown")
  return(as.character(fam))
}

# Classify higher-level species group to match the 3-color palette in the target PDF legend.
classify_species_group <- function(seq_id) {
  fam <- classify_family(seq_id)
  if (fam == "Cetacea") return("Cetacean")
  if (fam == "Hominidae") return("Other")
  if (fam %in% c("Bovidae", "Camelidae", "Cervidae", "Hippopotamidae", "Moschidae", "Suidae")) {
    return("Terrestrial_Artiodactyl")
  }
  return("Other")
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
  # Remove pipe separators (|) - tree uses pipes, PCA scores don't
  normalized <- gsub("\\|", "", normalized)
  # Normalize version numbers: .1 -> 1, .2 -> 2, etc. (remove dots from version numbers)
  normalized <- gsub("\\.([0-9]+)", "\\1", normalized)
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
# Try multiple possible variance file names
variance_file <- sub("_scores\\.csv$", "_eigenvalues.csv", pca_scores_file)
if (!file.exists(variance_file)) {
  variance_file <- sub("_scores\\.csv$", "_explained_variance.csv", pca_scores_file)
}
if (!file.exists(variance_file)) {
variance_file <- sub("_pca_scores\\.csv$", "_explained_variance.csv", pca_scores_file)
}
if (!file.exists(variance_file)) {
  variance_file <- sub("_ppca_scores\\.csv$", "_ppca_eigenvalues.csv", pca_scores_file)
}

if (file.exists(variance_file)) {
  variance_df <- read.csv(variance_file)
  # Check for different column name formats
  if ("Component" %in% colnames(variance_df) && "Variance_Explained" %in% colnames(variance_df)) {
    # pPCA eigenvalues format
    pc_x_idx <- which(variance_df$Component == pc_x)
    pc_y_idx <- which(variance_df$Component == pc_y)
    var_x <- if (length(pc_x_idx) > 0) variance_df$Variance_Explained[pc_x_idx] else NA
    var_y <- if (length(pc_y_idx) > 0) variance_df$Variance_Explained[pc_y_idx] else NA
  } else if ("PC" %in% colnames(variance_df) && "Explained_Variance_Ratio" %in% colnames(variance_df)) {
    # Regular PCA format
  pc_x_idx <- which(variance_df$PC == pc_x)
  pc_y_idx <- which(variance_df$PC == pc_y)
  var_x <- if (length(pc_x_idx) > 0) variance_df$Explained_Variance_Ratio[pc_x_idx] * 100 else NA
  var_y <- if (length(pc_y_idx) > 0) variance_df$Explained_Variance_Ratio[pc_y_idx] * 100 else NA
  } else {
    var_x <- NA
    var_y <- NA
  }
} else {
  var_x <- NA
  var_y <- NA
}

# Define base colors for the 3 species groups (Okabe-Ito palette, matching the target PDF legend).
group_colors <- c(
  "Cetacean" = "#0072B2",                  # Blue
  "Terrestrial_Artiodactyl" = "#E69F00", # Orange
  "Other" = "#F0E442"                    # Yellow
)

# Species group for each displayed family (drives marker fill colors in the plot and legend).
family_species_group <- c(
  "Bovidae" = "Terrestrial_Artiodactyl",
  "Camelidae" = "Terrestrial_Artiodactyl",
  "Cervidae" = "Terrestrial_Artiodactyl",
  "Cetacea" = "Cetacean",
  "Hippopotamidae" = "Terrestrial_Artiodactyl",
  "Hominidae" = "Other",
  "Moschidae" = "Terrestrial_Artiodactyl",
  "Suidae" = "Terrestrial_Artiodactyl"
)

# Family-specific marker shapes as in the target PDF legend.
# R pch mapping used by base graphics:
#   21 circle, 22 square, 23 diamond, 24 triangle up, 25 triangle down, 8 asterisk/star
family_shapes <- c(
  "Bovidae" = 24,        # triangle up
  "Camelidae" = 21,     # filled circle
  "Cervidae" = 23,      # filled diamond
  "Cetacea" = 21,       # filled circle
  "Hippopotamidae" = 8, # asterisk
  "Hominidae" = 21,     # filled circle
  "Moschidae" = 25,     # triangle down
  "Suidae" = 22          # square
)

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
     main = paste("Phylomorphospace (Colored by Species Group):", pc_x, "vs", pc_y),
     xaxt = "n", yaxt = "n")  # Suppress default axes to add custom tick marks

# Add custom tick marks with more intervals for better readability
# X-axis: Calculate nice round numbers for tick marks
x_range_span <- diff(x_range)
x_tick_interval <- if (x_range_span > 200) 50 else if (x_range_span > 100) 25 else if (x_range_span > 50) 10 else 5
x_tick_start <- floor(x_range[1] / x_tick_interval) * x_tick_interval
x_tick_end <- ceiling(x_range[2] / x_tick_interval) * x_tick_interval
x_ticks <- seq(x_tick_start, x_tick_end, by = x_tick_interval)
x_ticks <- x_ticks[x_ticks >= xlim[1] & x_ticks <= xlim[2]]  # Keep ticks within plot limits

# Y-axis: Calculate nice round numbers for tick marks
y_range_span <- diff(y_range)
y_tick_interval <- if (y_range_span > 200) 50 else if (y_range_span > 100) 25 else if (y_range_span > 50) 10 else 5
y_tick_start <- floor(y_range[1] / y_tick_interval) * y_tick_interval
y_tick_end <- ceiling(y_range[2] / y_tick_interval) * y_tick_interval
y_ticks <- seq(y_tick_start, y_tick_end, by = y_tick_interval)
y_ticks <- y_ticks[y_ticks >= ylim[1] & y_ticks <= ylim[2]]  # Keep ticks within plot limits

# Add axes with custom tick marks
axis(1, at = x_ticks, labels = x_ticks, cex.axis = 0.9, tck = -0.02)
axis(2, at = y_ticks, labels = y_ticks, cex.axis = 0.9, tck = -0.02, las = 1)

# Add minor tick marks for finer resolution
minor_x_interval <- x_tick_interval / 2
minor_x_ticks <- seq(x_tick_start, x_tick_end, by = minor_x_interval)
minor_x_ticks <- minor_x_ticks[minor_x_ticks >= xlim[1] & minor_x_ticks <= xlim[2]]
minor_x_ticks <- minor_x_ticks[!minor_x_ticks %in% x_ticks]  # Exclude major ticks

minor_y_interval <- y_tick_interval / 2
minor_y_ticks <- seq(y_tick_start, y_tick_end, by = minor_y_interval)
minor_y_ticks <- minor_y_ticks[minor_y_ticks >= ylim[1] & minor_y_ticks <= ylim[2]]
minor_y_ticks <- minor_y_ticks[!minor_y_ticks %in% y_ticks]  # Exclude major ticks

# Add minor tick marks (smaller, no labels)
axis(1, at = minor_x_ticks, labels = FALSE, tck = -0.01)
axis(2, at = minor_y_ticks, labels = FALSE, tck = -0.01)

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
# Strategy: Detect overlapping points and make them transparent
# Keep cetacean and human points on top, but make overlapping ones transparent

# Calculate pairwise distances to detect overlaps
# Use a threshold based on point size (in data units)
point_size_data_units <- diff(x_range) * 0.02  # ~2% of x-range as overlap threshold
overlap_threshold <- point_size_data_units

# Create distance matrix
n_points <- nrow(pc_matrix)
overlap_matrix <- matrix(FALSE, nrow = n_points, ncol = n_points)
rownames(overlap_matrix) <- rownames(pc_matrix)
colnames(overlap_matrix) <- rownames(pc_matrix)

for (i in 1:(n_points - 1)) {
  for (j in (i + 1):n_points) {
    dist_ij <- sqrt((pc_matrix[i, 1] - pc_matrix[j, 1])^2 + 
                    (pc_matrix[i, 2] - pc_matrix[j, 2])^2)
    if (dist_ij < overlap_threshold) {
      overlap_matrix[i, j] <- TRUE
      overlap_matrix[j, i] <- TRUE
    }
  }
}

# Identify which points overlap with others
has_overlap <- rowSums(overlap_matrix) > 0

# Define base colors and properties
base_point_size <- 1.4
base_border_lwd <- 1.2
overlap_alpha <- 0.7  # Transparency for overlapping points (increased from 0.4)
non_overlap_alpha <- 1.0  # Fully opaque for non-overlapping points

# Function to plot a set of points with appropriate transparency.
plot_points_with_transparency <- function(point_mask, family_name, is_priority = FALSE) {
  if (sum(point_mask) == 0) return()

  base_color <- group_colors[family_species_group[family_name]]
  if (is.na(base_color)) base_color <- "gray50"

  shape_pch <- family_shapes[family_name]
  if (is.na(shape_pch)) shape_pch <- 21L

  mask_overlap <- point_mask & has_overlap
  mask_non_overlap <- point_mask & !has_overlap

  # For pch 21-25, `bg` controls fill and `col` controls border.
  if (shape_pch %in% 21:25) {
    if (sum(mask_non_overlap) > 0) {
      points(pc_matrix[mask_non_overlap, 1], pc_matrix[mask_non_overlap, 2],
             pch = shape_pch,
             bg = base_color,
             cex = base_point_size,
             col = "black",
             lwd = if (is_priority) base_border_lwd * 1.3 else base_border_lwd)
    }
    if (sum(mask_overlap) > 0) {
      overlap_fill <- adjustcolor(base_color, alpha.f = overlap_alpha)
      points(pc_matrix[mask_overlap, 1], pc_matrix[mask_overlap, 2],
             pch = shape_pch,
             bg = overlap_fill,
             cex = base_point_size,
             col = "black",
             lwd = if (is_priority) base_border_lwd * 1.3 else base_border_lwd)
    }
  } else {
    # For pch 8 (asterisk), only `col` matters.
    if (sum(mask_non_overlap) > 0) {
      points(pc_matrix[mask_non_overlap, 1], pc_matrix[mask_non_overlap, 2],
             pch = shape_pch,
             cex = base_point_size,
             col = base_color,
             lwd = if (is_priority) base_border_lwd * 1.3 else base_border_lwd)
    }
    if (sum(mask_overlap) > 0) {
      overlap_col <- adjustcolor(base_color, alpha.f = overlap_alpha)
      points(pc_matrix[mask_overlap, 1], pc_matrix[mask_overlap, 2],
             pch = shape_pch,
             cex = base_point_size,
             col = overlap_col,
             lwd = if (is_priority) base_border_lwd * 1.3 else base_border_lwd)
    }
  }
}

# Plot all families; plot Cetacea and Hominidae last so they appear on top
unique_families <- unique(tip_families)
unique_families <- unique_families[!is.na(unique_families) & unique_families != "Unknown"]

priority_families <- intersect(unique_families, c("Cetacea", "Hominidae"))
other_families <- unique_families[!unique_families %in% priority_families]

# Track total points plotted for verification
total_points_plotted <- 0

for (fam in other_families) {
  family_mask <- tip_families == fam
  n_points_fam <- sum(family_mask)
  total_points_plotted <- total_points_plotted + n_points_fam
  plot_points_with_transparency(family_mask, fam, is_priority = FALSE)
}

for (fam in priority_families) {
  family_mask <- tip_families == fam
  n_points_fam <- sum(family_mask)
  total_points_plotted <- total_points_plotted + n_points_fam
  plot_points_with_transparency(family_mask, fam, is_priority = TRUE)
}

# Verify all points were plotted
message("  Total points plotted: ", total_points_plotted, " / ", nrow(pc_matrix))
if (total_points_plotted != nrow(pc_matrix)) {
  warning("  WARNING: Not all points were plotted! Missing: ", nrow(pc_matrix) - total_points_plotted)
}

# Add legend by family (matching the target PDF)
legend_families <- intersect(
  unique(tip_families[!is.na(tip_families) & tip_families != "Unknown"]),
  c("Bovidae", "Camelidae", "Cervidae", "Cetacea", "Hippopotamidae", "Hominidae", "Moschidae", "Suidae")
)

# Legend order as shown in the target PDF
legend_families <- legend_families[match(
  legend_families,
  c("Bovidae", "Camelidae", "Cervidae", "Cetacea", "Hippopotamidae", "Hominidae", "Moschidae", "Suidae")
)]

legend_labels <- legend_families
legend_pch <- as.integer(family_shapes[legend_families])
legend_fill <- group_colors[family_species_group[legend_families]]

# For pch 21-25, `col` is border and `pt.bg` is fill; for pch 8, only `col` matters.
legend_cols <- ifelse(legend_pch %in% 21:25, "black", legend_fill)
legend_bg <- legend_fill

legend("bottomright",
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

